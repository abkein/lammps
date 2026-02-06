/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// TODO: NUCC FILE

#include "compute_cluster_pe.h"
#include "compute_cluster_size_ext.h"
#include "nucc_cspan.hpp"

#include "comm.h"
#include "error.h"
#include "modify.h"
#include "update.h"

#include <cstring>
#include <unordered_map>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeClusterPE::ComputeClusterPE(LAMMPS* lmp, int narg, char** arg) : Compute(lmp, narg, arg)
{
  vector_flag     = 1;
  size_vector     = 0;
  extvector       = 0;
  local_flag      = 1;
  size_local_rows = 0;
  size_local_cols = 0;

  if (narg < 4) { utils::missing_cmd_args(FLERR, "compute pe/cluster", error); }

  // Parse arguments //

  // Get cluster/size compute
  compute_cluster_size = dynamic_cast<ComputeClusterSizeExt*>(lmp->modify->get_compute_by_id(arg[3]));
  if (compute_cluster_size == nullptr) { error->all(FLERR, "compute {}: Cannot find compute with style 'size/cluster' with id: {}", style, arg[3]); }

  size_cutoff = compute_cluster_size->get_size_cutoff();

  int iarg    = 4;
  while (iarg < narg) {
    if (::strcmp(arg[iarg], "cut") == 0) {
      const int t_size_cutoff = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (t_size_cutoff < 1) { error->all(FLERR, "size_cutoff for compute {} must be greater than 0", style); }
      if (t_size_cutoff > size_cutoff) {
        error->all(FLERR,
                   "size_cutoff for compute {} cannot be greater than it of "
                   "compute cluster/size",
                   style);
      }
      size_cutoff = std::min(size_cutoff, t_size_cutoff);
      iarg += 2;
    } else {
      error->all(FLERR, "{}: Unsupported argument: {}", style, arg[iarg]);
    }
  }

  // Get pe/atom compute
  auto computes = lmp->modify->get_compute_by_style("pe/atom");
  if (computes.empty()) { error->all(FLERR, "compute {}: Cannot find compute with style 'pe/atom'", style); }
  compute_pe_atom = computes[0];

  size_local_rows = size_cutoff + 1;
  size_vector  = size_cutoff + 1;
}

/* ---------------------------------------------------------------------- */

ComputeClusterPE::~ComputeClusterPE() noexcept(true)
{
  local_pes.destroy(memory);
  pes.destroy(memory);
}

/* ---------------------------------------------------------------------- */

void ComputeClusterPE::init()
{
  if ((modify->get_compute_by_style(style).size() > 1) && (comm->me == 0)) { error->warning(FLERR, "More than one compute {}", style); }

  local_pes.create(memory, size_local_rows, "compute:pe/cluster:local_pes");
  vector_local = local_pes.data();

  pes.create(memory, size_vector, "compute:pe/cluster:pes");
  vector = pes.data();
}

/* ---------------------------------------------------------------------- */

void ComputeClusterPE::compute_vector()
{
  if (invoked_vector == update->ntimestep) { return; }
  invoked_vector = update->ntimestep;

  compute_local();

  pes.reset();
  ::MPI_Allreduce(local_pes.data(), pes.data(), size_vector, MPI_DOUBLE, MPI_SUM, world);

  const double* const dist = compute_cluster_size->vector;
  for (int i = 0; i < size_vector; ++i) { pes[i] /= dist[i]; }
}

/* ---------------------------------------------------------------------- */

void ComputeClusterPE::compute_local()
{
  if (invoked_local == update->ntimestep) { return; }
  invoked_local = update->ntimestep;

  if (compute_cluster_size->invoked_vector != update->ntimestep) { compute_cluster_size->compute_vector(); }

  if (compute_pe_atom->invoked_peratom != update->ntimestep) { compute_pe_atom->compute_peratom(); }

  const double* const peratompes = compute_pe_atom->vector_atom;
  local_pes.reset();

  const int nclusters = dynamic_cast<ComputeClusterSizeExt*>(compute_cluster_size)->get_cluster_map().size();
  const auto& clusters = dynamic_cast<ComputeClusterSizeExt*>(compute_cluster_size)->get_clusters();
  for (int i = 0; i < nclusters; ++i) {
    const auto& clstr = clusters[i];
    const auto& atoms = clstr.atoms();
    if (clstr.g_size < size_cutoff) {
      for (int j = 0; j < clstr.l_size; ++j) {
        local_pes[clstr.g_size] += peratompes[atoms[j]];
      }
    }
  }

}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeClusterPE::memory_usage()
{
  return static_cast<double>(pes.memory_usage() + local_pes.memory_usage());
}

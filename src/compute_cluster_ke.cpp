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

#include "compute_cluster_ke.h"
#include "compute_cluster_size_ext.h"
#include "nucc_cspan.hpp"
#include "nucc_defs.hpp"

#include "comm.h"
#include "error.h"
#include "modify.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeClusterKE::ComputeClusterKE(LAMMPS* lmp, int narg, char** arg) : Compute(lmp, narg, arg)
{
  vector_flag     = 1;
  size_vector     = 0;
  extvector       = 0;
  local_flag      = 1;
  size_local_rows = 0;
  size_local_cols = 0;

  if (narg < 4) { utils::missing_cmd_args(FLERR, "compute ke/cluster", error); }

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

  // Get ke/atom compute
  auto computes = lmp->modify->get_compute_by_style("ke/atom");
  if (computes.empty()) { error->all(FLERR, "compute {}: Cannot find compute with style 'ke/atom'", style); }
  compute_ke_atom = computes[0];

  size_local_rows = size_cutoff + 1;
  size_vector     = size_cutoff + 1;
}

/* ---------------------------------------------------------------------- */

ComputeClusterKE::~ComputeClusterKE() noexcept(true)
{
  local_kes.destroy(memory);
  kes.destroy(memory);
}

/* ---------------------------------------------------------------------- */

void ComputeClusterKE::init()
{
  if ((modify->get_compute_by_style(style).size() > 1) && (comm->me == 0)) { error->warning(FLERR, "More than one compute {}", style); }

  local_kes.create(memory, size_local_rows, "compute:ke/cluster:local_kes");
  vector_local = local_kes.data();

  kes.create(memory, size_vector, "compute:ke/cluster:kes");
  vector = kes.data();
}

/* ---------------------------------------------------------------------- */

void ComputeClusterKE::compute_vector()
{
  invoked_vector = update->ntimestep;

  compute_local();

  kes.reset();
  ::MPI_Allreduce(local_kes.data(), kes.data(), size_vector, MPI_DOUBLE, MPI_SUM, world);

  const double* const dist = compute_cluster_size->vector;
  for (int i = 0; i < size_vector; ++i) {
    if (dist[i] > 0) { kes[i] /= dist[i]; }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeClusterKE::compute_local()
{
  invoked_local = update->ntimestep;

  if (compute_cluster_size->invoked_vector != update->ntimestep) { compute_cluster_size->compute_vector(); }

  if (compute_ke_atom->invoked_peratom != update->ntimestep) { compute_ke_atom->compute_peratom(); }

  const double* const peratomkes = compute_ke_atom->vector_atom;
  local_kes.reset();

  const int nclusters  = dynamic_cast<ComputeClusterSizeExt*>(compute_cluster_size)->get_cluster_map().size();
  const auto& clusters = dynamic_cast<ComputeClusterSizeExt*>(compute_cluster_size)->get_clusters();
  for (int i = 0; i < nclusters; ++i) {
    const auto& clstr = clusters[i];
    const auto& atoms = clstr.atoms();
    if (clstr.g_size < size_cutoff) {
      for (int i = 0; i < clstr.l_size; ++i) {
#ifdef __NUCC_ALGO_CHECK
        if (atoms[i] > atom->nlocal) {
          if (comm->me == 0) {
            utils::logmesg(lmp, "Cluster: {}\n", clid);
            utils::logmesg(lmp, "l_size: {}, g_size: {}\n", clstr.l_size, clstr.g_size);
            utils::logmesg(lmp, "Atoms:\n");
            for (int j = 0; j < clstr.l_size; ++j) { utils::logmesg(lmp, "{} ", atoms[j]); }
            utils::logmesg(lmp, "\n");
          }
          error->one(FLERR, "{}: {}: Atom indice exceeds nlocal", style, update->ntimestep);
        }
#endif    // __NUCC_ALGO_CHECK
        local_kes[clstr.g_size] += peratomkes[atoms[i]];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeClusterKE::memory_usage()
{
  return static_cast<double>(kes.memory_usage() + local_kes.memory_usage());
}

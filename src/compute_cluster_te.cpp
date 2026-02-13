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

#include "compute_cluster_te.h"
#include "compute_cluster_ke.h"
#include "compute_cluster_pe.h"
#include "compute_cluster_size_ext.h"
#include "nucc_cspan.hpp"

#include "comm.h"
#include "error.h"
#include "modify.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeClusterTE::ComputeClusterTE(LAMMPS* lmp, int narg, char** arg) : Compute(lmp, narg, arg)
{
  vector_flag     = 1;
  size_vector     = 0;
  extvector       = 0;
  local_flag      = 1;
  size_local_rows = 0;
  size_local_cols = 0;

  if (narg < 6) { utils::missing_cmd_args(FLERR, "compute te/cluster", error); }

  // Parse arguments //

  // Get cluster/size compute
  compute_cluster_size = dynamic_cast<ComputeClusterSizeExt*>(lmp->modify->get_compute_by_id(arg[3]));
  if (compute_cluster_size == nullptr) { error->all(FLERR, "compute {}: Cannot find compute with style 'size/cluster' with id: {}", style, arg[3]); }

  // Get source computes
  compute_cluster_pe = dynamic_cast<ComputeClusterPE*>(lmp->modify->get_compute_by_id(arg[4]));
  compute_cluster_ke = dynamic_cast<ComputeClusterKE*>(lmp->modify->get_compute_by_id(arg[5]));
  if (compute_cluster_ke == nullptr) { error->all(FLERR, "compute {}: Cannot find compute with style 'ke/cluster' with id: {}", style, arg[4]); }
  if (compute_cluster_pe == nullptr) { error->all(FLERR, "compute {}: Cannot find compute with style 'pe/cluster' with id: {}", style, arg[5]); }

  // Get the critical size
  size_cutoff = compute_cluster_size->get_size_cutoff();
  if ((narg >= 6) && (::strcmp(arg[6], "inherit") != 0)) {
    const int t_size_cutoff = utils::inumeric(FLERR, arg[6], true, lmp);
    if (t_size_cutoff < 1) { error->all(FLERR, "size_cutoff for compute {} must be greater than 0", style); }
    if (t_size_cutoff > size_cutoff) {
      error->all(FLERR,
                 "size_cutoff for compute {} cannot be greater than it of "
                 "compute size/cluster",
                 style);
    }
  }

  size_local_rows = size_cutoff + 1;
  size_vector     = size_cutoff + 1;
}

/* ---------------------------------------------------------------------- */

ComputeClusterTE::~ComputeClusterTE() noexcept(true)
{
  local_tes.destroy(memory);
  tes.destroy(memory);
}

/* ---------------------------------------------------------------------- */

void ComputeClusterTE::init()
{
  if ((modify->get_compute_by_style(style).size() > 1) && (comm->me == 0)) { error->warning(FLERR, "More than one compute {}", style); }

  local_tes.grow(memory, size_local_rows, "compute:te/cluster:local_tes");
  vector_local = local_tes.data();

  tes.grow(memory, size_vector, "compute:te/cluste:tes");
  vector = tes.data();
}

/* ---------------------------------------------------------------------- */

void ComputeClusterTE::compute_vector()
{
  if (invoked_vector == update->ntimestep) { return; }
  invoked_vector = update->ntimestep;

  compute_local();

  tes.reset();
  ::MPI_Allreduce(local_tes.data(), tes.data(), size_vector, MPI_DOUBLE, MPI_SUM, world);

  const double* const dist = compute_cluster_size->vector;
  for (int i = 0; i < size_vector; ++i) { tes[i] /= dist[i]; }
}

/* ---------------------------------------------------------------------- */

void ComputeClusterTE::compute_local()
{
  if (invoked_local == update->ntimestep) { return; }
  invoked_local = update->ntimestep;

  if (compute_cluster_ke->invoked_vector != update->ntimestep) { compute_cluster_ke->compute_vector(); }
  if (compute_cluster_pe->invoked_vector != update->ntimestep) { compute_cluster_pe->compute_vector(); }

  const double* const per_cluster_pes = compute_cluster_ke->vector_local;
  const double* const per_cluster_kes = compute_cluster_pe->vector_local;
  local_tes.reset();

  for (int i = 0; i < size_vector; i++) { local_tes[i] = per_cluster_pes[i] + per_cluster_kes[i]; }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeClusterTE::memory_usage()
{
  return static_cast<double>(tes.memory_usage() + local_tes.memory_usage());
}

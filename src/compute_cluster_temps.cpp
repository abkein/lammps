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

#include "compute_cluster_temps.h"
#include "compute_cluster_size_ext.h"

#include "comm.h"
#include "domain.h"
#include "error.h"
#include "modify.h"
#include "update.h"

#include <algorithm>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeClusterTemp::ComputeClusterTemp(LAMMPS *lmp, int narg, char **arg) : Compute(lmp, narg, arg)
{
  vector_flag = 1;
  size_vector = 0;
  extvector = 0;

  if (narg < 5) { utils::missing_cmd_args(FLERR, "compute temp/cluster", error); }

  // Parse arguments //

  // Get cluster/size compute
  compute_cluster_size = dynamic_cast<ComputeClusterSizeExt *>(lmp->modify->get_compute_by_id(arg[3]));
  if (compute_cluster_size == nullptr) {
    error->all(FLERR, "{}: Cannot find compute with style 'size/cluster' with id: {}", style,
               arg[3]);
  }

  compute_cluster_ke = lmp->modify->get_compute_by_id(arg[4]);
  if (compute_cluster_ke == nullptr) {
    error->all(FLERR, "{}: Cannot find compute with style 'ke/cluster' with id: {}", style,
               arg[4]);
  }

  // Get the critical size
  size_cutoff = compute_cluster_size->get_size_cutoff();
  if ((narg >= 5) && (::strcmp(arg[5], "inherit") != 0)) {
    const int t_size_cutoff = utils::inumeric(FLERR, arg[5], true, lmp);
    if (t_size_cutoff < 1) {
      error->all(FLERR, "size_cutoff for {} must be greater than 0: {}", style, arg[5]);
    }
    size_cutoff = std::min(size_cutoff, t_size_cutoff);
  }

  size_vector = size_cutoff + 1;
  temp.create(memory, size_vector, "temp/cluster:temp");
  vector = temp.data();
}

/* ---------------------------------------------------------------------- */

ComputeClusterTemp::~ComputeClusterTemp() noexcept(true)
{
  temp.destroy(memory);
}

/* ---------------------------------------------------------------------- */

void ComputeClusterTemp::init()
{
  if ((modify->get_compute_by_style(style).size() > 1) && (comm->me == 0)) {
    error->warning(FLERR, "More than one compute {}", style);
  }
}

/* ---------------------------------------------------------------------- */

void ComputeClusterTemp::compute_vector()
{
  if (invoked_vector == update->ntimestep) { return; }
  invoked_vector = update->ntimestep;

  if (compute_cluster_size->invoked_vector != update->ntimestep) {
    compute_cluster_size->compute_vector();
  }

  if (compute_cluster_ke->invoked_vector != update->ntimestep) {
    compute_cluster_ke->compute_vector();
  }

  temp.reset();
  const double* const kes = compute_cluster_ke->vector;
  const double* const dist = compute_cluster_size->vector;
  for (int i = 0; i < size_cutoff; ++i) {
    if (dist[i] > 0) { temp[i] = 2 * kes[i] / i / domain->dimension; }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeClusterTemp::memory_usage()
{
  return static_cast<double>(temp.memory_usage());
}

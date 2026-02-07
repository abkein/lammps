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

ComputeClusterTemp::ComputeClusterTemp(LAMMPS* lmp, int narg, char** arg) : Compute(lmp, narg, arg)
{
  vector_flag = 1;
  extvector   = 0;
  scalar_flag = 1;
  extscalar   = 0;

  if (narg < 5) { utils::missing_cmd_args(FLERR, "compute temp/cluster", error); }

  // Parse arguments //

  // Get cluster/size compute
  compute_cluster_size = dynamic_cast<ComputeClusterSizeExt*>(lmp->modify->get_compute_by_id(arg[3]));
  if (compute_cluster_size == nullptr) { error->all(FLERR, "{}: Cannot find compute with style 'size/cluster' with id: {}", style, arg[3]); }
  size_cutoff        = compute_cluster_size->get_size_cutoff();

  compute_cluster_ke = lmp->modify->get_compute_by_id(arg[4]);
  if (compute_cluster_ke == nullptr) { error->all(FLERR, "{}: Cannot find compute with style 'ke/cluster' with id: {}", style, arg[4]); }

  int iarg = 5;
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

  size_vector = size_cutoff + 1;
}

/* ---------------------------------------------------------------------- */

ComputeClusterTemp::~ComputeClusterTemp() noexcept(true)
{
  temp.destroy(memory);
}

/* ---------------------------------------------------------------------- */

void ComputeClusterTemp::init()
{
  if ((modify->get_compute_by_style(style).size() > 1) && (comm->me == 0)) { error->warning(FLERR, "More than one compute {}", style); }

  temp.create(memory, size_vector, "temp/cluster:temp");
  vector = temp.data();
}

/* ---------------------------------------------------------------------- */

void ComputeClusterTemp::compute_vector()
{
  invoked_vector = update->ntimestep;

  if (compute_cluster_ke->invoked_vector != update->ntimestep) { compute_cluster_ke->compute_vector(); }

  temp.reset();
  const double* const kes = compute_cluster_ke->vector;
  for (int i = 0; i < size_cutoff; ++i) { temp[i] = 2 * kes[i] / i / domain->dimension; }
}

/* ---------------------------------------------------------------------- */

double ComputeClusterTemp::compute_scalar()
{
  invoked_scalar = update->ntimestep;
  double sum     = 0;
  for (int i = 0; i < size_cutoff; ++i) { sum += temp[i]; }
  return sum / size_cutoff;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeClusterTemp::memory_usage()
{
  return static_cast<double>(temp.memory_usage());
}

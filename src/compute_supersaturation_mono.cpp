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

#include "compute_supersaturation_mono.h"
#include "compute_cluster_temps.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "modify.h"
#include "region.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeSupersaturationMono::ComputeSupersaturationMono(LAMMPS* lmp, int narg, char** arg) : Compute(lmp, narg, arg)
{

  scalar_flag = 1;
  extscalar   = 0;
  local_flag  = 1;

  if (narg < 8) { utils::missing_cmd_args(FLERR, "compute supersaturation/mono", error); }

  // Parse arguments //

  // Target region
  region = domain->get_region_by_id(arg[3]);
  if (region == nullptr) { error->all(FLERR, "{}: Cannot find target region {}", style, arg[3]); }

  // Get neighs compute
  compute_neighs = lmp->modify->get_compute_by_id(arg[4]);
  if (compute_neighs == nullptr) { error->all(FLERR, "{}: Cannot find compute with style 'coord/atom' with id: {}", style, arg[4]); }

  // Arrhenius coeffs
  coeffs[0] = utils::numeric(FLERR, arg[5], true, lmp);
  coeffs[1] = utils::numeric(FLERR, arg[6], true, lmp);
  coeffs[2] = utils::numeric(FLERR, arg[7], true, lmp);

  if (narg == 9) {
    if (::strcmp(arg[8], "uset1") == 0) {
      use_t1 = true;
    } else {
      error->all(FLERR, "{}: Uknown option {}", style, arg[8]);
    }
  }

  int iarg = 8;
  while (iarg < narg) {
    if (::strcmp(arg[iarg], "uset1") == 0) {
      use_t1 = true;
      iarg += 1;
    } else {
      error->all(FLERR, "{}: Uknown option {}", style, arg[iarg]);
    }
  }

  const auto& temp_computes = lmp->modify->get_compute_by_style("temp");
  if (temp_computes.empty()) { error->all(FLERR, "{}: Cannot find compute with style 'temp'.", style); }
  compute_temp = temp_computes[0];

  if (use_t1) {
    const auto& cl_temp_computes = lmp->modify->get_compute_by_style("cluster/temp");
    if (cl_temp_computes.empty()) { error->all(FLERR, "{}: Cannot find compute with style 'cluster/temp'.", style); }
    compute_cltemp = dynamic_cast<ComputeClusterTemp*>(cl_temp_computes[0]);
  }
}

/* ---------------------------------------------------------------------- */

ComputeSupersaturationMono::~ComputeSupersaturationMono() noexcept(true)
{
  mono_idx.destroy(memory);
}

/* ---------------------------------------------------------------------- */

void ComputeSupersaturationMono::init()
{
  if ((modify->get_compute_by_style(style).size() > 1) && (comm->me == 0)) { error->warning(FLERR, "More than one compute {}", style); }

  nloc = atom->nlocal;
  mono_idx.grow(memory, nloc, "compute supersaturation/mono:mono_idx");
}

/* ---------------------------------------------------------------------- */

double ComputeSupersaturationMono::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  if (invoked_local != update->ntimestep) { compute_local(); }

  bigint _local_monomers = local_monomers;
  ::MPI_Allreduce(&_local_monomers, &global_monomers, 1, MPI_LMP_BIGINT, MPI_SUM, world);

  scalar = static_cast<double>(global_monomers) / domain->volume() / execute_func();
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeSupersaturationMono::compute_local()
{
  invoked_local = update->ntimestep;

  if (nloc < atom->nlocal) {
    nloc = atom->nlocal;
    mono_idx.grow(memory, nloc, "compute supersaturation/mono:mono_idx");
  }

  mono_idx.reset();

  region->prematch();

  local_monomers = 0;
  if (compute_neighs->invoked_peratom != update->ntimestep) { compute_neighs->compute_peratom(); }
  if (use_t1) {
    if (compute_temp->invoked_vector != update->ntimestep) { compute_temp->compute_vector(); }
  } else {
    if (compute_temp->invoked_scalar != update->ntimestep) { compute_temp->compute_scalar(); }
  }
  for (int i = 0; i < atom->nlocal; ++i) {
    if (((atom->mask[i] & groupbit) != 0) && (compute_neighs->vector_atom[i] == 0) &&
        ((region->match(atom->x[i][0], atom->x[i][1], atom->x[i][2])) != 0)) {
      ++local_monomers;
      mono_idx[local_monomers] = i;
    }
  }

  local_scalar = static_cast<double>(local_monomers) / domain->subvolume() / execute_func();
}

/* ---------------------------------------------------------------------- */

double ComputeSupersaturationMono::execute_func() const
{
  double result = 0;
  if (use_t1) {
    result = coeffs[0] * ::exp(coeffs[1] - coeffs[2] / compute_cltemp->vector[1]);
  } else {
    result = coeffs[0] * ::exp(coeffs[1] - coeffs[2] / compute_temp->scalar);
  }
  return result;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeSupersaturationMono::memory_usage()
{
  return static_cast<double>(nloc * sizeof(int));
}

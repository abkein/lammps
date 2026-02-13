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

#include "compute_supersaturation_density.h"

#include "comm.h"
#include "domain.h"
#include "error.h"
#include "modify.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeSupersaturationDensity::ComputeSupersaturationDensity(LAMMPS* lmp, int narg, char** arg) : Compute(lmp, narg, arg)
{

  scalar_flag = 1;
  extscalar   = 0;

  if (narg < 8) { utils::missing_cmd_args(FLERR, "compute supersaturation/density", error); }

  // Parse arguments //

  // Get cluster/size compute
  compute_cluster_size = lmp->modify->get_compute_by_id(arg[3]);
  if (compute_cluster_size == nullptr) {
    error->all(FLERR,
               "{}: Cannot find compute with style 'cluster/size' "
               "with id: {}",
               style, arg[3]);
  }

  // Get kmax
  kmax = utils::inumeric(FLERR, arg[4], true, lmp);
  if (kmax < 1) { error->all(FLERR, "{}: kmax cannot be less than 1", style); }

  // Arrhenius coeffs
  coeffs[0]                 = utils::numeric(FLERR, arg[5], true, lmp);
  coeffs[1]                 = utils::numeric(FLERR, arg[6], true, lmp);
  coeffs[2]                 = utils::numeric(FLERR, arg[7], true, lmp);

  const auto& temp_computes = lmp->modify->get_compute_by_style("temp");
  if (temp_computes.empty()) { error->all(FLERR, "{}: Cannot find compute with style 'temp'.", style); }
  compute_temp = temp_computes[0];
}

/* ---------------------------------------------------------------------- */

ComputeSupersaturationDensity::~ComputeSupersaturationDensity() noexcept(true) = default;

/* ---------------------------------------------------------------------- */

void ComputeSupersaturationDensity::init()
{
  if ((modify->get_compute_by_style(style).size() > 1) && (comm->me == 0)) { error->warning(FLERR, "More than one compute {}", style); }
}

/* ---------------------------------------------------------------------- */

double ComputeSupersaturationDensity::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  if (compute_cluster_size->invoked_vector != update->ntimestep) { compute_cluster_size->compute_vector(); }
  if (compute_temp->invoked_scalar != update->ntimestep) { compute_temp->compute_scalar(); }

  const double* const dist = compute_cluster_size->vector;
  double sum               = 0;
  for (int size = 1; size <= kmax; ++size) { sum += size * dist[size]; }

  scalar = sum / domain->volume() / execute_func();
  return scalar;
}

/* ---------------------------------------------------------------------- */

double ComputeSupersaturationDensity::execute_func() const
{
  return coeffs[0] * ::exp(coeffs[1] - coeffs[2] / compute_temp->scalar);
}

/* ---------------------------------------------------------------------- */

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

#include "compute_supersaturation_mono_ext.h"
#include "compute_cluster_temps.h"

#include "comm.h"
#include "domain.h"
#include "error.h"
#include "modify.h"
#include "region.h"
#include "update.h"

#include <cmath>
#include <cstring>
#include <format>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeSupersaturationMonoExt::ComputeSupersaturationMonoExt(LAMMPS* lmp, int narg, char** arg) : Compute(lmp, narg, arg)
{

  scalar_flag     = 1;
  extscalar       = 0;
  local_flag      = 1;
  size_local_cols = 0;
  size_local_rows = 4;

  if (narg < 7) { utils::missing_cmd_args(FLERR, "compute supersaturation/mono", error); }

  // Parse arguments //

  compute_cluster_size = lmp->modify->get_compute_by_id(arg[3]);
  if (compute_cluster_size == nullptr) { error->all(FLERR, "{}: Cannot find compute with style 'cluster/size/ext' with id: {}", style, arg[3]); }

  // Arrhenius coeffs
  coeffs[0] = utils::numeric(FLERR, arg[4], true, lmp);
  coeffs[1] = utils::numeric(FLERR, arg[5], true, lmp);
  coeffs[2] = utils::numeric(FLERR, arg[6], true, lmp);

  int iarg  = 7;

  while (iarg < narg) {
    if (::strcmp(arg[iarg], "temp") == 0) {
      if (iarg + 2 > narg) { utils::missing_cmd_args(FLERR, std::format("{}: temp", style), error); }
      temp_avg = ::strcmp(arg[iarg + 1], "avg") == 0;
      iarg += 2;
    } else if (::strcmp(arg[iarg], "temp_compute") == 0) {
      if (iarg + 2 > narg) { utils::missing_cmd_args(FLERR, std::format("{}: temp_compute", style), error); }
      compute_temp = lmp->modify->get_compute_by_id(arg[iarg + 1]);
      if (compute_temp == nullptr) { error->all(FLERR, "{}: Cannot find compute with id {}.", style, arg[iarg + 1]); }
      iarg += 2;
    } else {
      error->all(FLERR, "{}: Uknown option {}", style, arg[iarg]);
    }
  }

  if (temp_avg) {
    const auto& temp_computes = lmp->modify->get_compute_by_style("temp");
    if (temp_computes.empty()) { error->all(FLERR, "{}: Cannot find compute with style 'temp'.", style); }
    compute_temp = temp_computes[0];
  } else {
    if (compute_temp == nullptr) { error->all(FLERR, "{}: Compute cluster/temp is not set.", style); }
  }

  vector_local = data.data();
}

/* ---------------------------------------------------------------------- */

ComputeSupersaturationMonoExt::~ComputeSupersaturationMonoExt() noexcept(true) = default;

/* ---------------------------------------------------------------------- */

void ComputeSupersaturationMonoExt::init()
{
  if ((modify->get_compute_by_style(style).size() > 1) && (comm->me == 0)) { error->warning(FLERR, "More than one compute {}", style); }
}

/* ---------------------------------------------------------------------- */

double ComputeSupersaturationMonoExt::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  if (invoked_local != update->ntimestep) { compute_local(); }

  const bigint _local_monomers = local_monomers;
  ::MPI_Allreduce(&_local_monomers, &global_monomers, 1, MPI_LMP_BIGINT, MPI_SUM, world);

  const double mult = domain->volume() * execute_func();
  data[0]           = static_cast<double>(global_monomers);
  data[1]           = mult;
  scalar            = static_cast<double>(global_monomers) / mult;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeSupersaturationMonoExt::compute_local()
{
  invoked_local  = update->ntimestep;

  local_monomers = 0;
  if (compute_cluster_size->invoked_local != update->ntimestep) { compute_cluster_size->compute_local(); }
  if (temp_avg) {
    if (compute_temp->invoked_scalar != update->ntimestep) { compute_temp->compute_scalar(); }
  } else {
    if (compute_temp->invoked_vector != update->ntimestep) { compute_temp->compute_vector(); }
  }
  local_monomers    = compute_cluster_size->size_local_rows;

  const double mult = domain->subvolume() * execute_func();
  data[2]           = static_cast<double>(local_monomers);
  data[3]           = mult;
  local_scalar      = static_cast<double>(local_monomers) / mult;
}

/* ---------------------------------------------------------------------- */

double ComputeSupersaturationMonoExt::execute_func() const
{
  double result = 0;
  if (temp_avg) {
    result = coeffs[0] * ::exp(coeffs[1] - coeffs[2] / compute_temp->scalar);
  } else {
    result = coeffs[0] * ::exp(coeffs[1] - coeffs[2] / compute_temp->vector[1]);
  }
  return result;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeSupersaturationMonoExt::memory_usage()
{
  return static_cast<double>(0);
}

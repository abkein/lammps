/* -*- c++ -*- ----------------------------------------------------------
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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(supersaturation/mono/ext,ComputeSupersaturationMonoExt);
// clang-format on
#else

#ifndef LMP_COMPUTE_SUPERSATURATION_MONO_EXT_H
#define LMP_COMPUTE_SUPERSATURATION_MONO_EXT_H

#include "compute.h"

#include <array>

namespace LAMMPS_NS {

class ComputeSupersaturationMonoExt : public Compute {
 public:
  ComputeSupersaturationMonoExt(class LAMMPS* lmp, int narg, char** arg);
  ~ComputeSupersaturationMonoExt() noexcept(true) override;
  void init() override;
  double compute_scalar() override;
  void compute_local() override;
  double memory_usage() override;

 private:
  Compute* compute_cluster_size = nullptr;
  Compute* compute_temp         = nullptr;

  std::array<double, 3> coeffs{0, 0, 0};    // [user-defined] arrhenius coefficients
  bool temp_avg          = true;            // [user-defined] use monomer temperature instead of average

  double local_scalar    = 0;    // local supersaturation
  int local_monomers     = 0;    // number of local monomers
  bigint global_monomers = 0;    // number of global monomers

  std::array<double, 4> data{0.0, 0.0, 0.0, 0.0};

  [[nodiscard]] double execute_func() const;    // monomer number density at saturation curve
};

}    // namespace LAMMPS_NS

#endif
#endif

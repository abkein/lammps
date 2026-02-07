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
ComputeStyle(supersaturation/density,ComputeSupersaturationDensity);
// clang-format on
#else

#ifndef LMP_COMPUTE_SUPERSATURATION_DENSITY_H
#define LMP_COMPUTE_SUPERSATURATION_DENSITY_H

#include "compute.h"

#include <array>

namespace LAMMPS_NS {

class ComputeSupersaturationDensity : public Compute {
 public:
  ComputeSupersaturationDensity(class LAMMPS* lmp, int narg, char** arg);
  ~ComputeSupersaturationDensity() noexcept(true) override;
  void init() override;
  double compute_scalar() override;

 private:
  Compute* compute_cluster_size = nullptr;
  Compute* compute_temp         = nullptr;

  std::array<double, 3> coeffs{0,0,0};
  int kmax = 0;    // max cluster size considered a vapor

  [[nodiscard]] double execute_func() const;    // number density at saturation curve
};

}    // namespace LAMMPS_NS

#endif
#endif

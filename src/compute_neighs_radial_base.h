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
// // clang-format off
// ComputeStyle(neighs/radial,ComputeNeighsRadial);
// // clang-format on
#else

#ifndef COMPUTE_NEIGHS_RADIAL_BASE_H
#define COMPUTE_NEIGHS_RADIAL_BASE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeNeighsRadialBase : public Compute {
 public:
  ComputeNeighsRadialBase(LAMMPS* lmp, int narg, char** arg);
  ~ComputeNeighsRadialBase() override;
  void init() override;
  void init_list(int, class NeighList*) override;

  inline constexpr const int get_nbins() const noexcept { return nbins; }
  inline constexpr const int get_cutoff() const noexcept { return cutoff; }
  inline constexpr const int get_delta_r() const noexcept { return delta; }

 protected:
  int nmax              = 0;    // previous `atom->nmax`, number of rows allocated for `rdf`
  class NeighList* list = nullptr;
  double cutoff;    // max radius compute the correlation function to; defined[user]
  double cutsq;     // cutoff squared; computed[user]

  double delta;    // width of a sigle bin, delta r; defined[user]
  int nbins;       // number of bins; computed[user]

  double** rdf = nullptr;    // computed number of neighbors
};

}    // namespace LAMMPS_NS

#endif
#endif

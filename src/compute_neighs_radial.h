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
ComputeStyle(neighs/radial,ComputeNeighsRadial);
// clang-format on
#else

#ifndef COMPUTE_NEIGHS_RADIAL_H
#define COMPUTE_NEIGHS_RADIAL_H

#include "compute.h"
#include "nucc_defs.hpp"
#ifdef __NUCC_NEIGHS_RADIAL_PRECOMPUTE_NORM
#include "nucc_cspan.hpp"
#endif    // #ifdef __NUCC_NEIGHS_RADIAL_PRECOMPUTE_NORM

namespace LAMMPS_NS {

class ComputeClusterNeighsRadial : public Compute {
 public:
  ComputeClusterNeighsRadial(class LAMMPS*, int, char**);
  ~ComputeClusterNeighsRadial() override;
  void init() override;
  void init_list(int, class NeighList*) override;
  void compute_peratom() override;
  double memory_usage() override;

  inline constexpr const int get_nbin() const noexcept { return nbins; }

 private:
  int nmax              = 0;
  class NeighList* list = nullptr;
  double cutoff;    // max radius compute the correlation function to; defined[user]
  double cutsq;     // cutoff squared; computed[user]

  double sigma;    // width of a sigle bin, delta r; defined[user]
  int nbins;       // number of bins; computed[user]
  double norm;     // normalization constant; computed[user]

  double** rdf = nullptr;    // computed number of neighbors, semi-normalized
#ifdef __NUCC_NEIGHS_RADIAL_PRECOMPUTE_NORM
  NUCC::cspan<double> norms;    // precomputed normalization constants for each bin
#endif    // __NUCC_NEIGHS_RADIAL_PRECOMPUTE_NORM
};

}    // namespace LAMMPS_NS

#endif
#endif

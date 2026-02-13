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
ComputeStyle(neighs/radial/cluster/size,ComputeNeighsRadialClusterSize);
// clang-format on
#else

#ifndef COMPUTE_NEIGHS_RADIAL_CLUSTER_SIZE_H
#define COMPUTE_NEIGHS_RADIAL_CLUSTER_SIZE_H

#include "compute.h"
#include "nucc_cspan.hpp"

namespace LAMMPS_NS {

class ComputeNeighsRadialClusterSize : public Compute {
 public:
  ComputeNeighsRadialClusterSize(class LAMMPS*, int, char**);
  ~ComputeNeighsRadialClusterSize() override;
  void   init() override;
  void   compute_local() override;
  void   compute_array() override;
  double memory_usage() override;

 private:
  int                            do_smooth     = 0;          // whether do the smoothing or not
  int                            cutoff        = 0;          // max size; defined[internal]
  int                            nbins         = 0;          // number of bins; defined[internal]
  double                         delta         = 0;          // width of a sigle bin, delta r; defined[internal]
  double                         sigma         = 0;          // smoothing kernel width; defined[user]
  double**                       weights       = nullptr;    // weights used for smoothing the distribution
  double**                       counts        = nullptr;    // unsmoothed counts
  double**                       counts2       = nullptr;    // unsmoothed counts
  double**                       counts_global = nullptr;    // global counts
  int                            max_neigh_bin = 0;          // max neighbor bin to count contribution to smoothing from
  double                         norm          = 0;          // normalization constant; computed[internal]
  NUCC::cspan<int>               atom_counts_by_size;        // number of local atoms belonging to the size

  class ComputeClusterSizeExt*   compute_cluster_size  = nullptr;
  class ComputeNeighsRadialBase* compute_neighs_radial = nullptr;

  // this should be made conditional
  // for example via template, depending on NUCC::Defines::NEIGHS_RADIAL_PRECOMPUTE_NORM
  NUCC::cspan<double>            norms;    // precomputed normalization constants for each bin
};

}    // namespace LAMMPS_NS

#endif
#endif

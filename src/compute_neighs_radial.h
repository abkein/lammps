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

#include "compute_neighs_radial_base.h"

namespace LAMMPS_NS {

class ComputeNeighsRadial : public ComputeNeighsRadialBase {
 public:
  ComputeNeighsRadial(class LAMMPS*, int, char**);
  void compute_peratom() override;
  double memory_usage() override;
};

}    // namespace LAMMPS_NS

#endif
#endif

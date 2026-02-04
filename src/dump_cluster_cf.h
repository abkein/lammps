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

#ifdef DUMP_CLASS
// clang-format off
DumpStyle(cf/cluster,DumpClusterCF);
// clang-format on
#else

#ifndef LMP_DUMP_CLUSTER_CF_H
#define LMP_DUMP_CLUSTER_CF_H

#include "dump.h"
#include <cstdio>

namespace LAMMPS_NS {

class DumpClusterCF : public Dump {
 public:
  DumpClusterCF(LAMMPS *, int, char **);
  ~DumpClusterCF() override;

 protected:
  class Compute* compute_cluster_size = nullptr;
  class Compute *compute_cf = nullptr;

  FILE *file = nullptr;

  void init_style() override;
  void pack(tagint *ids) override;
  void write_header(bigint ndump) override {}
  void write() override;
  void write_data(int n, double *mybuf) override {};
};

}    // namespace LAMMPS_NS

#endif
#endif

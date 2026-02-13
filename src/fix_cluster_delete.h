/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org
------------------------------------------------------------------------- */

// TODO: NUCC FILE

#ifdef FIX_CLASS
// clang-format off
FixStyle(cluster/delete,FixClusterDelete);
// clang-format on
#else

#ifndef LAMMPS_FIX_CLUSTER_DELETE_H
#define LAMMPS_FIX_CLUSTER_DELETE_H

#include "fix.h"
#include "nucc_cspan.hpp"

#include <cstdio>

namespace LAMMPS_NS {
class FixClusterDelete : public Fix {
 public:
  FixClusterDelete(class LAMMPS* lmp, int narg, char** arg);
  ~FixClusterDelete() noexcept(true) override;
  void init() override;
  int setmask() override;
  void pre_exchange() override;

 protected:
  class ComputeClusterSizeExt* compute_cluster_size = nullptr;

  FILE* fp                                          = nullptr;    // file write diagnostics to
  bigint next_step                                  = 0;          // next timestep wake up at

  int nloc                                          = 0;    // number of elements allocated in arrays, ~atom->nlocal
  NUCC::cspan<int> ids_a2m;                                 // local ids of atoms to delete

  // user-defined parameters
  int screenflag = 0;    // [user-defined] whether to output info to screen
  int fileflag   = 0;    // [user-defined] whether to output info into file
  int kmax       = 0;    // [user-defined] max size of clusters

  void deleteAtoms(const int to_delete) const noexcept(true);
};

}    // namespace LAMMPS_NS

#endif
#endif

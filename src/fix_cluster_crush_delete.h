/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org
------------------------------------------------------------------------- */

// TODO: NUCC FILE

#ifdef FIX_CLASS
// clang-format off
FixStyle(cluster/crush/delete,FixClusterCrushDelete);
// clang-format on
#else

#ifndef LAMMPS_FIX_CLUSTER_CRUSH_DELETE_H
#define LAMMPS_FIX_CLUSTER_CRUSH_DELETE_H

#include "fix.h"
#include "nucc_cspan.hpp"

#include <array>
#include <cstdint>
#include <cstdio>

enum class DIST : uint8_t { DIST_UNIFORM, DIST_GAUSSIAN };

namespace LAMMPS_NS {
class FixClusterCrushDelete : public Fix {
 public:
  FixClusterCrushDelete(class LAMMPS* lmp, int narg, char** arg);
  ~FixClusterCrushDelete() noexcept(true) override;
  void init() override;
  int setmask() override;
  void pre_exchange() override;

 protected:
  // necessary things for computation

  class Region* region                              = nullptr;
  class ComputeClusterSizeExt* compute_cluster_size = nullptr;
  class Compute* compute_temp                       = nullptr;
  class Compute* compute_cluster_temp               = nullptr;
  class Compute* compute_ss_mono                    = nullptr;

  FILE* fp                                          = nullptr;    // file write diagnostics to
  bigint next_step                                  = 0;          // next timestep wake up at

  int nloc                                          = 0;    // number of elements allocated in arrays, ~atom->nlocal
  NUCC::cspan<int> ids_a2m;                                 // local ids of atoms to move
  NUCC::cspan<int> count_a2m;                               // number of atoms to move per rank [comm->nprocs]
  NUCC::cspan<int> count_c2c;                               // number of clusters to crush per rank [comm->nprocs]
  std::array<double, 6> sbonds{};
  std::array<double, 6> vels{};    // [user-defined] velocities to assign to created atoms
  std::array<double, 3> xmid{};
  int balance             = 0;

  // user-defined parameters
  int screenflag          = 0;        // [user-defined] whether to output info to screen
  int fileflag            = 1;        // [user-defined] whether to output info into file
  int scaleflag           = 0;        // [user-defined]
  int kmax                = 0;        // [user-defined] max size of clusters
  double overlap          = 0;        // [user-defined] minimum distance to other atoms from the place atom teleports to
  double overlapsq        = 0;        //
  int maxtry              = 1000;     // [user-defined] max attempts to search for a new suitable location
  int ntype               = 0;        // [user-defined] type of atoms to create
  int groupid             = 0;        // [user-defined]
  bool keep_ss            = false;    // [user-defined] whether to keep mono supersaturation
  double supersaturation  = 0;        // [user-defined] desired mono supersaturation
  int insertion_rate      = 0;        // [user-defined] allow this many insertions per period (`nevery`). Zero if unlimited.

  // velocity and coordinates
  bool assign_temperature = false;                  // [user-defined] whether temperature of created atoms should be assigned at creation
  bool temp_fix           = true;                   // [user-defined] whether temperature of created atoms is fixed or tracked
  int temp_size           = 0;                      // [user-defined] 0 if average, otherwise size
  double atom_temperature = 0;                      // [user-defined] temperature of created atoms
  class RanPark* vrandom  = nullptr;                // random generator for velocities
  double vsigma           = 0;                      // MSD for velocities
  DIST vdist              = DIST::DIST_GAUSSIAN;    // distribution type for velocities
  class RanPark* xrandom  = nullptr;                // random generator for coordinates
  double xsigma           = 0;                      // MSD for coordinates
  DIST xdist              = DIST::DIST_UNIFORM;     // distribution type for coordinates

  int varflag             = 0;
  std::array<int, 4> vars{};

  void deleteAtoms(const int atoms2move_local) const noexcept(true);
  void postDelete() const noexcept(true);

  [[nodiscard]] int add(const int to_insert) const;
  void gen_pos(std::array<double, 3>& coord) const noexcept;
  [[nodiscard]] bool vartest(const std::array<double, 3>& coord) const noexcept;
  [[nodiscard]] int check_overlap(const std::array<double, 3>& coord) const noexcept;
  void create_atom(const std::array<double, 3>& coord, bigint tag) const noexcept;
  [[nodiscard]] int placement_check_me(const std::array<double, 3>& newcoord, const double* const sublo, const double* const subhi) const;
};

}    // namespace LAMMPS_NS

#endif
#endif

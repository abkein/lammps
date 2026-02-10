/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org
------------------------------------------------------------------------- */

// TODO: NUCC FILE

#ifdef FIX_CLASS
// clang-format off
FixStyle(supersaturation/mono/ext,FixSupersaturationMonoExt);
// clang-format on
#else

#ifndef LAMMPS_FIX_SUPERSATURATION_MONO_EXT_H
#define LAMMPS_FIX_SUPERSATURATION_MONO_EXT_H

#include "fix.h"
#include "nucc_cspan.hpp"

#include <array>
#include <cstdint>
#include <cstdio>


namespace LAMMPS_NS {
class FixSupersaturationMonoExt : public Fix {
 public:
  FixSupersaturationMonoExt(class LAMMPS* lmp, int narg, char** arg);
  ~FixSupersaturationMonoExt() noexcept(true) override;
  void init() override;
  int setmask() override;
  void pre_exchange() override;

 protected:
  enum class DIST : uint8_t { DIST_UNIFORM, DIST_GAUSSIAN };

  class Region* region                = nullptr;
  class Compute* compute_cluster_size = nullptr;
  class Compute* compute_temp         = nullptr;
  class Compute* compute_ss_mono      = nullptr;

  FILE* fp                            = nullptr;    // file write diagnostics to
  bigint next_step                    = 0;          // next timestep wake up at

  int nloc                            = 0;    // number of elements allocated in arrays, ~atom->nlocal
  NUCC::cspan<int> ids_a2m;                   // local ids of atoms to move
  std::array<double, 6> sbonds{};
  std::array<double, 6> vels{};    // [user-defined] velocities to assign to created atoms
  std::array<double, 3> xmid{};

  // user-defined parameters
  int screenflag          = 0;       // [user-defined] whether to output info to screen
  int fileflag            = 1;       // [user-defined] whether to output info into file
  int scaleflag           = 0;       // [user-defined]
  double overlap          = 0;       // [user-defined] minimum distance to other atoms from the place atom teleports to
  double overlapsq        = 0;       //
  int maxtry              = 1000;    // [user-defined] max attempts to search for a new suitable location
  int ntype               = 0;       // [user-defined] type of atoms to create
  int groupid             = 0;       // [user-defined]
  double supersaturation  = 0;       // [user-defined] desired mono supersaturation
  int rate_limit          = 0;       // [user-defined] allow this many deleteions/insertions per period (`nevery`). Zero if unlimited.

  // velocity and coordinates
  bool assign_temperature = false;                  // [user-defined] whether temperature of created atoms should be assigned at creation
  bool temp_fix           = true;                   // [user-defined] whether temperature of created atoms is fixed or tracked
  bool temp_use_mono      = false;                  // [user-defined] 0 if average, otherwise size
  double atom_temperature = 0;                      // [user-defined] temperature of created atoms
  class RanPark* vrandom  = nullptr;                // random generator for velocities
  double vsigma           = 0;                      // MSD for velocities
  DIST vdist              = DIST::DIST_GAUSSIAN;    // distribution type for velocities
  class RanPark* xrandom  = nullptr;                // random generator for coordinates
  double xsigma           = 0;                      // MSD for coordinates
  DIST xdist              = DIST::DIST_UNIFORM;     // distribution type for coordinates

  int varflag             = 0;
  std::array<int, 4> vars{};

  void deleteAtoms(const int to_delete) const noexcept(true);    // NOLINT(modernize-use-nodiscard)
  void postDelete() const noexcept(true);

  int add(const int to_insert) const;    // NOLINT(modernize-use-nodiscard)
  void gen_pos(std::array<double, 3>& coord) const noexcept;
  [[nodiscard]] bool vartest(const std::array<double, 3>& coord) const noexcept;
  [[nodiscard]] int check_overlap(const std::array<double, 3>& coord) const noexcept;
  void create_atom(const std::array<double, 3>& coord, bigint tag) const noexcept;
  [[nodiscard]] int placement_check_me(const std::array<double, 3>& newcoord, const double* const sublo, const double* const subhi) const;
};

}    // namespace LAMMPS_NS

#endif
#endif

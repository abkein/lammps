/*
 ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org
------------------------------------------------------------------------- */

// TODO: NUCC FILE

#include "fix_capture.h"
#include "compute.h"
#include "random_park.h"
#include "region.h"

#include "atom.h"
#include "atom_vec.h"
#include "atom_vec_body.h"
#include "atom_vec_ellipsoid.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixCapture::FixCapture(LAMMPS* lmp, int narg, char** arg) : Fix(lmp, narg, arg), action(ACTION::COUNT), screenflag(true), fileflag(false), allow(true)
{

  restart_pbc = 1;

  if (narg < 5) { utils::missing_cmd_args(FLERR, "fix capture", error); }

  // Parse arguments //

  // Target region
  region = domain->get_region_by_id(arg[3]);
  if (region == nullptr) { error->all(FLERR, "{}: Cannot find target region {}", style, arg[3]); }

  // Get number of sigmas
  nsigma   = utils::inumeric(FLERR, arg[4], true, lmp);

  int iarg = 5;
  fp       = nullptr;

  while (iarg < narg) {
    if (::strcmp(arg[iarg], "action") == 0) {

      if (::strcmp(arg[iarg + 1], "count") == 0) {

        action = ACTION::COUNT;
        iarg += 2;

      } else if (::strcmp(arg[iarg + 1], "delete") == 0) {

        action = ACTION::DELETE;
        iarg += 2;

      } else if (::strcmp(arg[iarg + 1], "cool") == 0) {

        action          = ACTION::SLOW;
        // Get the seed for velocity generator
        int const vseed = utils::inumeric(FLERR, arg[iarg + 2], true, lmp);
        vrandom         = new RanPark(lmp, vseed);
        iarg += 3;

      } else {
        error->all(FLERR, "{}: Unknown mode: {}", style, arg[iarg + 1]);
      }

    } else if (::strcmp(arg[iarg], "noscreen") == 0) {

      // Do not output to screen
      screenflag = false;
      iarg += 1;

    } else if (::strcmp(arg[iarg], "file") == 0) {

      // Write output to new file
      if (comm->me == 0) {
        fileflag = true;
        fp       = ::fopen(arg[iarg + 1], "w");
        if (fp == nullptr) { error->one(FLERR, "{}: Cannot open stats file {}: {}", style, arg[iarg + 1], utils::getsyserror()); }
      }
      iarg += 2;

    } else if (::strcmp(arg[iarg], "append") == 0) {

      // Append output to file
      if (comm->me == 0) {
        fileflag = true;
        fp       = ::fopen(arg[iarg + 1], "a");
        if (fp == nullptr) { error->one(FLERR, "{}: Cannot open stats file {}: {}", style, arg[iarg + 1], utils::getsyserror()); }
      }
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix capture command");
    }
  }

  // Get temp compute
  auto temp_computes = lmp->modify->get_compute_by_style("temp");
  if (temp_computes.empty()) { error->all(FLERR, "fix capture: Cannot find compute with style 'temp'."); }
  compute_temp = temp_computes[0];

  if (fileflag && (comm->me == 0)) {
    fmt::print(fp, "ntimestep,n\n");
    ::fflush(fp);
  }

  memory->create(rmins, atom->ntypes + 1, atom->ntypes + 1, "fix_capture:rmins");
  memory->create(vmax_coeffs, atom->ntypes + 1, "fix_capture:vmax_coeffs");
  memory->create(sigmas, atom->ntypes + 1, "fix_capture:sigmas");

}

/* ---------------------------------------------------------------------- */

FixCapture::~FixCapture() noexcept(true)
{
  delete vrandom;
  if (rmins != nullptr) { memory->destroy(rmins); }
  if (vmax_coeffs != nullptr) { memory->destroy(vmax_coeffs); }
  if (sigmas != nullptr) { memory->destroy(sigmas); }

  if ((fp != nullptr) && (comm->me == 0)) {
    ::fflush(fp);
    ::fclose(fp);
  }
}

/* ---------------------------------------------------------------------- */

int FixCapture::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCapture::init()
{
  if ((modify->get_fix_by_style(style).size() > 1) && (comm->me == 0)) { error->warning(FLERR, "More than one fix {}", style); }

  for (int i = 1; i <= atom->ntypes; ++i) {
    if (atom->mass_setflag[i] == 0) { error->all(FLERR, "fix capture: mass is not set for atom type {}.", i); }
  }

  constexpr long double eight_over_pi_sqrt = 1.5957691216057307117597842397375L;    // sqrt(8/pi)
  constexpr long double ssdd               = 0.6734396116428514837424685996751L;    // sqrt(3-8/pi)
  for (int i = 1; i <= atom->ntypes; ++i) {
    const double m = atom->mass[i];
    vmax_coeffs[i] = ::pow(eight_over_pi_sqrt * ::sqrt(1.0 / m) + nsigma * ssdd * ::sqrt(1.0 / m), 2);
  }
}

/* ---------------------------------------------------------------------- */

[[gnu::hot]] void FixCapture::pre_exchange()
{
  if (compute_temp->invoked_scalar != update->ntimestep) { compute_temp->compute_scalar(); }

  ::memset(Fcounts, 0, FIX_CAPTURE_FLAGS_COUNT * sizeof(bigint));

  for (int i = 1; i <= atom->ntypes; ++i) {
    sigmas[i] = ::sqrt(compute_temp->scalar / atom->mass[i]);
    for (int j = 0; j <= atom->ntypes; ++j) { rmins[i][j] = rmin(i, j > 0 ? j : i); }
  }

  group->vcm(igroup, group->mass(igroup), vmean);

  // region->prematch();

  for (int i = 0; i < atom->nlocal; ++i) {
    allow = true;
    if (((atom->mask[i] & groupbit) != 0)) {
      // && (region->match(atom->x[i][0], atom->x[i][1], atom->x[i][2]) != 0)
      test_xnonnum(i);
      test_vnonnum(i);
      test_superspeed<false>(&i);
      test_superspeed<true>(&i);
      test_overlap(i);
    }
  }

  if ((action == ACTION::DELETE) && (Fcounts[FIX_CAPTURE_OVERSPEED_FLAG] > 0)) { post_delete(); }

  bigint GFcounts[FIX_CAPTURE_FLAGS_COUNT]{};
  ::MPI_Allreduce(&Fcounts, &GFcounts, FIX_CAPTURE_FLAGS_COUNT, MPI_LMP_BIGINT, MPI_SUM, world);

  if (fileflag && (comm->me == 0)) {
    fmt::print(fp, "{}", update->ntimestep);
    for (const bigint GFcount : GFcounts) { fmt::print(fp, ",{}", GFcount); }
    fmt::print(fp, "\n");
    ::fflush(fp);
  }
}

/* ---------------------------------------------------------------------- */

[[gnu::hot]] inline long double FixCapture::rmin(const int i, const int j) noexcept(true)
{
  constexpr long double two_sqrt      = 1.4142135623730950488016887242096L;    // sqrt(2)
  constexpr long double two_one_third = 1.2599210498948731647672106072782L;    // 2^(1/3)
  constexpr long double one_over_six  = 0.1666666666666666666666666666666L;    // 1/6
  return two_one_third /
      ::pow((2 + two_sqrt * ::sqrt(2 + (atom->mass[i] * vmax_coeffs[i] + atom->mass[j] * vmax_coeffs[j]) * compute_temp->scalar)), one_over_six);
}

/* ---------------------------------------------------------------------- */

[[gnu::hot]] inline bool FixCapture::isnonnumeric(const double* const vec3) noexcept(true)
{
  return std::isnan(vec3[0]) || std::isnan(vec3[1]) || std::isnan(vec3[2]) || std::isinf(vec3[0]) || std::isinf(vec3[1]) || std::isinf(vec3[2]);
}

/* ---------------------------------------------------------------------- */

template <bool ISREL>
[[gnu::hot]] inline bool FixCapture::super_cond(const double* const v, const double* const vmeanx, const double sigma) const noexcept(true)
{
  if constexpr (ISREL) {
    return ::fabs(*v - *vmeanx) > nsigma * sigma;
  } else {
    return ::fabs(*v) > nsigma * sigma;
  }
}

template <bool ISREL>
[[gnu::hot]] void FixCapture::test_superspeed(int* atomid) noexcept(true)
{
  const int i        = *atomid;
  double* v          = atom->v[i];
  const double sigma = sigmas[atom->type[i]];
  bool flag          = false;
  double old_v[3]    = {v[0], v[1], v[2]};
  for (int j = 0; j < 3; ++j) {
    if (super_cond<ISREL>(v + j, vmean + j, sigma)) {
      flag = true;
      if (action == ACTION::SLOW) {
        if constexpr (ISREL) {
          v[j] = vrandom->gaussian() * sigma + vmean[j];
        } else {
          v[j] = vrandom->gaussian() * sigma;
        }
      } else if (action == ACTION::DELETE) {
        atom->avec->copy(atom->nlocal - 1, i, 1);
        --atom->nlocal;
        --(*atomid);
        break;
      } else {
        // for linter
      }
    }
  }
  if (flag) {
    if constexpr (ISREL) {
      ++Fcounts[FIX_CAPTURE_OVERSPEED_REL_FLAG];
    } else {
      ++Fcounts[FIX_CAPTURE_OVERSPEED_FLAG];
    }

    captured();
    utils::logmesg(lmp, "{}, {}####### SUPER ATOM VELOCITY ({}) ######{}\n", atom->tag[i], update->ntimestep, ISREL ? "REL" : "NOREL", comm->me);
    utils::logmesg(lmp, "{}           POS: {:.3f} {:.3f} {:.3f}\n", atom->tag[i], atom->x[i][0], atom->x[i][1], atom->x[i][2]);
    if constexpr (ISREL) {
      utils::logmesg(lmp, "{}      VELOCITY: {:.3f} {:.3f} {:.3f}\n", atom->tag[i], old_v[0], old_v[1], old_v[2]);
      old_v[0] -= vmean[0];
      old_v[1] -= vmean[1];
      old_v[2] -= vmean[2];
      double vtot = ::sqrt(old_v[0] * old_v[0] + old_v[1] * old_v[1] + old_v[2] * old_v[2]);
      utils::logmesg(lmp, "{}VELOCITY NOREL: {:.3f} {:.3f} {:.3f} ({:.3f} > {} * {:.3f})\n", atom->tag[i], old_v[0], old_v[1], old_v[2], vtot, nsigma,
                     sigma);
    } else {
      double vtot = ::sqrt(old_v[0] * old_v[0] + old_v[1] * old_v[1] + old_v[2] * old_v[2]);
      utils::logmesg(lmp, "{}      VELOCITY: {:.3f} {:.3f} {:.3f} ({:.3f} > {} * {:.3f})\n", atom->tag[i], old_v[0], old_v[1], old_v[2], vtot, nsigma,
                     sigma);
      old_v[0] -= vmean[0];
      old_v[1] -= vmean[1];
      old_v[2] -= vmean[2];
      utils::logmesg(lmp, "{}VELOCITY NOREL: {:.3f} {:.3f} {:.3f}\n", atom->tag[i], old_v[0], old_v[1], old_v[2]);
    }
    if (action == ACTION::DELETE) {
      utils::logmesg(lmp, "{} DELETED\n", atom->tag[i]);
    } else {
      utils::logmesg(lmp, "{}  NEW VELOCITY: {:.3f} {:.3f} {:.3f}\n", atom->tag[i], v[0], v[1], v[2]);
    }
  }
}

/* ---------------------------------------------------------------------- */

[[gnu::hot]] inline void FixCapture::captured() noexcept(true)
{
  if (allow) {
    ++Fcounts[FIX_CAPTURE_TOTAL_FLAG];
    allow = false;
  }
}

/* ---------------------------------------------------------------------- */

[[gnu::hot]] void FixCapture::test_xnonnum(int i) noexcept(true)
{
  if (isnonnumeric(atom->x[i])) {
    ++Fcounts[FIX_CAPTURE_xNaN_FLAG];
    captured();
    utils::logmesg(lmp, "{}, {}####### NON-NUMERIC ATOM COORDS ######{}\n", atom->tag[i], update->ntimestep, comm->me);
    utils::logmesg(lmp, "{}     POS: {:.3f} {:.3f} {:.3f}\n", atom->tag[i], atom->x[i][0], atom->x[i][1], atom->x[i][2]);
    utils::logmesg(lmp, "{}VELOCITY: {:.3f} {:.3f} {:.3f}\n", atom->tag[i], atom->v[i][0], atom->v[i][1], atom->v[i][2]);
  }
}

/* ---------------------------------------------------------------------- */

[[gnu::hot]] void FixCapture::test_vnonnum(int i) noexcept(true)
{
  if (isnonnumeric(atom->v[i])) {
    ++Fcounts[FIX_CAPTURE_vNaN_FLAG];
    captured();
    utils::logmesg(lmp, "{}, {}####### NON-NUMERIC ATOM VELOCITY ######{}\n", atom->tag[i], update->ntimestep, comm->me);
    utils::logmesg(lmp, "{}     POS: {:.3f} {:.3f} {:.3f}\n", atom->tag[i], atom->x[i][0], atom->x[i][1], atom->x[i][2]);
    utils::logmesg(lmp, "{}VELOCITY: {:.3f} {:.3f} {:.3f}\n", atom->tag[i], atom->v[i][0], atom->v[i][1], atom->v[i][2]);
  }
}

/* ---------------------------------------------------------------------- */

[[gnu::hot]] void FixCapture::test_overlap(int i) noexcept(true)
{
  double** x = atom->x;
  // for (int j = i + 1; j < atom->nmax; ++j) {
  for (int j = i + 1; j < atom->nlocal; ++j) {
    const double dx = x[i][0] - x[j][0];
    const double dy = x[i][1] - x[j][1];
    const double dz = x[i][2] - x[j][2];
    const double rm = rmins[atom->type[i]][j < atom->nlocal ? atom->type[j] : 0];
    if (dx * dx + dy * dy + dz * dz < rm * rm) {
      ++Fcounts[FIX_CAPTURE_OVERLAP_FLAG];
      captured();
      utils::logmesg(lmp, "{}, {}####### ATOM OVERLAP {} with {} ({})######{}\n", atom->tag[i], update->ntimestep, atom->tag[i], atom->tag[j],
                     j < atom->nlocal ? "owned" : "not owned", comm->me);
      utils::logmesg(lmp, "{} {}      POS: {:.3f} {:.3f} {:.3f}\n", atom->tag[i], atom->tag[i], atom->x[i][0], atom->x[i][1], atom->x[i][2]);
      utils::logmesg(lmp, "{} {} VELOCITY: {:.3f} {:.3f} {:.3f}\n", atom->tag[i], atom->tag[i], atom->v[i][0], atom->v[i][1], atom->v[i][2]);
      utils::logmesg(lmp, "{} {}      POS: {:.3f} {:.3f} {:.3f}\n", atom->tag[i], atom->tag[j], atom->x[j][0], atom->x[j][1], atom->x[j][2]);
      utils::logmesg(lmp, "{} {} VELOCITY: {:.3f} {:.3f} {:.3f}\n", atom->tag[i], atom->tag[j], atom->v[j][0], atom->v[j][1], atom->v[j][2]);
      utils::logmesg(lmp, "{}       DELTA: {:.3f} {:.3f} {:.3f}\n", atom->tag[i], dx, dy, dz);
      utils::logmesg(lmp, "{}    DISTANCE: {:.3f}  Rmin: {:.3f}\n", atom->tag[i], ::sqrt(dx * dx + dy * dy + dz * dz), rm);
      break;
    }
  }
}

/* ---------------------------------------------------------------------- */

[[gnu::cold]] void FixCapture::post_delete() noexcept(true)
{
  if (atom->molecular == Atom::ATOMIC) {
    tagint* tag      = atom->tag;
    int const nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; ++i) { tag[i] = 0; }
    atom->tag_extend();
  }

  // reset atom->natoms and also topology counts

  bigint nblocal = atom->nlocal;
  ::MPI_Allreduce(&nblocal, &atom->natoms, 1, MPI_LMP_BIGINT, MPI_SUM, world);

  // reset bonus data counts

  const auto* avec_ellipsoid = dynamic_cast<AtomVecEllipsoid*>(atom->style_match("ellipsoid"));
  const auto* avec_line      = dynamic_cast<AtomVecLine*>(atom->style_match("line"));
  const auto* avec_tri       = dynamic_cast<AtomVecTri*>(atom->style_match("tri"));
  const auto* avec_body      = dynamic_cast<AtomVecBody*>(atom->style_match("body"));
  bigint nlocal_bonus        = 0;

  if (atom->nellipsoids > 0) {
    nlocal_bonus = avec_ellipsoid->nlocal_bonus;
    ::MPI_Allreduce(&nlocal_bonus, &atom->nellipsoids, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  }
  if (atom->nlines > 0) {
    nlocal_bonus = avec_line->nlocal_bonus;
    ::MPI_Allreduce(&nlocal_bonus, &atom->nlines, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  }
  if (atom->ntris > 0) {
    nlocal_bonus = avec_tri->nlocal_bonus;
    ::MPI_Allreduce(&nlocal_bonus, &atom->ntris, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  }
  if (atom->nbodies > 0) {
    nlocal_bonus = avec_body->nlocal_bonus;
    ::MPI_Allreduce(&nlocal_bonus, &atom->nbodies, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  }

  // reset atom->map if it exists
  // set nghost to 0 so old ghosts of deleted atoms won't be mapped

  if (atom->map_style != Atom::MAP_NONE) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }
}

/* ---------------------------------------------------------------------- */

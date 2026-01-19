/*
 ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org
------------------------------------------------------------------------- */

// TODO: NUCC FILE

#include "fix_capture_vel.h"
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
#include "neigh_list.h"
#include "neighbor.h"

#include <cmath>
#include <cstring>
#include <unordered_set>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixCaptureVel::FixCaptureVel(LAMMPS* lmp, int narg, char** arg) : Fix(lmp, narg, arg), screenflag(true), fileflag(false)
{
  if (narg < 6) { utils::missing_cmd_args(FLERR, "fix capture/vel", error); }
  nevery = 1;

  // Parse arguments //

  // Target region
  region = domain->get_region_by_id(arg[3]);
  if (region == nullptr) { error->all(FLERR, "{}: Cannot find target region {}", style, arg[3]); }

  // Get number of sigmas
  nsigma   = utils::inumeric(FLERR, arg[4], true, lmp);
  nsigmasq = nsigma * nsigma;

  const int vseed = utils::inumeric(FLERR, arg[5], true, lmp);
  vrandom         = new RanPark(lmp, vseed);

  int iarg = 6;
  fp       = nullptr;

  while (iarg < narg) {
    if (::strcmp(arg[iarg], "noscreen") == 0) {

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

    } else if (::strcmp(arg[iarg], "delete_overlap") == 0) {
      delete_overlap = true;
      iarg += 1;
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

  ncaptured[0] = 0;
  ncaptured[1] = 0;

  memory->create(rmins, atom->ntypes + 1, atom->ntypes + 1, "fix_capture:rmins");
  memory->create(vmax_coeffs, atom->ntypes + 1, "fix_capture:vmax_coeffs");
  memory->create(sigmas, atom->ntypes + 1, "fix_capture:sigmas");
}

/* ---------------------------------------------------------------------- */

void FixCaptureVel::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

FixCaptureVel::~FixCaptureVel() noexcept(true)
{
  delete vrandom;
  memory->destroy(rmins);
  memory->destroy(vmax_coeffs);
  memory->destroy(sigmas);

  if ((fp != nullptr) && (comm->me == 0)) {
    ::fflush(fp);
    ::fclose(fp);
  }
}

/* ---------------------------------------------------------------------- */

int FixCaptureVel::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= PRE_FORCE;
  mask |= END_OF_STEP;
  // mask |= PRE_EXCHANGE;
  mask |= POST_NEIGHBOR;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCaptureVel::init()
{
  if ((modify->get_fix_by_style(style).size() > 1) && (comm->me == 0)) { error->warning(FLERR, "More than one fix {}", style); }

  for (int i = 1; i <= atom->ntypes; ++i) {
    if (atom->mass_setflag[i] == 0) { error->all(FLERR, "{}: mass is not set for atom type {}.", style, i); }
  }

  constexpr long double eight_over_pi_sqrt = 1.5957691216057307117597842397375L;    // sqrt(8/pi)
  constexpr long double ssdd               = 0.6734396116428514837424685996751L;    // sqrt(3-8/pi)
  for (int i = 1; i <= atom->ntypes; ++i) {
    const double m = atom->mass[i];
    vmax_coeffs[i] = ::pow(eight_over_pi_sqrt * ::sqrt(1.0 / m) + nsigma * ssdd * ::sqrt(1.0 / m), 2);
  }

  // Request a full neighbor list
  int list_flags = NeighConst::REQ_FULL;
  // need neighbors of the ghost atoms
  list_flags |= NeighConst::REQ_GHOST;
  neighbor->add_request(this, list_flags);

  to_delete.reserve(5);
}

/* ---------------------------------------------------------------------- */

void FixCaptureVel::pre_force(int vflag)
{
  initial_integrate(vflag);
}

/* ---------------------------------------------------------------------- */

// void FixCaptureVel::pre_exchange()
// {
//   std::unordered_map<bigint, bigint> ttm;
//   for (const auto[k, v] : flags) { if (update->ntimestep - v.second < 500) { ttm.emplace(k, v.second); } }
//   flags.clear();
//   // for (int i = 0; i < atom->nlocal; ++i) { flags.insert({atom->tag[i], {false, 0}}); }
//   for (int i = 0; i < atom->nlocal; ++i) {
//     if (ttm.count(atom->tag[i]) > 0) {
//       flags.insert({atom->tag[i], {true, ttm[atom->tag[i]]}});
//     } else {
//       flags.insert({atom->tag[i], {false, 0}});
//     }
//   }
// }

void FixCaptureVel::post_neighbor()
{
  std::unordered_map<bigint, bigint> ttm;
  for (const auto[k, v] : flags) { if (update->ntimestep - v.second < 500) { ttm.emplace(k, v.second); } }
  flags.clear();
  // for (int i = 0; i < atom->nlocal; ++i) { flags.insert({atom->tag[i], {false, 0}}); }
  for (int i = 0; i < atom->nlocal; ++i) {
    if (ttm.count(atom->tag[i]) > 0) {
      flags.insert({atom->tag[i], {true, ttm[atom->tag[i]]}});
    } else {
      flags.insert({atom->tag[i], {false, 0}});
    }
  }

  if (compute_temp->invoked_scalar != update->ntimestep) { compute_temp->compute_scalar(); }

  for (int i = 1; i <= atom->ntypes; ++i) {
    sigmas[i] = nsigmasq * compute_temp->scalar / atom->mass[i];
    for (int j = 0; j <= atom->ntypes; ++j) { rmins[i][j] = rminsq(i, j > 0 ? j : i); }
  }

  region->prematch();
  double **v = atom->v;
  double **x = atom->x;

  int   inum = list->inum + list->gnum;
  int*  ilist = list->ilist;
  int*  numneigh = list->numneigh;
  int** firstneigh = list->firstneigh;

  int *mask = atom->mask;

  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    if (i >= atom->nlocal) { continue; }
    if (((mask[i] & groupbit) != 0) && (region->match(atom->x[i][0], atom->x[i][1], atom->x[i][2]) != 0)) {
      const double vx = v[i][0];
      const double vy = v[i][1];
      const double vz = v[i][2];
      const double vm = vx*vx + vy*vy + vz*vz;
      if (vm > sigmas[atom->type[i]]) {
        ++ncaptured[0];
        flags[atom->tag[i]] = {true, update->ntimestep};
        const double sigma = ::sqrt(sigmas[atom->type[i]] / nsigmasq);
        v[i][0] *= sigma / v[i][0];
        v[i][1] *= sigma / v[i][1];
        v[i][2] *= sigma / v[i][2];
      }
      if (delete_overlap) {
        const double xtmp = x[i][0];
        const double ytmp = x[i][1];
        const double ztmp = x[i][2];
        int* jlist = firstneigh[i];
        int jnum = numneigh[i];
        for (int jj = 0; jj < jnum; jj++) {
          int j = jlist[jj];
          j &= NEIGHMASK;

          const double delx = xtmp - x[j][0];
          const double dely = ytmp - x[j][1];
          const double delz = ztmp - x[j][2];
          if (delx*delx + dely*dely + delz*delz < rmins[atom->type[i]][atom->type[j]]) {
            to_delete.emplace_back(i);
            break;
          }
        }
      }
    }
  }

  ncaptured[1] = to_delete.size();
  to_delete.clear();

  ncaptured_global[0] = 0;
  ncaptured_global[1] = 0;
  ::MPI_Allreduce(&ncaptured, &ncaptured_global, 2, MPI_LMP_BIGINT, MPI_SUM, world);

  if (fileflag && ((ncaptured_global[0] > 0) || (ncaptured_global[1] > 0)) && (comm->me == 0)) {
    fmt::print(fp, "{},{},{}\n", update->ntimestep, ncaptured_global[0], ncaptured_global[1]);
    ::fflush(fp);
  }

  ncaptured[0] = 0;
  ncaptured[1] = 0;
}

/* ---------------------------------------------------------------------- */

[[gnu::hot]] void FixCaptureVel::initial_integrate(int /*vflag*/)
{
  // flags.clear();
  // for (int i = 0; i < atom->nlocal; ++i) { flags.insert({atom->tag[i], false}); }

  if (compute_temp->invoked_scalar != update->ntimestep) { compute_temp->compute_scalar(); }

  for (int i = 1; i <= atom->ntypes; ++i) {
    sigmas[i] = nsigmasq * compute_temp->scalar / atom->mass[i];
    for (int j = 0; j <= atom->ntypes; ++j) { rmins[i][j] = rminsq(i, j > 0 ? j : i); }
  }

  region->prematch();
  double **v = atom->v;
  double **x = atom->x;

  int   inum = list->inum + list->gnum;
  int*  ilist = list->ilist;
  int*  numneigh = list->numneigh;
  int** firstneigh = list->firstneigh;

  int *mask = atom->mask;

  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    if (i >= atom->nlocal) { continue; }
    if (((mask[i] & groupbit) != 0) && (region->match(atom->x[i][0], atom->x[i][1], atom->x[i][2]) != 0)) {
      const double vx = v[i][0];
      const double vy = v[i][1];
      const double vz = v[i][2];
      const double vm = vx*vx + vy*vy + vz*vz;
      if (vm > sigmas[atom->type[i]]) {
        ++ncaptured[0];
        flags[atom->tag[i]] = {true, update->ntimestep};
        const double sigma = ::sqrt(sigmas[atom->type[i]] / nsigmasq);
        v[i][0] *= sigma / v[i][0];
        v[i][1] *= sigma / v[i][1];
        v[i][2] *= sigma / v[i][2];
      }
      if (delete_overlap) {
        const double xtmp = x[i][0];
        const double ytmp = x[i][1];
        const double ztmp = x[i][2];
        int* jlist = firstneigh[i];
        int jnum = numneigh[i];
        for (int jj = 0; jj < jnum; jj++) {
          int j = jlist[jj];
          j &= NEIGHMASK;

          const double delx = xtmp - x[j][0];
          const double dely = ytmp - x[j][1];
          const double delz = ztmp - x[j][2];
          if (delx*delx + dely*dely + delz*delz < rmins[atom->type[i]][atom->type[j]]) {
            to_delete.emplace_back(i);
            break;
          }
        }
      }
    }
  }

  if (to_delete.size() > 0) {
    for (int i : to_delete) {
      atom->avec->copy(atom->nlocal - 1, i, 1);
      --atom->nlocal;
    }
    ncaptured[1] += to_delete.size();
    to_delete.clear();
  }

  ncaptured_global[0] = 0;
  ncaptured_global[1] = 0;
  ::MPI_Allreduce(&ncaptured, &ncaptured_global, 2, MPI_LMP_BIGINT, MPI_SUM, world);

  if (ncaptured_global[1] > 0) { post_delete(); }

  if (fileflag && ((ncaptured_global[0] > 0) || (ncaptured_global[1] > 0)) && (comm->me == 0)) {
    fmt::print(fp, "{},{},{}\n", update->ntimestep, ncaptured_global[0], ncaptured_global[1]);
    ::fflush(fp);
  }

  ncaptured[0] = 0;
  ncaptured[1] = 0;
}

/* ---------------------------------------------------------------------- */

void FixCaptureVel::end_of_step() {
  bigint nloc = atom->nlocal;
  bigint bloc = 0;
  ::MPI_Allreduce(&nloc, &bloc, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  if (bloc < atom->natoms) {
    for (int i = 0; i < atom->nlocal; ++i) { flags.erase(atom->tag[i]); }
    for (const auto[k, v] : flags) {
      if (v.first) {
        utils::logmesg(lmp, "{}: Lost particle tag {}, flagged {} timesteps before\n", comm->me, k, update->ntimestep - v.second);
      } else {
        utils::logmesg(lmp, "{}: Lost particle tag {}, not flagged\n", comm->me, k);
      }
    }
    // int nrec_local = flags.size();
    // int* nrec;
    // int* displs;
    // memory->create(nrec, comm->nprocs, "capture/vel:counts");
    // memory->create(displs, comm->nprocs, "capture/vel:displs");
    // ::MPI_Allgather(&nrec_local, 1, MPI_INT, nrec, 1, MPI_INT, world);
    // bigint nrecv = 0;
    // for (int i = 0; i < comm->nprocs; ++i) { nrecv += nrec[i]; }
    // displs[0] = 0;
    // for (int i = 1; i < comm->nprocs; ++i) { displs[i] = displs[i - 1] + nrec[i]; }
    // bigint* to_send;
    // bigint* to_recv;
    // memory->create(to_send, nrec_local, "capture/vel:to_send");
    // memory->create(to_recv, nrecv, "capture/vel:to_recv");
    // int j = 0;
    // for (const auto[k, v] : flags) { to_send[j++] = k; }
    // ::MPI_Allgatherv(to_send, nrec_local, MPI_LMP_BIGINT, to_recv, nrec, displs, MPI_LMP_BIGINT, world);

    // if (comm->me) {}

    // memory->destroy(to_send);
    // memory->destroy(to_recv);
    // memory->destroy(displs);
    // memory->destroy(nrec);
  }
}

/* ---------------------------------------------------------------------- */

[[gnu::hot]] inline long double FixCaptureVel::rminsq(const int i, const int j) noexcept(true)
{
  // constexpr long double two_sqrt        = 1.4142135623730950488016887242096L;    // sqrt(2)
  // constexpr long double two_one_third   = 1.2599210498948731647672106072782L;    // 2^(1/3)
  // constexpr long double one_over_six    = 0.1666666666666666666666666666666L;    // 1/6
  // constexpr long double one_over_three  = 0.3333333333333333333333333333333L;    // 1/3

  constexpr double two_sqrt       = 1.41421356237310;
  constexpr double one_over_three = 0.33333333333333;
  constexpr double two_one_third  = 1.25992104989487;
  constexpr double two_one_third_sq  = two_one_third * two_one_third;
  return two_one_third_sq /
      ::pow((2 + two_sqrt * ::sqrt(2 + (atom->mass[i] * vmax_coeffs[i] + atom->mass[j] * vmax_coeffs[j]) * compute_temp->scalar)), one_over_three);
}

/* ---------------------------------------------------------------------- */

[[gnu::cold]] void FixCaptureVel::post_delete() noexcept(true)
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

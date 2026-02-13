/* ----------------------------------------------------------------------
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

#include "compute_cf_atom.h"
#include "nucc_defs.hpp"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

#ifdef __NUCC_CHECK_ACCESS
#include <cassert>
#endif // __NUCC_CHECK_ACCESS

// Number of bins above and below the central one that will be
// considered as affected by the gaussian kernel
// 3 seems a good compromise between speed and good mollification

constexpr std::array<int, 3> deltabins = {3, 3, 3};
constexpr int range_phi = 2 * deltabins[2] + 1;
constexpr double pi4over3 = (4./3.)*MY_PI;

/* ---------------------------------------------------------------------- */

ComputeCFAtom::ComputeCFAtom(LAMMPS *lmp, int narg, char **arg) : Compute(lmp, narg, arg), volume_loc(0)
{
  peratom_flag = 1;

  if (narg != 8) { error->all( FLERR, "Illegal compute cf/atom command; wrong number of arguments"); }

  // Arguments are: sigma cutoff avg yes/no cutoff2 local yes/no
  //   sigma is the gaussian width
  //   cutoff is the cutoff for the calculation of g(r)
  //   avg is optional and allows averaging the pair entropy over neighbors
  //   the next argument should be yes or no
  //   cutoff2 is the cutoff for the averaging
  //   local is optional and allows using the local density to normalize
  //     the g(r)

  sigmas[0] = utils::numeric(FLERR,arg[3],false,lmp);
  sigmas[1] = utils::numeric(FLERR,arg[4],false,lmp);
  sigmas[2] = utils::numeric(FLERR,arg[5],false,lmp);
  cutoff    = utils::numeric(FLERR,arg[6],false,lmp);
  if (sigmas[0] <= 0.0) { error->all(FLERR,"Illegal compute {} command; sigma r must be positive: {}",     style, arg[3]); }
  if (sigmas[1] <= 0.0) { error->all(FLERR,"Illegal compute {} command; sigma theta must be positive: {}", style, arg[4]); }
  if (sigmas[2] <= 0.0) { error->all(FLERR,"Illegal compute {} command; sigma phi must be positive: {}",   style, arg[5]); }
  if (cutoff    <= 0.0) { error->all(FLERR,"Illegal compute {} command; cutoff must be positive: {}",      style, arg[6]); }
  smooth = utils::logical(FLERR, arg[7], false, lmp);
  cutsq = cutoff*cutoff;

  // optional keywords
  // int iarg = 7;
  // while (iarg < narg) {
  //   if (strcmp(arg[iarg],"avg") == 0) {
  //     if (iarg+3 > narg) error->all(FLERR,"Illegal compute entropy/atom command");
  //     avg_flag = utils::logical(FLERR,arg[iarg+1],false,lmp);
  //     cutoff2 = utils::numeric(FLERR,arg[iarg+2],false,lmp);
  //     if (cutoff2 < 0.0) error->all(FLERR,"Illegal compute entropy/atom command; negative cutoff2");
  //     cutsq2 = cutoff2*cutoff2;
  //     iarg += 3;
  //   } else if (strcmp(arg[iarg],"local") == 0) {
  //     if (iarg+2 > narg) error->all(FLERR,"Illegal compute entropy/atom command");
  //     local_flag = utils::logical(FLERR,arg[iarg+1],false,lmp);
  //     iarg += 2;
  //   } else error->all(FLERR,"Illegal compute entropy/atom command");
  // }

  nbins[0] = static_cast<int>(cutoff / sigmas[0]) + 1;
  nbins[1] = static_cast<int>(2. / sigmas[1]) + 1;
  nbins[2] = static_cast<int>(MY_2PI / sigmas[2]) + 1;
  if (comm->me == 0) { utils::logmesg(lmp, "{}: {} r , {} thata, {} phi bins wiil be used\n", style, nbins[0], nbins[1], nbins[2]); }
  if (nbins[0] < deltabins[0]) { error->all(FLERR, "Insufficient number of r bins, decrease sigma_r"); }
  if (nbins[1] < deltabins[1]) { error->all(FLERR, "Insufficient number of theta bins, decrease sigma_theta"); }
  if (nbins[2] < deltabins[2]) { error->all(FLERR, "Insufficient number of phi bins, decrease sigma_phi"); }

  ddvol = sigmas[0] * sigmas[1] * sigmas[2];
  norm = 1. / (::pow(MY_2PI, 3./2.) * ddvol);
  norms[0] = 1. / (sigmas[0] * sigmas[0]);
  norms[1] = 1. / (sigmas[1] * sigmas[1]);
  norms[2] = 1. / (sigmas[2] * sigmas[2]);
  if (comm->me == 0) { utils::logmesg(lmp, "{}: Normalizations: {:.5f} global, {:.5f} r, {:.5f} thata, {:.5f} phi\n", style, norm, norms[0], norms[1], norms[2]); }

  size_peratom_cols = nbins[0] * nbins[1] * nbins[2] + 1;
  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeCFAtom::~ComputeCFAtom()
{
  memory->destroy(rdf);
  bins[0].destroy(memory);
  bins[1].destroy(memory);
  bins[2].destroy(memory);
}

/* ---------------------------------------------------------------------- */

void ComputeCFAtom::init()
{
  if ((modify->get_compute_by_style(style).size() > 1) && (comm->me == 0)) {
    error->warning(FLERR, "More than one compute {}", style);
  }

  if (force->pair == nullptr) { error->all(FLERR,"Compute {} requires a pair style be defined", style); }

  if ((cutoff) > (force->pair->cutforce  + neighbor->skin)) {
    error->all(FLERR,"Compute {} cutoff is longer than the"
                " pairwise cutoff. Increase the neighbor list skin"
                " distance.", style);
  }

  // Request a full neighbor list
  int list_flags = NeighConst::REQ_FULL;
  // need neighbors of the ghost atoms
  list_flags |= NeighConst::REQ_GHOST;
  neighbor->add_request(this, list_flags);

  bins[0].grow(memory, nbins[0], "r_bin");
  bins[1].grow(memory, nbins[1], "costheta_bin");
  bins[2].grow(memory, nbins[2], "phi_bin");

  for (int i = 0; i < nbins[0]; ++i) { bins[0][i] = (i + 0.5) * sigmas[0]; }
  for (int i = 0; i < nbins[1]; ++i) { bins[1][i] = (i + 0.5) * sigmas[1] - 1; }
  for (int i = 0; i < nbins[2]; ++i) { bins[2][i] = (i + 0.5) * sigmas[2]; }

  rbs.grow(memory, nbins[0], "r_bin_sq");
  for (int i = 0; i < nbins[0]; ++i) { rbs[i] = 1. / (bins[0][i] * bins[0][i] * ddvol); }

  if (comm->me == 0) { utils::logmesg(lmp, "{}: Created: {} r, {} thata, {} phi bins\n", style, nbins[0], nbins[1], nbins[2]); }

  array_atom = memory->create(rdf, static_cast<int>(atom->nmax*NUCC::Defines::ALLOC_COEFF), size_peratom_cols, "rdf");
  memory->create(rdf, static_cast<int>(atom->nmax*NUCC::Defines::ALLOC_COEFF), size_peratom_cols, "rdf");

  double neigh_cutoff = force->pair->cutforce  + neighbor->skin;
  double neigh_bin_vol = neigh_cutoff*neigh_cutoff*neigh_cutoff;
  double volume_loc = pi4over3*neigh_bin_vol;

  initialized_flag = 1;
}

/* ---------------------------------------------------------------------- */

void ComputeCFAtom::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeCFAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    array_atom = memory->grow(rdf, static_cast<int>(atom->nmax*NUCC::Defines::ALLOC_COEFF), size_peratom_cols, "rdf");
  }

  int   inum = list->inum + list->gnum;
  int*  ilist = list->ilist;
  int*  numneigh = list->numneigh;
  int** firstneigh = list->firstneigh;

  double **x = atom->x;
  int *mask = atom->mask;

  double gvolume = domain->xprd * domain->yprd * domain->zprd;
  double ginvdensity = gvolume / atom->natoms;



  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    if ((mask[i] & groupbit) != 0) {
      double xtmp = x[i][0];
      double ytmp = x[i][1];
      double ztmp = x[i][2];
      int* jlist = firstneigh[i];
      int jnum = numneigh[i];
      double* cfi = rdf[i];

      if (jnum < 3) { continue; }

      double invdensity = (jnum < 3) ? (ginvdensity) : (2 * volume_loc / (jnum - 1) / (jnum - 2));
      cfi[size_peratom_cols - 1] = invdensity;

      // loop over list of all neighbors within force cutoff

      // initialize cf
      ::memset(rdf[i], 0.0, size_peratom_cols * sizeof(double));

      for (int jj = 0; jj < jnum; jj++) {
        int j = jlist[jj];
        j &= NEIGHMASK;

        double delx = xtmp - x[j][0];
        double dely = ytmp - x[j][1];
        double delz = ztmp - x[j][2];
        double rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < cutsq) {
          // contribute to cf
          double r = ::sqrt(rsq);
          double cos_theta = delz / r;
          double phi = ::atan2(dely, delx);
          std::array<int, 3> ibins = {
            static_cast<int>((r - 0) / sigmas[0]),
            static_cast<int>((cos_theta + 1) / sigmas[1]),
            static_cast<int>(phi / sigmas[2])
          };
          ibins[2] %= nbins[2];
          if (ibins[2] < 0) { ibins[2] += nbins[2]; }

          if (smooth > 0) {
            int minbin_r = MIN(MAX(ibins[0] - deltabins[0], 0), nbins[0] - 1);
            int maxbin_r = MIN(MAX(ibins[0] + deltabins[0], 0), nbins[0] - 1);
            int minbin_theta = MIN(MAX(ibins[1] - deltabins[1], 0), nbins[1] - 1);
            int maxbin_theta = MIN(MAX(ibins[1] + deltabins[1], 0), nbins[1] - 1);
            int m = (ibins[2] - deltabins[2]) % nbins[2];
            if (m < 0) { m += nbins[2]; }
            // if (comm->me == 0) { utils::logmesg(lmp, "{}: _r: {} ({}-{}), _t: {}({}-{}), _p: {}\n", style,
            //   ibins[0], minbin_r, maxbin_r, ibins[1], minbin_theta, maxbin_theta, m
            // ); }

            for (int k = minbin_r; k < maxbin_r; ++k) {
              double __r = (bins[0][k] - r);
              double _r = __r * __r * norms[0];
              for (int l = minbin_theta; l < maxbin_theta; ++l) {
                double __t = (bins[1][l] - cos_theta);
                double _t = __t * __t * norms[1];
                for (int _m = 0; _m < range_phi; ++_m) {
                  double __p = bins[2][m] - phi;
                  double _p = __p + MY_PI;
                  int iquot = static_cast<int>(_p/MY_2PI);
                  _p -= iquot * MY_2PI + MY_PI;
                  _p *= _p * norms[2];

                  // if (comm->me == 0) { utils::logmesg(lmp, "{}: Sum: {}, dens: {}\n", style, asum, invdensity); }

                  int idx = k * (nbins[1] * nbins[2]) + l * nbins[2] + m;
                  #ifdef __NUCC_CHECK_ACCESS
                  assert(idx < size_array_cols);
                  #endif // __NUCC_CHECK_ACCESS

                  double fct = -0.5 * (_r + _t + _p);
                  cfi[idx] += ::exp(fct) * norm * invdensity;
                  // (k == ibins[0]) && (l == ibins[1]) && (m == ibins[2])
                  // if ((comm->me == 0) && (i==1)) {
                  //   // utils::logmesg(lmp, "{}: _r: {:.5f}({:.5f}), _t: {:.5f}({:.5f}), _p: {:.5f}({} - {:.5f}) -> {:.5f}\n", style, __r, _r, __t, _t, __p, iquot, _p, fct);
                  //   utils::logmesg(lmp, "{}: fct: {}, cf: {}\n", style, fct, rdf[i][idx]);
                  // }

                  if (++m >= nbins[2]) { m = 0; }
                }
              }
            }
          } else {
            const int idx = ibins[0] * (nbins[1] * nbins[2]) + ibins[1] * nbins[2] + ibins[2];
            #ifdef __NUCC_CHECK_ACCESS
            assert(idx < size_array_cols);
            #endif // __NUCC_CHECK_ACCESS
            cfi[idx] += 1;
          }
        }
      }
      // if (i == 1) {
      //   double asum = 0;
      //   for (int k = 0; k < size_peratom_cols; ++k) { asum += rdf[i][k]; }
      //   if (comm->me == 0) { utils::logmesg(lmp, "{}: Sum: {}, dens: {}\n", style, asum, invdensity); }
      // }
    }
  }

  if (smooth == 0) {
    for (int m = 0; m < atom->nmax; ++m) {
      double* cfi = rdf[m];
      double invdensity = cfi[size_peratom_cols-1];
      for (int i = 0; i < nbins[0]; ++i) {
        double dvol = rbs[i] * invdensity;
        int nx = i * (nbins[1] * nbins[2]);
        for (int j = 0; j < nbins[1]; ++j) {
          int ny = nx + j * nbins[1];
          for (int k = 0; k < nbins[2]; ++k) {
            cfi[ny + k] *= dvol;
          }
        }
      }
    }
  }

}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeCFAtom::memory_usage()
{
  std::size_t res = nmax * (size_peratom_cols * sizeof(double) + sizeof(double*));
  res += bins[0].memory_usage() + bins[1].memory_usage() + bins[2].memory_usage();
  return res;
}

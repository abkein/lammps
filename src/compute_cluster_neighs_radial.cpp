// clang-format off
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

#include "compute_cluster_neighs_radial.h"
#include "nucc_defs.hpp"

#include "atom.h"
#include "comm.h"
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

/* ---------------------------------------------------------------------- */

ComputeClusterNeighsRadial::ComputeClusterNeighsRadial(LAMMPS *lmp, int narg, char **arg) : Compute(lmp, narg, arg)
{
  peratom_flag = 1;

  if (narg != 6) { error->all( FLERR, "Illegal compute cf/atom command; wrong number of arguments"); }

  // Arguments are: sigma cutoff avg yes/no cutoff2 local yes/no
  //   sigma is the gaussian width
  //   cutoff is the cutoff for the calculation of g(r)
  //   avg is optional and allows averaging the pair entropy over neighbors
  //   the next argument should be yes or no
  //   cutoff2 is the cutoff for the averaging
  //   local is optional and allows using the local density to normalize
  //     the g(r)

  sigma = utils::numeric(FLERR,arg[3],false,lmp);
  cutoff    = utils::numeric(FLERR,arg[4],false,lmp);
  if (sigma <= 0.0) { error->all(FLERR,"Illegal compute {} command; sigma r must be positive: {}",     style, arg[3]); }
  if (cutoff    <= 0.0) { error->all(FLERR,"Illegal compute {} command; cutoff must be positive: {}",      style, arg[4]); }
  // Get cluster/atom compute
  compute_cluster_atom = lmp->modify->get_compute_by_id(arg[5]);
  if (compute_cluster_atom == nullptr) { error->all(FLERR, "{}: Cannot find compute with style 'cluster/atom' with given id: {}", style, arg[5]); }

  cutsq = cutoff*cutoff;
  norm = 1. / (MathConst::MY_4PI * sigma * sigma * sigma);

  nbins = static_cast<int>(::ceil(cutoff / sigma));
  if (comm->me == 0) { utils::logmesg(lmp, "{}: {} r bins will be used for each atom\n", style, nbins); }
  size_peratom_cols = nbins;
}

/* ---------------------------------------------------------------------- */

ComputeClusterNeighsRadial::~ComputeClusterNeighsRadial()
{
  memory->destroy(rdf);
#ifdef __NUCC_NEIGHS_RADIAL_PRECOMPUTE_NORM
  norms.destroy(memory);
#endif // __NUCC_NEIGHS_RADIAL_PRECOMPUTE_NORM

}

/* ---------------------------------------------------------------------- */

void ComputeClusterNeighsRadial::init()
{
  if ((modify->get_compute_by_style(style).size() > 1) && (comm->me == 0)) {
    error->warning(FLERR, "More than one compute {}", style);
  }

  if (force->pair == nullptr) { error->all(FLERR,"Compute {} requires a pair style be defined", style); }

  if ((cutoff) > (force->pair->cutforce  + neighbor->skin)) {
    error->all(FLERR,"Compute {} cutoff is longer than the"
                " pairwise cutoff+skin length. Increase the neighbor list skin"
                " distance.", style);// norms.create(memory, nbins, )
  }

#ifdef __NUCC_NEIGHS_RADIAL_PRECOMPUTE_NORM
  norms.create(memory, nbins, "compute:neighs/radial:norms");
  for (int k = 0; k < nbins; ++k) {
    const double tmp = 1./(0.5 + k);
    norms[k] = norm * tmp * tmp;
  }
#endif // __NUCC_NEIGHS_RADIAL_PRECOMPUTE_NORM

// Request a neighbor list
#ifdef __NUCC_NEIGHS_RADIAL_USE_HALF
  neighbor->add_request(this);
#else // __NUCC_NEIGHS_RADIAL_USE_HALF
  neighbor->add_request(this, NeighConst::REQ_FULL);
#endif // __NUCC_NEIGHS_RADIAL_USE_HALF

  nmax = atom->nmax;
  array_atom = memory->create(rdf, static_cast<int>(nmax*LMP_NUCC_ALLOC_COEFF), size_peratom_cols, "compute:neighs/radial:rdf");

  initialized_flag = 1;
}

/* ---------------------------------------------------------------------- */

void ComputeClusterNeighsRadial::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

#ifndef __NUCC_NEIGHS_RADIAL_USE_HALF
void ComputeClusterNeighsRadial::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  if (compute_cluster_atom->invoked_peratom != update->ntimestep) { compute_cluster_atom->compute_peratom(); }
  const double* const cluster_ids = compute_cluster_atom->vector_atom;

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    array_atom = memory->grow(rdf, static_cast<int>(nmax*LMP_NUCC_ALLOC_COEFF), size_peratom_cols, "compute:neighs/radial:rdf");
  }

  for (int i = 0; i < atom->nlocal; ++i) {
    ::memset(rdf[i], 0.0, size_peratom_cols * sizeof(double));
  }

  const int   inum = list->inum;
  const int*  ilist = list->ilist;
  const int*  numneigh = list->numneigh;
  int** firstneigh = list->firstneigh;

  double **x = atom->x;
  const int *mask = atom->mask;

  for (int ii = 0; ii < inum; ++ii) {
    const int i = ilist[ii];
    if ((mask[i] & groupbit) != 0) {
      const double xtmp = x[i][0];
      const double ytmp = x[i][1];
      const double ztmp = x[i][2];
      const int* jlist = firstneigh[i];
      const int jnum = numneigh[i];
      double* cfi = rdf[i];

      if (jnum == 0) { continue; }

      // loop over list of all neighbors within force cutoff

      // // initialize cf
      // ::memset(cfi, 0.0, size_peratom_cols * sizeof(double));

      for (int jj = 0; jj < jnum; ++jj) {
        const int j = jlist[jj] & NEIGHMASK;

        if (static_cast<int>(cluster_ids[i]) != static_cast<int>(cluster_ids[j])) { continue; }

        const double delx = xtmp - x[j][0];
        const double dely = ytmp - x[j][1];
        const double delz = ztmp - x[j][2];
        const double rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < cutsq) {
          // contribute to cf
          const double r = ::sqrt(rsq);
          const int ibin = static_cast<int>((r - 0) / sigma);

          #ifdef __NUCC_CHECK_ACCESS
          if (!(ibin < size_array_cols)) {
            error->one(FLERR, "{}@{}: cutoff: {}, sigma: {}, nbins: {}, cutsq: {}, rsq: {}, r: {}, ibin: {}. ibin > nbins", style, comm->me, cutoff, sigma, nbins, cutsq, rsq, r, ibin);
          }
          #endif // __NUCC_CHECK_ACCESS

          cfi[ibin] += 1;
        }
      }
      for (int k = 0; k < nbins; ++k) {
#ifdef __NUCC_NEIGHS_RADIAL_PRECOMPUTE_NORM
        cfi[k] *= norms[k];
#else // __NUCC_NEIGHS_RADIAL_PRECOMPUTE_NORM
        const double tmp = 1./(0.5 + k);
        cfi[k] *= norm * tmp * tmp;
#endif // __NUCC_NEIGHS_RADIAL_PRECOMPUTE_NORM
      }
    }
  }
}
#else // __NUCC_NEIGHS_RADIAL_USE_HALF

/* ---------------------------------------------------------------------- */

void ComputeNeighsRadial::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  if (compute_cluster_atom->invoked_peratom != update->ntimestep) { compute_cluster_atom->compute_peratom(); }
  const double* const cluster_ids = compute_cluster_atom->vector_atom;

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    array_atom = memory->grow(rdf, static_cast<int>(nmax*LMP_NUCC_ALLOC_COEFF), size_peratom_cols, "compute:neighs/radial:rdf");
  }

  for (int i = 0; i < atom->nlocal; ++i) {
    ::memset(rdf[i], 0.0, size_peratom_cols * sizeof(double));
  }

  const int   inum = list->inum;
  const int*  ilist = list->ilist;
  const int*  numneigh = list->numneigh;
  int** firstneigh = list->firstneigh;

  double **x = atom->x;
  const int *mask = atom->mask;

  for (int ii = 0; ii < inum; ++ii) {
    const int i = ilist[ii];
    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];
    const int* jlist = firstneigh[i];
    const int jnum = numneigh[i];

    if (jnum == 0) { continue; }

    // loop over list of all neighbors within force cutoff

    for (int jj = 0; jj < jnum; ++jj) {
      const int j = jlist[jj] & NEIGHMASK;

      if (((mask[i] | mask[j]) & groupbit) == 0) { continue; }
      if (static_cast<int>(cluster_ids[i]) != static_cast<int>(cluster_ids[j])) { continue; }

      const double delx = xtmp - x[j][0];
      const double dely = ytmp - x[j][1];
      const double delz = ztmp - x[j][2];
      const double rsq = delx*delx + dely*dely + delz*delz;
      if (rsq < cutsq) {
        // contribute to cf
        const double r = ::sqrt(rsq);
        const int ibin = static_cast<int>((r - 0) / sigma);

        #ifdef __NUCC_CHECK_ACCESS
        if (!(ibin < size_array_cols)) {
          error->one(FLERR, "{}: cutoff: {}, sigma: {}, nbins: {}, cutsq: {}, rsq: {}, r: {}, ibin: {}. ibin > nbins", style, cutoff, sigma, nbins, cutsq, rsq, r, ibin);
        }
        #endif // __NUCC_CHECK_ACCESS

        if ((mask[i] & groupbit) != 0) {
          rdf[i][ibin] += 1;
        }
        if ((j < atom->nlocal) && ((mask[j] & groupbit) != 0)) {
          rdf[j][ibin] += 1;
        }
      }
    }
  }
  for (int i = 0; i < atom->nlocal; ++i) {
    if ((mask[i] & groupbit) != 0) {
      for (int k = 0; k < nbins; ++k) {
#ifdef __NUCC_NEIGHS_RADIAL_PRECOMPUTE_NORM
        rdf[i][k] *= norms[k];
#else // __NUCC_NEIGHS_RADIAL_PRECOMPUTE_NORM
        const double tmp = 1./(0.5 + k);
        rdf[i][k] *= norm * tmp * tmp;
#endif // __NUCC_NEIGHS_RADIAL_PRECOMPUTE_NORM
      }
    }
  }
}
#endif // __NUCC_NEIGHS_RADIAL_USE_HALF


/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeClusterNeighsRadial::memory_usage()
{
  return nmax * (size_peratom_cols * sizeof(double) + sizeof(double*));
}

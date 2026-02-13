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

#include "compute_neighs_radial_cluster.h"
#include "nucc_defs.hpp"

#include "atom.h"
#include "comm.h"
#include "error.h"
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

ComputeNeighsRadialCluster::ComputeNeighsRadialCluster(LAMMPS* lmp, int narg, char** arg) : ComputeNeighsRadialBase(lmp, narg, arg)
{
  if (narg < 6) { error->all(FLERR, "Illegal compute neighs/raidal/cluster command; wrong number of arguments"); }

  // Get cluster/atom compute
  compute_cluster_atom = lmp->modify->get_compute_by_id(arg[5]);
  if (compute_cluster_atom == nullptr) { error->all(FLERR, "{}: Cannot find compute with style 'cluster/atom' with given id: {}", style, arg[5]); }
}

/* ---------------------------------------------------------------------- */

void ComputeNeighsRadialCluster::compute_peratom()
{
  if (invoked_peratom == update->ntimestep) { return; }
  invoked_peratom = update->ntimestep;

  if (compute_cluster_atom->invoked_peratom != update->ntimestep) { compute_cluster_atom->compute_peratom(); }
  const double* const cluster_ids = compute_cluster_atom->vector_atom;

  if (atom->nmax > nmax) {
    nmax       = atom->nmax;
    array_atom = memory->grow(rdf, static_cast<int>(nmax * NUCC::Defines::ALLOC_COEFF), size_peratom_cols, "compute:neighs/radial:rdf");
  }

  for (int i = 0; i < atom->nlocal; ++i) { ::memset(rdf[i], 0.0, size_peratom_cols * sizeof(double)); }

  const int                  inum       = list->inum;
  const int* const           ilist      = list->ilist;
  const int* const           numneigh   = list->numneigh;
  const int* const* const    firstneigh = list->firstneigh;

  const double* const* const x          = atom->x;
  const int* const           mask       = atom->mask;

  for (int ii = 0; ii < inum; ++ii) {
    const int i = ilist[ii];
    if ((mask[i] & groupbit) != 0) { continue; }
    const double     xtmp  = x[i][0];
    const double     ytmp  = x[i][1];
    const double     ztmp  = x[i][2];
    const int* const jlist = firstneigh[i];
    const int        jnum  = numneigh[i];

    // loop over list of all neighbors within force cutoff

    for (int jj = 0; jj < jnum; ++jj) {
      const int j = jlist[jj] & NEIGHMASK;

      if (static_cast<int>(cluster_ids[i]) != static_cast<int>(cluster_ids[j])) { continue; }

      const double delx = xtmp - x[j][0];
      const double dely = ytmp - x[j][1];
      const double delz = ztmp - x[j][2];
      const double rsq  = delx * delx + dely * dely + delz * delz;
      if (rsq < cutsq) {
        // contribute to cf
        const double r    = ::sqrt(rsq);
        const int    ibin = static_cast<int>((r - 0) / delta);

        if constexpr (NUCC::Defines::CHECK_ACCESS) {
          if (!(ibin < size_array_cols)) {
            error->one(FLERR, "{}@{}: cutoff: {}, delta: {}, nbins: {}, cutsq: {}, rsq: {}, r: {}, ibin: {}. ibin > nbins", style, comm->me, cutoff,
                       delta, nbins, cutsq, rsq, r, ibin);
          }
        }

        rdf[i][ibin] += 1;
        if constexpr (NUCC::Defines::NEIGHS_RADIAL_USE_HALF) {
          if ((j < atom->nlocal) && ((mask[j] & groupbit) != 0)) { rdf[j][ibin] += 1; }
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeNeighsRadialCluster::memory_usage()
{
  return static_cast<double>(nmax * (size_peratom_cols * sizeof(double) + sizeof(double*)));
}

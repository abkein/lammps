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

#include "compute_neighs_radial_base.h"

#include "nucc_defs.hpp"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeNeighsRadialBase::ComputeNeighsRadialBase(LAMMPS* lmp, int narg, char** arg) : Compute(lmp, narg, arg)
{
  peratom_flag = 1;

  if (narg < 5) { error->all(FLERR, "Illegal compute neighs/radial command; wrong number of arguments"); }

  delta  = utils::numeric(FLERR, arg[3], false, lmp);
  cutoff = utils::numeric(FLERR, arg[4], false, lmp);
  if (delta <= 0.0) { error->all(FLERR, "Illegal compute {} command; delta r must be positive: {}", style, arg[3]); }
  if (cutoff <= 0.0) { error->all(FLERR, "Illegal compute {} command; cutoff must be positive: {}", style, arg[4]); }
  cutsq = cutoff * cutoff;

  nbins = static_cast<int>(::ceil(cutoff / delta));
  if (comm->me == 0) { utils::logmesg(lmp, "{}: {} r bins will be used for each atom\n", style, nbins); }
  size_peratom_cols = nbins;
}

/* ---------------------------------------------------------------------- */

ComputeNeighsRadialBase::~ComputeNeighsRadialBase()
{
  memory->destroy(rdf);
}

/* ---------------------------------------------------------------------- */

void ComputeNeighsRadialBase::init()
{
  if ((modify->get_compute_by_style(style).size() > 1) && (comm->me == 0)) { error->warning(FLERR, "More than one compute {}", style); }

  if (force->pair == nullptr) { error->all(FLERR, "Compute {} requires a pair style be defined", style); }

  if ((cutoff) > (force->pair->cutforce + neighbor->skin)) {
    error->all(FLERR,
               "Compute {} cutoff is longer than the"
               " pairwise cutoff+skin length. Increase the neighbor list skin"
               " distance.",
               style);    // norms.create(memory, nbins, )
  }

  // Request a neighbor list
  if constexpr (NUCC::Defines::NEIGHS_RADIAL_USE_HALF) {
    neighbor->add_request(this);
  } else {
    neighbor->add_request(this, NeighConst::REQ_FULL);
  }

  nmax             = atom->nmax;
  array_atom       = memory->create(rdf, static_cast<int>(nmax * NUCC::Defines::ALLOC_COEFF), size_peratom_cols, "compute:neighs/radial:rdf");

  initialized_flag = 1;
}

/* ---------------------------------------------------------------------- */

void ComputeNeighsRadialBase::init_list(int /*id*/, NeighList* ptr)
{
  list = ptr;
}

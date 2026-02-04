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

#include "compute_neighs_radial_cluster_size.h"
#include "compute_cluster_size_ext.h"
#include "compute_neighs_radial_base.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "pair.h"
#include "update.h"

#include <algorithm>
#include <cmath>
#include <cstring>

constexpr int NEIGH_BIN_CUTOFF_COEFF = 2;

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeNeighsRadialClusterSize::ComputeNeighsRadialClusterSize(LAMMPS *lmp, int narg, char **arg) : Compute(lmp, narg, arg)
{
  local_flag = 1;
  array_flag = 1;
  extarray = 0;

  if (narg < 5) { error->all( FLERR, "Illegal compute cf/atom command; wrong number of arguments"); }

  // Get cluster/size compute
  compute_cluster_size = dynamic_cast<ComputeClusterSizeExt*>(lmp->modify->get_compute_by_id(arg[3]));
  if (compute_cluster_size == nullptr) { error->all(FLERR, "{}: Cannot find compute with style 'cluster/size' with given id: {}", style, arg[3]); }
  // Get neighs/radial compute
  compute_neighs_radial = dynamic_cast<ComputeNeighsRadialBase*>(lmp->modify->get_compute_by_id(arg[4]));
  if (compute_neighs_radial == nullptr) { error->all(FLERR, "{}: Cannot find compute with style 'neighs/radial' with given id: {}", style, arg[4]); }
  size_array_rows = size_local_rows = cutoff = compute_cluster_size->get_size_cutoff();
  size_array_cols = size_local_cols = nbins = compute_neighs_radial->get_nbins(); // compute_neighs_radial->size_peratom_cols
  delta = compute_neighs_radial->get_delta_r();

  if (narg > 5 && ::strcmp(arg[5], "smooth") == 0) {
    if (narg < 7) { error->all( FLERR, "Illegal compute cf/atom command; wrong number of arguments"); }
    do_smooth = 1;
    sigma = utils::numeric(FLERR,arg[6],false,lmp);
    if (sigma <= 0.0) { error->all(FLERR,"Illegal compute {} command; kernel width must be positive: {}",     style, arg[6]); }
    max_neigh_bin = static_cast<int>(::ceil(NEIGH_BIN_CUTOFF_COEFF*sigma/delta));
  }
  norm = 1. / (MathConst::MY_4PI * delta * delta * delta * comm->nprocs);
  if (comm->me == 0) {
    utils::logmesg(lmp, "{}: cutoff: {}, nbins: {}, delta: {}, norm: {:.8f}\n", style, cutoff, nbins, delta, norm);
  }
}

/* ---------------------------------------------------------------------- */

ComputeNeighsRadialClusterSize::~ComputeNeighsRadialClusterSize()
{
  memory->destroy(counts);
  memory->destroy(counts_global);
  if (do_smooth > 0) {
    memory->destroy(counts2);
    memory->destroy(weights);
  }
  atom_counts_by_size.destroy(memory);
  #ifdef __NUCC_NEIGHS_RADIAL_PRECOMPUTE_NORM
  norms.destroy(memory);
  #endif // __NUCC_NEIGHS_RADIAL_PRECOMPUTE_NORM
}

/* ---------------------------------------------------------------------- */

void ComputeNeighsRadialClusterSize::init()
{
  array_local = memory->create(counts, cutoff, nbins, "compute:neighs/radial/size:counts");
  array = memory->create(counts_global, cutoff, nbins, "compute:neighs/radial/size:counts_global");
  atom_counts_by_size.create(memory, cutoff, "compute:neighs/radial/size:atom_counts");

  if (do_smooth > 0) {
    array_local = memory->create(counts2, cutoff, nbins, "compute:neighs/radial/size:counts_raw");

    const int num_neighs = 2*max_neigh_bin + 1;
    memory->create(weights, nbins, num_neighs, "compute:neighs/radial/size:weights");
    for (int i = 0; i<nbins; ++i){
      ::memset(weights[i], 0.0, num_neighs * sizeof(double));
    }

    const double coeff = -delta*delta/2/(sigma*sigma);
    for (int i = 0; i<nbins; ++i){
      for (int j = -max_neigh_bin; j<max_neigh_bin+1; ++j){
        if ((i+j>=0) && (i+j<nbins)) {
          weights[i][max_neigh_bin+j] = ::exp(j*j*coeff);
        } else {
          weights[i][max_neigh_bin+j] = 0;
        }
      }
    }
    for (int i = 0; i<nbins; ++i){
      double sum = 0;
      for (int j = 0; j<num_neighs; ++j){
        sum += weights[i][j];
      }
      for (int j = 0; j<num_neighs; ++j){
        weights[i][j] /= sum;
      }
    }
  }

  #ifdef __NUCC_NEIGHS_RADIAL_PRECOMPUTE_NORM
  norms.create(memory, nbins, "compute:neighs/radial/size:norms");
  if (comm->me == 0) {utils::logmesg(lmp, "{}: norms:", style);}
  for (int nbin = 0; nbin < nbins; ++nbin) {
    const double tmp = 1./(0.5 + nbin);
    norms[nbin] = norm * tmp * tmp;
    if (comm->me == 0) {utils::logmesg(lmp, " {:.8f}", norms[nbin]);}
  }
  if (comm->me == 0) {utils::logmesg(lmp, "\n");}
  #endif // __NUCC_NEIGHS_RADIAL_PRECOMPUTE_NORM

  initialized_flag = 1;
}

/* ---------------------------------------------------------------------- */

void ComputeNeighsRadialClusterSize::compute_local()
{
  if (invoked_local == update->ntimestep) { return; }
  invoked_local = update->ntimestep;

  if (compute_cluster_size->invoked_peratom != update->ntimestep) { compute_cluster_size->compute_peratom(); }
  if (compute_neighs_radial->invoked_peratom != update->ntimestep) { compute_neighs_radial->compute_peratom(); }

  const double* const* const rdf = compute_neighs_radial->array_atom;
  const double* sizes = compute_cluster_size->vector_atom;

  const int *mask = atom->mask;

  for (int i = 0; i < cutoff; ++i) {
    ::memset(counts[i], 0.0, nbins * sizeof(double));
  }
  atom_counts_by_size.reset();

  for (int i = 0; i < atom->nlocal; ++i) {
    if ((mask[i] & groupbit) != 0) {
      const int size = static_cast<int>(sizes[i]);
      ++atom_counts_by_size[size];
      double* const cf_sum = counts[size];
      const double* const cfi = rdf[i];
      for (int nbin = 0; nbin < nbins; ++nbin){
        cf_sum[nbin] += cfi[nbin];
      }
    }
  }

  for (int size = 0; size < cutoff; ++size){
    double* const size_counts = counts[size];
    const int count = atom_counts_by_size[size];
    if (count == 0) { continue; }
    for (int nbin = 0; nbin < nbins; ++nbin){
      size_counts[nbin] /= count;
    }
  }

  if (do_smooth > 0){
    for (int i = 0; i < cutoff; ++i) {
      ::memset(counts2[i], 0.0, nbins * sizeof(double));
    }
    const int num_neighs = 2*max_neigh_bin + 1;
    for (int size = 0; size < cutoff; ++size) {
      const double* const size_counts = counts[size];
      double* const size_counts2 = counts2[size];
      for (int nbin = 0; nbin<nbins; ++nbin){
        double sum = 0;
        const double* const bin_weights = weights[nbin];

        const int ii = nbin-max_neigh_bin;
        for (int j = std::max(0, -ii); j<std::min(num_neighs, nbins-ii); ++j){
          sum += size_counts[ii+j]*bin_weights[j];
        }

        // more explicit version of the same loop
        // for (int j = -max_neigh_bin; j<max_neigh_bin+1; ++j){
        //   if ((nbin+j>=0) && (nbin+j<nbins)) {
        //     sum += size_counts[nbin+j]*bin_weights[max_neigh_bin+j];
        //   }
        // }

        size_counts2[nbin] = sum;
      }
    }
  }

  for (int size = 0; size<cutoff; ++size) {
    double* const bins_size = array_local[size];
    for (int nbin = 0; nbin < nbins; ++nbin) {
      #ifdef __NUCC_NEIGHS_RADIAL_PRECOMPUTE_NORM
      const double bin_norm = norms[nbin];
      #else // __NUCC_NEIGHS_RADIAL_PRECOMPUTE_NORM
      const double tmp = 1./(0.5 + nbin);
      const double bin_norm = norm * tmp * tmp;
      #endif // __NUCC_NEIGHS_RADIAL_PRECOMPUTE_NORM
      bins_size[nbin] *= bin_norm;
    }
  }

  // rearranged variant possibly having better performance
  // for (int nbin = 0; nbin < nbins; ++nbin) {
  //   #ifdef __NUCC_NEIGHS_RADIAL_PRECOMPUTE_NORM
  //   const double bin_norm = norms[nbin];
  //   #else // __NUCC_NEIGHS_RADIAL_PRECOMPUTE_NORM
  //   const double tmp = 1./(0.5 + nbin);
  //   const double bin_norm = norm * tmp * tmp;
  //   #endif // __NUCC_NEIGHS_RADIAL_PRECOMPUTE_NORM
  //   for (int size = 0; size<cutoff; ++size) {
  //     array_local[size][nbin] *= bin_norm;
  //   }
  // }
}

/* ---------------------------------------------------------------------- */

void ComputeNeighsRadialClusterSize::compute_array(){
  if (invoked_array == update->ntimestep) { return; }
  invoked_array = update->ntimestep;

  if (invoked_local != update->ntimestep) { compute_local(); }

  ::memset(counts_global[0], 0.0, cutoff * nbins * sizeof(double));

  ::MPI_Allreduce(array_local[0], counts_global[0], cutoff * nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeNeighsRadialClusterSize::memory_usage()
{
  return 2 * cutoff * (nbins + 1) * sizeof(double) + nbins * sizeof(double*);
}

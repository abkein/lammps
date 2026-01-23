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

#include "compute_cluster_cf.h"
#include "compute_cluster_size_ext.h"
#include "compute_cf_atom.h"

#include "comm.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

#include <cstring>
#include <unordered_map>

using namespace LAMMPS_NS;
// using NUCC::cspan;

/* ---------------------------------------------------------------------- */

ComputeClusterCF::ComputeClusterCF(LAMMPS* lmp, int narg, char** arg) : Compute(lmp, narg, arg)
{
  array_flag      = 1;
  size_array_cols = 0;
  size_array_rows = 0;
  extarray        = 0;
  local_flag      = 1;
  size_local_rows = 0;
  size_local_cols = 0;

  if (narg < 5) { utils::missing_cmd_args(FLERR, "compute cf/cluster/avg", error); }

  // Parse arguments //

  // Get cluster/size compute
  compute_cluster_size = dynamic_cast<ComputeClusterSizeExt*>(lmp->modify->get_compute_by_id(arg[3]));
  if (compute_cluster_size == nullptr) { error->all(FLERR, "compute {}: Cannot find compute with style 'size/cluster' with id: {}", style, arg[3]); }

  compute_rdf_atom = dynamic_cast<ComputeCFAtom*>(lmp->modify->get_compute_by_id(arg[4]));
  if (compute_rdf_atom == nullptr) { error->all(FLERR, "compute {}: Cannot find compute with style 'cf/atom' with id: {}", style, arg[4]); }

  size_cutoff = compute_cluster_size->get_size_cutoff();

  int iarg    = 5;
  while (iarg < narg) {
    if (iarg + 2 > narg) { error->all(FLERR, "Illegal command", style); }
    if (::strcmp(arg[iarg], "cut") == 0) {
      int t_size_cutoff = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (t_size_cutoff < 1) { error->all(FLERR, "size_cutoff for compute {} must be greater than 0", style); }
      if (t_size_cutoff > size_cutoff) {
        error->all(FLERR,
                   "size_cutoff for compute {} cannot be greater than it of "
                   "compute cluster/size",
                   style);
      }
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal {} command", style);
    }
  }
}

/* ---------------------------------------------------------------------- */

ComputeClusterCF::~ComputeClusterCF() noexcept(true)
{
  memory->destroy(cf);
  memory->destroy(cf_local);
}

/* ---------------------------------------------------------------------- */

void ComputeClusterCF::init()
{
  if ((modify->get_compute_by_style(style).size() > 1) && (comm->me == 0)) { error->warning(FLERR, "More than one compute {}", style); }

  size_local_rows = size_cutoff + 1;
  size_local_cols = compute_rdf_atom->size_peratom_cols;
  array_local     = memory->create(cf_local, size_local_rows, size_local_cols, "cf/cluster/avg:cf_local");

  size_array_rows = size_cutoff + 1;
  size_array_cols = compute_rdf_atom->size_peratom_cols;
  array           = memory->create(cf, size_array_rows, size_array_cols, "cf/cluster/avg:cf");
}

/* ---------------------------------------------------------------------- */

void ComputeClusterCF::compute_array()
{
  invoked_array = update->ntimestep;

  compute_local();

  for (int i = 1; i < size_array_rows; ++i) {
    ::memset(cf[i], 0.0, size_array_cols * sizeof(double));
    ::MPI_Allreduce(cf_local[i], cf[i], size_array_cols, MPI_DOUBLE, MPI_SUM, world);
  }

  const double* dist = compute_cluster_size->vector;
  for (int i = 1; i < size_array_rows; ++i) {
    if (dist[i] > 0) {
      double norm = dist[i] * i;
      for (int j = 0; j < size_array_cols; ++j) { cf[i][j] /= norm; }
    }
  }
  // if (comm->me == 0) {
  //   for (int i = 1; i < size_array_rows; ++i) {
  //     if (dist[i] > 0) {
  //       double gsum = 0;
  //       for (int j = 0; j < size_array_cols; ++j) { gsum += cf[i][j]; }
  //       utils::logmesg(lmp , "{}-{:.5f} ", i, gsum);
  //     }
  //   }
  //   utils::logmesg(lmp , "\n");
  // }
}

/* ---------------------------------------------------------------------- */

void ComputeClusterCF::compute_local()
{
  invoked_local = update->ntimestep;

  if (compute_cluster_size->invoked_vector != update->ntimestep) { compute_cluster_size->compute_vector(); }

  if (compute_rdf_atom->invoked_peratom != update->ntimestep) { compute_rdf_atom->compute_peratom(); }

  for (int i = 0; i < size_local_rows; ++i) { ::memset(cf_local[i], 0.0, size_local_cols * sizeof(double)); }

  double** const peratomcf = compute_rdf_atom->array_atom;
  int nclusters = compute_cluster_size->get_cluster_map().size();
  const auto& clusters = compute_cluster_size->get_clusters();

  for (int i = 0; i < nclusters; ++i) {
    const auto& clstr = clusters[i];
    if (clstr.g_size >= size_cutoff) { continue; }
    // if ((clstr.g_size >= size_cutoff) || (clstr.g_size >= size_local_rows)) { continue; }
    const auto atoms = clstr.atoms();
    for (int j = 0; j < clstr.l_size; ++j) {
      for (int k = 0; k < size_local_cols; ++k) {
        // if (k >= size_local_cols) { continue; }
        cf_local[clstr.g_size][k] += peratomcf[atoms[j]][k];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeClusterCF::memory_usage()
{
  std::size_t sum = size_array_rows * (size_array_cols * sizeof(double) + sizeof(double*));
  sum += size_local_rows * (size_local_cols * sizeof(double) + sizeof(double*));
  return static_cast<double>(sum);
}

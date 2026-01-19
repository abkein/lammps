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

#include "compute_cluster_size_ext.h"
#include "nucc_cspan.hpp"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

#include <cstddef>
#include <cstring>
#include <utility>
#include <algorithm>

using namespace LAMMPS_NS;
using namespace NUCC;

/* ---------------------------------------------------------------------- */

ComputeClusterSizeExt::ComputeClusterSizeExt(LAMMPS* lmp, int narg, char** arg) :
    Compute(lmp, narg, arg), nloc(0), nc_global(0), natom_loc(0), nloc_peratom(0)
{
  vector_flag = 1;
  extvector = 0;
  size_vector = 0;
  size_vector_variable = 1;

  peratom_flag = 1;
  size_peratom_cols = 0;

  if (comm->nprocs > LMP_NUCC_CLUSTER_MAX_OWNERS) {
    error->all(FLERR, "{}: Number of processor exceeds MAX_OWNER limit. Recompile with higher MAX_OWNER limit.", style);
  }

  if (narg < 4) { utils::missing_cmd_args(FLERR, "compute size/cluster", error); }

  // Parse arguments //

  // Get cluster/atom compute
  compute_cluster_atom = lmp->modify->get_compute_by_id(arg[3]);
  if (compute_cluster_atom == nullptr) { error->all(FLERR, "{}: Cannot find compute with style 'cluster/atom' with given id: {}", style, arg[3]); }

  // Get the critical size
  size_cutoff = utils::inumeric(FLERR, arg[4], true, lmp);
  if (size_cutoff < 1) { error->all(FLERR, "{}: size_cutoff must be greater than 0", style); }

  size_vector = size_cutoff + 1;

  // keeper1 = new MemoryKeeper(memory);
  // cluster_map_allocator = new MapAlloc_t<int, int>(keeper1);
  // cluster_map = new Map_t<int, int>(*cluster_map_allocator);

  // keeper2 = new MemoryKeeper(memory);
  // alloc_map_vec1 = new MapAlloc_t<int, Vec_t<int>>(keeper2);
  // cIDs_by_size = new Map_t<int, Vec_t<int>>(*alloc_map_vec1);

  // keeper3 = new MemoryKeeper(memory);
  // alloc_map_vec2 = new MapAlloc_t<int, Vec_t<int>>(keeper3);
  // cIDs_by_size_all = new Map_t<int, Vec_t<int>>(*alloc_map_vec2);

  counts_global.create(memory, comm->nprocs, "size/cluster/ext:counts_global");
  displs.create(memory, comm->nprocs, "size/cluster/ext:displs");

  dist_local.create(memory, size_vector, "size/cluster/ext:dist_local");
  dist.create(memory, size_vector, "size/cluster/ext:dist");
  vector = dist.data();

  cIDs_by_size.reserve(size_cutoff);
  cIDs_by_size_all.reserve(size_cutoff);


  // MPI_Datatype type[2] = {MPI_INT, MPI_INT};
  // int blocklen[2] = {1, 1};
  MPI_Aint disp[2];

  // Calculate displacements
  disp[0] = offsetof(cldata, id);
  disp[1] = offsetof(cldata, sz);

  // MPI_Type_create_struct(2, blocklen, disp, type, &MPI_CLDATA);
  // MPI_Type_commit(&MPI_CLDATA);
}

/* ---------------------------------------------------------------------- */

ComputeClusterSizeExt::~ComputeClusterSizeExt() noexcept(true)
{
  // MPI_Type_free(&MPI_CLDATA);

  dist.destroy(memory);
  dist_local.destroy(memory);
  counts_global.destroy(memory);
  displs.destroy(memory);
  clusters.destroy(memory);
  ns.destroy(memory);
  gathered.destroy(memory);
  monomers.destroy(memory);

  // delete cluster_map;
  // delete cluster_map_allocator;
  // delete keeper1;

  // delete cIDs_by_size;
  // delete alloc_map_vec1;
  // delete keeper2;

  // delete cIDs_by_size_all;
  // delete alloc_map_vec2;
  // delete keeper3;
}

/* ---------------------------------------------------------------------- */

void ComputeClusterSizeExt::init()
{
  if ((modify->get_compute_by_style(style).size() > 1) && (comm->me == 0)) { error->warning(FLERR, "More than one compute {}", style); }

  if ((nloc < atom->nlocal) || (ns.empty() && clusters.empty() && monomers.empty())) {
    nloc = static_cast<int>(atom->nlocal * LMP_NUCC_ALLOC_COEFF);
    cluster_map.reserve(nloc);

    // keeper1->pool_size<MapMember_t<int, int>>(nloc);
    // keeper2->pool_size<MapMember_t<int, Vec_t<int>>>(nloc);
    // keeper3->pool_size<MapMember_t<int, Vec_t<int>>>(nloc);

    clusters.grow(memory, nloc, "size/cluster/ext:clusters");
    ns.grow(memory, nloc, "size/cluster/ext:ns");
    monomers.grow(memory, nloc, "size/cluster/ext:monomers");
    monomers.reset();
    ns.reset_unsafe<int>(0);
    for (int i = 0; i < nloc; ++i) {
      clusters[i] = cluster_data();
    }
  }

  if (ns.empty() || clusters.empty() || monomers.empty()) { error->one(FLERR, "{}: Inconsistent arrays state", style); }

  if ((gathered.empty()) || (natom_loc < atom->natoms)) {
    natom_loc = static_cast<bigint>(static_cast<long double>(atom->natoms) * LMP_NUCC_ALLOC_COEFF);
    gathered.grow(memory, natom_loc, "size/cluster/ext:gathered");
    gathered.reset_unsafe<int>(0);
  }

  if ((peratom_size.empty()) || (nloc_peratom < atom->nlocal)) {
    nloc_peratom = static_cast<int>(atom->nlocal * LMP_NUCC_ALLOC_COEFF);
    peratom_size.grow(memory, nloc_peratom, "size/cluster/ext:peratom");
    peratom_size.reset();
    vector_atom = peratom_size.data();
  }

  initialized_flag = 1;
}

/* ---------------------------------------------------------------------- */

void ComputeClusterSizeExt::compute_vector()
{
  invoked_vector = update->ntimestep;

  auto& cmap = cluster_map;
  auto& cbs = cIDs_by_size;
  auto& cbs_all = cIDs_by_size_all;

  if (compute_cluster_atom->invoked_peratom != update->ntimestep) { compute_cluster_atom->compute_peratom(); }

  const double* const cluster_ids = compute_cluster_atom->vector_atom;
  cmap.clear();
  cbs.clear();
  cbs_all.clear();
  counts_global.reset();
  displs.reset();
  dist_local.reset();
  dist.reset();

  if (nloc_peratom < atom->nlocal) {
    nloc_peratom = static_cast<int>(atom->nlocal * LMP_NUCC_ALLOC_COEFF);
    peratom_size.grow(memory, nloc_peratom, "size/cluster/ext:peratom");
    peratom_size.reset();
    vector_atom = peratom_size.data();
  }

  if (nloc < atom->nlocal) {
    nloc = static_cast<int>(atom->nlocal * LMP_NUCC_ALLOC_COEFF);
    cmap.reserve(nloc);

    clusters.grow(memory, nloc, "size/cluster/ext:clusters");
    ns.grow(memory, nloc, "size/cluster/ext:ns");
    monomers.grow(memory, nloc, "size/cluster/ext:monomers");
    monomers.reset();
    ns.reset_unsafe<int>(0);
    for (int i = 0; i < nloc; ++i) {
      clusters[i] = cluster_data();
    }
  }

  // Sort atom IDs by cluster IDs
  for (int i = 0; i < atom->nlocal; ++i) {
    if ((atom->mask[i] & groupbit) != 0) {
      const int clid = static_cast<int>(cluster_ids[i]);
      if (cmap.count(clid) == 0) {
        const int clidx = cmap.size();
        cmap[clid] = clidx;
        clusters[clidx] = cluster_data(clid);
      }
      cluster_data& clstr = clusters[cmap[clid]];
      if (clstr.l_size < LMP_NUCC_CLUSTER_MAX_SIZE) {
        clstr.atoms<false>()[clstr.l_size++] = i;
      } else {
        error->warning(FLERR, "{}: Cluster size exceeds max size", style);
      }
    }
  }

  for (const auto& [clid, clidx] : cmap) {
    ns[clidx] = cldata(clid, clusters[clidx].l_size);
    // cluster_data& clstr = clusters[clidx];
    // const auto clatoms = clstr.atoms();
    // for (int i = 0; i < clstr.l_size; ++i) {
    //   if (clatoms[i] >= atom->nlocal) { error->one(FLERR, "{}/compute_vector_1:{}: particle index exceeds nlocal", style, comm->me); }
    // }
  }

  // // add ghost atoms
  // for (int i = atom->nlocal; i < atom->nmax; ++i) {
  //   if ((atom->mask[i] & groupbit) != 0) {
  //     const auto clid = static_cast<int>(cluster_ids[i]);
  //     if (cmap.count(clid) > 0) {
  //       cluster_data& clstr = clusters[cmap[clid]];
  //       // also possible segfault if number of ghost exceeds LMP_NUCC_CLUSTER_MAX_GHOST
  //       clstr.ghost_initial()[clstr.nghost++] = i;
  //     }
  //   }
  // }

  // for (const auto& [clid, clidx] : cmap) {
  //   cluster_data& clstr = clusters[clidx];
  //   const auto clatoms = clstr.atoms();
  //   for (int i = 0; i < ns[2 * clidx + 1]; ++i) {
  //     if (clatoms[i] >= atom->nlocal) { error->one(FLERR, "{}/compute_vector_2:{}: particle index exceeds nlocal", style, comm->me); }
  //   }
  // }

  // communicate about number of unique clusters
  {
    const int ncluster_local = 2 * cmap.size();
    ::MPI_Allgather(&ncluster_local, 1, MPI_INT, counts_global.data(), 1, MPI_INT, world);


    bigint tcon = counts_global[0];
    for (int i = 1; i < comm->nprocs; ++i) {
      tcon += counts_global[i];
      displs[i] = displs[i - 1] + counts_global[i - 1];
    }

    const bigint t_alloc = tcon / 2 + 1;
    if (t_alloc > natom_loc) {
      natom_loc = static_cast<bigint>(t_alloc * LMP_NUCC_ALLOC_COEFF);
      gathered.grow(memory, natom_loc, "gathered");
      gathered.reset_unsafe<int>(0);
    }

    // communicate about local cluster sizes
    ::MPI_Allgatherv(ns.data(), ncluster_local, MPI_INT, gathered.data(), counts_global.data(), displs.data(), MPI_INT, world);
  }

  // fill local data
  for (int i = 0; i < comm->nprocs; ++i) {
    for (int j = 0; j < counts_global[i] / 2; ++j) {
      // int const k = displs[i] + j;
      const cldata& cl = gathered[displs[i] / 2 + j];
      if (cmap.count(cl.id) > 0) {
        cluster_data& clstr = clusters[cmap[cl.id]];
        if (i != comm->me) { clstr.owners<false>()[clstr.nowners++] = i; }
        clstr.g_size += cl.sz;
        if (cl.sz > clstr.nhost) {
          clstr.host = i;
          clstr.nhost = cl.sz;
        }
      }
    }
  }

  // adjust local data and fill local size distribution
  nonexclusive = 0;
  nmono = 0;

  for (const auto& [clid, clidx] : cmap) {
    cluster_data& clstr = clusters[clidx];
    // const auto clatoms = clstr.atoms();
    // for (int i = 0; i < clstr.l_size; ++i) {
    //   if (clatoms[i] >= atom->nlocal) { error->one(FLERR, "{}/compute_vector_3:{}: particle index exceeds nlocal", style, comm->me); }
    // }
    // clstr.rearrange();
    // for (int i = 0; i < clstr.l_size; ++i) {
    //   if (clatoms[i] >= atom->nlocal) { error->one(FLERR, "{}/compute_vector_4:{}: particle index exceeds nlocal", style, comm->me); }
    //   peratom_size[clatoms[i]] = clstr.g_size;
    // }
    if ((clstr.g_size < size_cutoff) && (clstr.g_size > 1)) { cbs_all[clstr.g_size].push_back(clidx); }
    if (clstr.host == comm->me) {
      if (clstr.g_size < size_cutoff) { dist_local[clstr.g_size] += 1; }
      if (clstr.g_size == 1) {
        monomers[nmono++] = clidx;
      } else {
        cbs[clstr.g_size].push_back(clidx);
      }
    }
    if (clstr.nowners > 0) { ++nonexclusive; }
  }

  ::MPI_Allreduce(dist_local.data(), dist.data(), size_vector, MPI_DOUBLE, MPI_SUM, world);

  #ifdef __NUCC_CHECK_ACCESS
  for (const auto& [clid, clidx]: cluster_map) {
    const auto& clstr = clusters[clidx];
    if (clstr.clid != clid) { error->one(FLERR, "{}: Cluster ID does not equals to map ID", style); }
    if (clstr.l_size > LMP_NUCC_CLUSTER_MAX_SIZE) {
      if ((clstr.l_size > 0) && (clstr.l_size < atom->nlocal)) { error->one(FLERR, "{}: Local cluster size exceeds max, seems to be atom id", style); }
      error->one(FLERR, "{}: Local cluster size exceeds max", style);
    }
    const auto& atoms = clstr.atoms();
    for (int i = 0; i < clstr.l_size; ++i) {
      if (atoms[i] > atom->nlocal) { error->one(FLERR, "{}: Atom indice exceeds nlocal", style); }
    }
  }
  if (comm->me == 0) { utils::logmesg(lmp, "{}: {}: Check passed\n", style, update->ntimestep); }
  #endif // __NUCC_CHECK_ACCESS
}

/* ---------------------------------------------------------------------- */

void ComputeClusterSizeExt::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  if (invoked_vector != update->ntimestep) { compute_vector(); }
}

/* ----------------------------------------------------------------------
   memory usage of maps and dist
------------------------------------------------------------------------- */

double ComputeClusterSizeExt::memory_usage()
{
  std::size_t sum = dist.memory_usage() + dist_local.memory_usage();
  sum += counts_global.memory_usage() + displs.memory_usage();
  sum += clusters.memory_usage();
  sum += ns.memory_usage() + gathered.memory_usage();
  sum += monomers.memory_usage();
  // sum += keeper1->memory_usage();
  // sum += keeper2->memory_usage();
  // sum += keeper3->memory_usage();
  return static_cast<double>(sum);
}

/* ---------------------------------------------------------------------- */

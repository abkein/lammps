/*
 ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org
------------------------------------------------------------------------- */

// TODO: NUCC FILE

#include "fix_cluster_delete.h"
#include "compute_cluster_size_ext.h"
#include "compute_cluster_temps.h"
#include "fmt/core.h"
#include "nucc_cspan.hpp"
#include "nucc_defs.hpp"

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
#include "modify.h"
#include "region.h"
#include "update.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <format>
#include <functional>
#include <unordered_map>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixClusterDelete::FixClusterDelete(LAMMPS* lmp, int narg, char** arg) : Fix(lmp, narg, arg)
{
  restart_pbc = 1;
  nevery      = 1;

  if (narg < 9) { utils::missing_cmd_args(FLERR, "fix cluster/delete", error); }

  // Parse arguments //

  // Get the max size of clusters
  kmax = utils::inumeric(FLERR, arg[4], true, lmp);
  if (kmax < 2) { error->all(FLERR, "{}: kmax cannot be less than 2", style); }

  // Get cluster/size compute
  compute_cluster_size = dynamic_cast<ComputeClusterSizeExt*>(lmp->modify->get_compute_by_id(arg[5]));
  if (compute_cluster_size == nullptr) { error->all(FLERR, "{}: Cannot find compute of style 'cluster/size' with id: {}", style, arg[5]); }
  if (kmax > compute_cluster_size->get_size_cutoff()) {
    error->all(FLERR, "{}: kmax cannot be bigger than its value of compute size/cluster", style);
  }

  // Parse optional keywords
  int iarg = 10;

  while (iarg < narg) {
    if (::strcmp(arg[iarg], "noscreen") == 0) {
      // Do not output to screen
      screenflag = 0;
      iarg += 1;
    } else if (::strcmp(arg[iarg], "screen") == 0) {
      // Do output to screen
      screenflag = 1;
      iarg += 1;

    } else if (::strcmp(arg[iarg], "file") == 0) {
      if (iarg + 2 > narg) { utils::missing_cmd_args(FLERR, std::format("{}: file", style), error); }
      if (comm->me == 0) {
        // Write output to file
        if (fileflag != 0) { error->one(FLERR, "{}: Both file and append keywords are present.", style); }
        fp = ::fopen(arg[iarg + 1], "w");
        if (fp == nullptr) { error->one(FLERR, "{}: Cannot open output file {}: {}", style, arg[iarg + 1], utils::getsyserror()); }
        fileflag = 1;
      }
      iarg += 2;
    } else if (::strcmp(arg[iarg], "append") == 0) {
      if (iarg + 2 > narg) { utils::missing_cmd_args(FLERR, std::format("{}: append", style), error); }
      if (comm->me == 0) {
        // Append output to file
        if (fileflag != 0) { error->one(FLERR, "{}: Both file and append keywords are present.", style); }
        fp = ::fopen(arg[iarg + 1], "a");
        if (fp == nullptr) { error->one(FLERR, "{}: Cannot open output file {}: {}", style, arg[iarg + 1], utils::getsyserror()); }
        fileflag = 1;
      }
      iarg += 2;
    } else if (::strcmp(arg[iarg], "nevery") == 0) {
      if (iarg + 2 > narg) { utils::missing_cmd_args(FLERR, std::format("{}: nevery", style), error); }
      // Get execution period
      nevery = utils::inumeric(FLERR, arg[iarg + 1], true, lmp);
      iarg += 2;
    } else {
      error->all(FLERR, "{}: Illegal command option {}", style, arg[iarg]);
    }
  }

  // further setup and error check

  // error check and further setup for mode = MOLECULE

  if (atom->tag_enable == 0) { error->all(FLERR, "{}: Cannot work unless atoms have IDs", style); }

  if ((comm->me == 0) && (fileflag != 0)) {
    fmt::print(fp, "ntimestep,ntotal,clusters_deleted,atoms_deleted\n");
    ::fflush(fp);
  }
}

/* ---------------------------------------------------------------------- */

FixClusterDelete::~FixClusterDelete() noexcept(true)
{
  if ((fp != nullptr) && (comm->me == 0)) {
    ::fflush(fp);
    ::fclose(fp);
  }
  count_a2m.destroy(memory);
  count_c2c.destroy(memory);
  ids_a2m.destroy(memory);
}

/* ---------------------------------------------------------------------- */

void FixClusterDelete::init()
{
  if ((modify->get_fix_by_style(style).size() > 1) && (comm->me == 0)) { error->warning(FLERR, "More than one fix {}", style); }
  if (domain->dimension != 3) { error->all(FLERR, "{}: Can work only in 3D.", style); }
  if (atom->molecular != Atom::ATOMIC) { error->all(FLERR, "{}: Cannot use with molecular systems (atom deletion does not update topology)", style); }

  count_a2m.create(memory, comm->nprocs, "cluster/delete:count_a2m");
  count_c2c.create(memory, comm->nprocs, "cluster/delete:count_c2c");

  nloc = atom->nlocal;
  ids_a2m.grow(memory, nloc, "cluster/delete:ids_a2m");
}

/* ---------------------------------------------------------------------- */

int FixClusterDelete::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixClusterDelete::pre_exchange()
{
  if (update->ntimestep < next_step) { return; }
  next_step = update->ntimestep + nevery;

  if (compute_cluster_size->invoked_vector != update->ntimestep) { compute_cluster_size->compute_vector(); }
  const auto& cIDs_by_size = compute_cluster_size->get_clid_by_size();

  if (nloc < atom->nlocal) {
    nloc = atom->nlocal;
    ids_a2m.grow(memory, nloc, "cluster/delete:ids_a2m");
    ids_a2m.reset();
  }
  count_c2c.reset();
  count_a2m.reset();
  // ids_a2m.reset();   // we don't care of freeing this array because it's overwritten from the beginning and we keep track of its actual (used) size

  // Count amount of local clusters to delete
  int clusters2delete_local = 0;
  // Count amount of local atoms to delete
  int atoms2delete_local    = 0;

  const int nclusters       = compute_cluster_size->get_cluster_map().size();
  const auto& clusters      = compute_cluster_size->get_clusters();
  for (int i = 0; i < nclusters; ++i) {
    const auto& cluster = clusters[i];
    if (cluster.g_size > kmax) {
      ++clusters2delete_local;
#ifndef __NUCC_ALGO_CHECK

      std::copy(cluster.atoms().data(), cluster.atoms().offset(cluster.l_size), ids_a2m.offset(atoms2delete_local));
      atoms2delete_local += cluster.l_size;
#else
      const auto cluster_atoms = cluster.atoms();
      for (int j = 0; j < cluster.l_size; ++j) {
        if (cluster_atoms[j] >= atom->nlocal) { error->one(FLERR, "{}/pre_exchange:{}: particle index exceeds nlocal", style, comm->me); }
        ids_a2m[atoms2move_local++] = cluster_atoms[j];
      }
#endif    // !__NUCC_ALGO_CHECK
    }
  }

  // sort to delete atoms from end to lower the number of copy opretions
  std::sort(ids_a2m.data(), ids_a2m.data() + atoms2delete_local, std::greater<>());
  ::MPI_Allgather(&clusters2delete_local, 1, MPI_INT, count_c2c.data(), 1, MPI_INT, world);
  ::MPI_Allgather(&atoms2delete_local, 1, MPI_INT, count_a2m.data(), 1, MPI_INT, world);

  int atoms2delete_total    = 0;
  int clusters2delete_total = 0;
  for (int proc = 0; proc < comm->nprocs; ++proc) {
    atoms2delete_total += count_a2m[proc];
    clusters2delete_total += count_c2c[proc];
  }

  if (atoms2delete_total > 0) { deleteAtoms(atoms2delete_total); }

  bigint nblocal = atom->nlocal;
  ::MPI_Allreduce(&nblocal, &atom->natoms, 1, MPI_LMP_BIGINT, MPI_SUM, world);

  if (comm->me == 0) {
    // print status
    if (screenflag != 0) { utils::logmesg(lmp, "Deleted {} clusters -> deleted {} atoms.\n", clusters2delete_total, atoms2delete_total); }
    if (fileflag != 0) {
      utils::print(fp, "{},{},{},{}\n", update->ntimestep, atom->natoms, clusters2delete_total, atoms2delete_total);
      ::fflush(fp);
    }
  }

}    // void FixClusterCrush::pre_exchange()

/* ---------------------------------------------------------------------- */

void FixClusterDelete::deleteAtoms(const int to_delete) const noexcept(true)
{
  // delete local atoms
  // reset nlocal

  const int n_atom_local_prev = atom->nlocal;
  for (int i = 0; i < to_delete; i++) {
#ifdef __NUCC_ALGO_CHECK
    if (atom->nlocal < 0) { error->one(FLERR, "{}/deleteAtoms:{}: Negative nlocal", style, comm->me); }
    if (ids_a2m[i] < 0) { error->one(FLERR, "{}/deleteAtoms:{}: particle index less than 0", style, comm->me); }
    if (ids_a2m[i] >= atom->nlocal) { error->one(FLERR, "{}/deleteAtoms:{}: particle index exceeds nlocal", style, comm->me); }
#endif    // __NUCC_ALGO_CHECK
    atom->avec->copy(n_atom_local_prev - 1 - i, ids_a2m[i], 1);
  }
  atom->nlocal -= to_delete;

  if (atom->molecular == Atom::ATOMIC) {
    tagint* const tag = atom->tag;
    const int nlocal  = atom->nlocal;
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

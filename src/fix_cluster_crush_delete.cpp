/*
 ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org
------------------------------------------------------------------------- */

// TODO: NUCC FILE

#include "fix_cluster_crush_delete.h"
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
#include "group.h"
#include "input.h"
#include "lattice.h"
#include "modify.h"
#include "random_park.h"
#include "region.h"
#include "update.h"
#include "variable.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <format>
#include <functional>
#include <unordered_map>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixClusterCrushDelete::FixClusterCrushDelete(LAMMPS* lmp, int narg, char** arg) : Fix(lmp, narg, arg)
{
  restart_pbc = 1;
  nevery      = 1;

  if (narg < 9) { utils::missing_cmd_args(FLERR, "fix cluster/crush/delete", error); }

  // Parse arguments //

  // Target region
  region = domain->get_region_by_id(arg[3]);
  if (region == nullptr) { error->all(FLERR, "{}: Cannot find target region {}", style, arg[3]); }
  if (region->bboxflag == 0) { error->all(FLERR, "{}: region does not support a bounding box", style); }
  // if (region->dynamic_check() != 0) { error->all(FLERR,"{}: region cannot be dynamic", style); }

  // Get the max size of clusters
  kmax = utils::inumeric(FLERR, arg[4], true, lmp);
  if (kmax < 2) { error->all(FLERR, "{}: kmax cannot be less than 2", style); }

  // Get cluster/size compute
  compute_cluster_size = dynamic_cast<ComputeClusterSizeExt*>(lmp->modify->get_compute_by_id(arg[5]));
  if (compute_cluster_size == nullptr) { error->all(FLERR, "{}: Cannot find compute of style 'cluster/size' with id: {}", style, arg[5]); }
  if (kmax > compute_cluster_size->get_size_cutoff()) {
    error->all(FLERR, "{}: kmax cannot be bigger than its value of compute size/cluster", style);
  }

  // Minimum distance to other atoms from the place atom teleports to
  overlap = utils::numeric(FLERR, arg[6], true, lmp);
  if (overlap < 0) { error->all(FLERR, "{}: Minimum distance must be non-negative", style); }
  overlapsq       = overlap * overlap;

  // Get the seed for coordinate generator
  const int xseed = utils::inumeric(FLERR, arg[7], true, lmp);
  xrandom         = new RanPark(lmp, xseed);

  // Get the seed for velocity generator
  const int vseed = utils::inumeric(FLERR, arg[8], true, lmp);
  vrandom         = new RanPark(lmp, vseed);

  // Get the type of atoms to create
  ntype           = utils::inumeric(FLERR, arg[9], true, lmp);
  if (ntype <= 0) { error->all(FLERR, "{}: invalid atom type: {}<1", style, ntype); }
  if (ntype > atom->ntypes) { error->all(FLERR, "{}: invalid atom type: {} is bigger than number of types exist: {}", style, ntype, atom->ntypes); }

  // Parse optional keywords
  int iarg    = 10;
  int velsset = 0;
  const char *vstr{}, *xstr{}, *ystr{}, *zstr{};

  while (iarg < narg) {
    if (::strcmp(arg[iarg], "maxtry") == 0) {
      if (iarg + 2 > narg) { utils::missing_cmd_args(FLERR, std::format("{}: maxtry", style), error); }
      // Max attempts to search for a new suitable location
      maxtry = utils::inumeric(FLERR, arg[iarg + 1], true, lmp);
      if (maxtry < 1) { error->all(FLERR, "{}: maxtry cannot be less than 1", style); }
      iarg += 2;
    } else if (::strcmp(arg[iarg], "temp") == 0) {
      if (iarg + 3 > narg) { utils::missing_cmd_args(FLERR, std::format("{}: temp", style), error); }
      // Atom temperature
      assign_temperature = true;
      if (::strcmp(arg[iarg + 1], "fix") == 0) {
        temp_fix         = true;
        atom_temperature = utils::numeric(FLERR, arg[iarg + 2], true, lmp);
        if (atom_temperature < 0) { error->all(FLERR, "{}: Atom temperature cannot be negative", style); }
      } else if (::strcmp(arg[iarg + 1], "track") == 0) {
        temp_fix = false;
        if (::strcmp(arg[iarg + 2], "avg") == 0) {
          temp_size = 0;
        } else {
          temp_size = utils::inumeric(FLERR, arg[iarg + 2], true, lmp);
          if (iarg + 4 > narg) { utils::missing_cmd_args(FLERR, std::format("{}: temp track {}", style, temp_size), error); }
          compute_cluster_temp = dynamic_cast<ComputeClusterTemp*>(lmp->modify->get_compute_by_id(arg[iarg + 3]));
          iarg += 1;
        }
      } else {
        error->all(FLERR, "{}: Unrecognized style of `temp` keyword: {}. Possible values are `fix`, `track`.", style, arg[iarg + 1]);
      }
      vdist = DIST::DIST_GAUSSIAN;
      iarg += 3;

    } else if (::strcmp(arg[iarg], "noscreen") == 0) {
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
    } else if (::strcmp(arg[iarg], "rate") == 0) {
      if (iarg + 2 > narg) { utils::missing_cmd_args(FLERR, std::format("{}: rate", style), error); }
      insertion_rate = utils::inumeric(FLERR, arg[iarg + 1], true, lmp);
      iarg += 2;
    } else if (::strcmp(arg[iarg], "keep_ss") == 0) {
      if (iarg + 3 > narg) { utils::missing_cmd_args(FLERR, std::format("{}: keep_ss", style), error); }
      supersaturation = utils::numeric(FLERR, arg[iarg + 1], true, lmp);
      compute_ss_mono = lmp->modify->get_compute_by_id(arg[iarg + 2]);
      if (compute_ss_mono == nullptr) {
        error->all(FLERR, "{}: Cannot find compute of style 'supersaturation/mono' with id: {}", style, arg[iarg + 2]);
      }
      iarg += 3;
    } else if (::strcmp(arg[iarg], "nevery") == 0) {
      if (iarg + 2 > narg) { utils::missing_cmd_args(FLERR, std::format("{}: nevery", style), error); }
      // Get execution period
      nevery = utils::inumeric(FLERR, arg[iarg + 1], true, lmp);
      iarg += 2;

    } else if (::strcmp(arg[iarg], "units") == 0) {
      if (iarg + 2 > narg) { utils::missing_cmd_args(FLERR, std::format("{}: units", style), error); }
      // TODO: Not sure if this is handled properly
      if (::strcmp(arg[iarg + 1], "box") == 0) {
        scaleflag = 0;
      } else if (::strcmp(arg[iarg + 1], "lattice") == 0) {
        scaleflag = 1;
      } else {
        error->all(FLERR, "{}: Unknown units option {}", style, arg[iarg + 1]);
      }
      iarg += 2;

    } else if (strcmp(arg[iarg], "var") == 0) {
      if (iarg + 2 > narg) { utils::missing_cmd_args(FLERR, std::format("{}: var", style), error); }
      vstr    = utils::strdup(arg[iarg + 1]);
      varflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg], "set") == 0) {
      if (iarg + 3 > narg) { utils::missing_cmd_args(FLERR, std::format("{}: var", style), error); }
      if (strcmp(arg[iarg + 1], "x") == 0) {
        xstr = utils::strdup(arg[iarg + 2]);
      } else if (strcmp(arg[iarg + 1], "y") == 0) {
        ystr = utils::strdup(arg[iarg + 2]);
      } else if (strcmp(arg[iarg + 1], "z") == 0) {
        zstr = utils::strdup(arg[iarg + 2]);
      } else {
        error->all(FLERR, "{}: Unknown set option {}", style, arg[iarg + 2]);
      }
      iarg += 3;
    } else if (strcmp(arg[iarg], "group") == 0) {
      if (iarg + 2 > narg) { utils::missing_cmd_args(FLERR, std::format("{}: group", style), error); }
      groupid = group->find(arg[iarg + 1]);
      if (groupid <= 0) { error->all(FLERR, "Specified group not found or group all is used"); }
      iarg += 2;

    } else if (strcmp(arg[iarg], "vx") == 0) {
      if (iarg + 3 > narg) { utils::missing_cmd_args(FLERR, std::format("{}: vx", style), error); }
      vels[0] = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      vels[1] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      ++velsset;
      iarg += 3;
    } else if (strcmp(arg[iarg], "vy") == 0) {
      if (iarg + 3 > narg) { utils::missing_cmd_args(FLERR, std::format("{}: vy", style), error); }
      vels[2] = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      vels[3] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      ++velsset;
      iarg += 3;
    } else if (strcmp(arg[iarg], "vz") == 0) {
      if (iarg + 3 > narg) { utils::missing_cmd_args(FLERR, std::format("{}: vz", style), error); }
      vels[4] = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      vels[5] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      ++velsset;
      iarg += 3;
    } else if (strcmp(arg[iarg], "inpoint") == 0) {
      if (iarg + 5 > narg) { utils::missing_cmd_args(FLERR, std::format("{}: inpoint", style), error); }
      // TODO: Implement
      xmid[0] = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      xmid[1] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      xmid[2] = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      xsigma  = utils::numeric(FLERR, arg[iarg + 4], false, lmp);
      xdist   = DIST::DIST_GAUSSIAN;
      iarg += 5;
    } else {
      error->all(FLERR, "{}: Illegal command option {}", style, arg[iarg]);
    }
  }

  // further setup and error check

  if (atom->mass_setflag[ntype] == 0) { error->all(FLERR, "{}: Atom mass for atom type {} is not set!", style, ntype); }

  if (assign_temperature) {
    if (temp_fix) {
      vsigma = ::sqrt(atom_temperature / atom->mass[ntype]);
    } else {
      if (temp_size == 0) {
        // Get temp compute
        const auto& temp_computes = lmp->modify->get_compute_by_style("temp");
        if (temp_computes.empty()) { error->all(FLERR, "{}: Cannot find compute with style 'temp'.", style); }
        compute_temp = temp_computes[0];
      }
    }
  } else if (velsset != 3) {
    error->all(FLERR, "{}: Either velocities or temperature are not set.", style);
  } else {
    // No way
  }

  sbonds[0] = region->extent_xlo;
  sbonds[1] = region->extent_xhi;
  sbonds[2] = region->extent_ylo;
  sbonds[3] = region->extent_yhi;
  sbonds[4] = region->extent_zlo;
  sbonds[5] = region->extent_zhi;

  if (domain->triclinic == 0) {
    if (sbonds[0] < domain->boxlo[0] || sbonds[1] > domain->boxhi[0] || sbonds[2] < domain->boxlo[1] || sbonds[3] > domain->boxhi[1] ||
        sbonds[4] < domain->boxlo[2] || sbonds[5] > domain->boxhi[2]) {
      error->all(FLERR, "{}: Deposition region extends outside simulation box", style);
    }
  } else {
    if (sbonds[0] < domain->boxlo_bound[0] || sbonds[1] > domain->boxhi_bound[0] || sbonds[2] < domain->boxlo_bound[1] ||
        sbonds[3] > domain->boxhi_bound[1] || sbonds[4] < domain->boxlo_bound[2] || sbonds[5] > domain->boxhi_bound[2]) {
      error->all(FLERR, "{}: Deposition region extends outside simulation box", style);
    }
  }

  // error check and further setup for mode = MOLECULE

  if (atom->tag_enable == 0) { error->all(FLERR, "{}: Cannot work unless atoms have IDs", style); }

  // apply scaling factor for styles that use distance-dependent factors
  // setup scaling

  if (scaleflag != 0) {
    const double xscale = domain->lattice->xlattice;
    const double yscale = domain->lattice->ylattice;
    const double zscale = domain->lattice->zlattice;
    // apply scaling to all input parameters with dist/vel units

    overlap *= xscale;
    overlapsq *= xscale * xscale;
    vels[0] *= xscale;
    vels[1] *= xscale;
    vels[2] *= yscale;
    vels[3] *= yscale;
    vels[4] *= zscale;
    vels[5] *= zscale;
    xmid[0] *= xscale;
    xmid[1] *= yscale;
    xmid[2] *= zscale;
    xsigma *= ::pow(xscale * yscale * zscale, 1. / 3.);    // same as in region sphere
    vsigma *= ::pow(xscale * yscale * zscale, 1. / 3.);    // same as in region sphere
  }

  if ((vstr == nullptr) && ((xstr != nullptr) || (ystr != nullptr) || (zstr != nullptr))) {
    error->all(FLERR, "{}: Incomplete use of variables", style);
  }
  if ((vstr != nullptr) && ((xstr == nullptr) && (ystr == nullptr) && (zstr == nullptr))) {
    error->all(FLERR, "{}: Incomplete use of variables", style);
  }

  if (varflag != 0) {
    vars[3] = input->variable->find(vstr);
    if (vars[3] < 0) { error->all(FLERR, "{}: Variable {} does not exist", style, vstr); }
    if (input->variable->equalstyle(vars[3]) == 0) { error->all(FLERR, "{}: Variable {} is invalid style", style, vstr); }

    if (xstr != nullptr) {
      vars[0] = input->variable->find(xstr);
      if (vars[0] < 0) { error->all(FLERR, "{}: Variable {} does not exist", style, xstr); }
      if (input->variable->internalstyle(vars[0]) == 0) { error->all(FLERR, "{}: Variable {} is invalid style", style, xstr); }
    }
    if (ystr != nullptr) {
      vars[1] = input->variable->find(ystr);
      if (vars[1] < 0) { error->all(FLERR, "{}: Variable {} does not exist", style, ystr); }
      if (input->variable->internalstyle(vars[1]) == 0) { error->all(FLERR, "{}: Variable {} is invalid style", style, ystr); }
    }
    if (zstr != nullptr) {
      vars[2] = input->variable->find(zstr);
      if (vars[2] < 0) { error->all(FLERR, "{}: Variable {} does not exist", style, zstr); }
      if (input->variable->internalstyle(vars[2]) == 0) { error->all(FLERR, "{}: Variable {} is invalid style", style, zstr); }
    }
  }
  delete[] vstr;
  delete[] xstr;
  delete[] ystr;
  delete[] zstr;

  if ((comm->me == 0) && (fileflag != 0)) {
    fmt::print(fp, "ntimestep,ntotal,c2c,a2m,moved,a2mn\n");
    ::fflush(fp);
  }
}

/* ---------------------------------------------------------------------- */

FixClusterCrushDelete::~FixClusterCrushDelete() noexcept(true)
{
  if ((fp != nullptr) && (comm->me == 0)) {
    ::fflush(fp);
    ::fclose(fp);
  }
  count_a2m.destroy(memory);
  count_c2c.destroy(memory);
  ids_a2m.destroy(memory);

  delete xrandom;
  delete vrandom;
}

/* ---------------------------------------------------------------------- */

void FixClusterCrushDelete::init()
{
  if ((modify->get_fix_by_style(style).size() > 1) && (comm->me == 0)) { error->warning(FLERR, "More than one fix {}", style); }
  if (domain->dimension != 3) { error->all(FLERR, "{}: Can work only in 3D.", style); }
  if (atom->molecular != Atom::ATOMIC) { error->all(FLERR, "{}: Cannot use with molecular systems (atom deletion does not update topology)", style); }

  count_a2m.create(memory, comm->nprocs, "cluster/crush/delete:count_a2m");
  count_c2c.create(memory, comm->nprocs, "cluster/crush/delete:count_c2c");

  nloc = atom->nlocal;
  ids_a2m.grow(memory, nloc, "cluster/crush/delete:ids_a2m");
}

/* ---------------------------------------------------------------------- */

int FixClusterCrushDelete::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixClusterCrushDelete::pre_exchange()
{
  if (update->ntimestep < next_step) { return; }
  next_step = update->ntimestep + nevery;

  if (compute_cluster_size->invoked_vector != update->ntimestep) { compute_cluster_size->compute_vector(); }
  const auto& cIDs_by_size = compute_cluster_size->get_clid_by_size();

  if (nloc < atom->nlocal) {
    nloc = atom->nlocal;
    ids_a2m.grow(memory, nloc, "cluster/crush/delete:ids_a2m");
    ids_a2m.reset();
  }
  count_c2c.reset();
  count_a2m.reset();
  // ids_a2m.reset();   // we don't care of freeing this array because it's overwritten from the beginning and we keep track of its actual (used) size

  // Count amount of local clusters to crush
  int clusters2crush_local = 0;
  // Count amount of local atoms to move
  int atoms2move_local     = 0;

  const int nclusters      = compute_cluster_size->get_cluster_map().size();
  const auto& clusters     = compute_cluster_size->get_clusters();
  for (int i = 0; i < nclusters; ++i) {
    const auto& cluster = clusters[i];
    if (cluster.g_size > kmax) {
      ++clusters2crush_local;
#ifndef __NUCC_ALGO_CHECK
      std::copy(cluster.atoms().data(), cluster.atoms().offset(cluster.l_size), ids_a2m.offset(atoms2move_local));
      atoms2move_local += cluster.l_size;
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
  std::sort(ids_a2m.data(), ids_a2m.data() + atoms2move_local, std::greater<>());

  // count_c2c[comm->me] = clusters2crush_local;    // unnecessary
  ::MPI_Allgather(&clusters2crush_local, 1, MPI_INT, count_c2c.data(), 1, MPI_INT, world);

  // count_a2m[comm->me] = atoms2move_local;    // unnecessary
  ::MPI_Allgather(&atoms2move_local, 1, MPI_INT, count_a2m.data(), 1, MPI_INT, world);

  int atoms2move_total     = 0;
  int clusters2crush_total = 0;
  for (int proc = 0; proc < comm->nprocs; ++proc) {
    atoms2move_total += count_a2m[proc];
    clusters2crush_total += count_c2c[proc];
  }

  if (clusters2crush_total > 0) { deleteAtoms(atoms2move_local); }

  balance += atoms2move_total;
  const int to_insert_prev = balance;

  if (balance > 0) {
    if (assign_temperature && (!temp_fix)) {
      if (temp_size == 0) {
        const double ts_temp = (compute_temp->invoked_scalar != update->ntimestep) ? compute_temp->compute_scalar() : compute_temp->scalar;
        vsigma               = ::sqrt(ts_temp / atom->mass[ntype]);
      } else {
        if (compute_cluster_temp->invoked_vector != update->ntimestep) { compute_cluster_temp->compute_vector(); }
        vsigma = ::sqrt(compute_cluster_temp->vector[temp_size] / atom->mass[ntype]);
      }
    }

    balance -= add(insertion_rate > 0 ? insertion_rate : balance);
  }

  bigint nblocal = atom->nlocal;
  ::MPI_Allreduce(&nblocal, &atom->natoms, 1, MPI_LMP_BIGINT, MPI_SUM, world);

  if (comm->me == 0) {
    // print status
    if (screenflag != 0) { utils::logmesg(lmp, "Crushed {} clusters -> deleted {} atoms.\n", clusters2crush_total, atoms2move_total); }
    if (fileflag != 0) {
      utils::print(fp, "{},{},{},{},{},{}\n", update->ntimestep, atom->natoms, clusters2crush_total, atoms2move_total, to_insert_prev - balance,
                   balance);
      ::fflush(fp);
    }
  }

}    // void FixClusterCrush::pre_exchange()

/* ---------------------------------------------------------------------- */

void FixClusterCrushDelete::deleteAtoms(const int atoms2move_local) const noexcept(true)
{
  // delete local atoms
  // reset nlocal

  for (int i = 0; i < atoms2move_local; i++) {
#ifdef __NUCC_ALGO_CHECK
    if (atom->nlocal < 0) { error->one(FLERR, "{}/deleteAtoms:{}: Negative nlocal", style, comm->me); }
    if (ids_a2m[i] < 0) { error->one(FLERR, "{}/deleteAtoms:{}: particle index less than 0", style, comm->me); }
    if (ids_a2m[i] >= atom->nlocal) { error->one(FLERR, "{}/deleteAtoms:{}: particle index exceeds nlocal", style, comm->me); }
#endif    // __NUCC_ALGO_CHECK
    atom->avec->copy(atom->nlocal - 1 - i, ids_a2m[i], 1);
  }
  atom->nlocal -= atoms2move_local;

  postDelete();
}

/* ---------------------------------------------------------------------- */

int FixClusterCrushDelete::add(const int to_insert) const
{
  std::array<double, 3> coord = {0, 0, 0};

  // clear ghost count (and atom map) and any ghost bonus data
  //   internal to AtomVec
  // same logic as beginning of Comm::exchange()
  // do it now b/c inserting atoms will overwrite ghost atoms

  if (atom->map_style != Atom::MAP_NONE) { atom->map_clear(); }
  atom->nghost = 0;
  atom->avec->clear_bonus();

  const bool not_triclinic  = domain->triclinic == 0;
  const double* const boxlo = not_triclinic ? static_cast<double*>(domain->boxlo) : static_cast<double*>(domain->boxlo_lamda);
  const double* const boxhi = not_triclinic ? static_cast<double*>(domain->boxhi) : static_cast<double*>(domain->boxhi_lamda);

  const double* const sublo = not_triclinic ? static_cast<double*>(domain->sublo) : static_cast<double*>(domain->sublo_lamda);
  const double* const subhi = not_triclinic ? static_cast<double*>(domain->subhi) : static_cast<double*>(domain->subhi_lamda);

  // find maxid in case other fixes deleted/inserted atoms


  tagint maxtag_all         = 0;
  {
    tagint max = 0;
    const tagint* const tag   = atom->tag;
    for (int i = 0; i < atom->nlocal; ++i) { max = std::max(max, tag[i]); }
    ::MPI_Allreduce(&max, &maxtag_all, 1, MPI_LMP_TAGINT, MPI_MAX, world);
  }

  region->prematch();

  int ninserted = 0;
  for (int added = 0; added < to_insert; ++added) {
    // attempt an insertion until successful

    int success = 0;
    int attempt = 0;
    while (attempt < maxtry) {
      success = 0;
      ++attempt;

      // generate new position and write it to coord (automatic check against region)
      gen_pos(coord);

      // check against variable
      if ((varflag != 0) && vartest(coord)) { continue; }

      std::array<double, 3> lamda = {0, 0, 0};
      domain->x2lamda(coord.data(), lamda.data());
      const std::array<double, 3>& newcoord = domain->triclinic == 0 ? coord : lamda;

      // check against box
      const bool proceed = newcoord[0] >= boxlo[0] && newcoord[0] < boxhi[0] && newcoord[1] >= boxlo[1] && newcoord[1] < boxhi[1] &&
          newcoord[2] >= boxlo[2] && newcoord[2] < boxhi[2];
      if (!proceed) { continue; }

      // check for overlapping
      if (check_overlap(coord) != 0) { continue; }

      const int placement_flag = placement_check_me(newcoord, sublo, subhi);

      ::MPI_Allreduce(&placement_flag, &success, 1, MPI_INT, MPI_SUM, world);
      if (success > 1) { error->all(FLERR, "{}: Multiple procs ({} procs) tried to insert an atom (seems to be a fix bug)", style, success); }

      // if ok, create atom and generate velocity
      if (success != 0) {
        create_atom(coord, maxtag_all + 1);
        break;
      }
    }

    if (success != 0) {
      ++atom->natoms;
      ++maxtag_all;
      ++ninserted;
      if (atom->natoms < 0) { error->all(FLERR, "{}: Too many total atoms", style); }
      if (maxtag_all >= MAXTAGINT) { error->all(FLERR, "{}: New atom IDs exceed maximum allowed ID", style); }
    }
  }

  // warn if not there were unsuccessful insertion attempts
  const int diff = to_insert - ninserted;
  if ((diff > 0) && (comm->me == 0)) { error->warning(FLERR, "{}: {} particle depositions were unsuccessful", style, diff); }

  // rebuild atom map

  if (atom->map_style != Atom::MAP_NONE) {
    atom->map_init();
    atom->map_set();
  }

  return ninserted;
}

/* ---------------------------------------------------------------------- */

int FixClusterCrushDelete::placement_check_me(const std::array<double, 3>& newcoord, const double* const sublo, const double* const subhi) const
{
  bool ok = newcoord[0] >= sublo[0] && newcoord[0] < subhi[0] && newcoord[1] >= sublo[1] && newcoord[1] < subhi[1] && newcoord[2] >= sublo[2] &&
      newcoord[2] < subhi[2];

  if (!ok) {
    if (domain->dimension == 3 && newcoord[2] >= domain->boxhi[2]) {
      if (comm->layout != Comm::LAYOUT_TILED) {
        ok = comm->myloc[2] == comm->procgrid[2] - 1 && newcoord[0] >= sublo[0] && newcoord[0] < subhi[0] && newcoord[1] >= sublo[1] &&
            newcoord[1] < subhi[1];

      } else {
        ok = comm->mysplit[2][1] == 1.0 && newcoord[0] >= sublo[0] && newcoord[0] < subhi[0] && newcoord[1] >= sublo[1] && newcoord[1] < subhi[1];
      }
    } else if (domain->dimension == 2 && newcoord[1] >= domain->boxhi[1]) {
      if (comm->layout != Comm::LAYOUT_TILED) {
        ok = comm->myloc[1] == comm->procgrid[1] - 1 && newcoord[0] >= sublo[0] && newcoord[0] < subhi[0];
      } else {
        ok = comm->mysplit[1][1] == 1.0 && newcoord[0] >= sublo[0] && newcoord[0] < subhi[0];
      }
    }
  }

  return static_cast<int>(ok);
}

/* ---------------------------------------------------------------------- */

void FixClusterCrushDelete::create_atom(const std::array<double, 3>& coord, bigint tag) const noexcept
{
  atom->avec->create_atom(ntype, const_cast<double*>(coord.data()));    // aaagrrrh;
  const int n    = atom->nlocal - 1;
  atom->tag[n]   = tag;
  atom->mask[n]  = 1 | groupbit;
  atom->image[n] = (static_cast<imageint>(IMGMAX) << IMG2BITS) | (static_cast<imageint>(IMGMAX) << IMGBITS) | IMGMAX;
  if (assign_temperature) {
    atom->v[n][0] = vrandom->gaussian() * vsigma;
    atom->v[n][1] = vrandom->gaussian() * vsigma;
    atom->v[n][2] = vrandom->gaussian() * vsigma;
  } else {
    atom->v[n][0] = vels[0] + vrandom->uniform() * (vels[1] - vels[0]);
    atom->v[n][1] = vels[2] + vrandom->uniform() * (vels[3] - vels[2]);
    atom->v[n][2] = vels[4] + vrandom->uniform() * (vels[5] - vels[4]);
  }
  modify->create_attribute(n);

  if (groupid != 0) { atom->mask[n] |= (1 << groupid); }
}

/* ---------------------------------------------------------------------- */

void FixClusterCrushDelete::gen_pos(std::array<double, 3>& coord /*, int nparticle, int nattempt*/) const noexcept
{
  // choose random position for new particle within region
  if (xdist == DIST::DIST_UNIFORM) {
    do {
      coord[0] = sbonds[0] + xrandom->uniform() * (sbonds[1] - sbonds[0]);
      coord[1] = sbonds[2] + xrandom->uniform() * (sbonds[3] - sbonds[2]);
      coord[2] = sbonds[4] + xrandom->uniform() * (sbonds[5] - sbonds[4]);
    } while (region->match(coord[0], coord[1], coord[2]) == 0);
  } else if (xdist == DIST::DIST_GAUSSIAN) {
    do {
      coord[0] = xmid[0] + xrandom->gaussian() * xsigma;
      coord[1] = xmid[1] + xrandom->gaussian() * xsigma;
      coord[2] = xmid[2] + xrandom->gaussian() * xsigma;
    } while (region->match(coord[0], coord[1], coord[2]) == 0);
  } else {
    error->all(FLERR, "{}: Unknown particle distribution", style);
  }
}

/* ---------------------------------------------------------------------- */

int FixClusterCrushDelete::check_overlap(const std::array<double, 3>& coord) const noexcept
{
  const double* const* const x = atom->x;

  int flag                     = 0;
  for (int i = 0; i < atom->nlocal; i++) {
    double delx = coord[0] - x[i][0];
    double dely = coord[1] - x[i][1];
    double delz = coord[2] - x[i][2];
    domain->minimum_image(FLERR, delx, dely, delz);
    const double rsq = delx * delx + dely * dely + delz * delz;
    if (rsq < overlapsq) {
      flag = 1;
      break;
    }
  }
  int flagall = 0;
  ::MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_MAX, world);
  return flagall;
}

/* ---------------------------------------------------------------------- */

void FixClusterCrushDelete::postDelete() const noexcept(true)
{
  if (atom->molecular == Atom::ATOMIC) {
    ::memset(atom->tag, 0, atom->nlocal * sizeof(int));
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

/* ----------------------------------------------------------------------
   test a generated atom position against variable evaluation
   first set x,y,z values in internal variables
------------------------------------------------------------------------- */

bool FixClusterCrushDelete::vartest(const std::array<double, 3>& coord) const noexcept
{
  input->variable->internal_set(vars[0], coord[0]);
  input->variable->internal_set(vars[1], coord[1]);
  input->variable->internal_set(vars[2], coord[2]);

  return input->variable->compute_equal(vars[3]) == 0.0;
}

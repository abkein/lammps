/*
 ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org
------------------------------------------------------------------------- */

// TODO: NUCC FILE

#include "dump_cluster_cf.h"
#include "dump.h"

#include "comm.h"
#include "compute.h"
#include "error.h"
#include "fmt/core.h"
#include "modify.h"
#include "update.h"

#include <cstdio>

using namespace LAMMPS_NS;

#define DUMP_FLOAT_PRECISION ",{:.5f}"

/* ---------------------------------------------------------------------- */

DumpClusterCF::DumpClusterCF(LAMMPS* lmp, int narg, char** arg) : Dump(lmp, narg, arg)
{
  if (narg < 7) { error->all(FLERR, "Illegal dump {} command", style); }
  clearstep            = 1;
  first_flag           = 0;
  time_flag            = 1;
  has_id               = 0;
  sort_flag            = 0;
  filewriter           = static_cast<int>(comm->me == 0);

  compute_cluster_size = modify->get_compute_by_id(arg[5]);
  if (compute_cluster_size == nullptr) { error->all(FLERR, "{}: Cannot find compute size/cluster with id: {}", style, arg[5]); }

  compute_cf = modify->get_compute_by_id(arg[6]);
  if (compute_cf == nullptr) { error->all(FLERR, "{}: Cannot find compute cf/cluster with id: {}", style, arg[6]); }

  if (filewriter > 0) {
    file = ::fopen(filename, "a");
    if (file == nullptr) { error->one(FLERR, "{}: Cannot open file {}: {}", style, filename, utils::getsyserror()); }
  }
}

/* ---------------------------------------------------------------------- */

DumpClusterCF::~DumpClusterCF()
{
  if ((filewriter != 0) && (file != nullptr)) {
    ::fflush(file);
    ::fclose(file);
  }
}

/* ---------------------------------------------------------------------- */

void DumpClusterCF::init_style()
{
  // Initialize the vector data if needed
}

/* ---------------------------------------------------------------------- */

void DumpClusterCF::pack(tagint* ids) {}

/* ---------------------------------------------------------------------- */

void DumpClusterCF::write()
{
  if (update->ntimestep < 10) { return; }
  if (compute_cluster_size->invoked_vector != update->ntimestep) { compute_cluster_size->compute_vector(); }

  if (compute_cf->invoked_array != update->ntimestep) { compute_cf->compute_array(); }

  const double* const dist         = compute_cluster_size->vector;
  const double* const* const array = compute_cf->array;

  if (filewriter != 0) {
    // fmt::print(file, "{}\n", update->ntimestep);
    for (int i = 1; i < compute_cf->size_array_rows; ++i) {
      if (dist[i] > 0) {
        fmt::print(file, "{}", i);
        for (int j = 0; j < compute_cf->size_array_cols; ++j) { fmt::print(file, DUMP_FLOAT_PRECISION, array[i][j]); }
        fmt::print(file, "\n");
      }
    }
    fmt::print(file, "\n");
  }
}

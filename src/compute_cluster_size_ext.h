/* -*- c++ -*- ----------------------------------------------------------
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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(size/cluster/ext,ComputeClusterSizeExt);
// clang-format on
#else

#ifndef LMP_COMPUTE_CLUSTER_SIZE_ExT_H
#define LMP_COMPUTE_CLUSTER_SIZE_ExT_H

#include "compute.h"
#include "nucc_allocator.hpp"
#include "nucc_cluster_data.hpp"
#include "nucc_cspan.hpp"
#include "nucc_defs.hpp"

#include <unordered_map>
#include <vector>

namespace NUCC {
struct cldata {
  cldata() = default;
  cldata(int id, int sz) : id(id), sz(sz) {}
  int id = 0;
  int sz = 0;
};
}    // namespace NUCC

namespace LAMMPS_NS {

class ComputeClusterSizeExt : public Compute {
 public:
  ComputeClusterSizeExt(class LAMMPS* lmp, int narg, char** arg);
  ~ComputeClusterSizeExt() noexcept(true) override;
  void init() override;
  void compute_vector() override;
  void compute_peratom() override;
  void compute_local() override;
  double memory_usage() override;

  constexpr int get_size_cutoff() const noexcept(true) { return size_cutoff; }
  constexpr NUCC::cspan<const double> get_data() const noexcept { return NUCC::cspan<const double>(dist); }
  constexpr int get_nonexclusive() const noexcept(true) { return nonexclusive; }
  constexpr const std::unordered_map<int, int>& get_cluster_map() const noexcept(true) { return cluster_map; }
  //   inline constexpr const NUCC::Map_t<int, int> *get_cluster_map() const noexcept(true) { return cluster_map; }
  constexpr const std::unordered_map<int, std::vector<int>>& get_clid_by_size() const noexcept(true) { return clid_by_size; }
  //   inline constexpr const NUCC::Map_t<int, NUCC::Vec_t<int>> *get_cIDs_by_size_my() const noexcept { return cIDs_by_size; }
  constexpr const std::unordered_map<int, std::vector<int>>& get_clid_by_size_global() const noexcept(true) { return clid_by_size_global; }
  //   inline constexpr const NUCC::Map_t<int, NUCC::Vec_t<int>> *get_cIDs_by_size() const noexcept { return cIDs_by_size_all; }
  constexpr NUCC::cspan<const NUCC::cluster_data> get_clusters() const noexcept(true) { return NUCC::cspan<const NUCC::cluster_data>(clusters); }

 private:
  int size_cutoff;    // number of elements reserved in resulting distribution, i.e. the max size of cluster distribution will be built to

  //   NUCC::MemoryKeeper *keeper1;
  //   NUCC::MapAlloc_t<int, int> *cluster_map_allocator;
  //   NUCC::Map_t<int, int> *cluster_map;
  std::unordered_map<int, int> cluster_map;    // mapping from cluster id to its local indices (idx) in `clusters` and `ns`

  //   NUCC::MemoryKeeper *keeper2;
  //   NUCC::MapAlloc_t<int, NUCC::Vec_t<int>> *alloc_map_vec1;
  //   NUCC::Map_t<int, NUCC::Vec_t<int>> *cIDs_by_size;
  std::unordered_map<int, std::vector<int>> clid_by_size;    // mapping from cluster size to vector of local indices (vector(idx))

  //   NUCC::MemoryKeeper *keeper3;
  //   NUCC::MapAlloc_t<int, NUCC::Vec_t<int>> *alloc_map_vec2;
  //   NUCC::Map_t<int, NUCC::Vec_t<int>> *cIDs_by_size_all;
  std::unordered_map<int, std::vector<int>> clid_by_size_global;    // don't know if it's usable at all, maybe delete

  int nloc = 0;                                // possible count of clusters (nlocal*damp), size of allocated arrays
  NUCC::cspan<double> dist;                    // int, cluster size distribution (vector == dist)
  NUCC::cspan<double> dist_local;              // int, local cluster size distribution
  int nc_global = 0;                           // actual number of clusters (across all procs)
  NUCC::cspan<int> counts_global;              // counts for MPI communication
  NUCC::cspan<int> displs;                     // displacements for MPI communication
  NUCC::cspan<NUCC::cluster_data> clusters;    // cluster data
  NUCC::cspan<NUCC::cldata>
      ns;    // array of cluster ids and their sizes. same as `clusters`, but containing less data; used to communicate with other procs
  NUCC::cspan<NUCC::cldata> gathered;
  NUCC::cspan<double> peratom_size;    // int, size of cluster this atom belongs to
  bigint nloc_gather = 0;
  int nonexclusive   = 0;    // number of clusters not owned by this proc
  int nloc_peratom   = 0;    // number of elements allocated for peratom array

  int nmono          = 0;       // number of local monomers
  NUCC::cspan<double> monomers;    // int, local indices of monomers

  Compute* compute_cluster_atom = nullptr;

  // MPI_Datatype MPI_CLDATA;
};

}    // namespace LAMMPS_NS

#endif
#endif

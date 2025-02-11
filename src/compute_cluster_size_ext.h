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

#  ifndef LMP_COMPUTE_CLUSTER_SIZE_ExT_H
#    define LMP_COMPUTE_CLUSTER_SIZE_ExT_H

#    include "compute.h"
#    include "nucc_allocator.hpp"
#    include "nucc_defs.hpp"
#    include "nucc_cluster_data.hpp"
#    include "nucc_cspan.hpp"

#    include <array>
#    include <scoped_allocator>
#    include <span>
#    include <unordered_map>
#    include <vector>

namespace NUCC {
  struct cldata {
    cldata() = default;
    cldata(int id, int sz): id(id), sz(sz) {}
    int id = 0;
    int sz = 0;
  };
}

namespace LAMMPS_NS {

class ComputeClusterSizeExt : public Compute {
 public:
  ComputeClusterSizeExt(class LAMMPS *lmp, int narg, char **arg);
  ~ComputeClusterSizeExt() noexcept(true) override;
  void init() override;
  void compute_vector() override;
  void compute_peratom() override;
  double memory_usage() override;

  inline constexpr int get_size_cutoff() const noexcept(true) { return size_cutoff; }
  inline constexpr NUCC::cspan<const double> get_data() const noexcept { return dist; }
  inline constexpr int get_nonexclusive() const noexcept(true) { return nonexclusive; }
  inline constexpr const std::unordered_map<int, int> &get_cluster_map() const noexcept(true) { return cluster_map; }
  inline constexpr const std::unordered_map<int, std::vector<int>> &get_cIDs_by_size() const noexcept(true) { return cIDs_by_size; }
  inline constexpr const std::unordered_map<int, std::vector<int>> &get_cIDs_by_size_all() const noexcept(true) { return cIDs_by_size_all; }
//   inline constexpr const NUCC::Map_t<int, int> *get_cluster_map() const noexcept(true) { return cluster_map; }
//   inline constexpr const NUCC::Map_t<int, NUCC::Vec_t<int>> *get_cIDs_by_size_my() const noexcept { return cIDs_by_size; }
//   inline constexpr const NUCC::Map_t<int, NUCC::Vec_t<int>> *get_cIDs_by_size() const noexcept { return cIDs_by_size_all; }
  inline constexpr const NUCC::cspan<const NUCC::cluster_data> get_clusters() const noexcept(true) { return clusters; }

 private:
  int size_cutoff;    // number of elements reserved in dist

//   NUCC::MemoryKeeper *keeper1;
//   NUCC::MapAlloc_t<int, int> *cluster_map_allocator;
//   NUCC::Map_t<int, int> *cluster_map;
  std::unordered_map<int, int> cluster_map;    // clid -> idx

//   NUCC::MemoryKeeper *keeper2;
//   NUCC::MapAlloc_t<int, NUCC::Vec_t<int>> *alloc_map_vec1;
//   NUCC::Map_t<int, NUCC::Vec_t<int>> *cIDs_by_size;
  std::unordered_map<int, std::vector<int>> cIDs_by_size;    // size -> vector(idx)

//   NUCC::MemoryKeeper *keeper3;
//   NUCC::MapAlloc_t<int, NUCC::Vec_t<int>> *alloc_map_vec2;
//   NUCC::Map_t<int, NUCC::Vec_t<int>> *cIDs_by_size_all;
  std::unordered_map<int, std::vector<int>> cIDs_by_size_all;

  int nloc = 0;                      // number of reserved elements in atoms_by_cID and cIDs_by_size
  NUCC::cspan<double> dist;          // cluster size distribution (vector == dist)
  NUCC::cspan<double> dist_local;    // local cluster size distribution
  int nc_global = 0;                 // number of clusters total
  NUCC::cspan<int> counts_global;
  NUCC::cspan<int> displs;
  NUCC::cspan<NUCC::cluster_data> clusters;
  NUCC::cspan<NUCC::cldata> ns;
  NUCC::cspan<NUCC::cldata> gathered;
  NUCC::cspan<double> peratom_size;
  bigint natom_loc = 0;
  int nonexclusive = 0;
  int nloc_peratom = 0;

  NUCC::cspan<int> monomers;
  int nmono = 0;

  Compute *compute_cluster_atom = nullptr;

  // MPI_Datatype MPI_CLDATA;
};

}    // namespace LAMMPS_NS

#  endif
#endif

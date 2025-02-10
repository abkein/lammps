#ifndef __NUCC_CLUSTER_DATA_HPP
#define __NUCC_CLUSTER_DATA_HPP

#include "nucc_cspan.hpp"

namespace NUCC {

struct cluster_data {
  cluster_data(): clid(0), l_size(0), g_size(0), host(-1), nhost(0), nowners(0) {
    std::fill_n(_owners, LMP_NUCC_CLUSTER_MAX_OWNERS, 0);
    std::fill_n(_atoms, LMP_NUCC_CLUSTER_MAX_SIZE, 0);
  }

  explicit cluster_data(const int _clid): clid(_clid), l_size(0), g_size(0), host(-1), nhost(0), nowners(0) {
    std::fill_n(_owners, LMP_NUCC_CLUSTER_MAX_OWNERS, 0);
    std::fill_n(_atoms, LMP_NUCC_CLUSTER_MAX_SIZE, 0);
  }

  // void rearrange() noexcept { ::memcpy(_atoms + l_size, _ghost, nghost * sizeof(int)); }

  // NUCC::cspan<const int> atoms_all() const { return std::span<const int>(_atoms, l_size + nghost); }

  template <bool Protect = true>
    requires(!Protect)
  NUCC::cspan<int, LMP_NUCC_CLUSTER_MAX_SIZE> atoms()
  {
    return std::span<int, LMP_NUCC_CLUSTER_MAX_SIZE>(_atoms, LMP_NUCC_CLUSTER_MAX_SIZE);
  }

  NUCC::cspan<const int> atoms() const { return std::span<const int>(_atoms, l_size); }

  // NUCC::cspan<int, LMP_NUCC_CLUSTER_MAX_GHOST> ghost_initial() { return std::span<int, LMP_NUCC_CLUSTER_MAX_GHOST>(_ghost, nghost); }

  // NUCC::cspan<const int> ghost() const { return std::span<const int>(_atoms + l_size, nghost); }

  template <bool Protect = true>
    requires(!Protect)
  NUCC::cspan<int, LMP_NUCC_CLUSTER_MAX_OWNERS> owners()
  {
    return std::span<int, LMP_NUCC_CLUSTER_MAX_OWNERS>(_owners, LMP_NUCC_CLUSTER_MAX_OWNERS);
  }

  NUCC::cspan<const int> owners() const { return std::span<const int>(_owners, nowners); }

  int clid = 0;       // cluster ID
  int l_size = 0;     // local size
  int g_size = 0;     // global size
  int host = -1;      // host proc (me if <0)
  int nhost = 0;      // local cluster size of host proc
  int nowners = 0;    // number of owners
  // int nghost = 0;     // number of ghost atoms in cluster
 private:
  int _owners[LMP_NUCC_CLUSTER_MAX_OWNERS];    // procs owning some cluster's atoms
  int _atoms[LMP_NUCC_CLUSTER_MAX_SIZE];       // local ids of atoms
  // int _ghost[LMP_NUCC_CLUSTER_MAX_GHOST];      // local ids of ghost atoms
};

}    //  namespace NUCC

#endif    // !__NUCC_CLUSTER_DATA_HPP

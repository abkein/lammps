#ifndef NUCC_CLUSTER_DATA_HPP
#define NUCC_CLUSTER_DATA_HPP

#include "nucc_cspan.hpp"
#include "nucc_defs.hpp"

#include <algorithm>
#include <span>

namespace NUCC {

struct cluster_data {
  cluster_data() {
    std::fill_n(static_cast<int *>(_owners), LMP_NUCC_CLUSTER_MAX_OWNERS, 0);
    std::fill_n(static_cast<int *>(_atoms), LMP_NUCC_CLUSTER_MAX_SIZE, 0);
    std::fill_n(static_cast<int *>(_ghost), LMP_NUCC_CLUSTER_MAX_GHOST, 0);
  }

  explicit cluster_data(const int _clid): clid(_clid) {
    std::fill_n(static_cast<int *>(_owners), LMP_NUCC_CLUSTER_MAX_OWNERS, 0);
    std::fill_n(static_cast<int *>(_atoms), LMP_NUCC_CLUSTER_MAX_SIZE, 0);
    std::fill_n(static_cast<int *>(_ghost), LMP_NUCC_CLUSTER_MAX_GHOST, 0);
  }

  // void rearrange() noexcept { ::memcpy(_atoms + l_size, _ghost, (LMP_NUCC_CLUSTER_MAX_SIZE - l_size) * sizeof(int)); }

  // NUCC::cspan<const int> atoms_all() const { return std::span<const int>(_atoms, l_size + nghost); }

  // return a span rendering the underlying array, allowing changes
  template <bool Protect = true>
    requires(!Protect)
  [[nodiscard]] NUCC::cspan<int, LMP_NUCC_CLUSTER_MAX_SIZE> atoms()
  {
    return std::span<int, LMP_NUCC_CLUSTER_MAX_SIZE>(static_cast<int *>(_atoms), LMP_NUCC_CLUSTER_MAX_SIZE);
  }

  // return a span rendering the underlying array, const version
  [[nodiscard]] NUCC::cspan<const int> atoms() const { return std::span<const int>(static_cast<const int *>(_atoms), l_size); }

  // return a span rendering the underlying array, allowing changes
  template <bool Protect = true>
    requires(!Protect)
  [[nodiscard]] NUCC::cspan<int, LMP_NUCC_CLUSTER_MAX_SIZE> ghost()
  {
    return std::span<int, LMP_NUCC_CLUSTER_MAX_SIZE>(static_cast<int *>(_ghost), LMP_NUCC_CLUSTER_MAX_GHOST);
  }

  // return a span rendering the underlying array, const version
  [[nodiscard]] NUCC::cspan<const int> ghost() const { return std::span<const int>(static_cast<const int *>(_ghost), LMP_NUCC_CLUSTER_MAX_GHOST); }

  // return a span rendering the underlying array, allowing changes
  template <bool Protect = true>
    requires(!Protect)
  [[nodiscard]] NUCC::cspan<int, LMP_NUCC_CLUSTER_MAX_OWNERS> owners()
  {
    return std::span<int, LMP_NUCC_CLUSTER_MAX_OWNERS>(static_cast<int *>(_owners), LMP_NUCC_CLUSTER_MAX_OWNERS);
  }

  // return a span rendering the underlying array, const version
  [[nodiscard]] NUCC::cspan<const int> owners() const { return std::span<const int>(static_cast<const int *>(_owners), nowners); }

  int clid = 0;       // cluster ID
  int l_size = 0;     // local size
  int g_size = 0;     // global size
  int host = -1;      // host proc (me if <0)
  int nhost = 0;      // local cluster size of host proc
  int nowners = 0;    // number of owners
  int nghost = 0;     // number of ghost atoms in cluster

 private:
  int _owners[LMP_NUCC_CLUSTER_MAX_OWNERS];    // procs owning some cluster's atoms
  int _atoms[LMP_NUCC_CLUSTER_MAX_SIZE];       // local ids of atoms
  int _ghost[LMP_NUCC_CLUSTER_MAX_GHOST];      // local ids of ghost atoms
};

}    //  namespace NUCC

#endif    // !__NUCC_CLUSTER_DATA_HPP

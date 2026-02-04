#ifndef NUCC_CLUSTER_DATA_HPP
#define NUCC_CLUSTER_DATA_HPP

#include "nucc_cspan.hpp"
#include "nucc_defs.hpp"

#include <array>

namespace NUCC {

using atoms_t = std::array<int, LMP_NUCC_CLUSTER_MAX_SIZE>;
using owners_t = std::array<int, LMP_NUCC_CLUSTER_MAX_OWNERS>;
using ghost_t = std::array<int, LMP_NUCC_CLUSTER_MAX_GHOST>;

struct cluster_data {
  cluster_data() {
    _owners.fill(0);
    _atoms.fill(0);
    _ghost.fill(0);
  }

  explicit cluster_data(const int _clid): clid(_clid) {
    _owners.fill(0);
    _atoms.fill(0);
    _ghost.fill(0);
  }

  // void rearrange() noexcept { ::memcpy(_atoms + l_size, _ghost, (LMP_NUCC_CLUSTER_MAX_SIZE - l_size) * sizeof(int)); }

  // NUCC::cspan<const int> atoms_all() const { return std::span<const int>(_atoms, l_size + nghost); }

  // return a span rendering the underlying array, allowing changes
  template <bool Protect = true>
    requires(!Protect)
  [[nodiscard]] cspan<int, LMP_NUCC_CLUSTER_MAX_SIZE> atoms() {
    return cspan<int, LMP_NUCC_CLUSTER_MAX_SIZE>(_atoms.data(), LMP_NUCC_CLUSTER_MAX_SIZE);
  }

  // return a span rendering the underlying array, const version
  [[nodiscard]] cspan<const int> atoms() const {
    return cspan<const int>(static_cast<const int *>(_atoms.data()), l_size);
  }

  // return a span rendering the underlying array, allowing changes
  template <bool Protect = true>
    requires(!Protect)
  [[nodiscard]] cspan<int, LMP_NUCC_CLUSTER_MAX_SIZE> ghost()
  {
    return cspan<int, LMP_NUCC_CLUSTER_MAX_SIZE>(static_cast<int *>(_ghost.data()), LMP_NUCC_CLUSTER_MAX_GHOST);
  }

  // return a span rendering the underlying array, const version
  [[nodiscard]] cspan<const int> ghost() const {
    return cspan<const int>(static_cast<const int *>(_ghost.data()), nghost);
  }

  // return a span rendering the underlying array, allowing changes
  template <bool Protect = true>
    requires(!Protect)
  [[nodiscard]] cspan<int, LMP_NUCC_CLUSTER_MAX_OWNERS> owners()
  {
    return cspan<int, LMP_NUCC_CLUSTER_MAX_OWNERS>(static_cast<int *>(_owners.data()), LMP_NUCC_CLUSTER_MAX_OWNERS);
  }

  // return a span rendering the underlying array, const version
  [[nodiscard]] cspan<const int> owners() const {
    return cspan<const int>(static_cast<const int *>(_owners.data()), nowners);
  }

  int clid = 0;       // cluster ID
  int l_size = 0;     // local size
  int g_size = 0;     // global size
  int host = -1;      // host proc (me if <0)
  int nhost = 0;      // local cluster size of host proc
  int nowners = 0;    // number of owners
  int nghost = 0;     // number of ghost atoms in cluster

 private:
  atoms_t  _atoms;       // local ids of atoms
  owners_t _owners;    // procs owning some cluster's atoms
  ghost_t  _ghost;      // local ids of ghost atoms
};

}    //  namespace NUCC

#endif    // !__NUCC_CLUSTER_DATA_HPP

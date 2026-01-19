#ifndef __NUCC_CUSTOM_STL_ALLOCATOR_HPP
#define __NUCC_CUSTOM_STL_ALLOCATOR_HPP

#include "memory.h"

#include <cassert>
#include <cstddef>
#include <memory>
#include <type_traits>
#include <vector>
#include <numeric>

namespace NUCC {

class MemoryKeeper {
  //  private:
  //   struct PoolInfo {
  //     constexpr PoolInfo(char* ptr, const size_t size) noexcept : ptr(ptr), size(size) {}

  //     PoolInfo()  = delete;

  //     char* ptr   = nullptr;
  //     size_t size = 0;
  //   };

 public:
  MemoryKeeper()                               = delete;
  MemoryKeeper(const MemoryKeeper&)            = delete;
  MemoryKeeper(MemoryKeeper&&)                 = delete;
  MemoryKeeper& operator=(const MemoryKeeper&) = delete;
  MemoryKeeper& operator=(MemoryKeeper&&)      = delete;

  explicit MemoryKeeper(LAMMPS_NS::Memory* memory) noexcept: memory_(memory) {}
  ~MemoryKeeper() { clear(); }

  template <typename T>
  void store(T*& ptr, const size_t size) noexcept(noexcept(std::declval<std::vector<std::pair<char*, std::size_t>>>().emplace_back(ptr, size)))
  {
    infos.emplace_back(ptr, size);
  }

  void clear() noexcept(noexcept(std::declval<LAMMPS_NS::Memory>().destroy<void>(std::declval<void*&>())))
  {
    for (auto& pool : infos) { memory_->destroy(pool.first); }
  }

  inline constexpr std::size_t pool_size() const noexcept
  {
    assert(_pool_size > 0);
    return _pool_size;
  }

  template <typename T>
  inline constexpr void pool_size(std::size_t n) noexcept
  {
    _pool_size = n * sizeof(T);
  }

  inline constexpr void pool_size(std::size_t n) noexcept { _pool_size = n; }

  char* allocate(const std::size_t nbytes)
  {
    char* ptr = nullptr;
    if (nbytes > _pool_size) {
      // If requested size is larger than pool, allocate separately
      memory_->create<char>(ptr, nbytes, "CustomAllocator_Large");
      infos.emplace_back(std::make_pair(ptr, nbytes));
      return ptr;
    }
    if ((current == nullptr) || (left < nbytes)) {
      // Pool is full or not initialized, request a new pool
      memory_->create<char>(current, _pool_size, "CustomAllocator_Pool");
      infos.emplace_back(std::make_pair(current, _pool_size));
      left = _pool_size;
    }
    ptr = current;
    left -= nbytes;
    current = current + nbytes;
    return ptr;
  }

  std::size_t memory_usage() {
    return std::accumulate(
      infos.begin(), infos.end(), std::size_t{0},
      [](std::size_t acc, const auto& pair) { return acc + pair.second; }
    );
  }

 private:
  char* current                    = nullptr;
  std::size_t left                 = 0;
  LAMMPS_NS::Memory* const memory_ = nullptr;
  std::size_t _pool_size           = 0;
  std::vector<std::pair<char*, std::size_t>> infos;
};

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

template <typename T>
class CustomAllocator {
 public:
  using value_type = T;

  template <typename U>
  friend class CustomAllocator;

  CustomAllocator() = delete;

  constexpr explicit CustomAllocator(MemoryKeeper* const keeper) noexcept : keeper_(keeper) {}

  template <typename U>
  constexpr CustomAllocator(const CustomAllocator<U>& other) noexcept : keeper_(other.keeper_)
  {
  }

  inline T* allocate(std::size_t n) const { return reinterpret_cast<T*>(keeper_->allocate(n * sizeof(T))); }

  inline constexpr void deallocate(T* /*p*/, const std::size_t /*n*/) const noexcept
  {
    // Deallocation can be handled when the allocator is destroyed
    // For pool allocator, individual deallocations are often no-ops
  }

  // Equality operators
  template <typename U>
  inline constexpr bool operator==(const CustomAllocator<U>& other) const noexcept
  {
    return keeper_ == other.keeper_;
  }

  template <typename U>
  inline constexpr bool operator!=(const CustomAllocator<U>& other) const noexcept
  {
    return !(*this == other);
  }

 private:
  MemoryKeeper* const keeper_;
};

}    // namespace NUCC

#endif    // CUSTOM_ALLOCATOR_HPP

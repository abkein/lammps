#ifndef NUCC_DEFS_HPP
#define NUCC_DEFS_HPP

#include <cstddef>
#include <cstdlib>
#include <type_traits>
// #include <scoped_allocator>
// #include <unordered_map>
// #include <vector>

#define LMP_NUCC_HDR_CSPAN_CHECK

namespace NUCC {

#if defined(LMP_NUCC_HDR_CSPAN_CHECK)

#if defined(__clang__)
#define DBG_NOINLINE __attribute__((noinline))
#define DBG_OPTNONE __attribute__((optnone))
#elif defined(__GNUC__)
#define DBG_NOINLINE __attribute__((noinline))
#define DBG_OPTNONE __attribute__((optimize("O0")))
#elif defined(_MSC_VER)
#define DBG_NOINLINE __declspec(noinline)
#define DBG_OPTNONE
#else
#define DBG_NOINLINE
#define DBG_OPTNONE
#endif
#else
#define DBG_NOINLINE
#define DBG_OPTNONE
#endif

namespace Defines {
  inline constexpr double ALLOC_COEFF                 = 1.2;
  inline constexpr int CLUSTER_MAX_OWNERS             = 128;
  inline constexpr int CLUSTER_MAX_SIZE               = 300;
  inline constexpr int CLUSTER_MAX_GHOST              = 300;
  inline constexpr int NEIGH_BIN_CUTOFF_COEFF         = 2;
  inline constexpr bool CHECK_ACCESS                  = true;
  inline constexpr bool ALGO_CHECK                    = true;
  inline constexpr bool NEIGHS_RADIAL_USE_HALF        = false;
  inline constexpr bool NEIGHS_RADIAL_PRECOMPUTE_NORM = true;

  inline constexpr bool CSPAN_CHECK_ACCESS      = true;
  inline constexpr bool CSPAN_DEBUG_CALLS       = true;
  inline constexpr bool ESPAN_CHECK_ACCESS      = true;
  inline constexpr bool ESPAN_DEBUG_CALLS       = true;
  inline constexpr bool ESPAN_CHECK_DEREFERENCE = true;
}    // namespace Defines



DBG_NOINLINE DBG_OPTNONE static void debug_trap() noexcept
{
#if defined(_MSC_VER)
  __debugbreak();
#elif defined(__clang__)
  __builtin_debugtrap();
#else
  __builtin_trap();
#endif
  std::abort();
}

DBG_NOINLINE DBG_OPTNONE static void debug_check_index(std::size_t i, std::size_t n) noexcept
{
  if (i >= n) { debug_trap(); }
}

template <typename T>
DBG_NOINLINE DBG_OPTNONE static void debug_dereference(T* ptr) noexcept
{
  if (ptr == nullptr) { debug_trap(); }
}


template <typename T>
concept Zeroable =
    std::is_integral_v<T> || std::is_floating_point_v<T> || std::is_same_v<T, bool> || std::is_pointer_v<T> || std::is_same_v<T, std::nullptr_t>;

template <typename T>
  requires Zeroable<T> && (!std::is_pointer_v<T>)
T zero_value()
{
  return T{};
}

template <typename T>
  requires Zeroable<T> && std::is_pointer_v<T>
T zero_value()
{
  return nullptr;
}

// template <typename A>
// using VecAlloc_t = CustomAllocator<A>;

// template <typename A>
// using Vec_t = std::vector<A, std::scoped_allocator_adaptor<VecAlloc_t<A>>>;

// template <typename A, typename B>
// using MapMember_t = std::pair<const A, B>;

// template <typename A, typename B>
// using MapAlloc_t = CustomAllocator<MapMember_t<A, B>>;

// template <typename A, typename B>
// using Map_t = std::unordered_map<A, B, std::hash<A>, std::equal_to<A>, std::scoped_allocator_adaptor<MapAlloc_t<A, B>>>;

}    //  namespace NUCC

#endif    // !__NUCC_DEFS_HPP

#ifndef NUCC_DEFS_HPP
#define NUCC_DEFS_HPP

#include <scoped_allocator>
#include <cstddef>
#include <unordered_map>
#include <vector>

#    define LMP_NUCC_ALLOC_COEFF 1.2
#    define LMP_NUCC_CLUSTER_MAX_OWNERS 128
#    define LMP_NUCC_CLUSTER_MAX_SIZE 300
#    define LMP_NUCC_CLUSTER_MAX_GHOST 300

// #define __NUCC_CSPAN_CHECK_ACCESS
// #define __NUCC_CSPAN_DEBUG_CALLS
// #define __NUCC_CHECK_ACCESS
// #define __NUCC_ALGO_CHECK
// #define __NUCC_NEIGHS_RADIAL_USE_HALF
#define __NUCC_NEIGHS_RADIAL_PRECOMPUTE_NORM

namespace NUCC {

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

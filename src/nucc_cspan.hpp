#ifndef __NUCC_CUSTOM_CSPAN_HPP
#define __NUCC_CUSTOM_CSPAN_HPP

#include "memory.h"

#include <algorithm>
#include <span>
#include <type_traits>

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

#ifdef __NUCC_CSPAN_DEBUG_CALLS
#include <iostream>
// dumb thing used to catch calls that are inlined by default
struct Kallbeck {
  __attribute_noinline__ void call() const noexcept
  {
    int a = 5;
    std::cout << "dsds" << a << std::endl;
  }
};
#endif

template <typename T, std::size_t Extent = std::dynamic_extent>
class cspan {
 public:
  cspan() noexcept                        = default;
  cspan(const std::span<T, Extent>& span) = delete;
  cspan(const cspan<T, Extent>& other)    = delete;
  constexpr cspan(T* ptr, std::size_t n) noexcept : span_(std::span<T, Extent>(ptr, n)) {}
  constexpr cspan(T* begin, T* end) noexcept : span_(std::span<T, Extent>(begin, end)) {}
  constexpr cspan(std::span<T, Extent>&& span) noexcept : span_(span) {}
  constexpr cspan(cspan<T, Extent>&& other) noexcept : span_(std::move(other.span_)) {}
  template <typename U, std::size_t OtherExtent>
    requires(std::is_convertible_v<U (*)[], T (*)[]>)
  constexpr cspan(const cspan<U, OtherExtent>& other) noexcept : span_(other.data(), other.size())
  {
  }
  inline constexpr cspan<T, Extent>& operator=(std::span<T, Extent>&& span) noexcept
  {
    span_ = span;
    return *this;
  }
  inline constexpr cspan<T, Extent>& operator=(cspan<T, Extent>&& other) noexcept
  {
    span_ = std::move(other.span_);
    return *this;
  }

  inline constexpr T& at(std::size_t index)
#ifndef __NUCC_CSPAN_CHECK_ACCESS
      noexcept
#endif
  {
#ifdef __NUCC_CSPAN_CHECK_ACCESS
    if (index >= span_.size()) {
#ifdef __NUCC_CSPAN_DEBUG_CALLS
      beck.call();
#endif
      throw std::out_of_range("Index out of range");
    }
#endif
    return span_[index];
  }

  inline constexpr T* offset(std::size_t offset)
#ifndef __NUCC_CSPAN_CHECK_ACCESS
      noexcept
#endif
  {
#ifdef __NUCC_CSPAN_CHECK_ACCESS
    if (offset >= span_.size()) {
#ifdef __NUCC_CSPAN_DEBUG_CALLS
      beck.call();
#endif
      throw std::out_of_range("Index out of range");
    }
#endif
    return span_.data() + offset;
  }

  inline constexpr T& operator[](std::size_t index)
#ifndef __NUCC_CSPAN_CHECK_ACCESS
      noexcept
#endif
  {
#ifdef __NUCC_CSPAN_CHECK_ACCESS
    if (index >= span_.size()) {
#ifdef __NUCC_CSPAN_DEBUG_CALLS
      beck.call();
#endif
      throw std::out_of_range("Index out of range");
    }
#endif
    return span_[index];
  }

  inline constexpr const T& at(std::size_t index) const
#ifndef __NUCC_CSPAN_CHECK_ACCESS
      noexcept
#endif
  {
#ifdef __NUCC_CSPAN_CHECK_ACCESS
    if (index >= span_.size()) {
#ifdef __NUCC_CSPAN_DEBUG_CALLS
      beck.call();
#endif
      throw std::out_of_range("Index out of range");
    }
#endif
    return span_[index];
  }

  inline constexpr const T* offset(std::size_t offset) const
#ifndef __NUCC_CSPAN_CHECK_ACCESS
      noexcept
#endif
  {
#ifdef __NUCC_CSPAN_CHECK_ACCESS
    if (offset >= span_.size()) {
#ifdef __NUCC_CSPAN_DEBUG_CALLS
      beck.call();
#endif
      throw std::out_of_range("Index out of range");
    }
#endif
    return span_.data() + offset;
  }

  inline constexpr const T& operator[](std::size_t index) const
#ifndef __NUCC_CSPAN_CHECK_ACCESS
      noexcept
#endif
  {
#ifdef __NUCC_CSPAN_CHECK_ACCESS
    if (index >= span_.size()) {
#ifdef __NUCC_CSPAN_DEBUG_CALLS
      beck.call();
#endif
      throw std::out_of_range("Index out of range");
    }
#endif
    return span_[index];
  }

  inline constexpr T* data() noexcept { return span_.data(); }

  inline constexpr const T* data() const noexcept { return span_.data(); }

  inline constexpr std::size_t size() const noexcept { return span_.size(); }

  void destroy(LAMMPS_NS::Memory* memory) noexcept
  {
    T* ptr = span_.data();
    memory->destroy(ptr);
    span_ = std::span<T, Extent>();
  }

  void create(LAMMPS_NS::Memory* memory, std::size_t n, const char* name)
  {
    T* ptr;
    memory->create(ptr, n, name);
    span_ = std::span<T, Extent>(ptr, n);
  }

  void grow(LAMMPS_NS::Memory* memory, std::size_t n, const char* name)
  {
    T* ptr = span_.data();
    memory->grow(ptr, n, name);
    span_ = std::span<T, Extent>(ptr, n);
  }

  inline constexpr bool empty() const noexcept { return span_.empty(); }

  inline void reset()
    requires(!std::is_const_v<T>)
  {
    if (!span_.empty()) { std::fill_n(span_.data(), span_.size(), zero_value<T>()); }
  }

  template <typename U>
  inline void reset_unsafe(const U& value)
    requires(!std::is_const_v<T>)
  {
    static_assert(sizeof(T) % sizeof(U) == 0, "Object size is not multiple of target value size");
    if (!span_.empty()) { std::fill(reinterpret_cast<U*>(span_.data()), reinterpret_cast<U*>(span_.data() + span_.size()), value); }
  }

  template <typename U>
  inline void reset_unsafe(U&& value)
    requires(!std::is_const_v<T>)
  {
    static_assert(sizeof(T) % sizeof(U) == 0, "Object size is not multiple of target value size");
    if (!span_.empty()) {
      std::fill(reinterpret_cast<std::remove_const_t<U>*>(span_.data()), reinterpret_cast<std::remove_const_t<U>*>(span_.data() + span_.size()),
                value);
    }
  }

  inline constexpr std::size_t memory_usage() const noexcept { return size() * sizeof(T) + sizeof(std::size_t) + sizeof(T*); }

 private:
  std::span<T, Extent> span_;
#ifdef __NUCC_CSPAN_DEBUG_CALLS
  Kallbeck beck;
#endif
};

}    // namespace NUCC

#endif    // CUSTOM_CSPAN_HPP

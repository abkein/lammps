#ifndef NUCC_CUSTOM_CSPAN_HPP
#define NUCC_CUSTOM_CSPAN_HPP

#include "memory.h"
#include "nucc_defs.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <span>
#include <type_traits>

namespace NUCC {

template <typename T, std::size_t Extent = std::dynamic_extent>
class cspan {
 public:
  cspan() noexcept                        = default;
  cspan(const std::span<T, Extent>& span) = delete;
  cspan(const cspan<T, Extent>& other)    = delete;
  constexpr cspan(T* ptr, std::size_t n) noexcept : span_(std::span<T, Extent>(ptr, n)) {}
  constexpr cspan(T* begin, T* end) noexcept : span_(std::span<T, Extent>(begin, end)) {}
  constexpr cspan(std::span<T, Extent>&& span) noexcept : span_(span) {}    // NOLINT(hicpp-explicit-conversions)
  constexpr cspan(cspan<T, Extent>&& other) noexcept : span_(std::move(other.span_)) {}
  template <typename U, std::size_t OtherExtent>
    requires(std::is_convertible_v<U (*)[], T (*)[]>)    // NOLINT(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays,hicpp-avoid-c-arrays)
  constexpr cspan(const cspan<U, OtherExtent>& other) noexcept : span_(other.data(), other.size())    // NOLINT(hicpp-explicit-conversions)
  {
  }

  constexpr cspan<T, Extent>& operator=(std::span<T, Extent>&& span) noexcept
  {
    span_ = span;
    return *this;
  }
  constexpr cspan<T, Extent>& operator=(cspan<T, Extent>&& other) noexcept
  {
    span_ = std::move(other.span_);
    return *this;
  }

  [[nodiscard]] constexpr T& at(std::size_t index) noexcept
  {
    if constexpr (Defines::CSPAN_DEBUG_CALLS) { debug_check_index(index, span_.size()); }
    return span_[index];
  }

  [[nodiscard]] constexpr T* offset(std::size_t offset) noexcept
  {
    if constexpr (Defines::CSPAN_DEBUG_CALLS) { debug_check_index(offset, span_.size()); }
    return span_.data() + offset;
  }

  [[nodiscard]] constexpr T& operator[](std::size_t index) noexcept
  {
    if constexpr (Defines::CSPAN_DEBUG_CALLS) { debug_check_index(index, span_.size()); }
    return span_[index];
  }

  [[nodiscard]] constexpr const T& at(std::size_t index) const noexcept
  {
    if constexpr (Defines::CSPAN_DEBUG_CALLS) { debug_check_index(index, span_.size()); }
    return span_[index];
  }

  [[nodiscard]] constexpr const T* offset(std::size_t offset) const noexcept
  {
    if constexpr (Defines::CSPAN_DEBUG_CALLS) { debug_check_index(offset, span_.size()); }
    return span_.data() + offset;
  }

  [[nodiscard]] constexpr const T* begin() const noexcept { return span_.data(); }

  [[nodiscard]] constexpr const T* end() const noexcept { return span_.data() + span_.size(); }

  [[nodiscard]] constexpr T* begin() noexcept { return span_.data(); }

  [[nodiscard]] constexpr T* end() noexcept { return span_.data() + span_.size(); }

  [[nodiscard]] constexpr const T& operator[](std::size_t index) const noexcept
  {
    if constexpr (Defines::CSPAN_DEBUG_CALLS) { debug_check_index(index, span_.size()); }
    return span_[index];
  }

  [[nodiscard]] constexpr T* data() noexcept { return span_.data(); }

  [[nodiscard]] constexpr const T* data() const noexcept { return span_.data(); }

  [[nodiscard]] constexpr std::size_t size() const noexcept { return span_.size(); }

  void destroy(LAMMPS_NS::Memory* memory) noexcept
  {
    T* ptr = span_.data();
    memory->destroy(ptr);
    span_ = std::span<T, Extent>();
  }

  // void create(LAMMPS_NS::Memory* memory, std::size_t n, const char* name)
  // {
  //   T* ptr;
  //   memory->create(ptr, n, name);
  //   span_ = std::span<T, Extent>(ptr, n);
  // }

  void grow(LAMMPS_NS::Memory* memory, std::size_t n, const char* name)
  {
    T* ptr = span_.data();
    memory->grow(ptr, n, name);
    span_ = std::span<T, Extent>(ptr, n);
  }

  [[nodiscard]] constexpr bool empty() const noexcept { return span_.empty(); }

  void reset()
    requires(!std::is_const_v<T>)
  {
    if (!span_.empty()) { std::fill_n(span_.data(), span_.size(), zero_value<T>()); }
  }

  template <typename U>
  void reset_unsafe(const U& value)
    requires(!std::is_const_v<T>)
  {
    static_assert(sizeof(T) % sizeof(U) == 0, "Object size is not multiple of target value size");
    if (!span_.empty()) { std::fill(reinterpret_cast<U*>(span_.data()), reinterpret_cast<U*>(span_.data() + span_.size()), value); }
  }

  template <typename U>
  void reset_unsafe(U&& value)
    requires(!std::is_const_v<T>)
  {
    static_assert(sizeof(T) % sizeof(U) == 0, "Object size is not multiple of target value size");
    if (!span_.empty()) {
      std::fill(reinterpret_cast<std::remove_const_t<U>*>(span_.data()), reinterpret_cast<std::remove_const_t<U>*>(span_.data() + span_.size()),
                value);
    }
  }

  [[nodiscard]] constexpr std::size_t memory_usage() const noexcept { return size() * sizeof(T) + sizeof(std::size_t) + sizeof(T*); }

 private:
  std::span<T, Extent> span_;
};

}    // namespace NUCC

#endif    // CUSTOM_CSPAN_HPP

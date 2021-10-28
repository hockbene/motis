#pragma once

#include <array>

#include "motis/core/schedule/time.h"

namespace motis {
namespace csa {

constexpr duration MAX_TRAVEL_TIME = 1440;
constexpr duration MAX_TRANSFERS = 7;

struct journey_pointer {
  journey_pointer() = default;
  journey_pointer(csa_connection const* enter_con,
                  csa_connection const* exit_con, footpath const* footpath)
      : enter_con_(enter_con), exit_con_(exit_con), footpath_(footpath) {}

  bool valid() const {
    return enter_con_ != nullptr && exit_con_ != nullptr &&
           footpath_ != nullptr;
  }

  csa_connection const* enter_con_{nullptr};
  csa_connection const* exit_con_{nullptr};
  footpath const* footpath_{nullptr};
};

// https://stackoverflow.com/questions/49318316/initialize-all-elements-or-stdarray-with-the-same-constructor-arguments
template <typename T, std::size_t N, std::size_t Idx = N>
struct array_maker {
  template <typename... Ts>
  static std::array<T, N> make_array(const T& v, Ts... tail) {
    return array_maker<T, N, Idx - 1>::make_array(v, v, tail...);
  }
};

template <typename T, std::size_t N>
struct array_maker<T, N, 1> {
  template <typename... Ts>
  static std::array<T, N> make_array(const T& v, Ts... tail) {
    return std::array<T, N>{v, tail...};
  }
};

}  // namespace csa
}  // namespace motis

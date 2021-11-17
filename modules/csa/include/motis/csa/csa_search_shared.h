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

void add_journey_leg(search_dir Dir, const csa_timetable& tt,
                     csa_journey& journey, journey_pointer& jp) {
  if (jp.footpath_->from_station_ != jp.footpath_->to_station_) {
    if (Dir == search_dir::FWD) {
      journey.edges_.emplace_back(
          &tt.stations_[jp.footpath_->from_station_],
          &tt.stations_[jp.footpath_->to_station_],
          jp.enter_con_->departure_ - jp.footpath_->duration_,
          jp.enter_con_->departure_, -1);
    } else {
      journey.edges_.emplace_back(
          &tt.stations_[jp.footpath_->from_station_],
          &tt.stations_[jp.footpath_->to_station_], jp.exit_con_->arrival_,
          jp.exit_con_->arrival_ + jp.footpath_->duration_, -1);
    }
  }

  assert(jp.enter_con_->trip_ == jp.exit_con_->trip_);
  auto const& trip_connections = tt.trip_to_connections_[jp.exit_con_->trip_];
  auto const add_trip_edge = [&](csa_connection const* con) {
    auto const enter = con == jp.enter_con_;
    auto const exit = con == jp.exit_con_;
    utl::verify(con->light_con_ != nullptr, "invalid light connection");
    journey.edges_.emplace_back(con->light_con_,
                                &tt.stations_[con->from_station_],
                                &tt.stations_[con->to_station_], enter, exit,
                                con->departure_, con->arrival_);
  };
  auto in_trip = false;
  if (Dir == search_dir::FWD) {
    for (int i = 0; i < static_cast<int>(trip_connections.size()); ++i) {
      auto const con = trip_connections[i];
      if (con == jp.enter_con_) {
        in_trip = true;
      }
      if (in_trip) {
        add_trip_edge(con);
      }
      if (con == jp.exit_con_) {
        break;
      }
    }
  } else {
    for (int i = static_cast<int>(trip_connections.size()) - 1; i >= 0; --i) {
      auto const con = trip_connections[i];
      if (con == jp.exit_con_) {
        in_trip = true;
      }
      if (in_trip) {
        add_trip_edge(con);
      }
      if (con == jp.enter_con_) {
        break;
      }
    }
  }
}

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

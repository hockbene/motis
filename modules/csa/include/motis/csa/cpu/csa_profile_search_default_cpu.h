#pragma once

#include <algorithm>
#include <array>
#include <limits>
#include <set>
#include <unordered_set>
#include <vector>

#include "utl/verify.h"

#include "motis/core/schedule/interval.h"

#include "motis/csa/csa_journey.h"
#include "motis/csa/csa_profile_reconstruction.h"
#include "motis/csa/csa_search_shared.h"
#include "motis/csa/csa_statistics.h"
#include "motis/csa/csa_timetable.h"
#include "motis/csa/error.h"

namespace motis::csa::cpu {

template <search_dir Dir>
struct csa_profile_search {
  static constexpr auto INVALID = Dir == search_dir::FWD
                                      ? std::numeric_limits<time>::max()
                                      : std::numeric_limits<time>::min();

  using arrival_times = std::array<time, MAX_TRANSFERS + 1>;

  csa_profile_search(csa_timetable const& tt, interval const& search_interval,
                     csa_statistics& stats)
      : tt_{tt},
        search_interval_{search_interval},
        arrival_time_(
            tt.stations_.size(),
            {std::make_pair(
                INVALID,
                array_maker<time, MAX_TRANSFERS + 1>::make_array(INVALID))}),
        trip_reachable_(
            tt.trip_count_,
            array_maker<time, MAX_TRANSFERS + 1>::make_array(INVALID)),
        final_footpaths_(tt.stations_.size(), INVALID),
        is_trip_reachable_(tt.trip_count_, true),
        targets_(),
        starts_(),
        stats_{stats} {}

  void add_destination(csa_station const& station, time /*initial_duration*/) {
    // For now initial_duration is always 0
    targets_.insert(station.id_);
    // Initialize final_footpaths_
    if (Dir == search_dir::FWD) {
      for (auto const& fp : station.incoming_footpaths_) {
        final_footpaths_[fp.from_station_] =
            std::min(fp.duration_, final_footpaths_[fp.from_station_]);
      }
    } else {
      for (auto const& fp : station.footpaths_) {
        final_footpaths_[fp.to_station_] =
            std::min(fp.duration_, final_footpaths_[fp.to_station_]);
      }
    }

    stats_.destination_count_++;
  }

  bool dominates_profile(std::pair<time, arrival_times> const& a,
                         std::pair<time, arrival_times> const& b) {
    bool dominates = false;
    auto const dominates_pair = [&](time x1, time x2, time y1, time y2) {
      return Dir == search_dir::FWD
                 ? (x1 <= y1 && x2 <= y2) && (x1 < y1 || x2 < y2)
                 : (x1 >= y1 && x2 >= y2) && (x1 > y1 || x2 > y2);
    };
    for (auto i = 0; i < a.second.size(); ++i) {
      time a1 = a.first, a2 = a.second[i];
      time b1 = b.first, b2 = b.second[i];
      if (dominates_pair(b1, b2, a1, a2)) {
        return false;
      }
      dominates = dominates || dominates_pair(a1, a2, b1, b2);
    }
    return dominates;
  }

  arrival_times shift(arrival_times const& array) {
    auto result = array_maker<time, MAX_TRANSFERS + 1>::make_array(INVALID);
    for (auto i = 1; i < result.size(); ++i) {
      result[i] = array[i - 1];
    }
    if (array.size() >= 2) {
      result[result.size() - 1] =
          std::min(array[array.size() - 2], array[array.size() - 1]);
    }
    return result;
  }

  arrival_times cw_min(arrival_times a, arrival_times b) {
    auto c = array_maker<time, MAX_TRANSFERS + 1>::make_array(INVALID);
    if (Dir == search_dir::FWD) {
      std::transform(a.begin(), a.end(), b.begin(), c.begin(),
                     [](time x, time y) { return std::min(x, y); });
    } else {
      std::transform(a.begin(), a.end(), b.begin(), c.begin(),
                     [](time x, time y) { return std::max(x, y); });
    }
    return c;
  }

  arrival_times cw_min(arrival_times a, arrival_times b, arrival_times c) {
    return cw_min(cw_min(a, b), c);
  }

  arrival_times get_time_walking(const csa_connection& con) {
    auto x = Dir == search_dir::FWD ? final_footpaths_[con.to_station_]
                                    : final_footpaths_[con.from_station_];
    if (x != INVALID) {
      if (Dir == search_dir::FWD) {
        x = con.arrival_ + x;
      } else {
        x = con.departure_ - x;
      }
    }
    return array_maker<time, MAX_TRANSFERS + 1>::make_array(x);
  }

  arrival_times get_time_trip(const csa_connection& con) {
    return trip_reachable_[con.trip_];
  }

  arrival_times get_time_transfer(const csa_connection& con) {
    auto const station =
        Dir == search_dir::FWD ? con.to_station_ : con.from_station_;
    auto const timestamp =
        Dir == search_dir::FWD ? con.arrival_ : con.departure_;
    auto p_it = std::lower_bound(
        arrival_time_[station].begin(), arrival_time_[station].end(), timestamp,
        [](auto const& pair, time t) {
          return Dir == search_dir::FWD ? pair.first < t : pair.first > t;
        });
    return shift(p_it->second);
  }

  void search() {
    auto const& connections =
        Dir == search_dir::FWD ? tt_.fwd_connections_ : tt_.bwd_connections_;

    if (DEBUG_OUTPUT) {
      LOG(motis::logging::debug)
          << "[*] search interval (" << search_interval_.begin_ << ", "
          << search_interval_.end_ << ")";
      for (auto t : targets_) {
        LOG(motis::logging::debug) << "[*] target " << t;
      }
      for (auto my_station : tt_.stations_) {
        LOG(motis::logging::debug)
            << "[S]" << my_station.id_ << "_{"
            << "\"id\":" << my_station.id_ << ", "
            << "\"name\":"
            << "\"" << my_station.station_ptr_->name_ << "\", "
            << "\"lat\":" << my_station.station_ptr_->lat() << ", "
            << "\"long\":" << my_station.station_ptr_->lng() << ", "
            << "\"transfer\":" << my_station.transfer_time_ << ", "
            << "\"incoming_connections\":"
            << my_station.incoming_connections_.size() << ", "
            << "\"outgoing_connections\":"
            << my_station.outgoing_connections_.size() << ", "
            << "}";
        for (auto fp : my_station.footpaths_) {
          LOG(motis::logging::debug) << "[F]{"
                                     << "\"from\":" << fp.from_station_ << ", "
                                     << "\"to\":" << fp.to_station_ << ", "
                                     << "\"duration\":" << fp.duration_ << "}";
        }
      }
      for (auto con : connections) {
        LOG(motis::logging::debug)
            << "[C]{"
            << "\"dep_time\":" << con.departure_ << ","
            << "\"dep_stop\":" << con.from_station_ << ","
            << "\"arr_time\":" << con.arrival_ << ","
            << "\"arr_stop\":" << con.to_station_ << ","
            << "\"trip\":" << con.trip_ << "}";
      }
    }

    csa_connection const earliest_con{search_interval_.begin_};
    auto const last_connection = connections.rbegin();
    auto const first_connection = std::upper_bound(
        connections.rbegin(), connections.rend(), earliest_con,
        [&](auto const& a, auto const& b) {
          return Dir == search_dir::FWD ? a.departure_ > b.departure_
                                        : a.arrival_ < b.arrival_;
        });

    for (auto it = last_connection; it != first_connection; ++it) {
      auto const& con = *it;

      if (!is_trip_reachable_[con.trip_]) {
        continue;
      }

      auto const time_walking = get_time_walking(con);
      auto const time_trip = get_time_trip(con);
      auto const time_transfer = get_time_transfer(con);
      auto const min_arrival_time =
          cw_min(time_walking, time_trip, time_transfer);

      auto const curr_best_arrival_time =
          std::lower_bound(
              arrival_time_[Dir == search_dir::FWD ? con.from_station_
                                                   : con.to_station_]
                  .begin(),
              arrival_time_[Dir == search_dir::FWD ? con.from_station_
                                                   : con.to_station_]
                  .end(),
              Dir == search_dir::FWD ? con.departure_ : con.arrival_,
              [](auto const& pair, time t) {
                return Dir == search_dir::FWD ? pair.first < t : pair.first > t;
              })
              ->second;

      if (curr_best_arrival_time != min_arrival_time &&
          !is_dominated_in(std::make_pair(con.departure_, min_arrival_time),
                           arrival_time_[con.from_station_])) {
        auto const arrival_profile =
            cw_min(curr_best_arrival_time, min_arrival_time);
        auto const& footpaths =
            Dir == search_dir::FWD
                ? tt_.stations_[con.from_station_].incoming_footpaths_
                : tt_.stations_[con.from_station_].footpaths_;
        for (auto const& fp : footpaths) {
          auto const departure_time = Dir == search_dir::FWD
                                          ? con.departure_ - fp.duration_
                                          : con.arrival_ + fp.duration_;
          auto const footpath_start =
              Dir == search_dir::FWD ? fp.from_station_ : fp.to_station_;
          auto const profile_pair = make_pair(departure_time, arrival_profile);
          add_to_profile(profile_pair, footpath_start);
        }
      }

      trip_reachable_[con.trip_] = min_arrival_time;
      if (DEBUG_OUTPUT) {
        std::string min_arrival_time_string;
        for (auto t : min_arrival_time) {
          min_arrival_time_string += std::to_string(t);
          min_arrival_time_string += ", ";
        }
        LOG(motis::logging::debug)
            << "[T]" << con.trip_ << "_(" << min_arrival_time_string << ")";
      }
      stats_.connections_scanned_++;
    }
  }

  bool is_source_dominated(const std::pair<time, arrival_times>& profile_pair) {
    for (auto start_idx : starts_) {
      if (dominates_profile(arrival_time_[start_idx].front(), profile_pair)) {
        return true;
      }
    }
    return false;
  }

  bool is_dominated_in(
      const std::pair<time, arrival_times>& new_pair,
      const std::list<std::pair<time, arrival_times>>& profile) {
    return std::any_of(profile.begin(), profile.end(),
                       [&](auto const& profile_pair) {
                         return dominates_profile(profile_pair, new_pair);
                       });
  }

  void add_to_profile(const std::pair<time, arrival_times>& profile_pair,
                      const uint32_t station_id) {
    /*
    if (starts_.find(station_id) == starts_.end()) {
      if (is_source_dominated(profile_pair)) {
        return;
      }
    }*/

    auto& arrival_time = arrival_time_[station_id];
    auto const insert_position =
        std::lower_bound(arrival_time.begin(), arrival_time.end(), profile_pair,
                         [](auto const& p1, auto const& p2) {
                           return Dir == search_dir::FWD ? p1.first < p2.first
                                                         : p1.first > p2.first;
                         });

    auto const inserted_at = arrival_time.insert(insert_position, profile_pair);

    if (DEBUG_OUTPUT) {
      std::string profile;
      for (auto p : profile_pair.second) {
        profile += std::to_string(p);
        profile += ", ";
      }
      LOG(motis::logging::debug)
          << "[P](" << profile_pair.first << ", (" << profile << "))_"
          << station_id << "_" << tt_.stations_[station_id].station_ptr_->name_;
    }
    for (auto it = std::make_reverse_iterator(inserted_at);
         it != arrival_time.rend();) {
      if (dominates_profile(profile_pair, *it)) {
        it = std::make_reverse_iterator(
            arrival_time.erase(std::next(it).base()));
      } else {
        ++it;
      }
    }
  }

  std::vector<csa_journey> get_results(csa_station const& station,
                                       bool include_equivalent) {
    utl::verify_ex(!include_equivalent,
                   std::system_error{error::include_equivalent_not_supported});

    std::vector<csa_journey> journeys;
    auto journey_reconstruction = csa_profile_reconstruction<
        Dir, decltype(arrival_time_), decltype(final_footpaths_),
        decltype(trip_reachable_), decltype(targets_)>{
        tt_, arrival_time_, final_footpaths_, trip_reachable_, targets_};
    auto const& station_profile = arrival_time_[station.id_];
    for (auto it = station_profile.begin(); it != station_profile.end(); ++it) {
      auto const profile_pair = *it;
      auto const departure_time = profile_pair.first;
      if (departure_time > search_interval_.end_) {
        break;
      }
      for (auto transfers = 0; transfers <= MAX_TRANSFERS; ++transfers) {
        auto const arrival_time = profile_pair.second[transfers];
        auto const next_profile_arrival_time = next(it)->second[transfers];
        if (arrival_time != INVALID &&
            arrival_time != next_profile_arrival_time) {
          auto& journey = journeys.emplace_back(
              Dir, departure_time, arrival_time, transfers, nullptr);
          journey.start_station_ = &station;
          journey_reconstruction.extract_journey(journey);
        }
      }
    }
    return journeys;
  }

  csa_timetable const& tt_;
  interval search_interval_;
  std::vector<std::list<std::pair<time, arrival_times>>> arrival_time_;
  std::vector<arrival_times> trip_reachable_;
  std::vector<time> final_footpaths_;
  std::vector<bool> is_trip_reachable_;
  std::set<station_id> targets_;
  std::unordered_set<station_id> starts_;
  csa_statistics& stats_;
  bool DEBUG_OUTPUT = false;
};

}  // namespace motis::csa::cpu

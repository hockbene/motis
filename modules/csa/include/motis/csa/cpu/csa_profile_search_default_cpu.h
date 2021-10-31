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
    assert(a.second.size() == b.second.size());
    bool dominates = false;
    auto const dominates_triple = [&](time x1, time x2, time y1, time y2) {
      return Dir == search_dir::FWD
                 ? (x1 <= y1 && x2 <= y2) && (x1 < y1 || x2 < y2)
                 : (x1 >= y1 && x2 >= y2) && (x1 > y1 || x2 > y2);
    };
    for (auto i = 0; i < a.second.size(); ++i) {
      time a1 = a.first, a2 = a.second[i];
      time b1 = b.first, b2 = b.second[i];
      if (dominates_triple(b1, b2, a1, a2)) {
        return false;
      }
      dominates = dominates || dominates_triple(a1, a2, b1, b2);
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
    assert(a.size() == b.size());
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

  void search_bwd() {
    auto const& connections = tt_.bwd_connections_;
    csa_connection const start_at{search_interval_.end_};
    auto const first_connection =
        std::lower_bound(connections.begin(), connections.end(), start_at,
                         [&](csa_connection const& a, csa_connection const& b) {
                           return a.arrival_ > b.arrival_;
                         });
    if (first_connection == connections.begin()) {
      return;
    }

    for (auto it = first_connection; it-- != connections.begin();) {
      auto const& con = *it;

      if (!is_trip_reachable_[con.trip_]) {
        continue;
      }

      // tau_1 - erreichen per fußpfad
      auto x = final_footpaths_[con.from_station_];
      if (x != INVALID) {
        x = con.departure_ - x;
      }
      auto const tau_1 = array_maker<time, MAX_TRANSFERS + 1>::make_array(x);

      // tau_2 - erreichen per sitzenbleiben
      auto const tau_2 = trip_reachable_[con.trip_];

      // tau_3 - erreichen per umstieg
      /*
      auto p_it = arrival_time_[con.from_station].end();
      while ((*p_it).first > con.departure_) {
        p_it++;
      }
      */
      auto p_it = std::lower_bound(
          arrival_time_[con.from_station_].begin(),
          arrival_time_[con.from_station_].end(), con.departure_,
          [&](auto const pair, time t) { return pair.first > t; });

      auto const tau_3 = shift((*p_it).second);

      auto const tau_c = cw_min(cw_min(tau_1, tau_2), tau_3);
      auto const y_p = std::find_if(
          arrival_time_[con.to_station_].begin(),
          arrival_time_[con.to_station_].end(),
          [&](auto const pair) { return pair.first < con.arrival_; });
      auto const y = (*y_p).second;
      assert(y.size() == tau_c.size());
      auto const arrival_profile = cw_min(y, tau_c);

      if (y != tau_c && !is_dominated_in(std::make_pair(con.arrival_, tau_c),
                                         arrival_time_[con.to_station_])) {
        for (auto fp : tt_.stations_[con.to_station_].footpaths_) {
          auto const departure_time = con.arrival_ + fp.duration_;
          auto const profile_pair =
              std::make_pair(departure_time, arrival_profile);
          add_to_profile(profile_pair, fp.to_station_);
        }
      }
      trip_reachable_[con.trip_] = tau_c;
      stats_.connections_scanned_++;
    }
  }

  arrival_times get_tau_1(const csa_connection& con) {
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

  arrival_times get_tau_2(const csa_connection& con) {
    return trip_reachable_[con.trip_];
  }

  arrival_times get_tau_3(const csa_connection& con) {
    /*
    auto p_it = arrival_time_[con.to_station_].begin();
    while (p_it->first < con.arrival_) {
      p_it++;
    }*/
    auto const station =
        Dir == search_dir::FWD ? con.to_station_ : con.from_station_;
    auto const timestamp =
        Dir == search_dir::FWD ? con.arrival_ : con.departure_;
    auto p_it = std::lower_bound(
        arrival_time_[station].begin(), arrival_time_[station].end(), timestamp,
        [](auto const& pair, time t) {
          // TODO(root) may be <= / >= instead of < / >
          return Dir == search_dir::FWD ? pair.first < t : pair.first > t;
        });
    return shift(p_it->second);

    /*
     * auto p_it = std::lower_bound(
arrival_time_[con.from_station_].begin(),
arrival_time_[con.from_station_].end(), con.departure_,
[&](auto const pair, time t) { return pair.first > t; });

auto const tau_3 = shift((*p_it).second);
     */
  }

  void search() {
    if (Dir == search_dir::BWD) {
      search_bwd();
      return;
    }

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
    }
    if (DEBUG_OUTPUT) {
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
    int changes_in_s = 0;

    auto first_connection = connections.end();
    auto last_connection = connections.begin();

    // TODO(root) BWD
    // if (Dir == search_dir::FWD) {
    csa_connection const end_at{search_interval_.begin_};
    last_connection =
        std::lower_bound(begin(connections), end(connections), end_at,
                         [&](csa_connection const& a, csa_connection const& b) {
                           return a.departure_ < b.departure_;
                         });
    if (last_connection == connections.end()) {
      return;
    }
    // TODO(root) BWD
    /*} else {
      csa_connection const start_at{search_interval_.end_};
      first_connection = std::lower_bound(
          begin(connections), end(connections), start_at,
          [&](csa_connection const& a, csa_connection const& b) {
            return a.arrival_ > b.arrival_;
          });
      if (first_connection == connections.begin()) {
        return;
      }
    }
    */
    for (auto it = first_connection; it-- != last_connection;) {
      auto const& con = *it;

      if (!is_trip_reachable_[con.trip_]) {
        continue;
      }

      // tau's
      auto const tau_1 = get_tau_1(con);
      auto const tau_2 = get_tau_2(con);
      auto const tau_3 = get_tau_3(con);

      auto const tau_c = cw_min(cw_min(tau_1, tau_2), tau_3);

      // TODO(root) BWD
      auto const y_p = std::find_if(
          arrival_time_[con.from_station_].begin(),
          arrival_time_[con.from_station_].end(),
          [&](auto const pair) { return pair.first >= con.departure_; });
      auto const y = y_p->second;

      assert(y.size() == tau_c.size());

      auto const arrival_profile = cw_min(y, tau_c);

      // TODO(root) Limited Walking Optimization
      if (y != tau_c && !is_dominated_in(std::make_pair(con.departure_, tau_c),
                                         arrival_time_[con.from_station_])) {
        auto const& footpaths =
            Dir == search_dir::FWD
                ? tt_.stations_[con.from_station_].incoming_footpaths_
                : tt_.stations_[con.from_station_].footpaths_;
        for (auto const& fp : footpaths) {
          // TODO(root) BWD
          auto const departure_time = Dir == search_dir::FWD
                                          ? con.departure_ - fp.duration_
                                          : con.arrival_ + fp.duration_;
          auto const footpath_start =
              Dir == search_dir::FWD ? fp.from_station_ : fp.to_station_;
          auto const profile_pair = make_pair(departure_time, arrival_profile);
          add_to_profile(profile_pair, footpath_start);
        }
      }

      trip_reachable_[con.trip_] = tau_c;
      if (DEBUG_OUTPUT) {

        std::string tau_c_string;
        for (auto t : tau_c) {
          tau_c_string += std::to_string(t);
          tau_c_string += ", ";
        }
        LOG(motis::logging::debug)
            << "[T]" << con.trip_ << "_(" << tau_c_string << ")";
      }
      stats_.connections_scanned_++;
    }
    if (DEBUG_OUTPUT) {
      LOG(motis::logging::debug) << changes_in_s;
    }
    /*
    if (DEBUG_OUTPUT) {
      int total = 0;
      int d_fail = 0;
      int a_fail = 0;
      for (auto station_idx = 0; station_idx < arrival_time_.size();
    ++station_idx) { for (int i = 0; i <
    static_cast<int>(arrival_time_[station_idx].size()) - 1; ++i) { total++;
          auto d1 = arrival_time_[station_idx][i].first;
          auto a1 = arrival_time_[station_idx][i].second;
          auto d2 = arrival_time_[station_idx][i + 1].first;
          auto a2 = arrival_time_[station_idx][i + 1].second;
          if (!(d1 <= d2)) {
            d_fail++;
            LOG(motis::logging::debug)
                << "[*] d_fail at station " << station_idx;
          }
          for (int j = 0; j < MAX_TRANSFERS + 1; ++j) {
            if (!(a1[j] <= a2[j])) {
              a_fail++;
            }
          }
        }
      }
      LOG(motis::logging::debug) << "[*]  total " << total;
      LOG(motis::logging::debug) << "[*] d_fail " << d_fail;
      LOG(motis::logging::debug) << "[*] a_fail " << a_fail;
    }
    */
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

  // TODO(root) BWD
  void add_to_profile(const std::pair<time, arrival_times>& profile_pair,
                      const uint32_t station_id) {

    if (starts_.find(station_id) == starts_.end()) {
      if (is_source_dominated(profile_pair)) {
        return;
      }
    }

    auto it = arrival_time_[station_id].begin();
    if (Dir == search_dir::FWD) {
      while (it->first < profile_pair.first) {
        it++;
      }
    } else {
      while (it->first > profile_pair.first) {
        it++;
      }
    }

    if (it->first == profile_pair.first && it->second == profile_pair.second) {
      return;
    }

    if (std::any_of(it, arrival_time_[station_id].end(), [&](auto other_pair) {
          return !dominates_profile(profile_pair, other_pair);
        })) {
      return;
    }

    // it = arrival_time_[station_id].insert(arrival_time_[station_id].begin(),
    // profile_pair);
    it = arrival_time_[station_id].emplace(arrival_time_[station_id].begin(),
                                           profile_pair);

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

    while (it != arrival_time_[station_id].begin()) {
      if (dominates_profile((*it), profile_pair)) {
        it = arrival_time_[station_id].erase(it);
      } else {
        --it;
      }
    }

    /*
    auto entries_before_departure_time = std::vector<
        std::pair<time, arrival_times>>();
    while (arrival_time_[fp.from_station_][0].first < departure_time) {
      entries_before_departure_time.insert(
          entries_before_departure_time.begin(),
    arrival_time_[fp.from_station_][0]);
      arrival_time_[fp.from_station_].erase(arrival_time_[fp.from_station_].begin());
    }

    arrival_time_[fp.from_station_].insert(arrival_time_[fp.from_station_].begin(),
                                {departure_time, arrival_profile});

    while (!entries_before_departure_time.empty()) {
            auto const pair = entries_before_departure_time[0];
            entries_before_departure_time.erase(
                entries_before_departure_time.begin());
            //  pair.second = cw_min(pair.second, tau_c);
            arrival_time_[fp.from_station_].insert(arrival_time_[fp.from_station_].begin(),
    pair);
          }
          */
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

  csa_timetable const& tt_;  // Fahrplan
  interval search_interval_;  // Suchintervall
  std::vector<std::list<std::pair<time, arrival_times>>>
      arrival_time_;  // S - Profile
  std::vector<arrival_times> trip_reachable_;  // T - Trips
  std::vector<time> final_footpaths_;  // D - Fußpfade zum Ziel
  std::vector<bool> is_trip_reachable_;
  std::set<station_id> targets_;  // Zielstationen
  std::unordered_set<station_id> starts_;  // Startstationen
  csa_statistics& stats_;  // Stats
  bool DEBUG_OUTPUT = false;
};

}  // namespace motis::csa::cpu

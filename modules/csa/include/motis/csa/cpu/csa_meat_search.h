#pragma once

#include <limits>
#include <unordered_set>
#include <vector>

#include "utl/verify.h"

#include "motis/core/schedule/interval.h"

#include "motis/csa/csa_journey.h"
#include "motis/csa/csa_search_shared.h"
#include "motis/csa/csa_statistics.h"
#include "motis/csa/csa_timetable.h"
#include "motis/csa/error.h"

namespace motis::csa::cpu {

template <search_dir Dir>
struct csa_meat_search {
  static constexpr auto INVALID = Dir == search_dir::FWD
                                      ? std::numeric_limits<time>::max()
                                      : std::numeric_limits<time>::min();

  csa_meat_search(csa_timetable const& tt, interval const& search_interval,
                  csa_statistics& stats, time max_delay)
      : tt_{tt},
        search_interval_{search_interval},
        arrival_time_(tt.stations_.size(), {{INVALID, INVALID}}),
        trip_arrival_time_(tt.trip_count_, INVALID),
        final_footpaths_(tt.stations_.size(), INVALID),
        stats_{stats},
        MAX_DELAY(max_delay) {}

  void add_destination(csa_station const& station, time /*initial_duration*/) {
    // For now initial_duration is always 0
    targets_.emplace_back(station);
    // Initialize final_footpaths_
    for (auto const& fp : station.incoming_footpaths_) {
      final_footpaths_[fp.from_station_] =
          std::min(fp.duration_, final_footpaths_[fp.from_station_]);
    }
    stats_.destination_count_++;
  }

  time expected_value(time m, time d) {
    // Calculate E[X] with the CDF
    // E[X] = Integral from 0 to inf of (1-F(X))dx - Integral from -inf to 0 of
    // F(x) dx
    //      = (d*(11*ln(11*d)-11*ln(d)-10))/300 - (m*(4*ln(2*m)-4*ln(m)-5))/3
    // Simplifications reduce the term to two constants which can be
    // pre-calculated
    //      = d * (1/300 * (11 * log(11) - 10)) + m * (1/3 * (5 - 4 log(2)))
    //     ~= d * 0.05458 + m * 0.74247
    double const d_factor =
        0.0545894933359402532822712645253880743267959179777052964246808160;
    double const m_factor =
        0.7424704259200729207770238380557645758993331541863263278390933206;
    return std::ceil(d * d_factor + m * m_factor);
  }

  double f(double x, time m, time d) {
    if (x < 0) {
      return 0;
    } else if (0 <= x && x <= m) {
      return (2 * x) / (6 * m - 3 * x);
    } else if (m < x && x < m + d) {
      return (31 * (x - m) + 2 * d) / (30 * (x - m) + 3 * d);
    } else {
      return 1;
    }
  }

  time get_time_walking(const csa_connection& con) {
    auto time_walking = final_footpaths_[con.to_station_];
    if (time_walking != INVALID) {
      auto const expected_delay = expected_value(time_walking, MAX_DELAY);
      time_walking = con.arrival_ + expected_delay;
    }
    return time_walking;
  }

  time get_time_trip(const csa_connection& con) {
    return trip_arrival_time_[con.trip_];
  }

  time get_time_transfer(const csa_connection& con) {
    auto const transfer_time = tt_.stations_[con.to_station_].transfer_time_;

    auto const first_pair = std::lower_bound(
        arrival_time_[con.to_station_].begin(),
        arrival_time_[con.to_station_].end(), con.arrival_,
        [&](auto const& pair, time t) { return pair.first < t; });
    if (first_pair->first == INVALID) {
      return INVALID;
    }

    auto const safe_pair = std::upper_bound(
        arrival_time_[con.to_station_].begin(),
        arrival_time_[con.to_station_].end(), con.arrival_ + MAX_DELAY,
        [&](time t, auto const& pair) { return pair.first > t; });
    if (safe_pair->first == INVALID) {
      return INVALID;
    }

    std::list<std::pair<time, time>> relevant_pairs(first_pair,
                                                    std::next(safe_pair));
    relevant_pairs.insert(relevant_pairs.begin(), std::make_pair(0, 0));

    auto eat = 0;
    auto time_transfer = 0;
    for (auto prev_pair = relevant_pairs.begin();
         std::next(prev_pair) != relevant_pairs.end(); ++prev_pair) {
      auto const curr_pair = std::next(prev_pair);
      auto const probability =
          f(curr_pair->first - con.arrival_, transfer_time, MAX_DELAY) -
          f(prev_pair->first - con.arrival_, transfer_time, MAX_DELAY);

      if (curr_pair->second == INVALID) {
        time_transfer = INVALID;
      }
      eat += probability * curr_pair->second;
    }

    if (time_transfer != INVALID) {
      time_transfer = std::ceil(eat);
    }

    return time_transfer;
  }

  bool dominates(const std::pair<time, time>& a,
                 const std::pair<time, time>& b) {
    return (a.first <= b.first && a.second <= b.second) &&
           (a.first < b.first || a.second < b.second);
  }

  void search() {
    search([](auto const& con) {
      (void)con;
      return false;
    });
  }

  void search(std::function<bool(csa_connection const&)> skip_connection) {
    auto const& connections =
        Dir == search_dir::FWD ? tt_.fwd_connections_ : tt_.bwd_connections_;

    auto first_connection = connections.end();

    auto last_connection = connections.begin();
    csa_connection const end_at{search_interval_.begin_};
    last_connection =
        std::lower_bound(begin(connections), end(connections), end_at,
                         [&](csa_connection const& a, csa_connection const& b) {
                           return a.departure_ < b.departure_;
                         });
    if (last_connection == connections.end()) {
      return;
    }

    for (auto it = first_connection; it-- != last_connection;) {
      auto const& con = *it;
      if (skip_connection(con)) {
        continue;
      }

      auto const time_walking = get_time_walking(con);
      auto const time_trip = get_time_trip(con);
      auto const time_transfer = get_time_transfer(con);
      auto const min_arrival_time =
          std::min(time_walking, std::min(time_trip, time_transfer));

      if (min_arrival_time == INVALID) {
        continue;
      }
      auto const profile_pair =
          std::make_pair(con.departure_, min_arrival_time);

      bool is_non_dominated =
          std::none_of(arrival_time_[con.to_station_].begin(),
                       arrival_time_[con.to_station_].end(),
                       [&](std::pair<time, time> pair) {
                         return dominates(pair, profile_pair);
                       });

      if (is_non_dominated) {
        for (auto fp : tt_.stations_[con.from_station_].incoming_footpaths_) {
          add_to_profile(
              std::make_pair(con.departure_ - fp.duration_, min_arrival_time),
              fp.from_station_);
        }
      }

      trip_arrival_time_[con.trip_] = min_arrival_time;
    }
  }

  void add_to_profile(const std::pair<time, time>& profile_pair,
                      const uint32_t station_id) {

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

    // TODO(root) is this correct?
    if (std::any_of(it, arrival_time_[station_id].end(), [&](auto other_pair) {
          return !dominates(profile_pair, other_pair);
        })) {
      return;
    }

    it = arrival_time_[station_id].emplace(arrival_time_[station_id].begin(),
                                           profile_pair);
    for (auto r_it = std::make_reverse_iterator(it);
         r_it != arrival_time_[station_id].rend();) {
      if (dominates(*it, *r_it)) {
        r_it = std::make_reverse_iterator(
            arrival_time_[station_id].erase(std::next(r_it).base()));
      } else {
        ++it;
      }
    }
  }

  std::vector<csa_journey> get_results(csa_station const& station,
                                       bool include_equivalent) {
    utl::verify_ex(!include_equivalent,
                   std::system_error{error::include_equivalent_not_supported});

    (void)station;
    return {};
  }

  csa_timetable const& tt_;
  interval search_interval_;
  std::vector<std::list<std::pair<time, time>>> arrival_time_;
  std::vector<time> trip_arrival_time_;
  std::vector<csa_station> targets_;
  std::vector<time> final_footpaths_;
  csa_statistics& stats_;
  time MAX_DELAY;
  bool DEBUG_OUTPUT = false;
};

}  // namespace motis::csa::cpu

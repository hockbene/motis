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
        s_(tt.stations_.size(), {{INVALID, INVALID}}),
        t_(tt.trip_count_, INVALID),
        d_(tt.stations_.size(), INVALID),
        starts_(),
        stats_{stats},
        MAX_DELAY(max_delay) {}

  void add_destination(csa_station const& station, time /*initial_duration*/) {
    // For now initial_duration is always 0
    targets_.emplace_back(station);
    // Initialize d_
    for (auto const& fp : station.incoming_footpaths_) {
      d_[fp.from_station_] = std::min(fp.duration_, d_[fp.from_station_]);
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
    return std::floor(d * d_factor + m * m_factor);
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

  time max_D_c(const csa_connection& con) {
    return tt_.stations_[con.to_station_].transfer_time_ + MAX_DELAY;
  }

  time get_tau_1(const csa_connection& con, bool force_max_delay = false) {
    auto tau_1 = d_[con.to_station_];
    if (tau_1 != INVALID) {
      auto const expected_delay = force_max_delay
                                      ? tau_1 + MAX_DELAY
                                      : expected_value(tau_1, MAX_DELAY);
      tau_1 = con.arrival_ + expected_delay;
    }
    return tau_1;
  }

  time get_tau_2(const csa_connection& con, bool force_max_delay = false) {
    return t_[con.trip_];
  }

  time get_tau_3(const csa_connection& con, bool force_max_delay = false) {
    auto const transfer_time = tt_.stations_[con.to_station_].transfer_time_;

    auto const first_pair = std::lower_bound(
        s_[con.to_station_].begin(), s_[con.to_station_].end(), con.arrival_,
        [&](auto const& pair, time t) { return pair.first < t; });
    if (first_pair->first == INVALID) {
      return INVALID;
    }

    auto const safe_pair = std::upper_bound(
        s_[con.to_station_].begin(), s_[con.to_station_].end(),
        con.arrival_ + MAX_DELAY,
        [&](time t, auto const& pair) { return pair.first > t; });
    if (safe_pair->first == INVALID) {
      return INVALID;
    }

    if (force_max_delay) {
      return safe_pair->second;
    }

    std::list<std::pair<time, time>> relevant_pairs(first_pair,
                                                    std::next(safe_pair));
    relevant_pairs.insert(relevant_pairs.begin(), std::make_pair(0, 0));

    auto eat = 0;
    auto tau_3 = 0;
    for (auto it = relevant_pairs.begin();
         std::next(it) != relevant_pairs.end(); ++it) {
      auto const prev_p = *it;
      auto const curr_p = *std::next(it);
      auto const probability =
          f(curr_p.first - con.arrival_, transfer_time, MAX_DELAY) -
          f(prev_p.first - con.arrival_, transfer_time, MAX_DELAY);

      if (curr_p.second == INVALID) {
        tau_3 = INVALID;
      }
      eat += probability * curr_p.second;
    }

    if (tau_3 != INVALID) {
      tau_3 = std::floor(eat);
    }

    return tau_3;
  }

  bool dominates(const std::pair<time, time>& a,
                 const std::pair<time, time>& b) {
    return (a.first <= b.first && a.second <= b.second) &&
           (a.first < b.first || a.second < b.second);
  }

  void search(
      bool force_max_delay = false,
      std::function<bool(csa_connection const&)> skip_connection =
          [](auto const& con) { return false; }) {
    auto const& connections =
        Dir == search_dir::FWD ? tt_.fwd_connections_ : tt_.bwd_connections_;

    auto first_connection = connections.end();
    /*
    if (!force_max_delay) {
      csa_connection const start_at{search_interval_.end_};
      first_connection = std::upper_bound(
          begin(connections), end(connections), start_at,
          [&](csa_connection const& a, csa_connection const& b) {
            return a.departure_ > b.departure_;
          });
    }*/

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

      auto const tau_1 = get_tau_1(con, force_max_delay);
      auto const tau_2 = get_tau_2(con, force_max_delay);
      auto const tau_3 = get_tau_3(con, force_max_delay);
      auto const tau_c = std::min(tau_1, std::min(tau_2, tau_3));

      if (tau_c == INVALID) {
        continue;
      }
      auto const profile_pair = std::make_pair(con.departure_, tau_c);

      bool is_non_dominated =
          std::none_of(s_[con.to_station_].begin(), s_[con.to_station_].end(),
                       [&](std::pair<time, time> pair) {
                         return dominates(pair, profile_pair);
                       });

      if (is_non_dominated) {
        for (auto fp : tt_.stations_[con.from_station_].incoming_footpaths_) {
          add_to_profile(std::make_pair(con.departure_ - fp.duration_, tau_c),
                         fp.from_station_);
        }
      }

      t_[con.trip_] = tau_c;
    }

    /*
    int total = 0;
    int d_fail = 0;
    int a_fail = 0;
    for (auto station_idx = 0; station_idx < s_.size(); ++station_idx) {
      for (auto i = s_[station_idx].begin(); std::next(i) !=
    s_[station_idx].end(); ++i) { total++; auto d1 = i->first; auto a1 =
    i->second; auto d2 = std::next(i)->first; auto a2 = std::next(i)->second; if
    (!(d1 <= d2)) { d_fail++; LOG(motis::logging::debug) << "[*] d_fail at
    station " << station_idx;
        }
        if (!(a1 <= a2)) {
          a_fail++;
        }
      }
    }*/
  }

  void add_to_profile(const std::pair<time, time>& profile_pair,
                      const uint32_t station_id) {

    auto it = s_[station_id].begin();
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

    // TODO ist das richtig?
    if (std::any_of(it, s_[station_id].end(), [&](auto other_pair) {
          return !dominates(profile_pair, other_pair);
        })) {
      return;
    }

    // it = s_[station_id].insert(arrival_time_[station_id].begin(),
    // profile_pair);
    it = s_[station_id].emplace(s_[station_id].begin(), profile_pair);

    // TODO was wenn .begin() entfernt werden muss?
    /*
    while (it != s_[station_id].begin()) {
      if (dominates((*it), profile_pair)) {
        it = s_[station_id].erase(it);
      } else {
        --it;
      }
    }*/
    for (auto r_it = std::make_reverse_iterator(it);
         r_it != s_[station_id].rend();) {
      if (dominates(*it, *r_it)) {
        r_it = std::make_reverse_iterator(
            s_[station_id].erase(std::next(r_it).base()));
      } else {
        ++it;
      }
    }
  }

  std::vector<csa_journey> get_results(csa_station const& station,
                                       bool include_equivalent,
                                       bool force_max_delay = false) {
    utl::verify_ex(!include_equivalent,
                   std::system_error{error::include_equivalent_not_supported});

    std::vector<csa_journey> journeys;
    auto const& departure_time = s_[station.id_].begin()->first;
    auto const& arrival_time = s_[station.id_].begin()->second;
    csa_meat_reconstruction<Dir, decltype(s_), decltype(t_),
                            decltype(targets_)>{tt_, s_, t_, targets_}
        .extract_journey(
            journeys.emplace_back(Dir, departure_time, arrival_time, 0,
                                  &station, false),
            force_max_delay);
    return journeys;
  }

  csa_timetable const& tt_;  // Fahrplan
  interval search_interval_;  // Suchintervall
  std::vector<std::list<std::pair<time, time>>> s_;  // S - Profile
  std::vector<time> t_;  // T - Trips
  std::vector<csa_station> targets_;  // Zielstationen
  std::vector<time> d_;  // D - Fußpfade zum Ziel
  std::unordered_set<station_id> starts_;  // Startstationen
  csa_statistics& stats_;  // Stats
  time MAX_DELAY;
  bool DEBUG_OUTPUT = false;
};

}  // namespace motis::csa::cpu

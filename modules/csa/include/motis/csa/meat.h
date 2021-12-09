#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>
#include <unordered_map>

#include "motis/core/common/timing.h"
#include "motis/core/schedule/interval.h"
#include "motis/core/schedule/schedule.h"

#include "motis/csa/collect_start_times.h"
#include "motis/csa/csa_journey.h"
#include "motis/csa/csa_query.h"
#include "motis/csa/csa_statistics.h"
#include "motis/csa/csa_timetable.h"
#include "motis/csa/pareto_set.h"
#include "motis/csa/response.h"

#include "motis/csa/cpu/csa_meat_search.h"

namespace motis::csa {

template <typename CSAMeatSearch, typename CSAOnTripSearch>
struct meat {
  meat(schedule const& sched, csa_timetable const& tt, csa_query const& q,
       csa_statistics& stats)
      : sched_{sched},
        tt_{tt},
        q_{q},
        stats_{stats},
        schedule_begin_{SCHEDULE_OFFSET_MINUTES},
        schedule_end_{static_cast<motis::time>(
            (sched.schedule_end_ - sched.schedule_begin_) / 60)},
        search_interval_{q.search_interval_} {}

  time static constexpr MAX_DELAY = 30;
  double static constexpr ALPHA = 1.3;
  // TODO(root) Support for Relaxed Dominance
  // double static constexpr BETA = 0.0;

  response search() {
    csa_statistics stats;

    MOTIS_START_TIMING(total_timing);

    CSAMeatSearch eat_csa{tt_, search_interval_, stats, MAX_DELAY};
    for (auto const& dest_idx : q_.meta_dests_) {
      eat_csa.add_destination(tt_.stations_.at(dest_idx), 0);
    }

    CSAOnTripSearch ontrip_esat_csa{
        tt_, static_cast<time>(search_interval_.begin_), stats_};
    CSAOnTripSearch ontrip_eat_csa{
        tt_, static_cast<time>(search_interval_.begin_), stats_};
    for (auto const& start_idx : q_.meta_starts_) {
      ontrip_esat_csa.add_start(tt_.stations_.at(start_idx), 0);
      ontrip_eat_csa.add_start(tt_.stations_.at(start_idx), 0);
    }

    MOTIS_START_TIMING(search_timing);
    time tau_s = search_interval_.begin_;

    ontrip_esat_csa.search(
        MAX_DELAY, [&](auto const& con) { return con.departure_ < tau_s; });

    time esat_time = *std::min_element(
        ontrip_esat_csa.arrival_time_[q_.meta_dests_.at(0)].begin(),
        ontrip_esat_csa.arrival_time_[q_.meta_dests_.at(0)].end());
    time tau_last = std::ceil(tau_s + ALPHA * (esat_time - tau_s));

    ontrip_eat_csa.search(MAX_DELAY, [&](auto const& con) {
      return con.departure_ < tau_s || con.departure_ > tau_last;
    });

    std::vector<time> eat_times(tt_.stations_.size());
    std::transform(ontrip_eat_csa.arrival_time_.begin(),
                   ontrip_eat_csa.arrival_time_.end(), eat_times.begin(),
                   [&](auto const& arrival_times) {
                     return *std::min_element(arrival_times.begin(),
                                              arrival_times.end());
                   });

    eat_csa.search([&](auto const& con) {
      return con.arrival_ > tau_last ||
             eat_times[con.from_station_] > con.departure_;
    });
    MOTIS_STOP_TIMING(search_timing);

    MOTIS_START_TIMING(reconstruction_timing);
    std::vector<csa_journey> journeys;
    for (auto const& start_idx : q_.meta_starts_) {
      auto const& station_profile = eat_csa.arrival_time_[start_idx];
      auto const first_pair = std::lower_bound(
          station_profile.begin(), station_profile.end(), tau_s,
          [&](auto const& pair, time t) { return pair.first < t; });
      if (first_pair->first != INVALID_TIME) {
        extract_journeys(eat_csa.arrival_time_, journeys, {}, start_idx,
                         *first_pair);
      }
    }

    for (csa_journey& j : journeys) {
      results_.push_back(j);
    }

    std::unordered_map<station_id, std::vector<csa_journey::csa_edge>>
        graph_stations;
    std::vector<std::pair<csa_journey::csa_edge, csa_journey::csa_edge>> legs;

    for (csa_journey& j : journeys) {
      auto enter_con = csa_journey::csa_edge();
      auto exit_con = csa_journey::csa_edge();
      for (auto csa_edge : j.edges_) {
        if (csa_edge.enter_) {
          enter_con = csa_edge;
        }
        if (csa_edge.exit_) {
          exit_con = csa_edge;

          if (std::find_if(legs.begin(), legs.end(), [&](auto const& leg) {
                return leg.first.departure_ == enter_con.departure_ &&
                       leg.second.arrival_ == exit_con.arrival_;
              }) == legs.end()) {
            legs.emplace_back(enter_con, exit_con);
          }

          auto const enter_station = enter_con.from_->id_;
          auto const exit_station = exit_con.to_->id_;
          graph_stations.try_emplace(enter_station,
                                     std::vector<csa_journey::csa_edge>{});
          graph_stations.try_emplace(exit_station,
                                     std::vector<csa_journey::csa_edge>{});
          graph_stations.at(enter_station).emplace_back(enter_con);
          graph_stations.at(exit_station).emplace_back(exit_con);
        }
      }
    }

    if (!(graph_stations.empty() || legs.empty())) {
      std::string graph_string = generate_dot_graph(graph_stations, legs);
      std::ofstream graph_file;
      // TODO(root) have query id in filename
      std::string unique_id = std::to_string(q_.meta_starts_.at(0)) + "-" +
                              std::to_string(q_.meta_dests_.at(0)) + "-" +
                              std::to_string(search_interval_.begin_) + "-" +
                              std::to_string(search_interval_.end_);

      std::string filename = "responses-csa_meat-" + unique_id + ".dot";
      graph_file.open(filename);
      graph_file << graph_string;
      graph_file.close();

      // TODO(root) maybe already generate .png of graph?
      // command line call with
      // 'dot -Tpng -o responses-csa_meat.png responses-csa_meat.dot'
    }

    MOTIS_STOP_TIMING(reconstruction_timing);
    MOTIS_STOP_TIMING(total_timing);
    stats.search_duration_ += MOTIS_TIMING_MS(search_timing);
    stats.total_duration_ = MOTIS_TIMING_MS(total_timing);

    return {stats, std::move(results_.set_), search_interval_};
  }

  void extract_journeys(std::vector<std::list<std::pair<time, time>>> s_,
                        std::vector<csa_journey>& journeys,
                        csa_journey current_journey, station_id start_station,
                        std::pair<time, time> current_pair) {
    auto const transfer_time = tt_.stations_[start_station].transfer_time_;
    auto enter_con_it = std::find_if(
        tt_.stations_[start_station].outgoing_connections_.begin(),
        tt_.stations_[start_station].outgoing_connections_.end(),
        [&](auto const& con) {
          return con->departure_ == current_pair.first + transfer_time;
        });
    if (enter_con_it ==
        tt_.stations_[start_station].outgoing_connections_.end()) {
      return;
    }
    auto const& enter_con = *enter_con_it;
    auto const& trip_connections = tt_.trip_to_connections_[enter_con->trip_];
    auto const& first_connection =
        std::find_if(trip_connections.begin(), trip_connections.end(),
                     [enter_con](auto const& con) {
                       return con->departure_ == enter_con->departure_;
                     });
    for (auto it = trip_connections.rbegin();
         it != std::make_reverse_iterator(first_connection); ++it) {
      auto const& exit_con = *it;
      auto const& fp = tt_.stations_[exit_con->to_station_].footpaths_.at(0);
      auto jp = journey_pointer(enter_con, exit_con, &fp);

      if (is_destination(exit_con->to_station_)) {
        auto new_current_journey(current_journey);
        add_journey_leg(search_dir::FWD, tt_, new_current_journey, jp);

        new_current_journey.start_station_ =
            new_current_journey.edges_.front().from_;
        new_current_journey.start_time_ = new_current_journey.journey_begin();
        new_current_journey.destination_station_ =
            new_current_journey.edges_.back().to_;
        new_current_journey.arrival_time_ = new_current_journey.journey_end();
        new_current_journey.duration_ = new_current_journey.duration();
        journeys.emplace_back(new_current_journey);
        return;
      }

      auto const& station_profile = s_[exit_con->to_station_];
      auto const first_pair = std::lower_bound(
          station_profile.begin(), station_profile.end(), exit_con->arrival_,
          [&](auto const& pair, time t) { return pair.first < t; });
      if (first_pair->first == INVALID_TIME) {
        continue;
      }
      auto const safe_pair = std::upper_bound(
          station_profile.begin(), station_profile.end(),
          exit_con->arrival_ + MAX_DELAY,
          [&](time t, auto const& pair) { return pair.first > t; });
      if (safe_pair->first == INVALID_TIME) {
        continue;
      }
      bool found_exit = false;
      for (auto exit_pair = first_pair; exit_pair != std::next(safe_pair);
           ++exit_pair) {
        if (std::abs(exit_pair->second - current_pair.second) >= 1) {
          auto new_current_journey(current_journey);
          add_journey_leg(search_dir::FWD, tt_, new_current_journey, jp);
          new_current_journey.transfers_ += 1;
          extract_journeys(s_, journeys, new_current_journey,
                           exit_con->to_station_, *exit_pair);
          found_exit = true;
//          return;
        }
      }
      if (found_exit) {
        return;
      }
    }
  }

  //  TODO(root) support for compact representation
  std::string generate_dot_graph(
      std::unordered_map<station_id, std::vector<csa_journey::csa_edge>>&
          graph_stations,
      std::vector<std::pair<csa_journey::csa_edge, csa_journey::csa_edge>>&
          graph_legs) {
    std::ostringstream graph_str;
    graph_str << "digraph {\n"
              << "\trankdir=LR;\n"
              << "\tnode [ shape=record ];\n"
              << "\n";
    for (const auto& graph_station : graph_stations) {
      std::ostringstream station_str;
      station_str << "\t" << std::to_string(graph_station.first)
                  << " [ label = \""
                  << tt_.stations_.at(graph_station.first).station_ptr_->name_;
      std::set<time> connection_times;
      for (auto con : graph_station.second) {
        auto const timestamp = con.from_->id_ == graph_station.first
                                   ? con.departure_
                                   : con.arrival_;
        connection_times.insert(timestamp);
      }
      for (auto timestamp : connection_times) {
        station_str << "|<" + std::to_string(timestamp) + ">"
                    << format_time(timestamp).substr(0, 5);
      }
      station_str << "\"];\n";
      graph_str << station_str.str();
    }
    for (auto leg : graph_legs) {
      std::ostringstream leg_str;
      leg_str << "\t" << std::to_string(leg.first.from_->id_) << ":"
              << std::to_string(leg.first.departure_)
              //  << ":e"
              << " -> " << std::to_string(leg.second.to_->id_) << ":"
              << std::to_string(leg.second.arrival_);
      // leg_str << ":w"
      //  TODO(root) useful information about the trip (e.g. ICE 512)
      //   instead of MOTIS internal trip id
      leg_str << " [ label=\"" << std::to_string(leg.first.con_->trips_) << "\n"
              << "("
              << format_time(leg.second.arrival_ - leg.first.departure_)
                     .substr(0, 5)
              << "h)"
              << "\" ];\n";
      graph_str << leg_str.str();
    }
    graph_str << "}";
    return graph_str.str();
  }

private:
  bool is_start(station_id station_id) {
    return std::find_if(q_.meta_starts_.begin(), q_.meta_starts_.end(),
                        [&](auto const start_idx) {
                          return start_idx == station_id;
                        }) != q_.meta_starts_.end();
  }

  bool is_destination(station_id station_id) {
    return std::find_if(q_.meta_dests_.begin(), q_.meta_dests_.end(),
                        [&](auto const dest_idx) {
                          return dest_idx == station_id;
                        }) != q_.meta_dests_.end();
  }

  static bool dominates(csa_journey const& a, csa_journey const& b) {
    (void)a;
    (void)b;
    return false;
  }

  bool in_interval(csa_journey const& j) const {
    auto const t =
        query().dir_ == search_dir::FWD ? j.journey_begin() : j.journey_end();
    return t >= search_interval_.begin_ && t <= search_interval_.end_;
  }

  bool min_connection_count_reached() const {
    return std::count_if(begin(results_), end(results_),
                         [&](csa_journey const& j) {
                           return in_interval(j);
                         }) >= query().min_connection_count_;
  }

  bool max_interval_reached() const {
    return (!query().extend_interval_earlier_ ||
            search_interval_.begin_ == schedule_begin_) &&
           (!query().extend_interval_later_ ||  //
            search_interval_.end_ == schedule_end_);
  }

  motis::time map_to_interval(time const t) const {
    return std::min(schedule_end_, std::max(schedule_begin_, t));
  }

  csa_query const& query() const { return q_; }
  csa_statistics& stats() const { return stats_; }

  schedule const& sched_;
  csa_timetable const& tt_;
  csa_query const& q_;
  csa_statistics& stats_;
  motis::time schedule_begin_, schedule_end_;
  interval search_interval_{map_to_interval(q_.search_interval_.begin_),
                            map_to_interval(q_.search_interval_.end_)};
  pareto_set<csa_journey, decltype(&dominates)> results_{
      make_pareto_set<csa_journey>(&dominates)};
};

}  // namespace motis::csa

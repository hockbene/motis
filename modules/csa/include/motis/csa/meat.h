#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <queue>
#include <set>
#include <unordered_map>
#include <unordered_set>

#include "motis/core/common/timing.h"
#include "motis/core/schedule/interval.h"
#include "motis/core/schedule/schedule.h"

#include "motis/csa/collect_start_times.h"
#include "motis/csa/csa_query.h"
#include "motis/csa/csa_statistics.h"
#include "motis/csa/csa_timetable.h"
#include "motis/csa/pareto_set.h"
#include "motis/csa/response.h"

#include "motis/csa/cpu/csa_meat_search.h"

namespace motis::csa {

template <typename CSAMeatSearch>
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
  bool static constexpr FORCE_MAX_DELAY = true;
  double static constexpr ALPHA = 1.3;
  // TODO Support für Relaxed Dominance
  // relaxation tuning parameter beta
  // double static constexpr BETA = 0.0;

  response search() {
    // To solve the α-bounded (s, τ s , t )-MEAT problem, we perform the
    // following steps:
    csa_statistics stats;
    MOTIS_START_TIMING(total_timing);

    CSAMeatSearch esat_csa{tt_, search_interval_, stats, MAX_DELAY};

    for (auto const& dest_idx : q_.meta_dests_) {
      esat_csa.add_destination(tt_.stations_.at(dest_idx), 0);
    }

    MOTIS_START_TIMING(search_1_timing);
    // (1) Run a binary search on the connection set to determine the earliest
    //      connection c_first departing after τ_s .
    // (2) Run a one-to-one Connection Scan from s to t that assumes all
    // connections c are delayed by
    //      max_D_c to determine esat_time(s, τ_s, t).

    esat_csa.search(FORCE_MAX_DELAY);
    MOTIS_STOP_TIMING(search_1_timing);

    // (3) Let τ_last = τ_s + α · (esat_time(s, τ_s , t) − τ_s) and run a second
    // binary search on the
    //      connection set to find the last connection c_last departing before
    //      τ_last.
    // (4) Run a one-to-all Connection Scan from s restricted to the connections
    // from c_first to c_last
    //      to determine all eat(s, τ_s , ·).
    // TODO

    time tau_s = search_interval_.begin_;
    time esat_time =
        std::lower_bound(
            esat_csa.s_[q_.meta_starts_.at(0)].begin(),
            esat_csa.s_[q_.meta_starts_.at(0)].end(), tau_s,
            [&](auto const& pair, time t) { return pair.first < t; })
            ->second;
    time tau_last = tau_s + ALPHA * (esat_time - tau_s);

    CSAMeatSearch eat_csa{tt_, search_interval_, stats, MAX_DELAY};
    for (auto const& dest_idx : q_.meta_dests_) {
      eat_csa.add_destination(tt_.stations_.at(dest_idx), 0);
    }
    MOTIS_START_TIMING(search_2_timing);
    // (5) Run Phase 1 of the unbounded MEAT algorithm scanning the connections
    // from c_last to c_first,
    //      skipping connections c for which c_arr_time > τ_last or eat(s, τ_s ,
    //      c_dep_stop ) ≤ c_dep_time does not hold.
    /*
    auto skip_connections = [&](csa_connection const& con) {
      return con.arrival_ > tau_last ||
             std::any_of(q_.meta_starts_.begin(), q_.meta_starts_.end(),
                         [&](station_id const station) {
                           return eat(station, tau_s, con.from_station_) >
                                  con.departure_;
                         });
    };
    eat_csa.search(false, skip_connections);
     */
    MOTIS_STOP_TIMING(search_2_timing);

    // arrival time of an optimal safe earliest arrival journey
    auto const esat = [&](station_id start, time tau_s) {
      return find_earliest_profile_pair(esat_csa.s_[start], tau_s);
    };

    MOTIS_START_TIMING(reconstruction_timing);
    // (6) Finally, run Phase 2 of the unbounded MEAT algorithm; i.e., extract
    // the (s, τ s , t )-decision graph.

    // TODO Graph Extraction auslagern
    std::unordered_map<station_id, std::vector<csa_connection>> graph_stations;
    std::vector<std::pair<csa_connection, csa_connection>> legs;

    auto const queue_comparator = [](const csa_connection& a,
                                     const csa_connection& b) {
      return a.departure_ > b.departure_;
    };
    std::priority_queue<csa_connection, std::vector<csa_connection>,
                        decltype(queue_comparator)>
        queue(queue_comparator);

    for (auto const& start_idx : q_.meta_starts_) {
      auto const departure_time = esat(start_idx, tau_s).first;
      auto const initial_transfer_time =
          tt_.stations_[start_idx].transfer_time_;
      auto const earliest_connection_ptr = std::find_if(
          tt_.stations_[start_idx].outgoing_connections_.begin(),
          tt_.stations_[start_idx].outgoing_connections_.end(),
          [&](auto const& con) {
            return con->departure_ == departure_time + initial_transfer_time;
          });
      if (earliest_connection_ptr ==
          tt_.stations_[start_idx].outgoing_connections_.end()) {
        continue;
      }
      auto const earliest_connection = *earliest_connection_ptr;

      std::unordered_set<trip_id> used_trips;
      queue.push(*earliest_connection);

      while (!queue.empty()) {
        csa_connection const enter_con = queue.top();
        queue.pop();
        if (used_trips.find(enter_con.trip_) != used_trips.end()) {
          continue;
        }
        used_trips.insert(enter_con.trip_);
        csa_connection const exit_con = get_exit_con(enter_con);
        legs.emplace_back(enter_con, exit_con);
        if (!is_destination(exit_con.to_station_)) {
          auto const station_profile = esat_csa.s_[exit_con.to_station_];
          auto const exit_con_transfer_time =
              tt_.stations_[exit_con.to_station_].transfer_time_;
          auto const first_pair = std::lower_bound(
              station_profile.begin(), station_profile.end(), exit_con.arrival_,
              [&](auto const& pair, time t) { return pair.first < t; });
          auto const safe_pair = std::upper_bound(
              station_profile.begin(), station_profile.end(),
              exit_con.arrival_ + MAX_DELAY,
              [&](time t, auto const& pair) { return pair.first > t; });
          for (auto it = first_pair; it != std::next(safe_pair); ++it) {
            if (it->first == INVALID_TIME) {
              break;
            }
            auto const con = get_con_for_dep_time(
                it->first + exit_con_transfer_time, exit_con.to_station_);
            if (con.departure_ != INVALID_TIME) {
              queue.push(con);
            }
          }
        }
      }
    }

    // Extract all stations
    for (auto leg : legs) {
      auto const enter_station = leg.first.from_station_;
      auto const exit_station = leg.second.to_station_;
      graph_stations.try_emplace(enter_station, std::vector<csa_connection>{});
      graph_stations.try_emplace(exit_station, std::vector<csa_connection>{});
      graph_stations.at(enter_station).emplace_back(leg.first);
      graph_stations.at(exit_station).emplace_back(leg.second);
    }

    clean_graph(graph_stations, legs);

    // if (!(graph_stations.empty() || legs.empty())) {
    //  Generate .dot representation of graph
    //  TODO Support für kompakte Darstellung
    std::string output_string = generate_dot_graph(graph_stations, legs);
    std::ofstream graph_file;
    // TODO Query ID in den Dateinamen einbringen
    std::string unique_id = std::to_string(q_.meta_starts_.at(0)) + "-" +
                            std::to_string(q_.meta_dests_.at(0)) + "-" +
                            std::to_string(search_interval_.begin_) + "-" +
                            std::to_string(search_interval_.end_);

    std::string filename = "responses-csa_meat-" + unique_id + ".dot";
    graph_file.open(filename);
    graph_file << output_string;
    graph_file.close();

    // TODO Vielleicht schon .png generieren?
    // Aufruf der Kommandozeile mit
    // 'dot -Tpng -o responses-csa_meat.png responses-csa_meat.dot'
    //}

    MOTIS_STOP_TIMING(reconstruction_timing);
    MOTIS_STOP_TIMING(total_timing);
    stats.search_duration_ += MOTIS_TIMING_MS(search_1_timing);
    stats.search_duration_ += MOTIS_TIMING_MS(search_2_timing);
    stats.total_duration_ = MOTIS_TIMING_MS(total_timing);

    return {stats, std::move(results_.set_), search_interval_};
  }
  void clean_graph(
      std::unordered_map<station_id, std::vector<csa_connection>>&
          graph_stations,
      std::vector<std::pair<csa_connection, csa_connection>>& legs) {
    bool changes = true;
    while (changes) {
      changes = false;
      for (auto it = graph_stations.begin(); it != graph_stations.end();) {
        auto const graph_station = it->first;
        auto const connection_set = it->second;
        if (connection_set.empty() ||
            (!is_destination(graph_station) &&
             std::none_of(connection_set.begin(), connection_set.end(),
                          [&](auto const& con) {
                            return con.from_station_ == graph_station;
                          })) ||
            (!is_start(graph_station) &&
             std::none_of(connection_set.begin(), connection_set.end(),
                          [&](auto const& con) {
                            return con.to_station_ == graph_station;
                          }))) {
          it = graph_stations.erase(it);
          changes = true;
        } else {
          ++it;
        }
      }

      std::unordered_set<trip_id> seen_trips;
      for (auto it = legs.begin(); it != legs.end();) {
        auto const leg = *it;
        if (graph_stations.find(leg.first.from_station_) ==
                graph_stations.end() ||
            graph_stations.find(leg.second.to_station_) ==
                graph_stations.end() ||
            seen_trips.find(leg.first.trip_) != seen_trips.end()) {

          remove_leg_from_graph(graph_stations, leg);
          it = legs.erase(it);
          changes = true;
        } else {
          seen_trips.insert(leg.first.trip_);
          ++it;
        }
      }
    }
  }
  void remove_leg_from_graph(
      std::unordered_map<station_id, std::vector<csa_connection>>&
          graph_stations,
      const std::pair<csa_connection, csa_connection>& leg) const {
    auto const from_station = graph_stations.find(leg.first.from_station_);
    if (from_station != graph_stations.end()) {
      auto const con_it =
          std::find_if(from_station->second.begin(), from_station->second.end(),
                       [&](auto const& con) {
                         return leg.first.departure_ == con.departure_;
                       });
      if (con_it != from_station->second.end()) {
        from_station->second.erase(con_it);
      }
    }
    auto const to_station = graph_stations.find(leg.second.to_station_);
    if (to_station != graph_stations.end()) {
      auto const con_it = std::find_if(
          to_station->second.begin(), to_station->second.end(),
          [&](auto const& con) { return leg.second.arrival_ == con.arrival_; });
      if (con_it != to_station->second.end()) {
        to_station->second.erase(con_it);
      }
    }
  }

  std::string generate_dot_graph(
      std::unordered_map<station_id, std::vector<csa_connection>>&
          graph_stations,
      std::vector<std::pair<csa_connection, csa_connection>>& legs) {
    std::string output_string;
    output_string += "digraph {\n";
    output_string += "\trankdir=LR;\n";
    output_string += "\tnode [ shape=record ];\n";
    output_string += "\n";
    for (const auto& graph_station : graph_stations) {
      std::string station_str = "\t";
      station_str += std::to_string(graph_station.first);
      station_str += " [ label = \"";
      station_str += tt_.stations_.at(graph_station.first).station_ptr_->name_;
      std::set<time> connection_times;
      for (auto con : graph_station.second) {
        auto const timestamp = con.from_station_ == graph_station.first
                                   ? con.departure_
                                   : con.arrival_;
        connection_times.insert(timestamp);
      }
      for (auto timestamp : connection_times) {
        station_str += "|<" + std::to_string(timestamp) + ">";
        // TODO motis_time in lesbare Uhrzeit umwandeln (z.B. "%H:%M")
        station_str += std::to_string(timestamp);
      }
      station_str += "\"];\n";
      output_string += station_str;
    }
    for (auto leg : legs) {
      std::string leg_str = "\t";
      leg_str += std::to_string(leg.first.from_station_);
      leg_str += ":";
      leg_str += std::to_string(leg.first.departure_);
      // leg_str += ":e";
      leg_str += " -> ";
      leg_str += std::to_string(leg.second.to_station_);
      leg_str += ":";
      leg_str += std::to_string(leg.second.arrival_);
      // leg_str += ":w";
      //  TODO Sinnvolle Angabe des Trips (z.B. ICE 512) statt interner Trip-Id
      leg_str += " [ label=\"" + std::to_string(leg.first.trip_) + "\" ];\n";
      output_string += leg_str;
    }
    output_string += "}";
    return output_string;
  }

private:
  time eat(CSAMeatSearch csa, time tau_s, station_id start) {
    auto const& start_profile = csa.s_[start];
    return std::lower_bound(
               start_profile.begin(), start_profile.end(), tau_s,
               [&](auto const& pair, time t) { return pair.first < t; })
        ->second;
  }

  time eat(station_id start, time tau_s, station_id target) {
    CSAMeatSearch helper_csa{tt_, search_interval_, stats(), MAX_DELAY};
    helper_csa.add_destination(tt_.stations_.at(target), 0);
    helper_csa.search();
    auto const& start_profile = helper_csa.s_[start];
    return std::lower_bound(
               start_profile.begin(), start_profile.end(), tau_s,
               [&](auto const& pair, time t) { return pair.first < t; })
        ->second;
  }

  csa_connection get_con_for_dep_time(time departure_time,
                                      station_id station_id) {
    auto const earliest_connection =
        std::find_if(tt_.stations_[station_id].outgoing_connections_.begin(),
                     tt_.stations_[station_id].outgoing_connections_.end(),
                     [departure_time](auto const& con) {
                       return con->departure_ == departure_time;
                     });
    if (earliest_connection !=
        tt_.stations_[station_id].outgoing_connections_.end()) {
      return **earliest_connection;
    }
    return csa_connection(INVALID_TIME);
  }

  bool is_destination(station_id station_id) {
    return std::find_if(q_.meta_dests_.begin(), q_.meta_dests_.end(),
                        [&](auto const dest_idx) {
                          return dest_idx == station_id;
                        }) != q_.meta_dests_.end();
  }

  bool is_start(station_id station_id) {
    return std::find_if(q_.meta_starts_.begin(), q_.meta_starts_.end(),
                        [&](auto const start_idx) {
                          return start_idx == station_id;
                        }) != q_.meta_starts_.end();
  }

  csa_connection get_exit_con(csa_connection const& enter_con) {
    auto const& trip_connections = tt_.trip_to_connections_[enter_con.trip_];
    auto const& first_connection =
        std::find_if(trip_connections.begin(), trip_connections.end(),
                     [enter_con](auto const& con) {
                       return con->departure_ == enter_con.departure_;
                     });

    auto it = std::next(first_connection);
    for (; it != trip_connections.end(); ++it) {
      auto const exit_con = *it;
      auto const enter_con_eat =
          eat(enter_con.to_station_, enter_con.arrival_, q_.meta_dests_.at(0));
      auto const exit_con_eat =
          eat(exit_con->to_station_, exit_con->arrival_, q_.meta_dests_.at(0));
      if (enter_con_eat != exit_con_eat) {
        return *exit_con;
      }
    }
    return **std::prev(it);
  }

  std::pair<time, time> find_earliest_profile_pair(
      std::list<std::pair<time, time>> station_profile, time const tau_s) {
    return *std::lower_bound(
        station_profile.begin(), station_profile.end(), tau_s,
        [](auto const& pair, time t) { return pair.first < t; });
  }

  static bool dominates(csa_journey const& a, csa_journey const& b) {
    return a.journey_begin() >= b.journey_begin() &&
           a.journey_end() <= b.journey_end() && a.transfers_ <= b.transfers_;
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

  motis::time schedule_begin_, schedule_end_;
  interval search_interval_{map_to_interval(q_.search_interval_.begin_),
                            map_to_interval(q_.search_interval_.end_)};
  pareto_set<csa_journey, decltype(&dominates)> results_{
      make_pareto_set<csa_journey>(&dominates)};
  schedule const& sched_;
  csa_timetable const& tt_;
  csa_query const& q_;
  csa_statistics& stats_;
};

}  // namespace motis::csa

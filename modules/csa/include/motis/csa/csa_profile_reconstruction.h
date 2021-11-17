#pragma once

#include <unordered_map>

#include "utl/pipes.h"
#include "utl/verify.h"

#include "boost/filesystem.hpp"

#include "motis/core/common/logging.h"
#include "motis/core/access/service_access.h"
#include "motis/csa/csa_search_shared.h"
#include "motis/csa/csa_timetable.h"

namespace motis::csa {

template <search_dir Dir, typename ArrivalTimes, typename FinalFootpaths,
          typename TripReachable, typename Targets>
struct csa_profile_reconstruction {
  // TODO(root): fix getting type from ArrivalTimes
  /*
  using arrival_time_t = std::remove_reference_t<
      std::remove_cv_t<decltype(std::declval<ArrivalTimes>()[0])>>;
  */
  // for now hardcoded with motis::time
  static constexpr auto INVALID = Dir == search_dir::FWD
                                      ? std::numeric_limits<time>::max()
                                      : std::numeric_limits<time>::min();

  csa_profile_reconstruction(csa_timetable const& tt,
                             ArrivalTimes const& arrival_time,
                             FinalFootpaths const& final_footpaths,
                             TripReachable const& trip_reachable,
                             Targets const& targets)
      : tt_(tt),
        arrival_time_(arrival_time),
        final_footpaths_(final_footpaths),
        trip_reachable_(trip_reachable),
        targets_(targets) {}

  inline bool is_destination(station_id const station) {
    return targets_.find(station) != targets_.end();
  }

  inline bool is_destination(csa_station const& station) {
    return targets_.find(station.id_) != targets_.end();
  }

  // TODO(root) BWD
  // TODO(root) may be very inefficient.
  //  In the worst case, the profile gets scanned two times
  std::pair<time, std::array<time, MAX_TRANSFERS + 1>>
  find_earliest_profile_pair(csa_station const& station, int transfers,
                             time const tau_s) {
    auto const& profile = arrival_time_[station.id_];
    auto it = profile.begin();
    for (; it != profile.end(); ++it) {
      auto const& pair = *it;
      if (pair.first >= tau_s) {
        break;
      }
    }
    auto const& p = *it;

    auto const arrival_time = p.second[transfers];
    it = profile.end();
    for (; it-- != profile.begin();) {
      auto const& q = *it;
      if (q.second[transfers] == arrival_time) {
        return q;
      }
    }

    // return (INVALID, array_maker<time, MAX_TRANSFERS +
    // 1>::make_array(INVALID));
    return {};
  }
  // TODO(root) BWD
  void extract_journey(csa_journey& journey) {
    if (journey.is_reconstructed()) {
      return;
    }
    auto start = journey.start_station_;
    auto transfers = journey.transfers_;
    auto tau_s = journey.start_time_;
    auto const optimal_profile =
        find_earliest_profile_pair(*start, transfers, tau_s);

    tau_s = optimal_profile.first;
    journey.start_time_ = optimal_profile.first;

    for (; transfers != -1; --transfers) {

      auto jp = get_journey_pointer(*start, transfers, tau_s);

      if (jp.valid()) {
        assert(jp.enter_con_->trip_ == jp.exit_con_->trip_);
        add_journey_leg(journey, jp);
        start = &tt_.stations_[jp.exit_con_->to_station_];
        tau_s = jp.exit_con_->arrival_;
      }
    }
    journey.destination_station_ = start;
  }

  void add_journey_leg(csa_journey& journey, journey_pointer& jp) {
    if (jp.footpath_->from_station_ != jp.footpath_->to_station_) {
      if (Dir == search_dir::FWD) {
        journey.edges_.emplace_back(
            &tt_.stations_[jp.footpath_->from_station_],
            &tt_.stations_[jp.footpath_->to_station_],
            jp.enter_con_->departure_ - jp.footpath_->duration_,
            jp.enter_con_->departure_, -1);
      } else {
        journey.edges_.emplace_back(
            &tt_.stations_[jp.footpath_->from_station_],
            &tt_.stations_[jp.footpath_->to_station_], jp.exit_con_->arrival_,
            jp.exit_con_->arrival_ + jp.footpath_->duration_, -1);
      }
    }

    assert(jp.enter_con_->trip_ == jp.exit_con_->trip_);
    auto const& trip_connections =
        tt_.trip_to_connections_[jp.exit_con_->trip_];
    auto const add_trip_edge = [&](csa_connection const* con) {
      auto const enter = con == jp.enter_con_;
      auto const exit = con == jp.exit_con_;
      utl::verify(con->light_con_ != nullptr, "invalid light connection");
      journey.edges_.emplace_back(con->light_con_,
                                  &tt_.stations_[con->from_station_],
                                  &tt_.stations_[con->to_station_], enter, exit,
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

  void add_final_footpath(csa_journey& journey, const csa_station& from_station,
                          time arrival_time) {
    auto const duration = final_footpaths_[from_station.id_];
    auto const& footpaths = Dir == search_dir::FWD
                                ? from_station.footpaths_
                                : from_station.incoming_footpaths_;
    for (auto const& fp : footpaths) {
      auto const& to_station =
          Dir == search_dir::FWD ? fp.to_station_ : fp.from_station_;
      if (fp.duration_ == duration && is_destination(to_station)) {
        if (Dir == search_dir::FWD) {
          journey.edges_.emplace_back(
              from_station, journey.destination_station_, arrival_time,
              arrival_time + fp.duration_);
        } else {
          journey.edges_.emplace_back(journey.destination_station_,
                                      from_station, arrival_time - fp.duration_,
                                      arrival_time);
        }
        journey.destination_station_ = &tt_.stations_[to_station];
        return;
      }
    }

    LOG(motis::logging::warn) << "csa profile extract journey: could not find "
                                 "a matching final footpath to add";
  }

  journey_pointer get_journey_pointer(csa_station const& station, int transfers,
                                      time tau_s) {
    if (Dir == search_dir::FWD) {
      auto earliest_pair =
          find_earliest_profile_pair(station, transfers, tau_s);
      for (auto const& fp : station.footpaths_) {
        for (auto const& enter_con : get_enter_candidates(
                 tt_.stations_[fp.to_station_],
                 earliest_pair.first + fp.duration_, transfers)) {
          auto const& trip_connections =
              tt_.trip_to_connections_[enter_con->trip_];
          auto const& first_connection =
              std::find_if(trip_connections.begin(), trip_connections.end(),
                           [enter_con](auto const& con) {
                             return con->departure_ == enter_con->departure_;
                           });
          for (auto it = first_connection; it != trip_connections.end(); ++it) {
            auto const& exit_con = *it;
            auto const& arrival_station = tt_.stations_[exit_con->to_station_];
            auto const& exit_pair = find_earliest_profile_pair(
                arrival_station, transfers, exit_con->arrival_);
            // Exit here if
            // 1. The arrival time departing from the current stop (transfers -
            // 1) matches the arrival time when taking the current trip OR
            // 2. The connection ends at a target stop OR
            // 3. The connection ends at a stop from which the target can be
            // reached by foot
            if ((transfers > 0 && earliest_pair.second[transfers] ==
                                      exit_pair.second[transfers - 1]) ||
                (is_destination(exit_con->to_station_))) {
              return journey_pointer(enter_con, exit_con, &fp);
            }
          }
        }
      }
      /*
      // Otherwise, compute l_enter and l_exit
      auto earliest_pair =
          find_earliest_profile_pair(station, transfers, tau_s);
      std::unordered_map<csa_connection const*, footpath> candidates;

      // TODO(root) auslagern zu get_enter_candidates?
      for (auto const& fp : station.footpaths_) {
        for (auto const& candidate :
             tt_.stations_[fp.to_station_].outgoing_connections_) {
          if (candidate->departure_ == earliest_pair.first + fp.duration_) {
            candidates.emplace(candidate, fp);
          }
        }
      }
      // TODO(root)  Optionally prune candidate set
      for (auto const& candidate : candidates) {
        auto const& enter_con = candidate.first;
        auto const& trip_connections =
            tt_.trip_to_connections_[enter_con->trip_];
        auto const& first_connection =
            std::find_if(trip_connections.begin(), trip_connections.end(),
                         [enter_con](auto const& con) {
                           return con->departure_ == enter_con->departure_;
                         });
        for (auto it = first_connection; it != trip_connections.end(); ++it) {
          auto const& exit_con = *it;
          auto const& arrival_station = tt_.stations_[exit_con->to_station_];
          auto const& e_p = find_earliest_profile_pair(
              arrival_station, transfers, exit_con->arrival_);
          if ((transfers != 0 &&
               earliest_pair.second[transfers] == e_p.second[transfers - 1]) ||
              (std::find_if(targets_.begin(), targets_.end(),
                            [&](auto const& station) {
                              return station.id_ == arrival_station.id_;
                            }) != targets_.end())) {
              (is_destination(arrival_station))) {
                auto const fp = candidate.second;
                return journey_pointer(enter_con, exit_con, &fp);
              }
          }
        }*/
    } else {
      auto earliest_pair =
          find_earliest_profile_pair(station, transfers, tau_s);
      std::unordered_map<csa_connection const*, footpath> candidates;
      for (auto const& fp : station.incoming_footpaths_) {
        for (auto const& candidate :
             tt_.stations_[fp.from_station_].incoming_connections_) {
          if (candidate->arrival_ ==
              earliest_pair.second[transfers] - fp.duration_) {
            candidates.emplace(candidate, fp);
          }
        }
      }
      for (auto const& candidate : candidates) {
        auto const& enter_con = candidate.first;
        auto const& trip_connections =
            tt_.trip_to_connections_[enter_con->trip_];
        auto const& first_connection =
            std::find_if(trip_connections.begin(), trip_connections.end(),
                         [enter_con](auto const& con) {
                           return con->departure_ == enter_con->departure_;
                         });
        for (auto it = first_connection; it != trip_connections.end(); ++it) {
          auto const& exit_con = *it;
          auto const& arrival_station = tt_.stations_[exit_con->to_station_];
          auto const& e_p = find_earliest_profile_pair(
              arrival_station, transfers, exit_con->arrival_);
          if ((transfers != 0 && earliest_pair.first == e_p.first) ||
              (is_destination(arrival_station.id_))) {
            auto const fp = candidate.second;
            return journey_pointer(enter_con, exit_con, &fp);
          }
        }
      }
    }
    return {};
  }

  auto get_exit_candidates(csa_station const& arrival_station,
                           time arrival_time, int transfers) const {
    return utl::all(arrival_station.incoming_connections_)  //
           | utl::remove_if(
                 [this, arrival_time, transfers](csa_connection const* con) {
                   return con->arrival_ != arrival_time ||
                          // !trip_reachable_[con->trip_][transfers - 1] ||
                          // TODO(root) is this necessary?
                          !con->to_out_allowed_;
                 })  //
           | utl::iterable();
  }

  auto get_enter_candidates(csa_station const& departure_station,
                            time departure_time, int transfers) const {
    return utl::all(departure_station.outgoing_connections_)  //
           | utl::remove_if(
                 [this, departure_time, transfers](csa_connection const* con) {
                   return con->departure_ != departure_time ||
                          // !trip_reachable_[con->trip_][transfers - 1] ||
                          // TODO(root) is this necessary?
                          !con->from_in_allowed_;
                 })  //
           | utl::iterable();
  }

  csa_timetable const& tt_;
  ArrivalTimes const& arrival_time_;
  FinalFootpaths const& final_footpaths_;
  TripReachable const& trip_reachable_;
  Targets const& targets_;
};

}  // namespace motis::csa
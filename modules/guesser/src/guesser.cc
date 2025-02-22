#include "motis/guesser/guesser.h"

#include <algorithm>
#include <numeric>

#include "motis/hash_set.h"

#include "boost/program_options.hpp"

#include "utl/to_vec.h"

#include "motis/core/common/logging.h"
#include "motis/core/schedule/schedule.h"
#include "motis/protocol/Message_generated.h"

using namespace flatbuffers;
using namespace motis::module;
using motis::logging::info;

namespace motis::guesser {

std::string trim(std::string const& s) {
  auto first = s.find_first_not_of(' ');
  auto last = s.find_last_not_of(' ');
  if (first == last) {
    return "";
  } else {
    return s.substr(first, (last - first + 1));
  }
}

guesser::guesser() : module("Guesser Options", "guesser") {}

void guesser::init(motis::module::registry& reg) {
  update_stations();
  reg.register_op("/guesser", [this](msg_ptr const& m) { return guess(m); });
  reg.subscribe("/rt/update", [this](msg_ptr const& m) {
    using namespace motis::rt;
    auto const update = motis_content(RtUpdates, m);
    if (std::any_of(update->updates()->begin(), update->updates()->end(),
                    [](RtUpdate const* u) {
                      return u->content_type() == Content_RtStationAdded;
                    })) {
      update_stations();
    }
    return nullptr;
  });
}

void guesser::update_stations() {
  auto const& sched = get_sched();

  mcd::hash_set<std::string> station_names;
  station_indices_.clear();
  for (auto const& s : sched.stations_) {
    if (s->dummy_) {
      continue;
    }
    if (station_names.insert(s->name_.str()).second) {
      station_indices_.push_back(s->index_);
    }
  }

  auto stations = utl::to_vec(station_indices_, [&](auto const station_idx) {
    auto const& s = *sched.stations_[station_idx];
    float factor = 0;
    for (auto i = 0UL; i < s.dep_class_events_.size(); ++i) {
      factor +=
          std::pow(
              10,
              (static_cast<service_class_t>(service_class::NUM_CLASSES) - i) /
                  3) *
          s.dep_class_events_.at(i);
    }
    return std::make_pair(s.name_.str(), factor);
  });

  if (!stations.empty()) {
    auto max_importance =
        std::max(std::max_element(begin(stations), end(stations),
                                  [](std::pair<std::string, float> const& lhs,
                                     std::pair<std::string, float> const& rhs) {
                                    return lhs.second < rhs.second;
                                  })
                     ->second,
                 1.0F);
    for (auto& s : stations) {
      s.second = 1 + (s.second / max_importance) * 0.5;
    }
  } else {
    LOG(info) << "no stations found";
  }

  guesser_ = std::make_unique<guess::guesser>(stations);
}

msg_ptr guesser::guess(msg_ptr const& msg) {
  auto req = motis_content(StationGuesserRequest, msg);

  message_creator b;
  std::vector<Offset<Station>> guesses;
  for (auto const& match :
       guesser_->guess_match(trim(req->input()->str()), req->guess_count())) {
    auto const guess = match.index;
    auto const& station = *get_sched().stations_[station_indices_[guess]];
    auto const pos = Position(station.width_, station.length_);
    guesses.emplace_back(CreateStation(b, b.CreateString(station.eva_nr_),
                                       b.CreateString(station.name_), &pos));
  }

  b.create_and_finish(
      MsgContent_StationGuesserResponse,
      CreateStationGuesserResponse(b, b.CreateVector(guesses)).Union());

  return make_msg(b);
}

}  // namespace motis::guesser

#include "motis/ris/ribasis/ribasis_parser.h"

#include <cstdint>

#include "rapidjson/document.h"
#include "rapidjson/error/en.h"

#include "utl/get_or_create.h"
#include "utl/to_vec.h"
#include "utl/verify.h"

#include "motis/core/common/date_time_util.h"
#include "motis/core/common/logging.h"
#include "motis/core/common/unixtime.h"
#include "motis/core/schedule/event_type.h"

#include "motis/ris/ribasis/common.h"

using namespace motis::logging;
using namespace flatbuffers;

namespace motis::ris::ribasis {

Offset<StationInfo> parse_station(context& ctx, rapidjson::Value const& stop) {
  auto const eva = get_str(stop, "evanummer");
  auto const ds100 = get_optional_str(stop, "rl100");
  auto const name = get_str(stop, "bezeichnung");
  return utl::get_or_create(ctx.stations_, eva, [&]() {
    return CreateStationInfo(
        ctx.b_, ctx.b_.CreateSharedString(eva.data(), eva.size()),
        ctx.b_.CreateString(ds100), ctx.b_.CreateString(name));
  });
}

Offset<FullTripId> parse_trip_id(context& ctx, rapidjson::Value const& data) {
  auto const& rel = get_obj(data, "fahrtrelation");
  auto const& start_stop = get_obj(rel, "starthaltestelle");
  auto const& dest_stop = get_obj(rel, "zielhaltestelle");
  auto const train_nr =
      get_parsed_number<std::uint32_t>(rel, "startfahrtnummer");
  auto const line = get_optional_str(rel, "startlinie");
  auto const start_station_eva = get_str(start_stop, "evanummer");
  auto const dest_station_eva = get_str(dest_stop, "evanummer");

  auto const start_si = parse_station(ctx, start_stop);
  auto const dest_si = parse_station(ctx, dest_stop);

  auto const start_time = get_schedule_timestamp(ctx, rel, "startzeit");
  auto const dest_time = get_schedule_timestamp(ctx, rel, "zielzeit");

  return CreateFullTripId(
      ctx.b_,
      CreateTripId(ctx.b_,
                   ctx.b_.CreateSharedString(start_station_eva.data(),
                                             start_station_eva.size()),
                   train_nr, start_time,
                   ctx.b_.CreateSharedString(dest_station_eva.data(),
                                             dest_station_eva.size()),
                   dest_time, ctx.b_.CreateString(line)),
      start_si, dest_si);
}

void parse_categories(context& ctx, rapidjson::Value const& data) {
  for (auto const& g : get_array(data, "allGattung")) {
    utl::verify(g.IsObject(), "invalid allGattung entry");
    ctx.categories_.emplace(
        std::string{get_str(g, "gattungid")},
        CreateCategoryInfo(ctx.b_, ctx.b_.CreateString(get_str(g, "name")),
                           ctx.b_.CreateString(get_str(g, "code"))));
  }
}

void parse_lines(context& ctx, rapidjson::Value const& data) {
  for (auto const& l : get_array(data, "allLinie")) {
    utl::verify(l.IsObject(), "invalid allLinie entry");
    ctx.lines_.emplace(std::string{get_str(l, "linieid")},
                       ctx.b_.CreateString(get_str(l, "name")));
  }
}

void parse_providers(context& ctx, rapidjson::Value const& data) {
  for (auto const& p : get_array(data, "allVerwaltung")) {
    utl::verify(p.IsObject(), "invalid allVerwaltung entry");
    auto const& info = get_obj(p, "betreiber");
    auto const id = get_str(p, "verwaltungid");
    ctx.providers_.emplace(
        std::string{id},
        CreateProviderInfo(ctx.b_, ctx.b_.CreateString(id),
                           ctx.b_.CreateString(get_str(info, "name")),
                           ctx.b_.CreateString(get_str(info, "code"))));
  }
}

Offset<String> parse_track(context& ctx, rapidjson::Value const& ev,
                           char const* key) {
  auto const& obj = get_value(ev, key);
  if (obj.IsNull()) {
    return ctx.b_.CreateSharedString("");
  } else {
    auto const track = get_str(obj, "bezeichnung");
    return ctx.b_.CreateSharedString(track.data(), track.size());
  }
}

TimestampType parse_timestamp_type(rapidjson::Value const& ev,
                                   char const* key) {
  auto const s = get_str(ev, key);
  if (s == "FAHRPLAN") {
    return TimestampType_Schedule;
  } else if (s == "MELDUNG") {
    return TimestampType_Is;
  } else if (s == "PROGNOSE") {
    return TimestampType_Forecast;
  } else {
    // TODO(pablo): AUTOMAT
    return TimestampType_Unknown;
  }
}

Offset<TripEvent> parse_event(context& ctx, rapidjson::Value const& ev,
                              event_type const ev_type) {
  auto const station = parse_station(ctx, get_obj(ev, "haltestelle"));
  auto const schedule_time = get_schedule_timestamp(
      ctx, ev,
      ev_type == event_type::DEP ? "planabfahrtzeit" : "planankunftzeit");
  auto const current_time = get_schedule_timestamp(
      ctx, ev, ev_type == event_type::DEP ? "abfahrtzeit" : "ankunftzeit");
  auto const current_time_type = parse_timestamp_type(
      ev,
      ev_type == event_type::DEP ? "abfahrtzeitstatus" : "ankunftzeitstatus");
  auto const interchange_allowed = get_bool(ev, "fahrgastwechsel");
  auto const schedule_track = parse_track(
      ctx, ev,
      ev_type == event_type::DEP ? "planabfahrtort" : "planankunftort");
  auto const current_track = parse_track(
      ctx, ev, ev_type == event_type::DEP ? "abfahrtort" : "ankunftort");
  return CreateTripEvent(ctx.b_, station, schedule_time, current_time,
                         current_time_type, interchange_allowed, schedule_track,
                         current_track);
}

Offset<TripSection> parse_section(context& ctx, rapidjson::Value const& sec) {
  auto const train_nr = get_parsed_number<std::uint32_t>(sec, "fahrtnummer");
  auto const& category = ctx.categories_.at(get_str(sec, "gattungid"));
  auto const line_id = get_optional_str(sec, "linieid");
  auto const& line =
      !line_id.empty() ? ctx.lines_.at(line_id) : ctx.b_.CreateSharedString("");
  auto const& provider = ctx.providers_.at(get_str(sec, "verwaltungid"));
  return CreateTripSection(
      ctx.b_, train_nr, category, line, provider,
      parse_event(ctx, get_obj(sec, "abfahrt"), event_type::DEP),
      parse_event(ctx, get_obj(sec, "ankunft"), event_type::ARR));
}

Offset<Vector<Offset<TripSection>>> parse_sections(
    context& ctx, rapidjson::Value::ConstArray const& sections_data) {
  return ctx.b_.CreateVector(utl::to_vec(sections_data, [&](auto const& sec) {
    utl::verify(sec.IsObject(), "invalid allFahrtabschnitt entry");
    return parse_section(ctx, sec);
  }));
}

void ribasis_parser::to_ris_message(
    std::string_view s, const std::function<void(ris_message&&)>& cb) {
  rapidjson::Document doc;
  if (doc.Parse(s.data(), s.size()).HasParseError()) {
    doc.GetParseError();
    LOG(error) << "RI Basis: Bad JSON: "
               << rapidjson::GetParseError_En(doc.GetParseError())
               << " at offset " << doc.GetErrorOffset();
    return;
  }

  try {
    utl::verify(doc.IsObject(), "no root object");
    auto const& meta = get_obj(doc, "meta");
    auto const& data = get_obj(doc, "data");
    auto const created_at = get_timestamp(meta, "created", "%FT%H:%M:%S%Ez");
    auto ctx = context{created_at};
    auto const ext_trip_ref = get_str(data, "fahrtid");
    auto const trp_id = parse_trip_id(ctx, data);
    parse_categories(ctx, data);
    parse_lines(ctx, data);
    parse_providers(ctx, data);
    auto const sections_data = get_array(data, "allFahrtabschnitt");
    auto const sections = parse_sections(ctx, sections_data);
    auto const trip_msg = CreateFullTripMessage(
        ctx.b_, trp_id, ctx.b_.CreateString(ext_trip_ref), sections);
    ctx.b_.Finish(CreateMessage(ctx.b_, ctx.earliest_, ctx.latest_,
                                ctx.timestamp_, MessageUnion_FullTripMessage,
                                trip_msg.Union()));
    utl::verify(ctx.earliest_ != std::numeric_limits<unixtime>::max(),
                "earliest not set");
    utl::verify(ctx.latest_ != std::numeric_limits<unixtime>::min(),
                "latest not set");
    cb(ris_message{ctx.earliest_, ctx.latest_, ctx.timestamp_,
                   std::move(ctx.b_)});
  } catch (std::runtime_error const& e) {
    LOG(error) << "unable to parse RI Basis message: " << e.what();
  }
}

std::vector<ris_message> ribasis_parser::parse(std::string_view s) {
  std::vector<ris_message> msgs;
  to_ris_message(s, [&](ris_message&& m) { msgs.emplace_back(std::move(m)); });
  return msgs;
}

}  // namespace motis::ris::ribasis

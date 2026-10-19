#include "OscarReader.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

namespace {
constexpr int kMaxParticles = std::numeric_limits<unsigned short>::max();
constexpr int kElastic = 1;
constexpr int kDecay = 5;
constexpr int kWall = 6;
constexpr int kThermalization = 7;
constexpr int kHypersurfaceCrossing = 8;
constexpr int kFailedString = 47;

bool IsHadronicInteraction(int type) {
    return type != 0 && type != kWall && type != kHypersurfaceCrossing && type != kFailedString;
}

bool ReadCount(std::istringstream& fields, int& count) {
    return bool(fields >> count) && count >= 0 && count <= kMaxParticles;
}

bool AtEnd(std::istringstream& fields) {
    fields >> std::ws;
    return fields.eof();
}
} // namespace

bool OscarReader::Fail(const std::string& message) {
    error_ = "line " + std::to_string(line_number_) + ": " + message;
    return false;
}

bool OscarReader::ReadLine(std::string& line) {
    while (std::getline(input_, line)) {
        ++line_number_;
        const auto start = line.find_first_not_of(" \t\r");
        if (start == std::string::npos)
            continue;
        line.erase(0, start);
        return true;
    }
    if (input_.bad())
        Fail("input read failed");
    return false;
}

bool OscarReader::ReadParticles(int count, std::vector<Particle>& particles) {
    particles.clear();
    particles.reserve(static_cast<std::size_t>(count));
    std::unordered_set<int> ids;
    for (int i = 0; i < count; ++i) {
        std::string line;
        if (!ReadLine(line))
            return Fail("unexpected EOF in particle block");
        if (line.front() == '#')
            return Fail("particle block ended before its declared count");
        Particle p;
        std::istringstream fields(line);
        if (!(fields >> p.t >> p.x >> p.y >> p.z >> p.mass >> p.p0 >> p.px >> p.py >> p.pz >>
              p.pdg >> p.id >> p.charge)) {
            return Fail("invalid OSCAR2013 particle record");
        }
        for (double value : {p.t, p.x, p.y, p.z, p.mass, p.p0, p.px, p.py, p.pz}) {
            if (!std::isfinite(value))
                return Fail("non-finite particle value");
        }
        for (std::size_t column = 12; column < particle_columns_; ++column) {
            double value = 0.;
            if (!(fields >> value) || !std::isfinite(value))
                return Fail("invalid extended particle record");
        }
        if (!AtEnd(fields))
            return Fail("particle column count does not match format header");
        if (p.id < 0 || p.pdg == 0 || p.mass < 0. || p.p0 <= 0.) {
            return Fail("invalid particle ID, PDG, mass or energy");
        }
        if (!ids.insert(p.id).second)
            return Fail("duplicate particle ID in block");
        particles.push_back(p);
    }
    return true;
}

bool OscarReader::ReadNext(OscarEvent& event) {
    if (!error_.empty())
        return false;
    event = OscarEvent{};
    std::unordered_map<int, Particle> initial, live, last_interaction;
    std::unordered_set<int> all_ids, participants;
    bool started = false, final_seen = false, hadronic_seen = false;
    std::vector<Particle> block;
    std::string line;

    while (ReadLine(line)) {
        if (line.front() != '#')
            return Fail("particle outside a declared block");
        std::istringstream fields(line);
        std::string marker, kind;
        fields >> marker;
        if (marker.rfind("#!", 0) == 0) {
            std::string content;
            fields >> content;
            if (started || format_seen_)
                return Fail("unexpected format header");
            if ((marker != "#!OSCAR2013" && marker != "#!OSCAR2013Extended") ||
                content != "full_event_history") {
                return Fail("expected OSCAR2013 full_event_history with start/end lists");
            }
            for (const auto* expected :
                 {"t", "x", "y", "z", "mass", "p0", "px", "py", "pz", "pdg", "ID", "charge"}) {
                std::string column;
                if (!(fields >> column) || column != expected)
                    return Fail("unsupported OSCAR particle columns");
                ++particle_columns_;
            }
            std::string column;
            while (fields >> column)
                ++particle_columns_;
            format_seen_ = true;
            continue;
        }
        if (marker != "#")
            return Fail("invalid comment/header marker");
        fields >> kind;
        if (kind != "event" && kind != "interaction")
            continue;
        if (!format_seen_)
            return Fail("missing OSCAR2013 format header");

        if (kind == "interaction") {
            if (!started || final_seen)
                return Fail("interaction outside an active event");
            std::string key;
            int nin = -1, nout = -1, type = -1;
            if (!(fields >> key) || key != "in" || !ReadCount(fields, nin) || !(fields >> key) ||
                key != "out" || !ReadCount(fields, nout)) {
                return Fail("invalid interaction particle counts");
            }
            while (fields >> key) {
                if (key == "type") {
                    if (!(fields >> type) || type < 0 || !AtEnd(fields))
                        return Fail("invalid interaction type");
                    break;
                }
                double value = 0.;
                if (!(fields >> value) || !std::isfinite(value))
                    return Fail("invalid interaction metadata");
            }
            if (type < 0)
                return Fail("missing interaction type");
            int expected_out = -1;
            switch (type) {
            case 1:
                expected_out = 2;
                break;
            case 2:
                expected_out = 1;
                break;
            case 3:
                expected_out = 2;
                break;
            case 4:
                expected_out = 3;
                break;
            case 13:
                expected_out = 5;
                break;
            case 15:
                expected_out = 4;
                break;
            default:
                break;
            }
            if ((expected_out >= 0 && (nin != 2 || nout != expected_out)) ||
                (type == kDecay && (nin != 1 || nout < 1))) {
                return Fail("particle counts do not match interaction type");
            }
            std::vector<Particle> incoming, outgoing;
            if (!ReadParticles(nin, incoming) || !ReadParticles(nout, outgoing))
                return false;
            std::unordered_set<int> incoming_ids;
            const bool hadronic = IsHadronicInteraction(type);
            hadronic_seen = hadronic_seen || hadronic;
            for (const auto& p : incoming) {
                const auto found = live.find(p.id);
                if (found == live.end() || found->second.pdg != p.pdg) {
                    return Fail("incoming particle has no matching live history: ID " +
                                std::to_string(p.id));
                }
                incoming_ids.insert(p.id);
                const auto original = initial.find(p.id);
                if (hadronic && original != initial.end() && IsNucleon(original->second.pdg))
                    participants.insert(p.id);
            }
            for (const auto& p : outgoing) {
                if (all_ids.count(p.id) && !incoming_ids.count(p.id)) {
                    return Fail("outgoing particle reuses an unrelated ID: " +
                                std::to_string(p.id));
                }
            }
            for (const auto& p : incoming)
                live.erase(p.id);
            for (const auto& p : outgoing) {
                live[p.id] = p;
                all_ids.insert(p.id);
                if (hadronic)
                    last_interaction[p.id] = p;
            }
            if (hadronic && nin == 2 && type != kDecay && type != kThermalization) {
                if (event.binary_collisions == std::numeric_limits<int>::max())
                    return Fail("too many binary collisions");
                ++event.binary_collisions;
                if (IsNucleon(incoming[0].pdg) && IsNucleon(incoming[1].pdg)) {
                    if (type == kElastic)
                        ++event.nn_elastic_collisions;
                    else
                        ++event.nn_inelastic_collisions;
                }
            }
            continue;
        }

        int number = -1, ensemble = 0;
        std::string phase;
        if (!(fields >> number >> phase) || number < 0)
            return Fail("invalid event header");
        if (phase == "ensemble") {
            if (!(fields >> ensemble >> phase) || ensemble < 0)
                return Fail("invalid ensemble header");
        }
        if (phase == "in") {
            if (started)
                return Fail(
                    "new event before previous end (interleaved ensembles are unsupported)");
            int count = -1;
            if (!ReadCount(fields, count) || !AtEnd(fields))
                return Fail("invalid initial particle count");
            if (!ReadParticles(count, block))
                return false;
            event.number = number;
            event.ensemble = ensemble;
            for (const auto& p : block) {
                initial.emplace(p.id, p);
                live.emplace(p.id, p);
                all_ids.insert(p.id);
                if (IsNucleon(p.pdg))
                    ++event.initial_nucleons;
            }
            started = true;
            continue;
        }
        if (!started || number != event.number || ensemble != event.ensemble)
            return Fail("event/ensemble header does not match active event");
        if (phase == "out") {
            if (final_seen)
                return Fail("duplicate final particle block");
            int count = -1;
            if (!ReadCount(fields, count) || !AtEnd(fields))
                return Fail("invalid final particle count");
            if (!ReadParticles(count, block))
                return false;
            if (block.size() != live.size())
                return Fail("final list does not match complete interaction history");
            bool first = true;
            for (auto p : block) {
                const auto found = live.find(p.id);
                if (found == live.end() || found->second.pdg != p.pdg)
                    return Fail("final particle has no matching live history: ID " +
                                std::to_string(p.id));
                if (first) {
                    event.final_time = p.t;
                    first = false;
                } else if (std::abs(event.final_time - p.t) > 1e-5 * std::max(1., std::abs(p.t))) {
                    return Fail("final particle list contains different output times");
                }
                const auto original = initial.find(p.id);
                p.spectator = original != initial.end() && IsNucleon(original->second.pdg) &&
                              !participants.count(p.id);
                if (p.spectator)
                    ++event.spectators;
                if (!p.spectator) {
                    const auto last = last_interaction.find(p.id);
                    const Particle* emission = nullptr;
                    if (last != last_interaction.end())
                        emission = &last->second;
                    else if (original != initial.end())
                        emission = &original->second;
                    if (!emission)
                        return Fail("missing emission point for final particle");
                    // Final momenta are preserved; only emission coordinates change.
                    p.t = emission->t;
                    p.x = emission->x;
                    p.y = emission->y;
                    p.z = emission->z;
                }
                event.particles.push_back(p);
            }
            final_seen = true;
            continue;
        }
        if (phase != "end")
            return Fail("unsupported event phase: " + phase);
        if (!final_seen)
            return Fail("missing final particle block (enable Print_Start_End: Yes)");
        int zero = -1;
        std::string key;
        if (!(fields >> zero >> key >> event.impact_parameter) || zero != 0 || key != "impact" ||
            !std::isfinite(event.impact_parameter) || event.impact_parameter < 0.)
            return Fail("invalid event end/impact parameter");
        if (!(fields >> key)) {
            // Legacy fallback when no projectile/target flag was written.
            event.projectile_target_scattering = event.binary_collisions > 0;
        } else {
            std::string value;
            if ((key != "scattering_projectile_target" && key != "empty") || !(fields >> value) ||
                (value != "yes" && value != "no") || !AtEnd(fields))
                return Fail("invalid projectile/target scattering flag");
            event.projectile_target_scattering = key == "empty" ? value == "no" : value == "yes";
        }
        event.participants = static_cast<int>(participants.size());
        if (event.projectile_target_scattering && !hadronic_seen) {
            return Fail("scattering flag yes without interaction history");
        }
        std::sort(event.particles.begin(), event.particles.end(),
                  [](const Particle& a, const Particle& b) { return a.id < b.id; });
        return true;
    }
    if (!error_.empty())
        return false;
    if (started)
        return Fail("unexpected EOF before event end");
    if (!format_seen_)
        return Fail("missing OSCAR2013 format header");
    return false;
}

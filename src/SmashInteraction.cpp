#include "SmashInteraction.h"

#include <charconv>
#include <istream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

bool IsInelastic(int type) {
    // SMASH ProcessType / oscaroutput.cc. A failed string (47) falls back to
    // elastic scattering; wall/hypersurface crossings and artificial sources
    // do not establish inelastic participation.
    switch (type) {
    case 0:
    case 1:
    case 6:
    case 8:
    case 47:
    case 90:
        return false;
    case 2:
    case 3:
    case 4:
    case 5:
    case 7:
    case 9:
    case 10:
    case 11:
    case 12:
    case 13:
    case 14:
    case 15:
    case 41:
    case 42:
    case 43:
    case 44:
    case 45:
    case 46:
        return true;
    default:
        throw std::runtime_error("unsupported SMASH interaction type: " + std::to_string(type));
    }
}

int Integer(const std::string& token) {
    int value = -1;
    const auto result = std::from_chars(token.data(), token.data() + token.size(), value);
    if (result.ec != std::errc{} || result.ptr != token.data() + token.size() || value < 0)
        throw std::runtime_error("invalid interaction header integer: " + token);
    return value;
}

} // namespace

SmashInteraction ReadSmashInteraction(std::istream& fields) {
    SmashInteraction interaction;
    std::string key, value;
    if (!(fields >> key >> value) || key != "in")
        throw std::runtime_error("missing interaction incoming count");
    interaction.incoming = Integer(value);
    if (!(fields >> key >> value) || key != "out")
        throw std::runtime_error("missing interaction outgoing count");
    interaction.outgoing = Integer(value);
    constexpr int maximum = std::numeric_limits<unsigned short>::max();
    if (interaction.incoming > maximum || interaction.outgoing > maximum ||
        interaction.incoming + interaction.outgoing == 0)
        throw std::runtime_error("invalid interaction particle count");

    bool has_type = false;
    while (fields >> key) {
        if (!(fields >> value))
            throw std::runtime_error("missing interaction header value: " + key);
        if (key == "type") {
            if (has_type)
                throw std::runtime_error("duplicate interaction type");
            interaction.type = Integer(value);
            has_type = true;
        }
    }
    if (!has_type)
        throw std::runtime_error("missing interaction type");
    interaction.inelastic = IsInelastic(interaction.type);
    if (interaction.type == 47 && (interaction.incoming != 2 || interaction.outgoing != 2))
        throw std::runtime_error("failed string interaction must have two inputs and outputs");
    return interaction;
}

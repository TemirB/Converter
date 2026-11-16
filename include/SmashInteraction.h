#pragma once

#include <iosfwd>

struct SmashInteraction {
    int incoming = 0;
    int outgoing = 0;
    int type = 0;
    bool inelastic = false;
};

// Reads the part after "# interaction" in an OSCAR2013 collision header.
// Unknown process types and incomplete/invalid headers are errors.
SmashInteraction ReadSmashInteraction(std::istream& fields);

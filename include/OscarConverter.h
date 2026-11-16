#pragma once

#include <string>

enum class OscarVersion { Legacy = 0, Optimized = 1, InteractionTypes = 2 };

struct ConversionOptions {
    // Also retain events whose final particle count is exactly 394.
    bool keep_empty_events = false;
    // Validation only: retain OSCAR IDs alongside the McParticle array.
    bool record_source_ids = false;
    // 0: frozen legacy, 1: optimized legacy semantics, 2: inelastic-process spectators.
    OscarVersion version = OscarVersion::Optimized;
};

class OscarConverter {
  public:
    bool Convert(const std::string& inputFilename, const std::string& outputFilename,
                 const ConversionOptions& options = {});
};

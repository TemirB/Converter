#pragma once

#include <string>

struct ConversionOptions {
    // Also retain events whose final particle count is exactly 394.
    bool keep_empty_events = false;
    // Validation only: retain OSCAR IDs alongside the McParticle array.
    bool record_source_ids = false;
};

class OscarConverter {
  public:
    bool Convert(const std::string& inputFilename, const std::string& outputFilename,
                 const ConversionOptions& options = {});
};

#pragma once

#include <string>

struct ConversionOptions {
    // Also retain events whose final particle count is exactly 394.
    bool keep_empty_events = false;
};

class OscarConverter {
  public:
    bool Convert(const std::string& inputFilename, const std::string& outputFilename,
                 const ConversionOptions& options = {});
};

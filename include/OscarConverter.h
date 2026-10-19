#pragma once

#include <string>

struct ConversionOptions {
    bool keep_empty_events = false;
};

class OscarConverter {
  public:
    bool Convert(const std::string& inputFilename, const std::string& outputFilename,
                 const ConversionOptions& options = {});
};

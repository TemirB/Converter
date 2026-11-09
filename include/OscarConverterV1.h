#pragma once

#include <string>

// Frozen legacy algorithm supplied for comparison with OscarConverter.
class OscarConverterV1 {
  public:
    bool Convert(const std::string& inputFilename, const std::string& outputFilename,
                 bool recordSourceIds = false);
};

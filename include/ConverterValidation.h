#pragma once

#include <string>

enum class ValidationResult { Equal = 0, Error = 2, Different = 3 };

// Both McDst inputs must include the validation-only SourceId branch.
// Writes a CSV report and a ROOT histogram file with the same stem.
ValidationResult CompareMcDstVersions(const std::string& currentFilename,
                                      const std::string& v1Filename,
                                      const std::string& reportFilename);
ValidationResult ValidateOscarVersions(const std::string& inputFilename,
                                       const std::string& outputDirectory = "");

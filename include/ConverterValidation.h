#pragma once

#include <string>

enum class ValidationResult { Equal = 0, Error = 2, Different = 3 };

// Both McDst inputs must include the validation-only SourceId branch.
// Writes a CSV report and a ROOT histogram file with the same stem.
ValidationResult CompareMcDstVersions(const std::string& currentFilename,
                                      const std::string& v1Filename,
                                      const std::string& reportFilename,
                                      const std::string& currentLabel = "current",
                                      const std::string& referenceLabel = "v1");
ValidationResult ValidateOscarVersions(const std::string& inputFilename,
                                       const std::string& outputDirectory = "");
// Runs versions 0/1/2, writes versionN.McDst.root and all three pairwise diffs.
ValidationResult ValidateAllOscarVersions(const std::string& inputFilename,
                                          const std::string& outputDirectory = "");

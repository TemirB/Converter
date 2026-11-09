#pragma once

#include "OscarReader.h"

#include <memory>
#include <string>

// Writes beside the destination, then atomically replaces it after success.
class McDstWriter {
  public:
    McDstWriter(const std::string& filename, const std::string& definitions);
    ~McDstWriter();
    McDstWriter(const McDstWriter&) = delete;
    McDstWriter& operator=(const McDstWriter&) = delete;

    bool WriteEvent(const OscarEvent& event);
    bool Finish(long long input_events, long long skipped_events);
    const std::string& Error() const;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

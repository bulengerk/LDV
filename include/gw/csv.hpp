#pragma once
#include <string>
#include <vector>

namespace gw {

class CsvWriter {
public:
    void write(const std::string& path, double dt,
               const std::vector<double>& raw,
               const std::vector<double>& filt) const;
};

} // namespace gw

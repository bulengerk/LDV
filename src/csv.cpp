#include "gw/csv.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>

namespace gw {

void CsvWriter::write(const std::string& path, double dt,
                      const std::vector<double>& raw,
                      const std::vector<double>& filt) const
{
    std::ofstream f(path);
    f << "t,raw,filtered\n";
    const size_t n = std::min(raw.size(), filt.size());
    for (size_t i=0; i<n; ++i) {
        f << (i*dt) << "," << raw[i] << "," << filt[i] << "\n";
    }
    std::cout << "✓ CSV: " << path << "\n";
}

} // namespace gw

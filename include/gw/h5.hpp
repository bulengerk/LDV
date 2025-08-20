#pragma once
#include <string>
#include <vector>

namespace gw {

struct StrainSeries {
    std::vector<double> samples;
    double dt = 1.0 / 4096.0;
    double t0 = 0.0;
    std::string datasetPath;
};

class H5StrainReader {
public:
    StrainSeries read(const std::string& h5file) const;
};

} // namespace gw

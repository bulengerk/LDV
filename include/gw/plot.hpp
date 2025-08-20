#pragma once
#include <string>
#include <vector>

namespace gw {

class SvgPlotter {
public:
    void lineplot(const std::string& path,
                  double width, double height,
                  double t0, double dt,
                  const std::vector<double>& y) const;
};

} // namespace gw

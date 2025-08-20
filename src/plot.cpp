#include "gw/plot.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>

namespace gw {

void SvgPlotter::lineplot(const std::string& path,
                          double width, double height,
                          double t0, double dt,
                          const std::vector<double>& y) const
{
    if (y.empty()) return;

    const size_t N = y.size();
    const size_t maxp = 10000;
    const size_t step = std::max<size_t>(1, N / maxp);

    const auto [minIt, maxIt] = std::minmax_element(y.begin(), y.end());
    double ymin = *minIt, ymax = *maxIt;
    if (ymax <= ymin) ymax = ymin + 1.0;

    const double left=60, right=20, top=20, bottom=40;
    const double W = width, H = height;
    const double plotW = W - left - right;
    const double plotH = H - top - bottom;

    std::ofstream s(path);
    s << "<svg xmlns='http://www.w3.org/2000/svg' width='"<<W<<"' height='"<<H<<"' viewBox='0 0 "<<W<<" "<<H<<"'>\n";
    s << "<rect x='0' y='0' width='"<<W<<"' height='"<<H<<"' fill='white' stroke='none'/>\n";
    s << "<line x1='"<<left<<"' y1='"<<top<<"' x2='"<<left<<"' y2='"<<(H-bottom)
      <<"' stroke='black' stroke-width='1'/>\n";
    s << "<line x1='"<<left<<"' y1='"<<(H-bottom)<<"' x2='"<<(W-right)<<"' y2='"<<(H-bottom)
      <<"' stroke='black' stroke-width='1'/>\n";

    s << "<polyline fill='none' stroke='black' stroke-width='1' points='";
    for (size_t i=0; i<N; i+=step) {
        const double t = t0 + i*dt;
        const double x = left + (t - t0) / (dt*(N-1)) * plotW;
        const double yv = y[i];
        const double ymap = top + (ymax - yv) / (ymax - ymin) * plotH;
        s << x << "," << ymap << " ";
    }
    s << "'/>\n";
    s << "</svg>\n";

    std::cout << "✓ SVG: " << path << "\n";
}

} // namespace gw

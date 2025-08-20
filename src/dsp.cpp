#include "gw/dsp.hpp"
#include <cmath>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace gw::dsp {

double Biquad::process(double x) {
    const double y = b0*x + z1;
    z1 = b1*x - a1*y + z2;
    z2 = b2*x - a2*y;
    return y;
}

void Biquad::reset() { z1 = 0; z2 = 0; }

static inline Biquad make_peaking(double b0, double b1, double b2, double a0, double a1, double a2) {
    Biquad bi;
    bi.b0 = b0 / a0; bi.b1 = b1 / a0; bi.b2 = b2 / a0;
    bi.a1 = a1 / a0; bi.a2 = a2 / a0;
    return bi;
}

Biquad design_lpf(double fs, double fc, double Q) {
    const double w0 = 2.0 * M_PI * fc / fs;
    const double c = std::cos(w0), s = std::sin(w0);
    const double alpha = s / (2.0 * Q);
    const double b0 = (1.0 - c) / 2.0;
    const double b1 = 1.0 - c;
    const double b2 = (1.0 - c) / 2.0;
    const double a0 = 1.0 + alpha;
    const double a1 = -2.0 * c;
    const double a2 = 1.0 - alpha;
    return make_peaking(b0, b1, b2, a0, a1, a2);
}

Biquad design_hpf(double fs, double fc, double Q) {
    const double w0 = 2.0 * M_PI * fc / fs;
    const double c = std::cos(w0), s = std::sin(w0);
    const double alpha = s / (2.0 * Q);
    const double b0 =  (1.0 + c) / 2.0;
    const double b1 = -(1.0 + c);
    const double b2 =  (1.0 + c) / 2.0;
    const double a0 = 1.0 + alpha;
    const double a1 = -2.0 * c;
    const double a2 = 1.0 - alpha;
    return make_peaking(b0, b1, b2, a0, a1, a2);
}

void apply_cascade(std::vector<double>& x, const std::vector<Biquad>& sos) {
    for (auto& bi : sos) {
        for (auto& v : x) v = const_cast<Biquad&>(bi).process(v);
    }
}

void normalize_std(std::vector<double>& x) {
    double mean = 0.0, M2 = 0.0;
    size_t n = 0;
    for (double v : x) { ++n; double d=v-mean; mean+=d/n; M2+=d*(v-mean); }
    const double var = (n>1)? M2/(n-1): 1.0;
    double s = (var>0)? std::sqrt(var):1.0;
    if (s==0) s=1.0;
    for (auto& v : x) v /= s;
}

void bandpass_and_whiten(std::vector<double>& x, double fs, double hp, double lp) {
    std::vector<Biquad> sos;
    sos.push_back(design_hpf(fs, hp, 0.7071));
    sos.push_back(design_hpf(fs, hp, 0.7071));
    sos.push_back(design_lpf(fs, lp, 0.7071));
    sos.push_back(design_lpf(fs, lp, 0.7071));
    apply_cascade(x, sos);
    normalize_std(x);
}

} // namespace gw::dsp

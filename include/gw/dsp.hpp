#pragma once
#include <vector>

namespace gw::dsp {

class Biquad {
public:
    Biquad() = default;
    double process(double x);
    void reset();
    double b0{1}, b1{0}, b2{0}, a1{0}, a2{0};
private:
    double z1{0}, z2{0};
};

Biquad design_lpf(double fs, double fc, double Q = 0.70710678);
Biquad design_hpf(double fs, double fc, double Q = 0.70710678);
void apply_cascade(std::vector<double>& x, const std::vector<Biquad>& sos);
void normalize_std(std::vector<double>& x);
void bandpass_and_whiten(std::vector<double>& x, double fs, double hp, double lp);

} // namespace gw::dsp

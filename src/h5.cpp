#include "gw/h5.hpp"

#include <H5Cpp.h>
#include <hdf5.h>
#include <H5Epublic.h>
#include <H5Lpublic.h>

#include <algorithm>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace gw {

static bool dataset_exists(H5::H5File& f, const std::string& path) {
    H5E_auto2_t old_func;
    void* old_client;
    H5Eget_auto2(H5E_DEFAULT, &old_func, &old_client);
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
    const htri_t ex = H5Lexists(f.getId(), path.c_str(), H5P_DEFAULT);
    H5Eset_auto2(H5E_DEFAULT, old_func, old_client);
    return ex > 0;
}

static std::optional<double> try_attr_double(H5::DataSet& ds, const char* name) {
    if (!H5Aexists(ds.getId(), name)) return std::nullopt;
    auto a = ds.openAttribute(name);
    if (a.getTypeClass() != H5T_FLOAT) return std::nullopt;
    H5::FloatType t = a.getFloatType();
    const size_t sz = t.getSize();
    if (sz == sizeof(double)) { double v=0; a.read(H5::PredType::NATIVE_DOUBLE, &v); return v; }
    else { float v=0; a.read(H5::PredType::NATIVE_FLOAT, &v); return static_cast<double>(v); }
}

StrainSeries H5StrainReader::read(const std::string& h5file) const {
    H5::Exception::dontPrint();
    H5::H5File f(h5file, H5F_ACC_RDONLY);

    const std::vector<std::string> candidates = {
        "/strain/Strain", "strain/Strain",
        "/strain", "strain",
        "/Strain", "Strain",
        "H1:STRAIN", "L1:STRAIN", "V1:STRAIN"
    };

    H5::DataSet ds;
    std::string used;
    bool found = false;
    for (const auto& c : candidates) {
        if (dataset_exists(f, c)) {
            ds = f.openDataSet(c);
            used = c;
            found = true;
            break;
        }
    }
    if (!found) throw std::runtime_error("no strain dataset found in: " + h5file);

    H5::DataSpace sp = ds.getSpace();
    if (sp.getSimpleExtentNdims() != 1) {
        throw std::runtime_error("strain dataset is not 1D");
    }
    hsize_t npts = 0;
    sp.getSimpleExtentDims(&npts, nullptr);

    StrainSeries out;
    out.samples.resize(static_cast<size_t>(npts));
    ds.read(out.samples.data(), H5::PredType::NATIVE_DOUBLE);
    out.datasetPath = used;

    out.dt = 1.0 / 4096.0;
    out.t0 = 0.0;
    if (auto a = try_attr_double(ds, "Xspacing"); a && *a > 0)           out.dt = *a;
    if (auto a = try_attr_double(ds, "SampleRate"); a && *a > 1e-9)      out.dt = 1.0 / *a;
    if (auto a = try_attr_double(ds, "SamplingRate"); a && *a > 1e-9)    out.dt = 1.0 / *a;
    if (auto a = try_attr_double(ds, "Xstart"); a)                        out.t0 = *a;

    return out;
}

} // namespace gw

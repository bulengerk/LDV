#include <H5Cpp.h>
#include <curl/curl.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <vector>
#include <sys/stat.h>

using std::string;
using std::vector;

static bool is_http_url(const string& s) {
    return s.rfind("http://", 0) == 0 || s.rfind("https://", 0) == 0;
}

static size_t curl_writefile(void* ptr, size_t size, size_t nmemb, void* stream) {
    FILE* fp = static_cast<FILE*>(stream);
    return fwrite(ptr, size, nmemb, fp) * size;
}

static void ensure_dir(const string& path) {
#ifdef _WIN32
    _mkdir(path.c_str());
#else
    ::mkdir(path.c_str(), 0755);
#endif
}

static string download_to(const string& url, const string& outdir) {
    ensure_dir(outdir);
    // Spróbuj zachować nazwę pliku z URL-a.
    auto slash = url.find_last_of('/');
    string fname = (slash == string::npos) ? "data.hdf5" : url.substr(slash + 1);
    if (fname.find('.') == string::npos) fname += ".hdf5";
    string outpath = outdir + "/" + fname;

    CURL* curl = curl_easy_init();
    if (!curl) throw std::runtime_error("curl init failed");

    FILE* fp = std::fopen(outpath.c_str(), "wb");
    if (!fp) {
        curl_easy_cleanup(curl);
        throw std::runtime_error("cannot open output file: " + outpath);
    }

    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, curl_writefile);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, fp);
    curl_easy_setopt(curl, CURLOPT_USERAGENT, "gwosc-plotter/1.0");

    CURLcode res = curl_easy_perform(curl);
    std::fclose(fp);
    long code = 0;
    curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &code);
    curl_easy_cleanup(curl);

    if (res != CURLE_OK || code >= 400) {
        std::remove(outpath.c_str());
        std::ostringstream oss;
        oss << "download failed (HTTP " << code << "): " << curl_easy_strerror(res);
        throw std::runtime_error(oss.str());
    }

    std::cout << "✓ Downloaded: " << outpath << "\n";
    return outpath;
}

struct H5Data {
    vector<double> strain;
    double dt = 1.0 / 4096.0;   // domyślnie 4096 Hz
    double t0 = 0.0;            // opcjonalny start (UTC lub GPS), tu tylko informacyjnie
    string ds_path;
};

static bool dataset_exists(H5::H5File& f, const string& path) {
    H5E_auto2_t old_func;
    void* old_client;
    H5Eget_auto2(H5E_DEFAULT, &old_func, &old_client);
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr); // wycisz
    bool ok = H5Lexists(f.getId(), path.c_str(), H5P_DEFAULT) > 0;
    H5Eset_auto2(H5E_DEFAULT, old_func, old_client);
    return ok;
}

static double read_attr_double(const H5::Attribute& a) {
    // spróbuj jako double, jeśli nie – jako float
    H5T_class_t cls = a.getTypeClass();
    if (cls == H5T_FLOAT) {
        H5::FloatType t = a.getFloatType();
        size_t sz = t.getSize();
        if (sz == sizeof(double)) {
            double v = 0.0;
            a.read(H5::PredType::NATIVE_DOUBLE, &v);
            return v;
        } else {
            float v = 0.f;
            a.read(H5::PredType::NATIVE_FLOAT, &v);
            return static_cast<double>(v);
        }
    }
    // nie-liczbowe: rzuć
    throw std::runtime_error("attribute is not float");
}

static std::optional<double> try_read_attr_double(H5::DataSet& ds, const char* name) {
    if (!H5Aexists(ds.getId(), name)) return std::nullopt;
    auto a = ds.openAttribute(name);
    try { return read_attr_double(a); }
    catch (...) { return std::nullopt; }
}

static H5Data load_h5_strain(const string& path) {
    H5::Exception::dontPrint();
    H5::H5File f(path, H5F_ACC_RDONLY);

    // typowe ścieżki w GWOSC:
    const vector<string> candidates = {
        "/strain/Strain", "strain/Strain",
        "/strain", "strain",
        "/Strain", "Strain",
        "H1:STRAIN", "L1:STRAIN", "V1:STRAIN"
    };

    H5::DataSet ds;
    string used;
    bool found = false;
    for (const auto& c : candidates) {
        if (dataset_exists(f, c)) {
            ds = f.openDataSet(c);
            used = c;
            found = true;
            break;
        }
    }
    if (!found) throw std::runtime_error("no strain dataset found");

    H5::DataSpace sp = ds.getSpace();
    if (sp.getSimpleExtentNdims() != 1) throw std::runtime_error("dataset not 1D");
    hsize_t npts = 0;
    sp.getSimpleExtentDims(&npts, nullptr);

    vector<double> buf;
    buf.resize(static_cast<size_t>(npts));
    // czytaj jako double niezależnie od typu bazowego
    ds.read(buf.data(), H5::PredType::NATIVE_DOUBLE);

    H5Data out;
    out.strain = std::move(buf);
    out.ds_path = used;

    // atrybuty czasu: Xspacing ~ dt, Xstart ~ t0 (GPS/UTC zależnie od pliku)
    if (H5Aexists(ds.getId(), "Xspacing")) {
        try {
            double dx = *try_read_attr_double(ds, "Xspacing");
            if (dx > 0) out.dt = dx;
        } catch (...) {}
    }
    if (H5Aexists(ds.getId(), "SampleRate")) {
        try {
            double fs = *try_read_attr_double(ds, "SampleRate");
            if (fs > 1e-9) out.dt = 1.0 / fs;
        } catch (...) {}
    }
    if (H5Aexists(ds.getId(), "SamplingRate")) {
        try {
            double fs = *try_read_attr_double(ds, "SamplingRate");
            if (fs > 1e-9) out.dt = 1.0 / fs;
        } catch (...) {}
    }
    if (H5Aexists(ds.getId(), "Xstart")) {
        try { out.t0 = *try_read_attr_double(ds, "Xstart"); } catch (...) {}
    }
    return out;
}

// Prosty biquad (RBJ cookbook). Przetwarzanie w miejscu (in/out).
struct Biquad {
    double b0=1, b1=0, b2=0, a1=0, a2=0;
    double z1=0, z2=0;
    inline double process(double x) {
        double y = b0*x + z1;
        z1 = b1*x - a1*y + z2;
        z2 = b2*x - a2*y;
        return y;
    }
};

static Biquad design_lpf(double fs, double fc, double Q=0.70710678) {
    double w0 = 2.0 * M_PI * fc / fs;
    double c = std::cos(w0), s = std::sin(w0);
    double alpha = s / (2.0 * Q);
    double b0 = (1.0 - c) / 2.0;
    double b1 = 1.0 - c;
    double b2 = (1.0 - c) / 2.0;
    double a0 = 1.0 + alpha;
    double a1 = -2.0 * c;
    double a2 = 1.0 - alpha;
    Biquad bi;
    bi.b0 = b0 / a0; bi.b1 = b1 / a0; bi.b2 = b2 / a0;
    bi.a1 = a1 / a0; bi.a2 = a2 / a0;
    return bi;
}

static Biquad design_hpf(double fs, double fc, double Q=0.70710678) {
    double w0 = 2.0 * M_PI * fc / fs;
    double c = std::cos(w0), s = std::sin(w0);
    double alpha = s / (2.0 * Q);
    double b0 =  (1.0 + c) / 2.0;
    double b1 = -(1.0 + c);
    double b2 =  (1.0 + c) / 2.0;
    double a0 = 1.0 + alpha;
    double a1 = -2.0 * c;
    double a2 = 1.0 - alpha;
    Biquad bi;
    bi.b0 = b0 / a0; bi.b1 = b1 / a0; bi.b2 = b2 / a0;
    bi.a1 = a1 / a0; bi.a2 = a2 / a0;
    return bi;
}

static void apply_cascade(vector<double>& x, const vector<Biquad>& sos) {
    for (auto& bi : sos) {
        for (auto& v : x) v = const_cast<Biquad&>(bi).process(v);
    }
}

static void normalize_std(vector<double>& x) {
    double mean = 0.0, M2 = 0.0;
    size_t n = 0;
    for (double v : x) { n++; double d=v-mean; mean+=d/n; M2+=d*(v-mean); }
    double var = (n>1)? M2/(n-1): 1.0;
    double s = (var>0)? std::sqrt(var):1.0;
    if (s==0) s=1.0;
    for (auto& v : x) v /= s;
}

static void write_csv(const string& path, double dt, const vector<double>& raw, const vector<double>& filt) {
    std::ofstream f(path);
    f << "t,raw,filtered\n";
    size_t n = std::min(raw.size(), filt.size());
    for (size_t i=0;i<n;i++) {
        f << (i*dt) << "," << raw[i] << "," << filt[i] << "\n";
    }
    std::cout << "✓ CSV: " << path << "\n";
}

static void write_svg_lineplot(const string& path, double width, double height,
                               double t0, double dt, const vector<double>& y) {
    if (y.empty()) return;
    // decymacja do max ~10k punktów
    size_t N = y.size();
    size_t maxp = 10000;
    size_t step = std::max<size_t>(1, N / maxp);

    double ymin = *std::min_element(y.begin(), y.end());
    double ymax = *std::max_element(y.begin(), y.end());
    if (ymax <= ymin) { ymax = ymin + 1.0; }

    double left=60, right=20, top=20, bottom=40;
    double W = width, H = height;
    double plotW = W - left - right;
    double plotH = H - top - bottom;

    std::ofstream s(path);
    s << "<svg xmlns='http://www.w3.org/2000/svg' width='"<<W<<"' height='"<<H<<"' viewBox='0 0 "<<W<<" "<<H<<"'>\n";
    s << "<rect x='0' y='0' width='"<<W<<"' height='"<<H<<"' fill='white' stroke='none'/>\n";
    // osie
    s << "<line x1='"<<left<<"' y1='"<<top<<"' x2='"<<left<<"' y2='"<<(H-bottom)
      <<"' stroke='black' stroke-width='1'/>\n";
    s << "<line x1='"<<left<<"' y1='"<<(H-bottom)<<"' x2='"<<(W-right)<<"' y2='"<<(H-bottom)
      <<"' stroke='black' stroke-width='1'/>\n";

    // etykiety min/max
    s << "<text x='"<<5<<"' y='"<<top+10<<"' font-size='12'>"
      << "max " << ymax << "</text>\n";
    s << "<text x='"<<5<<"' y='"<<(H-bottom)<<"' font-size='12'>"
      << "min " << ymin << "</text>\n";

    // polyline
    s << "<polyline fill='none' stroke='black' stroke-width='1' points='";
    for (size_t i=0;i<N;i+=step) {
        double t = t0 + i*dt;
        double x = left + (t - t0) / (dt*(N-1)) * plotW;
        double yv = y[i];
        double ymap = top + (ymax - yv) / (ymax - ymin) * plotH;
        s << x << "," << ymap << " ";
    }
    s << "'/>\n";
    s << "</svg>\n";
    std::cout << "✓ SVG: " << path << "\n";
}

struct Args {
    string source;     // URL lub ścieżka pliku .hdf5
    string outdir = "out";
    double hp = 20.0;
    double lp = 300.0;
    bool no_download = false;
};

static void usage(const char* prog) {
    std::cerr <<
        "Usage:\n"
        "  " << prog << " <url_or_hdf5_path> [--out outdir] [--hp 20] [--lp 300] [--no-download]\n\n"
        "Examples:\n"
        "  " << prog << " https://www.gw-openscience.org/archive/links/GW150914/v3/H-H1_GWOSC_4KHZ_R1-1126259446-32.hdf5\n"
        "  " << prog << " H-H1_GWOSC_4KHZ_R1-1126259446-32.hdf5 --hp 30 --lp 350\n";
}

static Args parse_args(int argc, char** argv) {
    if (argc < 2) { usage(argv[0]); std::exit(1); }
    Args a;
    a.source = argv[1];
    for (int i=2;i<argc;i++) {
        string k = argv[i];
        if (k=="--out" && i+1<argc) { a.outdir = argv[++i]; }
        else if (k=="--hp" && i+1<argc) { a.hp = std::stod(argv[++i]); }
        else if (k=="--lp" && i+1<argc) { a.lp = std::stod(argv[++i]); }
        else if (k=="--no-download") { a.no_download = true; }
        else { std::cerr << "Unknown arg: " << k << "\n"; usage(argv[0]); std::exit(1); }
    }
    return a;
}

int main(int argc, char** argv) {
    try {
        auto args = parse_args(argc, argv);
        string h5path = args.source;

        if (is_http_url(args.source)) {
            if (args.no_download) {
                std::cerr << "URL provided with --no-download. Nothing to do.\n";
                return 1;
            }
            curl_global_init(CURL_GLOBAL_DEFAULT);
            h5path = download_to(args.source, args.outdir);
            curl_global_cleanup();
        }

        auto data = load_h5_strain(h5path);
        const size_t N = data.strain.size();
        const double fs = 1.0 / data.dt;
        std::cout << "Dataset: " << data.ds_path << "\n";
        std::cout << "Samples: " << N << ", fs ≈ " << fs << " Hz, dt=" << data.dt << " s\n";

        // Kopie: raw (do wykresu) i filtrowany
        vector<double> raw = data.strain;
        vector<double> y = data.strain;

        // Prosty bandpass: (HPF fc=hp) -> (LPF fc=lp), kaskada 2x żeby uzyskać ~4 rząd.
        vector<Biquad> sos;
        sos.push_back(design_hpf(fs, args.hp, 0.7071));
        sos.push_back(design_hpf(fs, args.hp, 0.7071));
        sos.push_back(design_lpf(fs, args.lp, 0.7071));
        sos.push_back(design_lpf(fs, args.lp, 0.7071));
        apply_cascade(y, sos);

        // Normalizacja (z grubsza "whitening")
        normalize_std(y);

        ensure_dir(args.outdir);
        write_csv(args.outdir + "/trace.csv", data.dt, raw, y);
        write_svg_lineplot(args.outdir + "/raw.svg", 1200, 400, 0.0, data.dt, raw);
        write_svg_lineplot(args.outdir + "/filtered.svg", 1200, 400, 0.0, data.dt, y);

        std::cout << "Done. See " << args.outdir << "/raw.svg and /filtered.svg\n";
        return 0;
    } catch (const H5::Exception& e) {
        std::cerr << "HDF5 error: " << e.getDetailMsg() << "\n";
        return 2;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 3;
    }
}


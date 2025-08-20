#include "gw/app.hpp"
#include "gw/util.hpp"
#include "gw/http.hpp"
#include "gw/h5.hpp"
#include "gw/dsp.hpp"
#include "gw/plot.hpp"
#include "gw/csv.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstdlib>

namespace gw {

void App::usage(const char* prog) {
    std::cerr <<
        "Usage:\n"
        "  " << prog << " <url_or_hdf5_path> [--out outdir] [--hp 20] [--lp 300] [--no-download]\n\n"
        "Examples:\n"
        "  " << prog << " https://<GWOSC_FILE>.hdf5\n"
        "  " << prog << " file.hdf5 --hp 30 --lp 350\n";
}

Args App::parse_args(int argc, char** argv) {
    if (argc < 2) { usage(argv[0]); std::exit(1); }
    Args a;
    a.source = argv[1];
    for (int i=2; i<argc; ++i) {
        std::string k = argv[i];
        if (k=="--out" && i+1<argc) { a.outdir = argv[++i]; }
        else if (k=="--hp" && i+1<argc) { a.hp = std::stod(argv[++i]); }
        else if (k=="--lp" && i+1<argc) { a.lp = std::stod(argv[++i]); }
        else if (k=="--no-download") { a.no_download = true; }
        else { std::cerr << "Unknown arg: " << k << "\n"; usage(argv[0]); std::exit(1); }
    }
    return a;
}

int App::run(int argc, char** argv) {
    try {
        Args args = parse_args(argc, argv);
        std::string h5path = args.source;

        if (is_http_url(args.source)) {
            if (args.no_download) {
                std::cerr << "URL provided with --no-download. Nothing to do.\n";
                return 1;
            }
            HttpClient http;
            h5path = http.download_to(args.source, args.outdir);
        }

        H5StrainReader reader;
        StrainSeries ss = reader.read(h5path);

        const double fs = 1.0 / ss.dt;
        std::cout << "Dataset: " << ss.datasetPath << "\n"
                  << "Samples: " << ss.samples.size()
                  << ", fs ≈ " << fs << " Hz, dt=" << ss.dt << " s\n";

        ensure_dir(args.outdir);

        std::vector<double> raw = ss.samples;
        std::vector<double> filtered = ss.samples;

        dsp::bandpass_and_whiten(filtered, fs, args.hp, args.lp);

        CsvWriter csv;
        csv.write(args.outdir + "/trace.csv", ss.dt, raw, filtered);

        SvgPlotter plot;
        plot.lineplot(args.outdir + "/raw.svg",      1200, 400, 0.0, ss.dt, raw);
        plot.lineplot(args.outdir + "/filtered.svg", 1200, 400, 0.0, ss.dt, filtered);

        std::cout << "Done. See " << args.outdir << "/raw.svg and /filtered.svg\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 2;
    }
}

} // namespace gw

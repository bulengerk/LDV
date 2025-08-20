#pragma once
#include <string>

namespace gw {

struct Args {
    std::string source;
    std::string outdir = "out";
    double hp = 20.0;
    double lp = 300.0;
    bool no_download = false;
};

class App {
public:
    int run(int argc, char** argv);
private:
    static void usage(const char* prog);
    static Args parse_args(int argc, char** argv);
};

} // namespace gw

#include "gw/util.hpp"
#include <sys/stat.h>
#ifdef _WIN32
  #include <direct.h>
#endif
#include <string>

namespace gw {

void ensure_dir(const std::string& path) {
#ifdef _WIN32
    _mkdir(path.c_str());
#else
    ::mkdir(path.c_str(), 0755);
#endif
}

bool is_http_url(const std::string& s) {
    return s.rfind("http://", 0) == 0 || s.rfind("https://", 0) == 0;
}

} // namespace gw

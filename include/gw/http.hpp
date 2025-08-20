#pragma once
#include <string>

namespace gw {

class HttpClient {
public:
    // Download URL to outdir; returns local filepath. Throws on HTTP error.
    std::string download_to(const std::string& url, const std::string& outdir) const;
};

} // namespace gw

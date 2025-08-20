#include "gw/http.hpp"
#include "gw/util.hpp"

#include <curl/curl.h>
#include <cstdio>
#include <stdexcept>
#include <string>

namespace gw {

static size_t curl_writefile(void* ptr, size_t size, size_t nmemb, void* stream) {
    FILE* fp = static_cast<FILE*>(stream);
    return fwrite(ptr, size, nmemb, fp) * size;
}

std::string HttpClient::download_to(const std::string& url, const std::string& outdir) const {
    ensure_dir(outdir);
    const auto slash = url.find_last_of('/');
    std::string fname = (slash == std::string::npos) ? "data.hdf5" : url.substr(slash + 1);
    if (fname.find('.') == std::string::npos) fname += ".hdf5";
    std::string outpath = outdir + "/" + fname;

    curl_global_init(CURL_GLOBAL_DEFAULT);
    CURL* curl = curl_easy_init();
    if (!curl) throw std::runtime_error("curl init failed");

    FILE* fp = std::fopen(outpath.c_str(), "wb");
    if (!fp) {
        curl_easy_cleanup(curl);
        curl_global_cleanup();
        throw std::runtime_error("cannot open output file: " + outpath);
    }

    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, curl_writefile);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, fp);
    curl_easy_setopt(curl, CURLOPT_USERAGENT, "gwosc-plotter/3.0");

    const CURLcode res = curl_easy_perform(curl);
    long code = 0;
    curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &code);
    std::fclose(fp);
    curl_easy_cleanup(curl);
    curl_global_cleanup();

    if (res != CURLE_OK || code >= 400) {
        std::remove(outpath.c_str());
        throw std::runtime_error("download failed (HTTP " + std::to_string(code) + "): " +
                                 curl_easy_strerror(res));
    }

    std::printf("✓ Downloaded: %s\n", outpath.c_str());
    return outpath;
}

} // namespace gw

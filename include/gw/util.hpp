#pragma once
#include <string>

namespace gw {

void ensure_dir(const std::string& path);
bool is_http_url(const std::string& s);

} // namespace gw

#pragma once

#include <filesystem>

namespace zone {

namespace fs = std::filesystem;
static fs::path current_path = fs::current_path();

}
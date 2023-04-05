#pragma once

#include <filesystem>
#include <iostream>

namespace zone {

namespace fs = std::filesystem;
static fs::path current_path = fs::current_path();

}
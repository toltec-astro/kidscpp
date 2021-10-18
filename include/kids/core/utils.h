#pragma once

#include "constants.h"
#include <tula/grppi.h>

namespace kids {

namespace utils {

template <typename Config> auto exmode(const Config &config) {
    namespace ck = kids::config_keys;
    return tula::grppi_utils::dyn_ex(config.get_str(ck::exmode));
}

} // namespace utils

} // namespace kids

#include "kids/core/io.h"

using kids::KidsDataIOError;

auto KidsDataIOError::NotImplemented(std::string_view name) -> KidsDataIOError {
    return KidsDataIOError{fmt::format("{} not implemented", name)};
}

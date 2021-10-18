#pragma once

#include "kidsdata.h"
#include <tula/logging.h>
#include <tula/meta.h>

namespace kids {

struct KidsDataIOError : public std::runtime_error {
    using std::runtime_error::runtime_error;
    static auto NotImplemented(std::string_view name) -> KidsDataIOError;
};

/// @brief Kids data IO base class.
template <typename Derived, typename IO>
struct KidsDataIO {

    static_assert(
        requires { IO::read; } || requires { IO::write; },
        "IO CLASS DO NOT HAVE R/W METHODS");

    static auto from_filepath(const std::string &filepath) {
        return IO(filepath);
    }

    template <typename... Args>
    static auto read(Args... args) {
        if constexpr (requires { IO::read; }) {
            return IO::read(std::forward<decltype(args)>(args)...);
        } else {
            throw KidsDataIOError::NotImplemented(
                fmt::format("{}::{}", Derived::format, "read"));
        }
    }
    template <typename... Args>
    static auto read_slice(Args... args) {
        if constexpr (requires { IO::read_slice; }) {
            return IO::read_slice(std::forward<decltype(args)>(args)...);
        } else {
            throw KidsDataIOError::NotImplemented(
                fmt::format("{}::{}", Derived::format, "read_slice"));
        }
    }

    template <typename... Args>
    static auto write(Args... args) {
        if constexpr (requires { IO::write; }) {
            return IO::write(FWD(args)...);
        } else {
            throw KidsDataIOError::NotImplemented(
                fmt::format("{}::{}", Derived::format, "write"));
        }
    }
};

} // namespace kids

#pragma once

#include <kids/core/kidsdata.h>

#include <filesystem>
#include <stdexcept>
#include <string_view>
#include <tula/container.h>

namespace kids::toltec {

inline constexpr std::string_view raw_timestream_spec{
    "toltec.raw-timestream.1"};

using RawTimeStream = KidsData<KidsDataKind::RawTimeStream>;
using RawTimeStreamMeta = RawTimeStream::meta_t;
using SampleSlice = tula::container_utils::IndexSlice;

class RawTimeStreamIOError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

[[nodiscard]] auto
get_raw_timestream_meta(const std::filesystem::path &filepath)
    -> RawTimeStreamMeta;

[[nodiscard]] auto read_raw_timestream_slice(
    const std::filesystem::path &filepath,
    SampleSlice slice = {}) -> RawTimeStream;

} // namespace kids::toltec

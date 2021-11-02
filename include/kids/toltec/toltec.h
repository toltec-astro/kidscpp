#pragma once

#include "../core/io.h"
#include "../core/kidsdata.h"
#include <regex>
#include <tula/config/flatconfig.h>
#include <tula/container.h>
#include <tula/enum.h>
#include <tula/filesystem.h>
#include <tula/nc.h>
#include <tula/nddata/cacheddata.h>

// data spec toltec.1

/// @brief TolTEC data spec
namespace kids::toltec {

inline namespace v1 {

inline constexpr std::string_view name{"toltec.1"};

// clang-format off
/// @brief Supported toltec data io formats
TULA_ENUM(DataFormat, int,
           NcFile
);
// clang-format on

namespace internal {

using kids::internal::RMatrixXd;
using kids::internal::RMatrixXi;

using meta_t = KidsData<>::meta_t;
using tone_axis_t = kids::ToneAxis;
using sweep_axis_t = kids::SweepAxis;
using time_axis_t = kids::TimeAxis;
using sweep_data_t = Eigen::MatrixXcd;
using timestream_data_t = std::pair<RMatrixXd, RMatrixXd>;

template <DataFormat format_>
struct IO;

/// @brief IO handler for NetCDF file format.
template <>
struct IO<DataFormat::NcFile> {
    using Self = IO<DataFormat::NcFile>;
    using IndexSlice = tula::container_utils::IndexSlice;
    static auto read(const std::string &filepath,
                     KidsDataKind kind = KidsDataKind::Any) -> KidsData<>;
    static auto read_slice(const std::string &filepath, IndexSlice slice,
                           KidsDataKind kind = KidsDataKind::Any) -> KidsData<>;
    static auto write(const KidsData<> &data, const std::string &filepath)
        -> bool;

    using NcFile = netCDF::NcFile;
    IO(const std::string &filepath);
    ~IO();
    IO(const IO &other) = delete;
    IO(IO &&other) = delete;
    auto operator=(const IO &other) = delete;
    auto operator=(IO &&other) = delete;
    auto ncfile() const -> const NcFile & { return m_ncfile; }

private:
    std::string m_source{};
    NcFile m_ncfile{};
    using data_reducer_segments_t =
        std::vector<std::pair<Eigen::Index, Eigen::Index>>;
    struct data_reducer_segments_evaluator {
        static auto evaluate(const IO &) -> data_reducer_segments_t;
    };
    tula::nddata::CachedData<data_reducer_segments_t, &data_reducer_segments_evaluator::evaluate>
        m_data_reducer_segments{};

#define io_lazy_getter(name, type)                                             \
                                                                               \
public:                                                                        \
    auto name() const->const type & { return m_##name(*this); }                \
                                                                               \
private:                                                                       \
    struct name##_evaluator {                                                  \
        static auto evaluate(const IO &) -> type;                              \
    };                                                                         \
    tula::nddata::CachedData<type, &name##_evaluator::evaluate>      \
        m_##name{};

    // clang-format off
    io_lazy_getter(kind, KidsDataKind)
    io_lazy_getter(meta, meta_t)
    io_lazy_getter(tone_axis, tone_axis_t)
    io_lazy_getter(sweep_axis, sweep_axis_t)
    io_lazy_getter(time_axis, time_axis_t)
    io_lazy_getter(sweep_data, sweep_data_t)
    io_lazy_getter(timestream_data, timestream_data_t)
    // clang-format  on
#undef io_lazy_getter

public:
    // read the data bypassing the cache
    auto sweep_data_nocache() const -> sweep_data_t { return sweep_data_evaluator::evaluate(*this); }
    auto timestream_data_nocache() const -> timestream_data_t { return timestream_data_evaluator::evaluate(*this); }
    auto timestream_data_slice(IndexSlice slice) const -> timestream_data_t;
    auto time_axis_slice(IndexSlice slice) const -> time_axis_t;
};

} // namespace internal

/// @brief TolTEC Data IO Class
template <DataFormat format_>
struct KidsDataIO
    : kids::KidsDataIO<KidsDataIO<format_>, internal::IO<format_>> {
    static constexpr auto format = format_;
};

/// @brief Load TolTEC data from path.
template <KidsDataKind kind = KidsDataKind::Any>
auto read_data(const std::string &filepath) -> KidsData<kind> {
    if (std::regex_match(filepath, std::regex(".+\\.[Nn][Cc]$"))) {
        return KidsDataIO<DataFormat::NcFile>::read(filepath, kind);
    }
    throw KidsDataIOError(fmt::format("invalid path: {}", filepath));
}

/// @brief Load TolTEC data from path, given sample slice.
template <KidsDataKind kind = KidsDataKind::Any>
auto read_data_slice(const std::string &filepath, tula::container_utils::IndexSlice slice) -> KidsData<kind> {
    if (std::filesystem::is_directory(filepath)) {
        throw KidsDataIOError(fmt::format("invalid path: {}", filepath));
    }
    using result_t = KidsData<kind>;
    if (std::regex_match(filepath, std::regex(".+\\.[Nn][Cc]$"))) {
        auto kidsdata = KidsDataIO<DataFormat::NcFile>::read_slice(filepath, slice, kind);
        if constexpr (result_t::is_variant) {
            return kidsdata;
        } else {
            return std::get<result_t>(std::move(kidsdata));
        }
    }
    throw KidsDataIOError(fmt::format("invalid path: {}", filepath));
}


/// @brief Get meta.
template <KidsDataKind kind = KidsDataKind::Any>
auto get_meta(const std::string &filepath) -> std::pair<KidsDataKind, internal::meta_t> {
    if (std::filesystem::is_directory(filepath)) {
        throw KidsDataIOError(fmt::format("invalid path: {}", filepath));
    }
    if (std::regex_match(filepath, std::regex(".+\\.[Nn][Cc]$"))) {
        using IO = KidsDataIO<DataFormat::NcFile>;
        auto io = IO::from_filepath(filepath);
        return std::make_pair(io.kind(), io.meta());
    }
    throw KidsDataIOError(fmt::format("invalid path: {}", filepath));
}


}  // namespace v1
} // namespace kids::toltec

#include <kids/toltec/timestream.h>

#include <netcdf>
#include <fmt/format.h>
#include <tula/nc.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <iterator>
#include <regex>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace kids::toltec {
namespace {

using netCDF::NcFile;
using netCDF::NcVar;

struct BoundedSlice {
    Eigen::Index start;
    Eigen::Index stop;
    Eigen::Index step;
    Eigen::Index size;
};

[[nodiscard]] auto require_var(
    const NcFile &file,
    const std::string &name) -> NcVar
{
    auto var = file.getVar(name);
    if (var.isNull()) {
        throw RawTimeStreamIOError{"missing NetCDF variable " + name};
    }
    return var;
}

[[nodiscard]] auto require_dim_size(
    const NcFile &file,
    const std::string &name) -> int
{
    const auto dim = file.getDim(name);
    if (dim.isNull()) {
        throw RawTimeStreamIOError{"missing NetCDF dimension " + name};
    }
    return static_cast<int>(dim.getSize());
}

template <typename T>
[[nodiscard]] auto scalar(
    const NcFile &file,
    const std::string &name) -> T
{
    return tula::nc_utils::getscalar<T>(require_var(file, name));
}

template <typename T>
void set_scalar_if_present(
    RawTimeStreamMeta &meta,
    const NcFile &file,
    const std::string &key,
    const std::string &name)
{
    const auto var = file.getVar(name);
    if (!var.isNull()) {
        meta.set(key, scalar<T>(file, name));
    }
}

[[nodiscard]] auto bounded_slice(
    const SampleSlice &slice,
    Eigen::Index length) -> BoundedSlice
{
    const auto &[start_value, stop_value, step_value] = slice;
    const auto step = step_value.value_or(1);
    if (step <= 0) {
        throw RawTimeStreamIOError{
            "sample slice step must be a positive integer"};
    }

    auto normalize = [length](Eigen::Index index) {
        return index < 0 ? length + index : index;
    };
    const auto start =
        std::clamp(normalize(start_value.value_or(0)), Eigen::Index{0}, length);
    const auto stop = std::clamp(
        normalize(stop_value.value_or(length)), Eigen::Index{0}, length);
    if (stop < start) {
        throw RawTimeStreamIOError{
            "sample slice stop must not precede its start"};
    }
    const auto size = stop == start ? 0 : 1 + (stop - start - 1) / step;
    return {start, stop, step, size};
}

[[nodiscard]] auto read_meta(
    const std::filesystem::path &filepath,
    const NcFile &file) -> RawTimeStreamMeta
{
    RawTimeStreamMeta meta;

    auto source = filepath.string();
    const auto source_var = file.getVar("Header.Toltec.Filename");
    if (!source_var.isNull()) {
        source = tula::nc_utils::getstr(source_var);
    }
    const auto filename = std::filesystem::path{source}.filename().string();
    meta.set("source", source);
    meta.set("filename", filename);
    meta.set("instru", "toltec");
    meta.set("roachid", -1);
    meta.set("obsid", -1);
    meta.set("subobsid", -1);
    meta.set("scanid", -1);
    meta.set("ut", "");
    meta.set("kindstr", "");
    meta.set("fileext", "");

    const std::regex filename_pattern{
        R"(^toltec(\d+)_(\d+)_(\d+)_(\d+)_(\d{4}_\d{2}_\d{2}_\d{2}_\d{2}_\d{2})(?:_([^./]+))?\.(.+)$)"};
    std::smatch match;
    if (std::regex_match(filename, match, filename_pattern)) {
        meta.set("roachid", std::stoi(match[1].str()));
        meta.set("obsid", std::stoi(match[2].str()));
        meta.set("subobsid", std::stoi(match[3].str()));
        meta.set("scanid", std::stoi(match[4].str()));
        meta.set("ut", match[5].str());
        if (match[6].matched) {
            meta.set("kindstr", match[6].str());
        }
        meta.set("fileext", match[7].str());
    }

    meta.set("ntones_design", require_dim_size(file, "loclen"));
    meta.set("ntimes_all", require_dim_size(file, "time"));
    meta.set("ntones", require_dim_size(file, "iqlen"));
    meta.set("ntones_", require_dim_size(file, "toneFreqLen"));
    meta.set("nsweeps_all", require_dim_size(file, "numSweeps"));
    meta.set(
        "ntonemodelparams", require_dim_size(file, "modelParamsNum"));

    for (const auto &[key, name] : std::array{
             std::pair{"ntones_max", "Header.Toltec.MaxNumTones"},
             std::pair{"kindvar", "Header.Toltec.ObsType"},
             std::pair{"roachid", "Header.Toltec.RoachIndex"},
             std::pair{"obsid", "Header.Toltec.ObsNum"},
             std::pair{"subobsid", "Header.Toltec.SubObsNum"},
             std::pair{"scanid", "Header.Toltec.ScanNum"},
             std::pair{"cal_roachid", "Header.Toltec.RoachIndex"},
             std::pair{"cal_obsid", "Header.Toltec.TargSweepObsNum"},
             std::pair{
                 "cal_subobsid", "Header.Toltec.TargSweepSubObsNum"},
             std::pair{"cal_scanid", "Header.Toltec.TargSweepScanNum"},
             std::pair{
                 "nreps", "Header.Toltec.NumSamplesPerSweepStep"},
             std::pair{"nsweepsteps", "Header.Toltec.NumSweepSteps"},
             std::pair{"accumlen", "Header.Toltec.AccumLen"},
             std::pair{"master", "Header.Toltec.Master"}}) {
        set_scalar_if_present<std::int32_t>(meta, file, key, name);
    }
    for (const auto &[key, name] : std::array{
             std::pair{"flo_center", "Header.Toltec.LoCenterFreq"},
             std::pair{"fsmp", "Header.Toltec.SampleFreq"},
             std::pair{"atten_sense", "Header.Toltec.SenseAtten"},
             std::pair{"atten_drive", "Header.Toltec.DriveAtten"}}) {
        set_scalar_if_present<double>(meta, file, key, name);
    }

    const auto kind = meta.get_typed<int>("kindvar", 1);
    if (kind != 1) {
        throw RawTimeStreamIOError{
            "expected TolTEC raw timestream ObsType=1, found " +
            std::to_string(kind)};
    }
    const auto roachid = meta.get_typed<int>("roachid", -1);
    meta.set("roachname", "toltec" + std::to_string(roachid));

    const auto cal_roachid = meta.get_typed<int>("cal_roachid", -1);
    const auto cal_obsid = meta.get_typed<int>("cal_obsid", -1);
    const auto cal_subobsid = meta.get_typed<int>("cal_subobsid", -1);
    const auto cal_scanid = meta.get_typed<int>("cal_scanid", -1);
    if (cal_obsid > 0) {
        meta.set(
            "cal_file",
            fmt::format(
                "toltec{}_{:06d}_{:03d}_{:04d}.+\\.txt",
                cal_roachid, cal_obsid, cal_subobsid, cal_scanid));
    }
    if (meta.has("nsweepsteps") && meta.has("nreps")) {
        const auto samples_per_sweep =
            meta.get_typed<int>("nsweepsteps") *
            meta.get_typed<int>("nreps");
        if (samples_per_sweep > 0) {
            meta.set(
                "nsweeps",
                static_cast<double>(meta.get_typed<int>("ntimes_all")) /
                    samples_per_sweep);
        }
    }

    const auto ntones = meta.get_typed<int>("ntones");
    if (ntones != meta.get_typed<int>("ntones_")) {
        throw RawTimeStreamIOError{
            "iqlen and toneFreqLen dimensions do not agree"};
    }
    for (const auto &name :
         {"Data.Toltec.Is", "Data.Toltec.Qs",
          "Header.Toltec.ToneFreq"}) {
        static_cast<void>(require_var(file, name));
    }
    if (file.getVar("Data.Toltec.Ts").isNull() &&
        file.getVar("Data.Toltec.Xs").isNull()) {
        throw RawTimeStreamIOError{
            "missing NetCDF time variable Data.Toltec.Ts"};
    }
    return meta;
}

[[nodiscard]] auto read_tone_axis(
    const NcFile &file,
    const RawTimeStreamMeta &meta) -> ToneAxis
{
    const auto ntones = meta.get_typed<int>("ntones");
    const auto nsweeps = meta.get_typed<int>("nsweeps_all");
    const auto flo_center = meta.get_typed<double>("flo_center");
    const auto tone_var = require_var(file, "Header.Toltec.ToneFreq");

    Eigen::VectorXd tone_frequencies{ntones};
    if (tone_var.getDimCount() == 1) {
        tone_var.getVar(tone_frequencies.data());
    } else if (tone_var.getDimCount() == 2 && nsweeps > 0) {
        // The v1 raw-timestream reader used the first tone/model block.
        // Preserve that selection even if the file carries extra blocks.
        const std::vector<std::size_t> start{
            0, 0};
        const std::vector<std::size_t> count{
            1, static_cast<std::size_t>(ntones)};
        tone_var.getVar(start, count, tone_frequencies.data());
    } else {
        throw RawTimeStreamIOError{
            "ToneFreq must be a one- or two-dimensional variable"};
    }
    tone_frequencies.array() += flo_center;

    std::vector<std::string> labels{"f_tone"};
    Eigen::MatrixXd tone_data{1, ntones};
    tone_data.row(0) = tone_frequencies.transpose();

    const auto header_var = file.getVar("Header.Toltec.ModelParamsHeader");
    const auto model_var = file.getVar("Header.Toltec.ModelParams");
    if (!header_var.isNull() || !model_var.isNull()) {
        if (header_var.isNull() || model_var.isNull()) {
            throw RawTimeStreamIOError{
                "tone model data and its header must both be present"};
        }
        auto model_labels = tula::nc_utils::getstrs(header_var);
        const auto nparams = meta.get_typed<int>("ntonemodelparams");
        if (static_cast<int>(model_labels.size()) != nparams) {
            throw RawTimeStreamIOError{
                "ModelParamsHeader length does not match modelParamsNum"};
        }
        internal::RMatrixXd model_data{nparams, ntones};
        const std::vector<std::size_t> start{
            0, 0, 0};
        const std::vector<std::size_t> count{
            1, static_cast<std::size_t>(nparams),
            static_cast<std::size_t>(ntones)};
        model_var.getVar(start, count, model_data.data());
        tone_data.conservativeResize(1 + nparams, ntones);
        tone_data.bottomRows(nparams) = model_data;
        labels.insert(
            labels.end(),
            std::make_move_iterator(model_labels.begin()),
            std::make_move_iterator(model_labels.end()));
    }
    return {std::move(tone_data), std::move(labels)};
}

[[nodiscard]] auto read_time_axis(
    const NcFile &file,
    const BoundedSlice &slice) -> TimeAxis
{
    static constexpr std::array labels{
        "time0", "pps_count", "clock_count",
        "packet_count", "clock_at_pps", "time1"};
    auto var = file.getVar("Data.Toltec.Ts");
    if (var.isNull()) {
        // Compatibility with early TolTEC files, matching the v1 reader.
        var = require_var(file, "Data.Toltec.Xs");
    }
    if (var.getDimCount() != 2 || var.getDim(1).getSize() != labels.size()) {
        throw RawTimeStreamIOError{
            "Data.Toltec.Ts must have shape (time, 6)"};
    }

    internal::RMatrixXi data{slice.size, std::ssize(labels)};
    const std::vector<std::size_t> start{
        static_cast<std::size_t>(slice.start), 0};
    const std::vector<std::size_t> count{
        static_cast<std::size_t>(slice.size), labels.size()};
    const std::vector<std::ptrdiff_t> stride{
        static_cast<std::ptrdiff_t>(slice.step), 1};
    if (slice.size > 0) {
        var.getVar(start, count, stride, data.data());
    }
    return {
        std::move(data),
        std::vector<std::string>{labels.begin(), labels.end()}};
}

[[nodiscard]] auto read_iq(
    const NcFile &file,
    const std::string &name,
    const BoundedSlice &slice,
    Eigen::Index ntones) -> internal::RMatrixXd
{
    const auto var = require_var(file, name);
    if (var.getDimCount() != 2 ||
        var.getDim(1).getSize() != static_cast<std::size_t>(ntones)) {
        throw RawTimeStreamIOError{name + " has an unexpected shape"};
    }
    internal::RMatrixXd data{slice.size, ntones};
    const std::vector<std::size_t> start{
        static_cast<std::size_t>(slice.start), 0};
    const std::vector<std::size_t> count{
        static_cast<std::size_t>(slice.size),
        static_cast<std::size_t>(ntones)};
    const std::vector<std::ptrdiff_t> stride{
        static_cast<std::ptrdiff_t>(slice.step), 1};
    if (slice.size > 0) {
        var.getVar(start, count, stride, data.data());
    }
    return data;
}

template <typename Function>
[[nodiscard]] auto with_netcdf_errors(
    const std::filesystem::path &filepath,
    Function &&function) -> std::invoke_result_t<Function>
{
    try {
        return std::forward<Function>(function)();
    } catch (const RawTimeStreamIOError &) {
        throw;
    } catch (const netCDF::exceptions::NcException &error) {
        throw RawTimeStreamIOError{
            "failed to read TolTEC NetCDF file " + filepath.string() +
            ": " + error.what()};
    }
}

} // namespace

auto get_raw_timestream_meta(const std::filesystem::path &filepath)
    -> RawTimeStreamMeta
{
    return with_netcdf_errors(filepath, [&] {
        const NcFile file{filepath.string(), NcFile::read};
        return read_meta(filepath, file);
    });
}

auto read_raw_timestream_slice(
    const std::filesystem::path &filepath,
    SampleSlice slice) -> RawTimeStream
{
    return with_netcdf_errors(filepath, [&] {
        const NcFile file{filepath.string(), NcFile::read};
        auto result = RawTimeStream{};
        result.meta = read_meta(filepath, file);
        const auto bounds = bounded_slice(
            slice, result.meta.get_typed<int>("ntimes_all"));
        const auto ntones = result.meta.get_typed<int>("ntones");

        result.wcs.tone_axis = read_tone_axis(file, result.meta);
        result.wcs.time_axis = read_time_axis(file, bounds);
        result.is.data =
            read_iq(file, "Data.Toltec.Is", bounds, ntones);
        result.qs.data =
            read_iq(file, "Data.Toltec.Qs", bounds, ntones);
        result.meta.set("sample_slice_start", static_cast<int>(bounds.start));
        result.meta.set("sample_slice_stop", static_cast<int>(bounds.stop));
        result.meta.set("sample_slice_step", static_cast<int>(bounds.step));
        result.meta.set("sample_slice_size", static_cast<int>(bounds.size));
        return result;
    });
}

} // namespace kids::toltec

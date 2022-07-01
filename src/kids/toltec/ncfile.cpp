#include "kids/toltec/toltec.h"

#include <cmath>
#include <re2/re2.h>
#include <tula/algorithm/ei_uniquefy.h>
#include <tula/filename.h>
#include <tula/formatter/enum.h>
#include <tula/nc.h>
#include <tula/switch_invoke.h>
#include <tula/algorithm/ei_iterclip.h>
#include <tula/algorithm/ei_stats.h>
//#include <tula/nddata/cached_>

using Eigen::Dynamic;
using Eigen::Index;
using Eigen::Matrix;
using Eigen::MatrixXcd;
using Eigen::MatrixXd;
using Eigen::RowMajor;
using Eigen::VectorXd;
using RMatrixXd = Matrix<double, Dynamic, Dynamic, RowMajor>;
using RMatrixXi = Matrix<int, Dynamic, Dynamic, RowMajor>;
using tula::nc_utils::NcNodeMapper;

namespace {

// some helpers
template <typename T, typename U, typename F>
auto check_kind(T kind, U expected, F &&on_check_fail) {
    if (kind & expected) {
        return kind;
    }
    return std::forward<decltype(on_check_fail)>(on_check_fail)(kind, expected);
};
template <typename T, typename U>
auto check_kind(T kind, U expected) {
    return check_kind(
        kind, expected, [](auto kind, auto expected) -> decltype(kind) {
            throw kids::KidsDataIOError(fmt::format(
                "expect kind={} but found {} instead", expected, kind));
        });
};

} // namespace

namespace fmt {

template <typename Char>
struct formatter<kids::ToneAxis, Char>
    : formatter<kids::wcs::FrameBase<kids::ToneAxis>, Char> {};
} // namespace fmt

namespace kids::toltec {
inline namespace v1 {
namespace internal {

IO<DataFormat::NcFile>::IO(const std::string &filepath) try
    : m_source(filepath), m_ncfile(m_source, NcFile::read) {
    SPDLOG_TRACE("nc file \"{}\" {}", m_source,
                 tula::nc_utils::pprint(m_ncfile));
} catch (netCDF::exceptions::NcException &e) {
    SPDLOG_ERROR("{}", e.what());
    throw KidsDataIOError(
        fmt::format("failed to load data from file path {}", filepath));
}
IO<DataFormat::NcFile>::~IO() {
    // the close is automatically called here by the ncfile destructor.
    // m_ncfile.close();
    SPDLOG_TRACE("closed nc file {}", m_source);
}

KidsDataKind IO<DataFormat::NcFile>::kind_evaluator::evaluate(
    const IO<DataFormat::NcFile> &io) {
    auto nc = NcNodeMapper(io.ncfile(), {{"kind", "Header.Toltec.ObsType"},
                                         {"is", "Data.Toltec.Is"},
                                         {"qs", "Data.Toltec.Qs"},
                                         {"sweeps", "Data.Toltec.SweepFreq"},
                                         {"flos", "Data.Toltec.LoFreq"},
                                         {"tones", "Header.Toltec.ToneFreq"}});
    // get hint of the data kind from the obs type header
    KidsDataKind hint{KidsDataKind::Any};
    if (nc.has_var("kind")) {
        auto kindvar = tula::nc_utils::getscalar<int>(nc.var("kind"));
        SPDLOG_TRACE("found kindvar={} from {}", kindvar, nc._["kind"]);
        switch (kindvar) {
        case 1: {
            hint = KidsDataKind::RawTimeStream;
            break;
        }
        case 2: {
            hint = KidsDataKind::VnaSweep;
            break;
        }
        case 3: {
            hint = KidsDataKind::TargetSweep;
            break;
        }
        case 4: {
            // tone file.
            // it contains multiple target sweeps, and we only load
            // the last one
            hint = KidsDataKind::TargetSweep;
            break;
        }
        default: {
            // fallback to timestream
            hint = KidsDataKind::RawTimeStream;
            SPDLOG_WARN("kindvar={} unrecognized, assume {}", kindvar, hint);
            break;
        }
        }
    }
    // further validate with data entries
    SPDLOG_TRACE("validate with kind hint={}", hint);
    if (!nc.has_var("is", "qs")) {
        // solved timestream does not have Is and Qs
        hint = check_kind(KidsDataKind::SolvedTimeStream, hint,
                          [&](auto kind, auto) { return kind; });
        SPDLOG_TRACE("validated kind hint={} according to missing data {} {}",
                     hint, nc._["is"], nc._["qs"]);
    }
    // 191206: All kinds of data now contain LO freqs per sample.
    if (nc.has_var("flos")) {
        return hint;
    }
    // Compat for older data:
    if (!nc.has_var("sweeps")) {
        // raw timestream does not have sweep freqs
        hint = check_kind(KidsDataKind::RawTimeStream, hint,
                          [&](auto kind, auto) { return kind; });
        SPDLOG_TRACE("validated kind hint={} according to missing data {}",
                     hint, nc._["sweeps"]);
    }
    return hint;
}

meta_t IO<DataFormat::NcFile>::meta_evaluator::evaluate(
    const IO<DataFormat::NcFile> &io) {

    auto nc = NcNodeMapper(
        io.ncfile(), {{"ntones_design", "loclen"},
                      {"ntones_max", "Header.Toltec.MaxNumTones"},
                      // settings
                      // 20200925 this LO freq get removed
                      // and is replace with flo_center
                      {"flo_center", "Header.Toltec.LoCenterFreq"},
                      // {"flo_offset", "Header.Toltec.LoOffset"},
                      {"flos", "Data.Toltec.LoFreq"},
                      {"fsmp", "Header.Toltec.SampleFreq"},
                      {"atten_sense", "Header.Toltec.SenseAtten"},
                      {"atten_drive", "Header.Toltec.DriveAtten"},
                      // paths
                      {"source", "Header.Toltec.Filename"},
                      {"master", "Header.Toltec.Master"},
                      // info
                      {"kindvar", "Header.Toltec.ObsType"},
                      {"roachid", "Header.Toltec.RoachIndex"},
                      {"obsid", "Header.Toltec.ObsNum"},
                      {"subobsid", "Header.Toltec.SubObsNum"},
                      {"scanid", "Header.Toltec.ScanNum"},
                      // association
                      {"cal_roachid", "Header.Toltec.RoachIndex"},
                      {"cal_obsid", "Header.Toltec.TargSweepObsNum"},
                      {"cal_subobsid", "Header.Toltec.TargSweepSubObsNum"},
                      {"cal_scanid", "Header.Toltec.TargSweepScanNum"},
                      // data shape
                      {"ntimes_all", "time"},
                      {"ntones", "iqlen"},
                      {"ntones_", "toneFreqLen"},
                      {"nreps", "Header.Toltec.NumSamplesPerSweepStep"},
                      {"nsweepsteps", "Header.Toltec.NumSweepSteps"},
                      {"nsweeps_all", "numSweeps"},
                      {"ntonemodelparams", "modelParamsNum"},
                      {"accumlen", "Header.Toltec.AccumLen"}});
    meta_t meta{};
    auto set_meta = [&](const auto &key, const auto &val) {
        if (meta.has(key)) {
            meta.set(key, val);
            if constexpr (!std::is_floating_point_v<
                              std::decay_t<decltype(val)>>) {
                auto old = meta.at(key);
                if (old != meta.at(key)) {
                    SPDLOG_WARN("update meta key={} val={} -> {}", key, old,
                                val);
                }
            } else {
                // ignore check for floating number
            }
            return;
        }
        meta.set(key, val);
    };

    {
        // parse source. the info here may be overwritten by the header entries
        std::string source;
        if (!nc.has_var("source")) {
            // some old test files do not have source loc in header
            source = io.m_source;
        } else {
            source = tula::nc_utils::getstr(nc.var("source"));
        }
        set_meta("source", source);
        namespace filesystem = std::filesystem;
        filesystem::path path{meta.get_str("source")};
        set_meta("filename", path.filename().string());

        re2::RE2 re_filename(
            "^(?P<instru>toltec)(?P<roachid>\\d+)_"
            "(?P<obsid>\\d+)_(?P<subobsid>\\d+)_(?P<scanid>\\d+)_"
            "(?P<ut>\\d{4}_\\d{2}_\\d{2}(?:_\\d{2}_\\d{2}_\\d{2}))(?:_"
            "(?P<kindstr>[^_\\/.]+))?\\.(?P<fileext>.+)$");
        assert(re_filename.ok());
        set_meta("instru", "");
        set_meta("roachid", -1);
        set_meta("obsid", -1);
        set_meta("subobsid", -1);
        set_meta("scanid", -1);
        set_meta("ut", "");
        set_meta("kindstr", "");
        set_meta("fileext", "");
        // SPDLOG_TRACE("meta dict before parsing filename: {}",
        // meta.pformat());
        auto nargs = re_filename.NumberOfCapturingGroups();
        std::vector<RE2::Arg> args(TULA_SIZET(nargs));
        for (const auto &[k, i] : re_filename.NamedCapturingGroups()) {
            SPDLOG_TRACE("capture group {} -> {}/{}", k, i, nargs);
            auto &v = meta.at(k);
            auto &a = args[TULA_SIZET(i) - 1];
            if (std::holds_alternative<int>(v)) {
                a = RE2::Arg(&std::get<int>(v));
            } else if (std::holds_alternative<std::string>(v)) {
                a = RE2::Arg(&std::get<std::string>(v));
            } else {
                throw std::runtime_error(
                    fmt::format("entry for named group {} not initialized", k));
            }
        }
        std::vector<RE2::Arg *> ptr_args(args.size());
        for (std::size_t i = 0; i < args.size(); ++i) {
            ptr_args[i] = &args[i];
        }
        if (!re2::RE2::FullMatchN(meta.get_str("filename"), re_filename,
                                  ptr_args.data(), nargs)) {
            SPDLOG_WARN("failed parse filename {}", meta.get_str("filename"));
        }
        // SPDLOG_TRACE("meta dict after parsing filename: {}", meta.pformat());
    }
    // parse headers

    // string
    //     for (const auto &key : {"master"}) {
    //         if (nc.has_var(key)) {
    //             set_meta(key, nc_utils::getstr(nc.var(key)));
    //         }
    //     }
    // int
    for (const auto &key :
         {"master", "ntones_design", "ntones_max", "kindvar", "roachid",
          "obsid", "subobsid", "scanid", "cal_roachid", "cal_obsid",
          "cal_subobsid", "cal_scanid"}) {
        if (nc.has_var(key)) {
            set_meta(key, tula::nc_utils::getscalar<int>(nc.var(key)));
        } else if (nc.has_dim(key)) {
            set_meta(key, tula::meta::size_cast<int>(nc.dim(key).getSize()));
        }
    }
    set_meta("roachname", fmt::format("{}{}", meta.get_str("instru"),
                                      meta.get_typed<int>("roachid")));
    // parse shape information
    for (const auto &key :
         {"ntimes_all", "ntones", "ntones_", "nreps", "nsweepsteps",
          "nsweeps_all", "ntonemodelparams", "accumlen"}) {
        if (nc.has_var(key)) {
            set_meta(key, tula::nc_utils::getscalar<int>(nc.var(key)));
        } else if (nc.has_dim(key)) {
            set_meta(key, tula::meta::size_cast<int>(nc.dim(key).getSize()));
        }
    }
    // some sanity check
    auto m = [&](const auto &key) { return meta.get_typed<int>(key); };
    if (meta.has("ntones") && meta.has("ntones_")) {
        assert(m("ntones") == m("ntones_"));
    }
    if (meta.has("nsweepsteps")) {
        if (meta.get_str("kindstr") == "tune" && m("nsweepsteps") == 491) {
            SPDLOG_WARN("inconsistency nsweepsteps {} with kindstr {}",
                        m("nsweepsteps"), meta.get_str("kindstr"));
            set_meta("nsweepsteps",
                     m("ntimes_all") / m("nsweeps_all") / m("nreps"));
        }
        if (meta.has("nreps")) {
            // 2019/08/06:
            //     if nreps is present, nsweepsteps is the actual
            //     number of sweep steps. ntimes = nsweepsteps * nreps.
            //     If nreps is not present, we assume nreps=10, as
            //     this was the only value used. and ntimes = nsweepsteps
            // assert(m("nsweepsteps") * m("nreps") == m("ntimes_all"));
        } else {
            int nreps_default{10};
            SPDLOG_WARN("assume nreps={} for legacy data", nreps_default);
            meta.set("nreps", nreps_default);
            // assert(m("nsweepsteps") % m("nreps") == 0);
            meta.set("nsweepsteps", m("nsweepsteps") / m("nreps"));
        }
        auto ntimespersweep = m("nsweepsteps") * m("nreps");
        // assert(m("ntimes_all") % ntimespersweep == 0);

        if (meta.has("nsweeps_all")) {
            // on repeated data, packets could be lost, so we have to
            // set this to double value.
            meta.set("nsweeps", double(m("ntimes_all")) / ntimespersweep);
        }
    }
    // double
    // 191206: For data has flos, they do not have flo and flo_offset header
    // because we want to represent the tones in absolute frequency
    // rather than relative to LO, we need to get a LO nevertheless.
    // 200925: get rid of flo_offset and flo. There key words should always be
    // present
    for (const auto &key :
         {"flo_center", "fsmp", "atten_sense", "atten_drive"}) {
        // SPDLOG_INFO("{} {} {}", key, nc.has_var(key), meta.pformat());
        if (nc.has_var(key)) {
            set_meta(key, tula::nc_utils::getscalar<double>(nc.var(key)));
        }
    }
    // calib id
    if (meta.has("cal_obsid") && m("cal_obsid") > 0) {
        auto pattern =
            fmt::format("{}{:d}_{:06d}_{:03d}_{:04d}.+\\.txt",
                        meta.get_str("instru"), m("cal_roachid"),
                        m("cal_obsid"), m("cal_subobsid"), m("cal_scanid"));
        meta.set("cal_file", pattern);
    }
    SPDLOG_TRACE("meta={}", meta.pformat());
    return meta;
}

tone_axis_t IO<DataFormat::NcFile>::tone_axis_evaluator::evaluate(
    const IO<DataFormat::NcFile> &io) {
    const auto &meta = io.meta();
    auto ntones = meta.get_typed<int>("ntones", -1);
    if (ntones <= 0) {
        throw std::runtime_error("no tones found");
    }

    auto nc = NcNodeMapper(
        io.ncfile(),
        {{"tones", "Header.Toltec.ToneFreq"},
         {"tonemodelparams", "Header.Toltec.ModelParams"},
         {"tonemodelparamsheader", "Header.Toltec.ModelParamsHeader"}});

    if (!nc.has_var("tones")) {
        throw std::runtime_error("no tone data found");
    }

    // tone param header
    std::vector<std::string> toneparamsheader{"f_tone"};
    if (nc.has_var("tonemodelparamsheader")) {
        // model header
        auto hdr = tula::nc_utils::getstrs(nc.var("tonemodelparamsheader"));
        SPDLOG_TRACE("tone model params header: {}", hdr);
        // populate the tone params table
        assert(TULA_SIZET(meta.get_typed<int>("ntonemodelparams")) == hdr.size());
        toneparamsheader.reserve(toneparamsheader.size() + hdr.size());
        toneparamsheader.insert(toneparamsheader.end(),
                                std::make_move_iterator(hdr.begin()),
                                std::make_move_iterator(hdr.end()));
    }
    SPDLOG_DEBUG("tone params header: {}", toneparamsheader);

    // tone param data
    Eigen::MatrixXd toneparams(1, ntones);
    // number of sweep blocks for sweep data.
    // 20200326: repeated data may have incomplete sweeps.
    // default to 1 for timestream
    int nsweeps = 1;
    auto nsweeps_all = meta.get_typed<int>("nsweeps_all");
    // check nsweeps for sweep
    {
        auto n = meta.get_typed<double>("nsweeps", 1);
        assert(n <= 2.001); // currently nsweeps can only be below 2.
        if (n > 1.9) {
            nsweeps = 2;
            SPDLOG_WARN("incomplete number of sweeps {}, round to {}", n,
                        nsweeps);
        }
    }
    SPDLOG_TRACE("read nsweeps={} in n_sweeps_all={}", nsweeps, nsweeps_all);
    {
        const auto &v_tones = nc.var("tones");
        VectorXd tfs(ntones);
        if (v_tones.getDimCount() == 1) {
            assert(nsweeps == 1);
            // single set of tones compat with old data
            v_tones.getVar(&tfs.coeffRef(0));
        } else {
            // must be sweep data
            assert(nsweeps > 0);
            RMatrixXd tmp(nsweeps_all, ntones);
            v_tones.getVar(&tmp.coeffRef(0));
            // multiple set of tones, read in the last set
            // auto i = std::vector<std::size_t>{SIZET(nsweeps - 1), 0};
            // auto s = std::vector<std::size_t>{1, SIZET(ntones)};
            // SPDLOG_TRACE("read tfs for sweep {} with i={} s={}", nsweeps - 1,
            // i,
            //             s);
            // v_tones.getVar(i, s, &tfs.coeffRef(0));
            tfs = tmp.row(nsweeps - 1);
            // log the difference between the rows
            SPDLOG_TRACE("d_tfs{}", tmp.row(nsweeps - 1) - tmp.row(0));
        }
        tfs.array() += meta.get_typed<double>("flo_center");
        SPDLOG_TRACE("tfs{}", tfs);
        toneparams.row(0) = std::move(tfs);
    }
    // model params data
    if (nc.has_var("tonemodelparams")) {
        auto ntonemodelparams = meta.get_typed<int>("ntonemodelparams");
        RMatrixXd tps(ntonemodelparams, ntones);
        const auto &v_tonemodelparams = nc.var("tonemodelparams");
        assert(nsweeps > 0 && v_tonemodelparams.getDimCount() == 3);
        auto i = std::vector<std::size_t>{TULA_SIZET(nsweeps - 1), 0, 0};
        auto s = std::vector<std::size_t>{1, TULA_SIZET(ntonemodelparams),
                                          TULA_SIZET(ntones)};
        SPDLOG_TRACE("read tone model params for sweep {} with i={} s={}",
                     nsweeps - 1, i, s);
        v_tonemodelparams.getVar(i, s, &tps.coeffRef(0));
        SPDLOG_TRACE("tps{}", tps);
        toneparams.conservativeResize(toneparams.rows() + ntonemodelparams,
                                      ntones);
        toneparams.bottomRows(ntonemodelparams) = std::move(tps);
    }
    SPDLOG_TRACE("tone params{}", toneparams);
    assert(TULA_SIZET(toneparams.rows()) == toneparamsheader.size());
    return {std::move(toneparams), std::move(toneparamsheader)};
}

IO<DataFormat::NcFile>::data_reducer_segments_t
IO<DataFormat::NcFile>::data_reducer_segments_evaluator::evaluate(
    const IO<DataFormat::NcFile> &io) {
    const auto &kind = io.kind();
    const auto &meta = io.meta();
    data_reducer_segments_t segments;
    if (kind & KidsDataKind::Sweep) {
        int nsweeps = 1;
        {
            auto n = meta.get_typed<double>("nsweeps");
            // 20200326: fix repeated missing data
            if (n > 1.9) {
                nsweeps = 2;
            }
        }
        auto nreps = meta.get_typed<int>("nreps");
        auto nsweepsteps_all = meta.get_typed<int>("nsweepsteps");

        auto ntimes_all = meta.get_typed<int>("ntimes_all");
        // load the last sweep, and round the samples to nrep
        auto nsweepsteps = (ntimes_all / nreps * nreps) % nsweepsteps_all;
        if (nsweepsteps == 0) {
            nsweepsteps = nsweepsteps_all;
        }
        auto ntimes = nsweepsteps * nreps;
        auto i0 = (nsweeps - 1) * nsweepsteps_all * nreps;
        assert(ntimes + i0 <= ntimes_all);
        SPDLOG_TRACE(
            "build reducer segments for sweep {} using {} of {} sweepsteps",
            nsweeps - 1, nsweepsteps, nsweepsteps_all);
        segments.reserve(TULA_SIZET(nsweepsteps));
        for (Index i = 0; i < nsweepsteps; ++i) {
            segments.emplace_back(i0 + i * nreps, i0 + (i + 1) * nreps);
        }
        SPDLOG_TRACE("sweep reducer segments{}",
                     tula::eigen_utils::to_eigen(segments));
    } else if (kind & KidsDataKind::TimeStream) {
        SPDLOG_TRACE("timestream reducer segments{}", segments);
    }
    return segments;
}

sweep_axis_t IO<DataFormat::NcFile>::sweep_axis_evaluator::evaluate(
    const IO<DataFormat::NcFile> &io) {
    check_kind(io.kind(), KidsDataKind::Sweep);
    const auto &meta = io.meta();

    auto nc = NcNodeMapper(io.ncfile(), {{"sweeps", "Data.Toltec.SweepFreq"},
                                         {"flos", "Data.Toltec.LoFreq"}});
    // 191206: flos now contains all sweeps info
    if (!nc.has_var("flos")) {
        // compat for old data
        // 20200925: just error out with this old behavior.
        // if (!nc.has_var("sweeps")) {
        throw std::runtime_error("no sweep data found");
        // }
    }
    auto get_sweep_var = [&nc]() { return nc.var("flos"); };
    VectorXd sfs;
    {
        const auto &sweep_reducer_segment = io.m_data_reducer_segments(io);
        assert(sweep_reducer_segment.size() > 0);
        auto nsweepsteps =
            tula::meta::size_cast<Index>(sweep_reducer_segment.size());
        SPDLOG_TRACE("read {} of {} sweepsteps", nsweepsteps,
                     meta.get_typed<int>("nsweepsteps"));
        const auto &v_sweeps = get_sweep_var();
        auto i0 = sweep_reducer_segment.front().first;
        auto ntimes = sweep_reducer_segment.back().second - i0;
        auto i = std::vector<std::size_t>{TULA_SIZET(i0), 0};
        auto s = std::vector<std::size_t>{TULA_SIZET(ntimes), 1};
        SPDLOG_TRACE("read sfs with i={} s={}", i, s);
        sfs.resize(ntimes);
        v_sweeps.getVar(i, s, &sfs.coeffRef(0));
        // reduce the sweep frequencies
        for (Index i = 0; i < nsweepsteps; ++i) {
            // s, e are w.r.t 0
            auto [s, e] = sweep_reducer_segment.at(TULA_SIZET(i));
            sfs.coeffRef(i) = sfs.segment(s - i0, e - s).maxCoeff();
            //  some data have this bug
            if (sfs.coeff(i) == 0.) {
                SPDLOG_WARN("corrupted sweep data at {}", i);
            }
        }
        sfs.conservativeResize(nsweepsteps);
        sfs.array() -= meta.get_typed<double>("flo_center");
    }
    SPDLOG_TRACE("sfs{}", sfs);
    return {std::move(sfs)};
}

time_axis_t IO<DataFormat::NcFile>::time_axis_evaluator::evaluate(
    const IO<DataFormat::NcFile> &io) {
    check_kind(io.kind(), KidsDataKind::TimeStream);
    const auto &meta = io.meta();

    // COMPAT: some old data has Xs for time variables
    auto nc = NcNodeMapper(io.ncfile(), {{"times", "Data.Toltec.Ts"},
                                         {"times0", "Data.Toltec.Xs"}});
    if ((!nc.has_var("times")) && nc.has_var("times0")) {
        SPDLOG_WARN("file contains old time variable {}", nc._["times0"]);
        nc._["times"] = nc._["times0"];
    }
    if (!nc.has_var("times")) {
        throw std::runtime_error("no time data found");
    }
    // ClockTime (sec), PpsCount (pps ticks), ClockCount (clock ticks),
    // PacketCount (packet ticks), PpsTime (clock ticks), ClockTimeNanoSec
    // (nsec)
    std::vector<std::string> timesheader{"time0",        "pps_count",
                                         "clock_count",  "packet_count",
                                         "clock_at_pps", "time1"};
    RMatrixXi times;
    {
        const auto &time_reducer_segment = io.m_data_reducer_segments(io);
        if (time_reducer_segment.size() > 0) {
            // reduce time
            throw std::runtime_error("reducing time data is not implemented");
        } else {
            tula::logging::scoped_timeit _0("read time axis all");
            // load all times
            const auto &v_times = nc.var("times");
            auto ntimes_all = meta.get_typed<int>("ntimes_all");
            auto ctimes = tula::meta::size_cast<Index>(timesheader.size());
            assert(TULA_SIZET(ctimes) == v_times.getDim(1).getSize());
            RMatrixXi buf(ntimes_all, ctimes);
            v_times.getVar(buf.data());
            times = std::move(buf);
            /*
            for (Index c = 0; c < ctimes; ++c) {
                auto i = std::vector<std::size_t>{0, SIZET(c)};
                auto s = std::vector<std::size_t>{SIZET(ntimes_all), 1};
                auto ss = std::vector<std::ptrdiff_t>{1, ctimes};
                v_times.getVar(i, s, ss, times.col(c).data());
            }
            */
            SPDLOG_TRACE("timesheader{}", timesheader);
            SPDLOG_TRACE("times{}", times);
        }
    }
    return {std::move(times), std::move(timesheader)};
}

time_axis_t IO<DataFormat::NcFile>::time_axis_slice(IndexSlice slice) const {
    auto &io = *this;
    check_kind(io.kind(), KidsDataKind::TimeStream);
    const auto &meta = io.meta();

    // COMPAT: some old data has Xs for time variables
    auto nc = NcNodeMapper(io.ncfile(), {{"times", "Data.Toltec.Ts"},
                                         {"times0", "Data.Toltec.Xs"}});
    if ((!nc.has_var("times")) && nc.has_var("times0")) {
        SPDLOG_WARN("file contains old time variable {}", nc._["times0"]);
        nc._["times"] = nc._["times0"];
    }
    if (!nc.has_var("times")) {
        throw std::runtime_error("no time data found");
    }
    // ClockTime (sec), PpsCount (pps ticks), ClockCount (clock ticks),
    // PacketCount (packet ticks), PpsTime (clock ticks), ClockTimeNanoSec
    // (nsec)
    std::vector<std::string> timesheader{"time0",        "pps_count",
                                         "clock_count",  "packet_count",
                                         "clock_at_pps", "time1"};
    RMatrixXi times;
    {
        const auto &time_reducer_segment = io.m_data_reducer_segments(io);
        if (time_reducer_segment.size() > 0) {
            // reduce time
            throw std::runtime_error("reducing time data is not implemented");
        } else {
            tula::logging::scoped_timeit _0("read time axis slice");
            // load the time slice
            const auto &v_times = nc.var("times");
            auto ntimes_all = meta.get_typed<int>("ntimes_all");
            auto ctimes = tula::meta::size_cast<Index>(timesheader.size());
            assert(TULA_SIZET(ctimes) == v_times.getDim(1).getSize());
            auto _slice = tula::container_utils::to_indices(slice, ntimes_all);
            const auto &[start, stop, step, size] = _slice;
            SPDLOG_TRACE("read time with index range {}", _slice);
            std::vector<std::size_t> i{TULA_SIZET(start), 0};
            std::vector<std::size_t> s{TULA_SIZET(size), TULA_SIZET(ctimes)};
            std::vector<std::ptrdiff_t> d{
                tula::meta::size_cast<std::ptrdiff_t>(step), 1};
            RMatrixXi buf(size, ctimes);
            v_times.getVar(i, s, d, buf.data());
            times = std::move(buf);
            /*
            for (Index c = 0; c < ctimes; ++c) {
                auto i = std::vector<std::size_t>{0, SIZET(c)};
                auto s = std::vector<std::size_t>{SIZET(ntimes_all), 1};
                auto ss = std::vector<std::ptrdiff_t>{1, ctimes};
                v_times.getVar(i, s, ss, times.col(c).data());
            }
            */
            SPDLOG_TRACE("timesheader{}", timesheader);
            SPDLOG_TRACE("times{}", times);
        }
    }
    return {std::move(times), std::move(timesheader)};
}

sweep_data_t IO<DataFormat::NcFile>::sweep_data_evaluator::evaluate(
    const IO<DataFormat::NcFile> &io) {
    check_kind(io.kind(), KidsDataKind::Sweep);

    auto nc = NcNodeMapper(
        io.ncfile(), {{"is", "Data.Toltec.Is"}, {"qs", "Data.Toltec.Qs"}});
    if (!nc.has_var("is", "qs")) {
        throw std::runtime_error("no iqs data found");
    }
    const auto &meta = io.meta();
    const auto &tone_axis = io.tone_axis();
    auto ntones = tone_axis.array_shape();

    sweep_data_t iqs;
    {
        const auto &sweep_reducer_segment = io.m_data_reducer_segments(io);
        assert(sweep_reducer_segment.size() > 0);
        auto nsweepsteps =
            tula::meta::size_cast<Index>(sweep_reducer_segment.size());
        SPDLOG_TRACE("get iqs from {} of {} sweepsteps", nsweepsteps,
                     meta.get_typed<int>("nsweepsteps"));

        auto i0 = sweep_reducer_segment.front().first;
        auto ntimes = sweep_reducer_segment.back().second - i0;
        auto i = std::vector<std::size_t>{TULA_SIZET(i0), 0};
        auto s =
            std::vector<std::size_t>{TULA_SIZET(ntimes), TULA_SIZET(ntones)};
        
        double reduce_clip_sigma = 2;
        int reduce_clip_maxiter = 5;
         auto iterclip = tula::alg::iterclip(
                    // use median-MAD for robust statistics
                    [](const auto &v) { return tula::alg::nanmedmad(v); },
                    // the snr cut, we only cut the positive side
                    [n = reduce_clip_sigma](const auto &v, const auto &c,
                                           const auto &s) { return (v <= c + n * s) && (v >= c - n * s); },
                    reduce_clip_maxiter);

        auto read = [&](const auto &key) {
            const auto &v_data = nc.var(key);
            RMatrixXd data(ntimes, ntones);
            SPDLOG_TRACE("read {} with i={} s={}", key, i, s);
            v_data.getVar(i, s, data.data());
            // reduce the data
            for (Index i = 0; i < nsweepsteps; ++i) {
                // s, e are w.r.t 0
                auto [s, e] = sweep_reducer_segment.at(TULA_SIZET(i));
                data.row(i) =
                    data.block(s - i0, 0, e - s, ntones).colwise().mean();
                // iter clip the data
                /*
                for (Index j = 0; j < ntones; ++j) {
                    std::tie(
                        std::ignore,
                         std::ignore,
                         data.coeffRef(i, j),
                         std::ignore) = iterclip(
                        data.col(j).segment(s-i0, e - s));
                }
                */
            }
            data.conservativeResize(nsweepsteps, ntones);
            return data;
        };
        iqs.resize(nsweepsteps, ntones);
        iqs.real() = read("is");
        iqs.imag() = read("qs");
    }
    return iqs;
};

timestream_data_t IO<DataFormat::NcFile>::timestream_data_evaluator::evaluate(
    const IO<DataFormat::NcFile> &io) {
    check_kind(io.kind(), KidsDataKind::TimeStream);

    auto nc = NcNodeMapper(io.ncfile(), {{"is", "Data.Toltec.Is"},
                                         {"qs", "Data.Toltec.Qs"},
                                         {"xs", "Kidsproc.Xs"},
                                         {"rs", "Kidsproc.Rs"}});
    if ((!nc.has_var("is", "qs")) && (!nc.has_var("xs", "rs"))) {
        throw std::runtime_error("no data found");
    }
    const auto &meta = io.meta();
    const auto &tone_axis = io.tone_axis();
    auto ntones = tone_axis.array_shape();

    timestream_data_t iqs;
    {
        const auto &time_reducer_segment = io.m_data_reducer_segments(io);
        if (time_reducer_segment.size() > 0) {
            // reduce time
            throw std::runtime_error(
                "pre-reducing timestream data is not implemented");
        } else {
            tula::logging::scoped_timeit _0("read data");
            // load all
            auto ntimes_all = meta.get_typed<int>("ntimes_all");
            auto i = std::vector<std::size_t>{0, 0};
            auto s = std::vector<std::size_t>{TULA_SIZET(ntimes_all),
                                              TULA_SIZET(ntones)};
            auto read = [&](const auto &key, auto &&data) {
                const auto &v_data = nc.var(key);
                SPDLOG_TRACE("read {} with i={} s={}", key, i, s);
                data.resize(ntimes_all, ntones);
                v_data.getVar(i, s, data.data());
            };
            if (nc.has_var("is", "qs")) {
                read("is", iqs.first);
                read("qs", iqs.second);
            } else if (nc.has_var("xs", "rs")) {
                read("xs", iqs.first);
                read("rs", iqs.second);
            }
        }
    }
    return iqs;
};

timestream_data_t
IO<DataFormat::NcFile>::timestream_data_slice(IndexSlice slice) const {
    auto &io = *this;
    check_kind(io.kind(), KidsDataKind::TimeStream);

    auto nc = NcNodeMapper(io.ncfile(), {{"is", "Data.Toltec.Is"},
                                         {"qs", "Data.Toltec.Qs"},
                                         {"xs", "Kidsproc.Xs"},
                                         {"rs", "Kidsproc.Rs"}});
    if ((!nc.has_var("is", "qs")) && (!nc.has_var("xs", "rs"))) {
        throw std::runtime_error("no data found");
    }
    const auto &meta = io.meta();
    const auto &tone_axis = io.tone_axis();
    auto ntones = tone_axis.array_shape();

    timestream_data_t iqs;
    {
        const auto &time_reducer_segment = io.m_data_reducer_segments(io);
        if (time_reducer_segment.size() > 0) {
            // reduce time
            throw std::runtime_error(
                "pre-reducing timestream data is not implemented");
        } else {
            tula::logging::scoped_timeit _0("read data");

            // load the time slice
            auto ntimes_all = meta.get_typed<int>("ntimes_all");
            auto _slice = tula::container_utils::to_indices(slice, ntimes_all);
            const auto &[start, stop, step, size] = _slice;
            std::vector<std::size_t> i{TULA_SIZET(start), 0};
            std::vector<std::size_t> s{TULA_SIZET(size), TULA_SIZET(ntones)};
            std::vector<std::ptrdiff_t> d{
                tula::meta::size_cast<std::ptrdiff_t>(step), 1};
            auto read = [&, size = size](const auto &key, auto &&data) {
                const auto &v_data = nc.var(key);
                SPDLOG_TRACE("read {} with i={} s={} d={}", key, i, s, d);
                data.resize(size, ntones);
                v_data.getVar(i, s, d, data.data());
            };
            if (nc.has_var("is", "qs")) {
                read("is", iqs.first);
                read("qs", iqs.second);
            } else if (nc.has_var("xs", "rs")) {
                read("xs", iqs.first);
                read("rs", iqs.second);
            }
        }
    }
    SPDLOG_TRACE("timestream data{}", iqs);
    return iqs;
};

KidsData<> IO<DataFormat::NcFile>::read(const std::string &filepath,
                                        KidsDataKind expected_kind) {
    auto io = Self(filepath);
    const auto &kind = io.kind();
    check_kind(kind, expected_kind);
    SPDLOG_TRACE("read nc file {} kind={} meta={}", filepath, io.kind(),
                 io.meta().pformat());
    // dispatch
    auto dispatch = [&](auto kind_) {
        constexpr auto kind = std::decay_t<decltype(kind_)>::value;
        KidsData<kind> result;
        result.meta = io.meta();
        result.wcs.tone_axis = io.tone_axis();
        if constexpr (kind & KidsDataKind::Sweep) {
            result.wcs.sweep_axis = io.sweep_axis();
            SPDLOG_TRACE("tone_axis={} sweep_axis={}", io.tone_axis(),
                         io.sweep_axis());
            result.iqs() = io.sweep_data_nocache();
        } else if constexpr (kind & KidsDataKind::TimeStream) {
            result.wcs.time_axis = io.time_axis();
            SPDLOG_TRACE("tone_axis={} time_axis={}", result.wcs.tone_axis,
                         result.wcs.time_axis);
            if constexpr (kind == KidsDataKind::RawTimeStream) {
                std::tie(result.is.data, result.qs.data) =
                    io.timestream_data_nocache();
            } else if constexpr (kind == KidsDataKind::SolvedTimeStream) {
                std::tie(result.xs.data, result.rs.data) =
                    io.timestream_data_nocache();
            } else {
                throw KidsDataIOError(fmt::format("unexpect kind={}", kind));
            }
        } else {
            throw KidsDataIOError(fmt::format("unexpect kind={}", kind));
            // static_assert(meta::always_false<KidsDataKind>::value,
            //              "UNEXPECTED KIND");
        }
        return result;
    };
    switch (kind) {
    case KidsDataKind::VnaSweep: {
        return {dispatch(tula::meta::scalar_t<KidsDataKind::VnaSweep>{})};
    }
    case KidsDataKind::TargetSweep: {
        return {dispatch(tula::meta::scalar_t<KidsDataKind::TargetSweep>{})};
    }
    case KidsDataKind::RawTimeStream: {
        return {dispatch(tula::meta::scalar_t<KidsDataKind::RawTimeStream>{})};
    }
    case KidsDataKind::SolvedTimeStream: {
        return {
            dispatch(tula::meta::scalar_t<KidsDataKind::SolvedTimeStream>{})};
    }
    default:
        throw KidsDataIOError(fmt::format("unexpected kind={}", kind));
    }
}

KidsData<> IO<DataFormat::NcFile>::read_slice(const std::string &filepath,
                                              IndexSlice slice,
                                              KidsDataKind expected_kind) {
    auto io = Self(filepath);
    const auto &kind = io.kind();
    check_kind(kind, expected_kind);
    SPDLOG_TRACE("read nc data slice={} kind={} meta={}", slice, io.kind(),
                 io.meta().pformat());
    // dispatch
    auto dispatch = [&](auto kind_) {
        constexpr auto kind = std::decay_t<decltype(kind_)>::value;
        KidsData<kind> result;
        result.meta = io.meta();
        result.wcs.tone_axis = io.tone_axis();
        if constexpr (kind & KidsDataKind::TimeStream) {
            result.wcs.time_axis = io.time_axis_slice(slice);
            SPDLOG_TRACE("tone_axis={} time_axis={}", result.wcs.tone_axis,
                         result.wcs.time_axis);
            // update meta with slice
            const auto &[start, stop, step, size] =
                tula::container_utils::to_indices(
                    slice, result.meta.template get_typed<int>("ntimes_all"));
            result.meta.set("sample_slice_start",
                            tula::meta::size_cast<int>(start));
            result.meta.set("sample_slice_stop",
                            tula::meta::size_cast<int>(stop));
            result.meta.set("sample_slice_step",
                            tula::meta::size_cast<int>(step));
            result.meta.set("sample_slice_size",
                            tula::meta::size_cast<int>(size));
            if constexpr (kind == KidsDataKind::RawTimeStream) {
                std::tie(result.is.data, result.qs.data) =
                    io.timestream_data_slice(slice);
            } else if constexpr (kind == KidsDataKind::SolvedTimeStream) {
                std::tie(result.xs.data, result.rs.data) =
                    io.timestream_data_slice(slice);
            } else {
                throw KidsDataIOError(fmt::format("unexpect kind={}", kind));
            }
        } else {
            throw KidsDataIOError(fmt::format("unexpect kind={}", kind));
            // static_assert(meta::always_false<KidsDataKind>::value,
            //              "UNEXPECTED KIND");
        }
        return result;
    };
    switch (kind) {
    case KidsDataKind::VnaSweep: {
        return {dispatch(tula::meta::scalar_t<KidsDataKind::VnaSweep>{})};
    }
    case KidsDataKind::TargetSweep: {
        return {dispatch(tula::meta::scalar_t<KidsDataKind::TargetSweep>{})};
    }
    case KidsDataKind::RawTimeStream: {
        return {dispatch(tula::meta::scalar_t<KidsDataKind::RawTimeStream>{})};
    }
    case KidsDataKind::SolvedTimeStream: {
        return {
            dispatch(tula::meta::scalar_t<KidsDataKind::SolvedTimeStream>{})};
    }
    default:
        throw KidsDataIOError(fmt::format("unexpected kind={}", kind));
    }
}

bool IO<DataFormat::NcFile>::write(
    [[maybe_unused]] const KidsData<> &data,
    [[maybe_unused]] const std::string &filepath) {
    SPDLOG_TRACE("write {:s} to ncfile {}", data, filepath);
    throw KidsDataIOError("writing of NcFile is not yet implemneted");
    return true;
}

} // namespace internal
} // namespace v1
} // namespace kids::toltec

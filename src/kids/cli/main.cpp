#include "kids/core/kidsdata.h"
// #include "kids/proc/proc.h"
#include "kids/sweep/finder.h"
#include "kids/sweep/fitter.h"
#include "kids/timestream/solver.h"
#include "kids/toltec/toltec.h"
#include "tula/filename.h"
#include <tula/cli.h>
#include <tula/algorithm/ei_stats.h>
#include <tula/algorithm/ei_iterclip.h>
#include <tula/config/flatconfig.h>
#include <tula/formatter/container.h>
#include <tula/grppi.h>
#include <tula/logging.h>

#include <cstdlib>
#include <kidscpp_config/config.h>
#include <kidscpp_config/gitversion.h>

// Some implementation typedefs
using keymap_t = std::unordered_map<std::string, std::string>;
using config_t = tula::config::FlatConfig;

/// @brief Helper to construct Finder from CLI.
struct finder_options {

    using Finder = kids::SweepKidsFinder;

    constexpr static auto cli = [](auto &r) {
        using namespace tula::cli::clipp_builder;
        // clang-format off
        return "finder configs" % g(
        r(p("finder_threshold"     ), "Detection threshold",
                                     10., doub()),
        r(p("finder_resample"      ), "Frequency grid to use for the resampling,"
                                     " specified in '[start]:[end][:step]'",
                                     ":", str()),
        r(p("finder_smooth_size"   ), "Window size to smooth D21",
                                     5, int_()),
        r(p("finder_stats_clip_sigma"   ), "N sigmas to clip for compute stddev",
                                     3., doub()),
        r(p("finder_use_savgol_deriv"), "Use SavGol deriv."),
        r(p("finder_plot_fit"), "Plot fit result when --plot is set."),
        r(p("output_d21"          ), "Save the computed D21",
                                     "{stem}{suffix}.{ext}", opt_str("dest")),
        r(p("output_processed"    ), "Save the reduced data",
                                     "{stem}{suffix}.{ext}", opt_str("dest"))
        ); // clang-format on
    };
    constexpr static auto finder = [](auto &rc) {
        Finder::Config conf{};
        for (auto &[confkey, rckey] : keymap_t{
                 {"output_processed", "output_processed"},
                 {"output_d21", "output_d21"},
                 {"detect_min_thresh", "finder_threshold"},
                 {"resample", "finder_resample"},
                 {"stats_clip_sigma", "finder_stats_clip_sigma"},
                 {"data_smooth_size", "finder_smooth_size"},
                 {"use_savgol_deriv", "finder_use_savgol_deriv"},
                 {"fitter_weight_window_type", "fitter_weight_window_type"},
                 {"fitter_weight_window_Qr", "fitter_weight_window_Qr"},
                 {"fitter_lim_Qr_min", "fitter_lim_Qr_min"},
                 {"fitter_lim_Qr_max", "fitter_lim_Qr_max"},
                 {"fitter_lim_gain_min", "fitter_lim_gain_min"},
                 {"plot_fit", "finder_plot_fit"},
                 //
                 // {"fitter_weight_window_fwhm", "fitter_weight_window_fwhm"},
                 {"exmode", "grppiex"}}) {
            if (rc.is_set(rckey)) {
                conf.set(confkey, rc.at(rckey));
            }
        }
        return Finder{std::move(conf)};
    };
};

/// @brief Helper to construct Fitter from CLI.
struct fitter_options {

    using Fitter = kids::SweepFitter;

    constexpr static auto cli = [](auto &r) {
        using namespace tula::cli::clipp_builder;
        // clang-format off
        return "fitter configs" % g(
        r(p("fitter_weight_window_type"), "Fit with weight window of this"
                                         " type applied",
                                         Fitter::WeightOption::lorentz,
                                         list(Fitter::WeightOption{})),
        //_(rc, p("fitter_weight_window_fwhm"), "The FWHM of the weight window in Hz",
        //                               1.5e4, doub()),
        r(p("fitter_weight_window_Qr"  ), "The Qr of the detectors to compute weight window",
                                         2e4, doub()),
        r(p("fitter_lim_Qr_min"  ), "The min Qr",
                                         500., doub()),
        r(p("fitter_lim_Qr_max"  ), "The max Qr",
                                         50000., doub()),
        r(p("fitter_lim_gain_min"  ), "The min gain",
                                         1., doub()),
        r(p("fitter_modelspec"         ), "The spec of S21 model to use",
                                         Fitter::ModelSpec::gainlintrend,
                                         list(Fitter::ModelSpec{})),
        r(p("fitter_auto_global_shift"), "Do fitting with two-passes to account for global shift."),
        r(p("output_processed"     ), "Save the reduced data",
                                         "{stem}{suffix}.{ext}", opt_str("dest"))
        ); // clang-format on
    };
    constexpr static auto fitter = [](auto &rc) {
        Fitter::Config conf{};
        for (auto &[confkey, rckey] :
             keymap_t{{"output_processed", "output_processed"},
                      {"weight_window_type", "fitter_weight_window_type"},
                      // {"weight_window_fwhm", "fitter_weight_window_fwhm"},
                      {"weight_window_Qr", "fitter_weight_window_Qr"},
                      {"lim_Qr_min", "fitter_lim_Qr_min"},
                      {"lim_Qr_max", "fitter_lim_Qr_max"},
                      {"lim_gain_min", "fitter_lim_gain_min"},
                      {"modelspec", "fitter_modelspec"},
                     {"finder_smooth_size", "finder_smooth_size"},
                     {"finder_use_savgol_deriv", "finder_use_savgol_deriv"},
                      {"exmode", "grppiex"}}) {
            if (rc.is_set(rckey)) {
                conf.set(confkey, rc.at(rckey));
            }
        }
        return Fitter{std::move(conf)};
    };
};

/// @brief Helper to construct Solver from CLI.
struct solver_options {

    using Solver = kids::TimeStreamSolver;

    constexpr static auto cli = [](auto &r) {
        using namespace tula::cli::clipp_builder;
        // clang-format off
        return "solver configs" % g(
        r(p("solver_fitreportdir"      ), "Look for fitreport file in this directory",
                                         ".", str("dir")),
        r(p("solver_fitreportfile"      ), "Use this fitreport file",
                                         undef{}, str("file")),
        r(p("solver_extra_output"     ), "Compute extra output"),
        r(p("solver_chunk_size"       ), "Solve timestream by chunk",
                                         undef{}, opt_int()),
        r(p("solver_sample_slice"     ), "Use this range of samples.",
                                         ":", str("slice"))
        ); // clang-format on
    };
    constexpr static auto solver = [](auto &rc) {
        Solver::Config conf{};
        for (auto &[confkey, rckey] : keymap_t{
                 {"fitreportdir", "solver_fitreportdir"},
                 {"fitreportfile", "solver_fitreportfile"},
                 {"extra_output", "solver_extra_output"},
                 {"exmode", "grppiex"},
             }) {
            if (rc.is_set(rckey)) {
                conf.set(confkey, rc.at(rckey));
            }
        }
        return Solver{std::move(conf)};
    };
};

auto parse_args(int argc, char *argv[]) {
    // disable logger before parse
    spdlog::set_level(spdlog::level::off);
    using namespace tula::cli::clipp_builder;

    // some of the option specs
    auto ver_str =
        fmt::format("{} ({})", KIDSCPP_GIT_VERSION, KIDSCPP_BUILD_TIMESTAMP);
    constexpr auto level_names = tula::logging::active_level_names;
    auto default_level_name = []() {
        auto v = spdlog::level::debug;
        if (v < tula::logging::active_level) {
            v = tula::logging::active_level;
        }
        return tula::logging::get_level_name(v);
    }();
    using ex_config = tula::grppi_utils::ex_config;
    // clang-format off
    auto parse = config_parser<config_t, config_t>{};
    auto screen = tula::cli::screen{
    // =======================================================================
                        "kids",  KIDSCPP_PROJECT_NAME, ver_str,
                                 KIDSCPP_PROJECT_DESCRIPTION};
    auto [cli, rc, cc] = parse([&](auto &r, auto &c) { return (
    // rc -- runtime config
    // cc -- cli config
    // =======================================================================
    c(p(           "h", "help"), "Print help information and exit"),
    c(p(             "version"), "Print version information and exit"),
    // =======================================================================
                ("server mode" % g(
    r(p(           "p", "port"), "The port to use", 55437, opt_int()))
    // -----------------------------------------------------------------------
                  | "cmd mode" % g(
    r(                "source" , "The path or uri of input data", str())),
    // =======================================================================
              "common options" % g(
    c(p(      "l", "log_level"), "Set the log level.",
                                 default_level_name, list(level_names)),
    r(p(         "o", "output"), "Output dest",
                                 "{stem}{suffix}.{ext}", opt_str("dest")),
    r(p(                "plot"), "Make diagnostic plot"),
    r(p(        "plot_backend"), "Matplotlib backend to use",
                                 "default", str()),
    r(p(         "plot_output"), "Plot output dest",
                                 "{stem}.png", opt_str("dest")),
    r(p(             "grppiex"), "GRPPI executioon policy",
                                 ex_config::default_mode(),
                                 list(ex_config::mode_names_supported()))),
                 finder_options::cli(r),
                 fitter_options::cli(r),
                solver_options::cli(r))
    // =======================================================================
    );}, screen, argc, argv);
    // clang-format on
    // SPDLOG_TRACE("cc: {}", cc.pformat());
    if (cc.get_typed<bool>("help")) {
        screen.manpage(cli);
        std::exit(EXIT_SUCCESS);
    } else if (cc.get_typed<bool>("version")) {
        screen.version();
        std::exit(EXIT_SUCCESS);
    }
    {
        auto log_level_str = cc.get_str("log_level");
        auto log_level = spdlog::level::from_str(log_level_str);
        spdlog::set_level(log_level);
        SPDLOG_INFO("reconfigure logger to level={}", log_level_str);
    }
    // pass on the runtime config
    return std::move(rc);
}

auto run_server(int port) -> int {
    SPDLOG_INFO("start server on port {}", port);
    SPDLOG_WARN("this functionality is not implemented yet");
    SPDLOG_DEBUG("this is a test debug message");
    SPDLOG_TRACE("this is a test trace message");
    std::cin.get();
    SPDLOG_INFO("server shutdown");
    return EXIT_SUCCESS;
}

/*
auto run_cmdproc(const config_t &rc) -> int {
    tula::logging::scoped_timeit TULA_X{"run_cmdproc"};
    using kids::KidsData;
    using kids::KidsDataKind;
    using tula::logging::timeit;
    // IO spec
    namespace io_spec = kids::toltec;
    SPDLOG_INFO("use data spec: {}", io_spec::name);

    try {
        // get metadata
        auto input_filepath = rc.get_str("source");
        auto [kind, meta] = io_spec::get_meta<>(input_filepath);
        SPDLOG_INFO("read kidsdata kind={} meta={}", kind, meta.pformat());
        using proc_t = kids::KidsDataProc;
        using index_t = proc_t::index_t;
        // support solve by chunk for timestreams
        if (kind & KidsDataKind::TimeStream) {
            // check solver range and pass that to the reader
            auto sample_slice = tula::container_utils::parse_slice<index_t>(
                rc.get_str("solver_sample_slice"));
            auto sample_range =
                proc_t::resolve_sample_slice(sample_slice, meta);
            const auto &[start, stop, step, size] = sample_range;
            // check size, if chunk size is larger than data size, ignore chunk
            bool solve_by_chunk{false};
            auto opt_chunk_size =
                rc.get_typed<std::optional<int>>("solver_chunk_size");
            if (opt_chunk_size.has_value()) {
                auto chunk_size = opt_chunk_size.value();
                solve_by_chunk = (chunk_size < size);
                if (solve_by_chunk) {
                    SPDLOG_INFO("solve with chunk_size={} data_size={}",
                                chunk_size, size);
                } else {
                    SPDLOG_INFO("solver chunk size ({}) ignored because it is "
                                "greater than the data size ({})",
                                chunk_size, size);
                }
            }
            if (solve_by_chunk) {
                // solve by chunk conflicts with extra output with computes psds
                if (rc.get_typed<bool>("solver_extra_output")) {
                    SPDLOG_WARN(
                        "solve by chunk cannot produce extra output, ignored");
                }
                // solve by chunk requires an output file
                if (!rc.is_set("output")) {
                    throw std::runtime_error("solve by chunk requires output");
                }
                auto output_filepath = rc.get_str("output");
                auto solver = solver_options::solver(rc);
                proc_t::solve(solver, input_filepath, output_filepath,
                              sample_slice, opt_chunk_size.value());
                return EXIT_SUCCESS;
            }
        }

    } catch (const std::runtime_error &re) {
        SPDLOG_ERROR("runtime error: {}. abort", re.what());
        return EXIT_FAILURE;
    } catch (const std::exception &e) {
        SPDLOG_ERROR("unexpected error: {}. abort", e.what());
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
*/
int run_cmdproc(const config_t &rc) {
    using kids::KidsData;
    using kids::KidsDataKind;
    using tula::logging::timeit;
    // IO spec
    namespace spec = kids::toltec;
    SPDLOG_INFO("use data spec: {}", spec::name);
    try {
        // read data
        auto [kind, meta] = spec::get_meta<>(rc.get_str("source"));
        KidsData<> kidsdata;
        if (kind & KidsDataKind::TimeStream) {
            using index_t = Eigen::Index;
            // check solver range and pass that to the reader
            auto sample_slice = tula::container_utils::parse_slice<index_t>(
                rc.get_str("solver_sample_slice"));
            auto ntimes = meta.template get_typed<int>("ntimes_all");
            auto sample_range =
                tula::container_utils::to_indices(sample_slice, ntimes);
            SPDLOG_INFO("solve range {} out of {}", sample_range, ntimes);
            using slice_t = std::decay_t<decltype(sample_slice)>;
            using range_t = std::decay_t<decltype(sample_range)>;
            auto &[start, stop, step, size] = sample_range;

            // check size, if chunk size is larger than data size, ignore chunk
            bool solve_by_chunk{false};
            if (rc.is_set("solver_chunk_size")) {
                auto chunksize = rc.get_typed<int>("solver_chunk_size");
                SPDLOG_INFO("solver chunk size: {}", chunksize);
                solve_by_chunk = (chunksize < size);
            }
            SPDLOG_INFO("solve by chunk: {}", solve_by_chunk ? "yes" : "no");
            // make chunks
            if (solve_by_chunk) {
                if (rc.get_typed<bool>("solver_extra_output")) {
                    SPDLOG_WARN(
                        "solve by chunk cannot produce extra output, ignored");
                }
                if (!rc.is_set("output")) {
                    throw std::runtime_error("solve by chunk requires output");
                }
                std::vector<range_t> chunks{};
                auto chunksize = rc.get_typed<int>("solver_chunk_size");
                if (chunksize > size) {
                    // one chunk
                    chunks.push_back(sample_range);
                } else {
                    // split to multiple chunks
                    index_t c_start = start;
                    index_t c_stop = start;
                    while (true) {
                        // make t_slice the chunk size
                        c_stop = c_start + step * chunksize;
                        if (c_stop > stop) {
                            c_stop = stop;
                        }
                        chunks.emplace_back(tula::container_utils::to_indices(
                            slice_t{c_start, c_stop, step}, ntimes));
                        c_start = c_stop;
                        assert((std::get<3>(chunks.back()) == chunksize) ||
                               (c_stop == stop));
                        if (c_start == stop) {
                            break;
                        }
                    }
                }
                auto nchunks = chunks.size();
                SPDLOG_INFO("solve by chunks size={} nchunks={}", chunksize,
                            nchunks);
                auto ex = tula::grppi_utils::dyn_ex(rc.get_str("grppiex"));
                {
                    tula::logging::scoped_timeit l0("solve by chunk");
                    tula::logging::progressbar pb0(
                        [](const auto &msg) { SPDLOG_INFO("{}", msg); }, 60,
                        "solve by chunk ");
                    using rts_t =
                        kids::KidsData<kids::KidsDataKind::RawTimeStream>;
                    using netCDF::NcFile;
                    auto filepath = tula::filename_utils::parse_pattern(
                        rc.get_str("output"), meta.get_str("source"),
                        fmt::arg("ext", "nc"),
                        fmt::arg("suffix", "_processed"));
                    SPDLOG_INFO("save result to {}", filepath);
                    auto solver = solver_options::solver(rc);
                    kids::TimeStreamSolverResult::NcFileIO io{filepath};
                    std::size_t i = 0;
                    std::mutex io_mutex;
                    grppi::pipeline(
                        ex,
                        [&]() mutable -> std::optional<rts_t> {
                            std::scoped_lock lock(io_mutex);
                            // tula::logging::scoped_loglevel<spdlog::level::debug> l0{};
                            if (i >= chunks.size()) {
                                return std::nullopt;
                            }
                            auto [start, stop, step, size] = chunks.at(i);
                            ++i;
                            slice_t slice{start, stop, step};
                            SPDLOG_INFO("read slice {}", slice);
                            return timeit(
                                "read kids data slice",
                                spec::read_data_slice<
                                    kids::KidsDataKind::RawTimeStream>,
                                rc.get_str("source"), slice);
                        },
                        [&](auto kidsdata) {
                            SPDLOG_DEBUG("kidsdata ntimes={}",
                                         kidsdata.meta.template get_typed<int>(
                                             "sample_slice_size"));
                            auto result = timeit(
                                "solve xs", solver, kidsdata,
                                config_t{{"extra_output", {false}}});
                            {
                                std::scoped_lock lock(io_mutex);
                                result.append_to_nc(io);
                            }
                            pb0.count(nchunks, nchunks / 10);
                            return EXIT_SUCCESS;
                        });
                }
                return EXIT_SUCCESS;
            } else {
                kidsdata =
                    timeit("read kids data slice", spec::read_data_slice<>,
                           rc.get_str("source"), sample_slice);
            }
        } else {
            kidsdata = timeit("read kids data all", spec::read_data<>,
                              rc.get_str("source"));
        }
        SPDLOG_TRACE("kids data: {}", kidsdata);
        // logging::scoped_loglevel<spdlog::level::info> _0{};
        // process result callback
        auto handle_result = [&rc](const auto &result) {
            if (rc.is_set("output")) {
                // logging::scoped_loglevel<spdlog::level::trace> _1{};
                result.save(rc.get_str("output"));
            }
            if (rc.is_set("plot") and rc.get_typed<bool>("plot")) {
                result.plot();
            }
            return EXIT_SUCCESS;
        };
        // dispatch and do process
        auto exitcode = std::visit(
            [&rc, &handle_result](const auto &data) {
                using Data = std::decay_t<decltype(data)>;
                constexpr auto kind = Data::kind();
                using kids::KidsDataKind;
                if constexpr (kind == KidsDataKind::VnaSweep) {
                    auto finder = finder_options::finder(rc);
                    auto result = timeit("detect kids", finder, data);
                    if (rc.is_set("output_d21")) {
                        result.save_d21(rc.get_str("output_d21"));
                    }
                    if (rc.is_set("output_processed")) {
                        result.save_processed(rc.get_str("output_processed"));
                    }
                    return handle_result(result);
                } else if constexpr (kind == KidsDataKind::TargetSweep) {
                    auto fitter = fitter_options::fitter(rc);
                    namespace eiu = tula::eigen_utils;
                    std::vector<double> fs_init = eiu::to_stdvec(data.tones());
                    if (rc.is_set("fitter_auto_global_shift") and rc.get_typed<bool>("fitter_auto_global_shift")) {
                        // run pre-fit without weight window to derive intial shift
                        SPDLOG_INFO("run prefit with config to get global shift"); 
                        auto result = timeit(
                            "prefit", fitter, data, fs_init
                            //, config_t{{"weight_window_type", {"none"}}}
                            );
                        // calcuate global shift
                        auto fin = data.tones();
                        auto ntones = fin.size();
                        auto fout = result.output.col(0);
                        auto flag_ = result.output.col(2);
                        std::vector<double> f_shifts;
                        using Flag = kids::SweepFitResult::Flag;
                        const auto accept_flag =
                            Flag::D21LargeOffset | Flag::D21NotConverged | Flag::D21QrOutOfRange | Flag::D21OutOfRange | Flag::D21FitsBetter; // | Flag::LargeOffset
                        for (Eigen::Index i = 0; i < ntones; ++i) {
                            auto flag = static_cast<Flag>(int(flag_.coeff(i)));
                            if ((flag  | accept_flag) == accept_flag ) {
                                f_shifts.push_back(fout.coeff(i) - fin.coeff(i));
                            }
                        }
                        auto iterclip = tula::alg::iterclip(
                            // use median-MAD for robust statistics
                            [](const auto &v) { return tula::alg::nanmedmad(v); },
                            [n = 5](const auto &v, const auto &c,
                                                       const auto &s) { return (v <= c + n * s) && (v >= c - n * s); },
                            5);
                        double f_shift{0.};
                        if (!f_shifts.empty()) {
                            auto [s, c, center, std]  = iterclip(tula::eigen_utils::as_eigen(f_shifts));
                            SPDLOG_DEBUG("iterclip f_shifts converged={} center={} std={}", c, center, std);
                            if (!c) {
                                SPDLOG_INFO("f_shifts iterclip not converged, set global shift to 0.");
                                f_shift = 0.;
                            } else {
                                f_shift = center;
                            }
                        }
                        fs_init = eiu::to_stdvec(fin.array() + f_shift);
                        SPDLOG_INFO("apply global shift={} computed from n_tones={} to initial fs: {}", f_shift, f_shifts.size(), fs_init);
                    }
                    auto result = timeit("fit s21 model", fitter, data, fs_init);
                    if (rc.is_set("output_processed")) {
                        result.save_processed(rc.get_str("output_processed"));
                    }
                    return handle_result(result);
                    // return EXIT_SUCCESS;
                } else if constexpr (kind == KidsDataKind::RawTimeStream) {
                    auto solver = solver_options::solver(rc);
                    auto result = timeit("solve xs", solver, data);
                    return handle_result(result);
                } else {
                    SPDLOG_TRACE(
                        "processing of data kind {} is not yet implemented",
                        kind);
                    return EXIT_FAILURE;
                }
            },
            kidsdata);
        return exitcode;
    } catch (const std::runtime_error &re) {
        SPDLOG_ERROR("{}. abort", re.what());
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
int main(int argc, char *argv[]) {
    tula::logging::init();
    // spdlog::set_error_handler([](const std::string& msg) {
    //     throw std::runtime_error(msg);
    // });
    auto rc = parse_args(argc, argv);
    SPDLOG_TRACE("rc {}", rc.pformat());
    // server mode
    if (rc.is_set("port")) {
        return run_server(rc.get_typed<int>("port"));
    }
    // cmd mode
    if (rc.is_set("source")) {
        return run_cmdproc(rc);
    }
    fmt::print("Invalid argument. Type --help for usage.\n");
}

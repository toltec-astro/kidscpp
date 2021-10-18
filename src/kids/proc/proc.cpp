#include "kids/proc/proc.h"
#include <tula/formatter/container.h>
#include <tula/logging.h>

auto kids::KidsDataProc::resolve_sample_slice(const slice_t &sample_slice,
                                              const meta_t &meta) -> range_t {
    auto ntimes = meta.template get_typed<int>("ntimes_all");
    auto sample_range = tula::container_utils::to_indices(sample_slice, ntimes);
    SPDLOG_INFO("resolved sample range {} out of {}", sample_range, ntimes);
    return sample_range;
}

void kids::KidsDataProc::solve(const TimeStreamSolver &solver,
                               const std::string &input_filepath,
                               const std::string &output_filepath,
                               const std::optional<slice_t> &slice,
                               const std::optional<index_t> &chunk_size) {
    SPDLOG_INFO("solve input timestream file {}", input_filepath);
    SPDLOG_INFO("output {}", output_filepath);

    auto [kind, meta] = io_spec::get_meta<>(input_filepath);
    SPDLOG_INFO("read kidsdata kind={} meta={}", kind, meta.pformat());


    std::vector<range_t> chunks{};
                if (chunk_size > size) {
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
                        chunks.emplace_back(container_utils::to_indices(
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
                    logging::scoped_timeit l0("solve by chunk");
                    logging::progressbar pb0(
                        [](const auto &msg) { SPDLOG_INFO("{}", msg); }, 60,
                        "solve by chunk ");
                    using rts_t =
                        kids::KidsData<kids::KidsDataKind::RawTimeStream>;
                    using netCDF::NcFile;
                    auto filepath = filename_utils::parse_pattern(
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
                            logging::scoped_loglevel<spdlog::level::debug> l0{};
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
                                config::Config{{"extra_output", {false}}});
                            {
                                std::scoped_lock lock(io_mutex);
                                result.append_to_nc(io);
                            }
                            pb0.count(nchunks, nchunks / 10);
                            return EXIT_SUCCESS;
                        });
                }

}

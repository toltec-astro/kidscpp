#pragma once

#include "kids/core/kidsdata.h"
#include "kids/sweep/finder.h"
#include "kids/sweep/fitter.h"
#include "kids/timestream/solver.h"
#include <tula/container.h>

namespace kids {

struct KidsDataProc {

    using index_t = Eigen::Index;
    using slice_t = tula::container_utils::IndexSlice;
    using range_t = tula::container_utils::BoundedSlice<index_t>;
    using meta_t = KidsData<>::meta_t;
    constexpr static auto chunk_size_default = 10000;

    /// @brief Resolve slice to range with metadata.
    static auto resolve_sample_slice(const slice_t &sample_slice,
                                     const meta_t &meta) -> range_t;

    /// @brief Resolve slice to range with metadata.
    // static auto make_chunk_ranges(const range_t& range, std::chunk_size) ->
    // range_t;

    /// @brief Solve timestream file and write to disk, optionally in chunks.
    static void
    solve(const TimeStreamSolver &solver, const std::string &input_filepath,
          const std::string &output_filepath,
          const std::optional<slice_t> &slice = std::nullopt,
          const std::optional<index_t> &chunk_size = chunk_size_default);
};

} // namespace kids

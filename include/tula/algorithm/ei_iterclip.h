#pragma once

#include <tula/meta.h>
#include <tula/eigen.h>
// #include "../eigeniter.h"
#include <tula/logging.h>

namespace tula::alg {

/**
 * @brief Algorithm to find elements from iterative thresholding.
 * @param statsfunc Function that computes [mean, dev] pair of a series.
 * @param selectfunc Function that returns true for elements to be kept. It
 * takes [elem, center, std] as the input, and return.
 * @param max_iter The max number of iterations for iterative clipping.
 * @return A callable that returns a vector of indexes for input data such that
 * data[i] satisfies \p select function.
 */
template <typename F1, typename F2>
auto iterclip(F1 &&statsfunc, F2 &&selectfunc, int max_iter = 20) {
    return [max_iter = max_iter,
            fargs = TULA_FWD_CAPTURE(statsfunc, selectfunc)](const auto &data) {
        auto &&[statsfunc, selectfunc] = fargs;
        namespace eiu = tula::eigen_utils;
        using Eigen::Index;
        // copy the data
        auto clipped = eiu::to_stdvec(data);
        // SPDLOG_TRACE("size before clip: {}", clipped.size());
        double center = 0;
        double std = 0;
        bool converged = false;
        bool selectingcenter =
            false; // check with selectfunc to determin the clip direction
        for (int i = 0; i < max_iter; ++i) {
            auto old_size = clipped.size();
            std::tie(center, std) = std::forward<decltype(statsfunc)>(statsfunc)(eiu::as_eigen(clipped));
            if (std::forward<decltype(selectfunc)>(selectfunc)(center, center, std)) {
                // if center is whithin the selecting range,
                // the clip will remove outside the selection range (selectfunc
                // is keepfunc) otherwise the clip will remove the selection
                // range (selectfunc is clipfunc)
                selectingcenter = true;
            }
            // check each element in clipped to erase or keep
            clipped.erase(std::remove_if(
                              clipped.begin(), clipped.end(),
                              [&center, &std, &selectingcenter,
                               fargs = TULA_FWD_CAPTURE(selectfunc)](const auto &v) {
                                  auto &&[selectfunc] = fargs;
                                  auto select = std::forward<decltype(selectfunc)>(selectfunc)(v, center, std);
                                  if (selectingcenter) {
                                      // clip outside of selected
                                      return !select;
                                  }
                                  // clip selected
                                  return select;
                              }),
                          clipped.end());
            if (clipped.size() == old_size) {
                // SPDLOG_TRACE("clipping coverged after {} iters", i + 1);
                converged = true;
                break; // converged
            }
        }
        // reaches max_iter
        // if (!converged) {
        //     SPDLOG_DEBUG("clip fails to converge after {} iterations",
        //                  max_iter);
        // }
        // SPDLOG_TRACE("size after clip: {}", clipped.size());
        // get the original selected indexes
        std::vector<Index> selectindex;
        selectindex.reserve(data.size());
        // auto [begin, end] = eigeniter::iters(data);
        auto begin = data.begin();
        auto end = data.end();
        auto index = 0;
        for (auto it = begin; it != end; ++it) {
            auto select = std::forward<decltype(selectfunc)>(selectfunc)(*it, center, std);
            // SPDLOG_TRACE("select {} v={}, m={} s={}", select, *it, center, std);
            if (select) {
                selectindex.push_back(index);
            }
            ++index;
        }
        return std::make_tuple(std::move(selectindex), converged, center, std);
    };
}

}  // namespace tula::alg

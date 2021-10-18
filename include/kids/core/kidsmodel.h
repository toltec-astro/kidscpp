#pragma once

#include <Eigen/Core>
#include <tula/enum.h>
#include <tula/meta.h>

namespace kids {

/// @brief Kids model base class.
template <typename Derived, typename ModelImpl>
struct KidsModel {
    static_assert(
        requires { ModelImpl::solve; } || requires { ModelImpl::probe; },
        "MODEL CLASS DO NOT HAVE SOLVE/PROBE METHODS");

    template <typename... Args>
    static auto solve(Args... args) {
        if constexpr (requires { ModelImpl::solve; }) {
            return ModelImpl::solve(std::forward<decltype(args)>(args)...);
        } else {
            throw KidsModelError::NotImplemented(
                fmt::format("{}::{}", Derived::name, "solve"));
        }
    }
    template <typename... Args>
    static auto probe(Args... args) {
        if constexpr (requires { ModelImpl::probe; }) {
            return ModelImpl::probe(std::forward<decltype(args)>(args)...);
        } else {
            throw KidsModelError::NotImplemented(
                fmt::format("{}::{}", Derived::name, "probe"));
        }
    }
};

/// @brief Kids model impl classes
namespace models {

using Eigen::Dynamic;
using Eigen::Index;

template <auto param, Index param_offset,
          REQUIRES_(std::is_enum<decltype(param)>)>
struct Parameter {
    constexpr static auto offset = param_offset;
    constexpr static auto name = enum_utils::name(param);
    template <typename T>
    static auto value(const T *const params) -> decltype(auto) {
        return params[offset];
    }
};

template <typename... Parameters>
struct Model {
    constexpr static auto n_params = sizeof...(Parameters);
}

/// @brief Model to transform between (I, Q) and (r, x)
template <>
struct ResonanceCircleTransform2 : Model<> {
    constexpr static auto n_inputs = 2;
    constexpr static auto n_outputs = 2;

    /// @tparam T Used here instead double for autodiff.
    /// @param params The parameter values.
    /// @param result The output data. Interleaved real and imaginary.
    /// @param data The input data. Interleaved real and imaginary.
    template <Index n_data, typename T = double>
    static auto evaluate(const T *const params, T *result, double const *data) {
        constexpr Index I = 0;
        constexpr Index Q = 1;
        constexpr Index r = 0;
        constexpr Index x = 1;
        // (I + jQ) = 0.5 / (r + jx)
        // the inverse is the same
        auto factor = 0.5 / (data[r] * data[r] + data[x] * data[x]);
        result[I] = data[r] * factor;
        result[Q] = -data[x] * factor;
        return;
    }
};

/// @brief Model to transform bwtween Qr and r
template <>
struct Qr2r : Model<> {
    constexpr static auto n_inputs = 1;
    constexpr static auto n_outputs = 1;

    /// @tparam T Used here instead double for autodiff.
    /// @param params The parameter values.
    /// @param result The output data.
    /// @param data The input data.
    template <Index n_data, typename T = double>
    static auto evaluate(const T *const params, T *result, double const *data) {
        constexpr Index Qr = 0;
        constexpr Index r = 0;
        // r = 0.5 / Qr
        // the inverse is the same
        result[r] = 0.5 / data[Qr];
        return;
    }
};

/// @brief Model to compute x from freqs.
template <typename fr>
struct InstrumentalDetune : Model<fr> {
    constexpr static auto n_inputs = 1;  // fp
    constexpr static auto n_outputs = 1; // x

    /// @tparam T Used here instead double for autodiff.
    /// @param params The parameter values.
    /// @param result The output data.
    /// @param data The input data.
    template <Index n_data, typename T = double>
    static auto evaluate(const T *const params, T *result, double const *data) {
        constexpr Index fp = 0;
        constexpr Index x = 0;
        // x = (fp / fr) - 1
        result[x] = data[fp] / fr::value(params) - 1.;
        return;
    }
};

/// @brief Model to transform resonance circle (I, Q) composed of a gain factor
/// and a linear baseline.
template <typename g0, typename g1, typename f0, typename k0, typename k1,
          typename m0, typename m1>
struct ReadoutGainWithLinTrend
    : Model<typename g0, typename g1, typename f0, typename k0, typename k1,
            typename m0, typename m1> {
    constexpr static auto n_inputs = 1;  // fp
    constexpr static auto n_outputs = 2; // I, Q

    /// @tparam T Used here instead double for autodiff.
    /// @param params The parameter values.
    /// @param result The output data.
    /// @param data The input data.
    template <Index n_data, typename T = double>
    static auto evaluate(const T *const params, T *result, double const *data) {
        constexpr Index fp = 0;
        constexpr Index I = 0;
        constexpr Index Q = 1;
        // complex gain
        auto I1 = g0::value(params) * result[I] - g1::value(params) * data[Q];
        auto Q1 = g0::value(params) * result[Q] + g1::value(params) * data[I];
        // add linear term
        result[I] = I1 + (data[fp] - f0::value(params)) * result[i] return;
    }

f, g0, g1, f0, k0, k1, m0, m1):
        gg = g0 + 1.j * g1
        kk = k0 + 1.j * k1
        mm = m0 + 1.j * m1
        return gg * S + kk * (f - f0) + mm
};

} // namespace models

} // namespace kids

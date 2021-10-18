#pragma once

#include <tula/algorithm/ei_ceresfitter.h>
#include <tula/enum.h>

namespace kids {

// clang-format off
TULA_ENUM(SweepModel, int,
          D21,
          S21Basic,
          S21WithGain, S21WithGainLinTrend, S21WithTrans, S21WithTransLinTrend);
// clang-format on

namespace internal {

template <SweepModel model>
struct sweepmodel_t {
    using Self = sweepmodel_t<model>;
    // model implemetation
    using Index = Eigen::Index;
    // define data
    Index nx = 0;                  // input_size
    Index ny = 0;                  // output_size
    double const *xdata = nullptr; //
    double const *ydata = nullptr; // interleaved complex
    double const *yerr = nullptr;  // same as ydata

    constexpr static auto param_names = []() {
        constexpr auto names = [](auto &&...ns) {
            return std::array<std::string_view, sizeof...(ns)>{ns...};
        };
        if constexpr (model == SweepModel::D21) {
            return names("fp", "Qr", "fr", "normI_d21", "normQ_d21", "slopeI",
                         "slopeQ");
            // return names("fp", "Qr", "Qc" "fr", "A");
        } else if constexpr (model == SweepModel::S21Basic) {
            return names("fp", "Qr", "Qc", "fr", "A");
        } else if constexpr (model == SweepModel::S21WithGain) {
            return names("fp", "Qr", "Qc", "fr", "A", "normI", "normQ");
        } else if constexpr (model == SweepModel::S21WithGainLinTrend) {
            return names("fp", "Qr", "Qc", "fr", "A", "normI", "normQ",
                         "slopeI", "slopeQ", "interceptI", "interceptQ");
        } else if constexpr (model == SweepModel::S21WithTrans) {
            return names("fp", "Qr", "Qc", "fr", "A", "shiftI", "shiftQ",
                         "rot");
        } else if constexpr (model == SweepModel::S21WithTransLinTrend) {
            return names("fp", "Qr", "Qc", "fr", "A", "shiftI", "shiftQ", "rot",
                         "slopeI", "slopeQ");
        }
    }();

    using ParamSetting = tula::alg::ceresfit::ParamSetting<double>;
    inline static tula::alg::ceresfit::ParamSettings<> param_settings{
        {"rot", ParamSetting::getbounded(static_cast<double>(-EIGEN_PI),
                                         static_cast<double>(EIGEN_PI))},
        {"normI_d21", ParamSetting{1.}},
        {"normQ_d21", ParamSetting{0.}},
        {"normI", ParamSetting{1.}},
        {"normQ", ParamSetting{0.}},
        {"shiftI", ParamSetting{0.}},
        {"shiftQ", ParamSetting{0.}},
    };

    template <typename T, auto model_ = model>
              requires (model_ == SweepModel::D21)
    inline static void _eval(const T *const p, T *r, Index nx,
                             double const *xdata) {
        // 0   1  2   3     4     5  6
        // fp Qr fr normI normQ  slopeI slopeQ
        for (Index i = 0; i < nx; ++i) {
            auto f = xdata[i];          // f
            auto x = (f - p[2]) / p[2]; // x
            // denom = 1.+2.j*Qr*x
            // den_r = 1.
            auto den_i = 2. * p[1] * x;
            // s21 = 1 / denom
            auto den2 = 1. + den_i * den_i;
            auto mdl_r0 = 1. / den2;
            auto mdl_i0 = -den_i / den2;
            // d21 = g * s21^2
            auto mdl_r1 = mdl_r0 * mdl_r0 - mdl_i0 * mdl_i0;
            auto mdl_i1 = 2. * mdl_r0 * mdl_i0;
            r[2 * i] = p[3] * mdl_r1 - p[4] * mdl_i1 + p[5];
            r[2 * i + 1] = p[3] * mdl_i1 + p[4] * mdl_r1 + p[6];
        }
    }

    // use T for auto diff residual
    // here r is is layouted as [real, imag, real, imag, ...]
    template <typename T, auto model_ = model> requires (model_ != SweepModel::D21)
    inline static void _eval(const T *const p, T *r, Index nx,
                             double const *xdata) {
        // p: 0    1    2    3    4
        //    fp   Qr   Qc   fr   A
        // it seems that we have to separate the real and imag part for
        // autodiff to work
        // Eq 3.60
        // http://mnemosyne.phys.columbia.edu/~bjohnson/web/files/thesis/flanigan_thesis_2018.pdf
        // numerator = Qr/Qc * (1.+1.j*A)
        auto num_r = p[1] / p[2];
        auto num_i = p[1] / p[2] * p[4];
        for (Index i = 0; i < nx; ++i) {
            auto f = xdata[i];          // f
            auto x = (f - p[3]) / p[3]; // x
            // denom = 1.+2.j*Qr*x
            // den_r = 1.
            auto den_i = 2. * p[1] * x;
            // s21 = numerator / denom
            auto den2 = 1. + den_i * den_i;
            auto mdl_r0 = (num_r + num_i * den_i) / den2;
            auto mdl_i0 = (num_i - num_r * den_i) / den2;
            if constexpr (model == SweepModel::S21Basic) {
                // no extra step
                r[2 * i] = mdl_r0;
                r[2 * i + 1] = mdl_i0;
                continue;
            } else if constexpr ((model == SweepModel::S21WithGain) ||
                                 (model == SweepModel::S21WithGainLinTrend)) {
                // complex gain * s21
                // 5      6      7      8      9        10
                // normI normQ slopeI slopeQ interceptI interceptQ
                // auto mdl_r1 = p[5] * (1. - mdl_r0) - p[6] * (-mdl_i0);
                // auto mdl_i1 = p[5] * (-mdl_i0) + p[6] * (1. - mdl_r0);
                auto mdl_r1 = p[5] * mdl_r0 - p[6] * mdl_i0;
                auto mdl_i1 = p[5] * mdl_i0 + p[6] * mdl_r0;
                if constexpr (model == SweepModel::S21WithGainLinTrend) {
                    // lintrend slope * (f - fp) + intercept
                    r[2 * i] = mdl_r1 + p[7] * (f - p[0]) + p[9];
                    r[2 * i + 1] = mdl_i1 + p[8] * (f - p[0]) + p[10];
                    continue;
                } else {
                    r[2 * i] = mdl_r1;
                    r[2 * i + 1] = mdl_i1;
                    continue;
                }
            } else if constexpr ((model == SweepModel::S21WithTrans) ||
                                 (model == SweepModel::S21WithTransLinTrend)) {
                // model = shift - rot (s21)
                // 5       6      7
                // shiftI shiftQ rot
                // auto rot = atan(p[7]) * 2.;
                auto mdl_r1 =
                    p[5] + (cos(p[7]) * (mdl_r0)-sin(p[7]) * (mdl_i0));
                auto mdl_i1 =
                    p[6] + (sin(p[7]) * (mdl_r0) + cos(p[7]) * (mdl_i0));
                if constexpr (model == SweepModel::S21WithTransLinTrend) {
                    // lintrend slope * (f - fp)
                    // 8       9
                    // slopeI slopeQ
                    r[2 * i] = mdl_r1 + p[8] * (f - p[0]);
                    r[2 * i + 1] = mdl_i1 + p[9] * (f - p[0]);
                    continue;
                } else {
                    r[2 * i] = mdl_r1;
                    r[2 * i + 1] = mdl_i1;
                    continue;
                }
            }
        }
    }
    template <typename T>
    inline void eval(const T *const p, T *r) const {
        Self::_eval(p, r, this->nx, this->xdata);
    }
};

template <SweepModel model, typename Model = sweepmodel_t<model>,
          typename Fitter =
              tula::alg::ceresfit::Fitter<Model::param_names.size(), 1, 1, double>>
struct fitter_t : public Model, Fitter {
    using Index = Eigen::Index;
    // residual
    template <typename T>
    inline bool operator()(const T *const p, T *r) const {
        // compute model
        Model::_eval(p, r, this->nx, this->xdata);
        // compute residule
        for (Index i = 0; i < this->nx; ++i) {
            // compute residual
            auto y_r = this->ydata[2 * i];     // S21_r
            auto y_i = this->ydata[2 * i + 1]; // S21_i
            auto e_r = this->yerr[2 * i];      // e_S21_r
            auto e_i = this->yerr[2 * i + 1];  // e_S21_i
            r[2 * i] = (y_r - r[2 * i]) / e_r;
            r[2 * i + 1] = (y_i - r[2 * i + 1]) / e_i;
        }
        return true;
    }
};

} // namespace internal

/// @brief Model class for \p model.
template <SweepModel model_>
struct Model : internal::fitter_t<model_> {
    static constexpr auto model = model_;
};

} // namespace kids

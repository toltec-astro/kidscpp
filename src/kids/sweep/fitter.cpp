#include "kids/sweep/fitter.h"
#include "kids/sweep/model.h"
#include <tula/algorithm/ei_ceresfitter.h>
#include <tula/algorithm/ei_numdiff.h>
#include <tula/algorithm/ei_stats.h>
#include <tula/algorithm/ei_convolve.h>
#include <tula/container.h>
#include <tula/eigen.h>
#include <tula/formatter/enum.h>
#include <tula/grppi.h>
#include <tula/switch_invoke.h>
#include <gram_savitzky_golay/gram_savitzky_golay.h>

using Eigen::Index;
using kids::SweepFitResult;
using kids::SweepFitter;

SweepFitter::SweepFitter(const Config &config_) : config(config_) {}

SweepFitResult SweepFitter::operator()(const TargetSweepData &data,
                                       const std::vector<double> &fs_init,
                                       const Config &config_) {
    // merge the config object
    Config config{this->config};

    config.update(config_);

    SPDLOG_DEBUG("run kids fitter with config: {}", config.pformat());
    auto ex = tula::grppi_utils::dyn_ex(config.get_str("exmode"));

    // data geometry
    const auto &sfs = data.sweeps();
    Index nsweeps{0};
    Index ntones{0};
    double sweepstep{0};
    Eigen::MatrixXd toneranges;
    {
        std::tie(nsweeps, ntones) = tula::alg::shape(data.iqs());
        assert(nsweeps == data.sweeps().size());
        assert(ntones == data.tones().size());
        [[maybe_unused]] auto sweepspan = tula::alg::span(sfs);
        sweepstep = tula::alg::step(sfs);
        [[maybe_unused]] auto tonespan = tula::alg::span(data.tones());

        toneranges.resize(2, ntones);
        toneranges.row(0) = data.fs().colwise().minCoeff();
        toneranges.row(1) = data.fs().colwise().maxCoeff();
        SPDLOG_TRACE(
            "toneranges={} tonespan={:g} sweepspan={:g} sweepstep={:g}",
            toneranges, tonespan, sweepspan, sweepstep);
        SPDLOG_TRACE("fit iqs{}", data.iqs());
    }
    auto toneindices = tula::container_utils::index(ntones);

    // initial values
    Eigen::VectorXd fs0;
    if (fs_init.empty()) {
        fs0 = data.tones();
    } else {
        fs0 = tula::eigen_utils::as_eigen(fs_init);
    }
    // intial flag
    using Flag = SweepFitResult::Flag;
    const auto flag0 = bitmask::bitmask{Flag::Good};
    // flags container
    auto flag = std::vector(TULA_SIZET(ntones), flag0);
    // helper to update flag
    auto update_flag = [&flag](auto f, auto c) {
        flag[c] |= f;
        SPDLOG_TRACE("add flag {} to tone #{}", f, c);
    };

    // handle weight
    WeightOption weight_option{WeightOption::lorentz};
    if (auto type =
            WeightOption_meta::from_name(config.get_str("weight_window_type"));
        type.has_value()) {
        weight_option = type.value().value;
    } else {
        SPDLOG_WARN("unknown weight window type {}, use {} instead",
                    config.get_str("weight_window_type"), weight_option);
    }
    using WeightTypes = tula::meta::cases<WeightOption::boxcar, WeightOption::gauss,
                                    WeightOption::lorentz>;
    SPDLOG_DEBUG("use window type {}", weight_option);
    // auto window_width = config.get_typed<double>("weight_window_fwhm");
    auto window_Qr = config.get_typed<double>("weight_window_Qr");
    auto window_width_at_500MHz = 5e8 / window_Qr;
    SPDLOG_DEBUG("use window width from Qr={} ({} Hz at 500MHz)", window_Qr,
                 window_width_at_500MHz);
    // TODO: make it a config entry?
    // auto large_offset_at_500MHz = 6 * window_width_at_500MHz;
    // double large_offset_at_500MHz = 50e3; // 2x fwhm with Qr=20k
    double large_offset_at_500MHz = window_width_at_500MHz * 6;
    SPDLOG_DEBUG("use large offset limit {} Hz at 500MHz",
                 large_offset_at_500MHz);
    auto lim_Qr_min = config.get_typed<double>("lim_Qr_min");
    auto lim_Qr_max = config.get_typed<double>("lim_Qr_max");
    SPDLOG_DEBUG("use Qr lims {} {}", lim_Qr_min, lim_Qr_max);
    auto lim_gain_min = config.get_typed<double>("lim_gain_min");
    SPDLOG_DEBUG("use gain min={}", lim_gain_min);
    Eigen::MatrixXcd uncertainty(nsweeps, ntones);
    uncertainty.setConstant(std::complex<double>{1., 1.});
    auto update_uncertainty = [&](auto weight_type_, auto c, auto f_c) {
        constexpr auto weight_type = std::decay_t<decltype(weight_type_)>::value;
        auto window_width = f_c / window_Qr;
        auto f_l = f_c - window_width / 2.;
        auto f_u = f_c + window_width / 2.;
        if constexpr (weight_type == WeightOption::boxcar) {
            uncertainty.col(c) =
                ((data.fs().col(c).array() >= f_l) &&
                 (data.fs().col(c).array() < f_u))
                    .select(
                        uncertainty.col(c),
                        std::complex(std::numeric_limits<double>::infinity(),
                                     std::numeric_limits<double>::infinity()));

            return;
        }
        auto x = 2. / window_width * (data.fs().col(c).array() - f_c);
        SPDLOG_TRACE("x{}", x);
        if constexpr (weight_type == WeightOption::gauss) {
            uncertainty.col(c).real().array() = (x.square() / 2.).exp();
        } else if constexpr (weight_type == WeightOption::lorentz) {
            uncertainty.col(c).real().array() = (1. + x.square()).sqrt();
        }
        uncertainty.col(c).imag() = uncertainty.col(c).real();
    };
    // fit params
    using ParamSetting = tula::alg::ceresfit::ParamSetting<double>;

    // setup inital parameters for s21
    double qr = window_Qr;
    double qc = 5e4;
    using S21 = kids::Model<kids::SweepModel::S21Basic>;
    Eigen::MatrixXd s21ps(S21::NP, ntones);
    s21ps.setZero(); // initialize all other params to zero except the below
    s21ps.row(0) = data.tones();  // fp, these are always fixed
    s21ps.row(1).setConstant(qr); // Qr
    s21ps.row(2).setConstant(qc); // Qc
    s21ps.row(3) = fs0;           // initial freq position
    s21ps.row(4).setConstant(0.); // A is degenerated, fix at 0
    SPDLOG_TRACE("initial s21 params{}", s21ps);
    // update all uncertainties centered at the initial position
    tula::meta::switch_invoke<WeightTypes>(
        [&](auto type) {
            grppi::map(ex, toneindices, toneindices, [&](auto c) {
                update_uncertainty(type, c, s21ps.coeff(3, c));
                return c;
            });
        },
        weight_option);
    // helper to update cost
    auto update_cost = [](const auto &costs_, const auto &summary) {
        auto &costs = const_cast<std::decay_t<decltype(costs_)> &>(costs_);
        costs.coeffRef(0) = summary.initial_cost;
        costs.coeffRef(1) = summary.final_cost;
    };
    // create d21 and fit d21 to update s21ps
    Eigen::MatrixXcd diqs(nsweeps, ntones);
    diqs.setConstant(std::complex(0., 0.));

    using D21 = kids::Model<kids::SweepModel::D21>;
    Eigen::MatrixXd d21costs(2, ntones);
    Eigen::MatrixXd d21ps(D21::NP, ntones);
    d21ps.setZero(); // initialize all other params to zero except the below
    d21ps.row(0) = s21ps.row(0); // fp, these are always fixed
    d21ps.row(1) = s21ps.row(1); // Qr
    d21ps.row(2) = s21ps.row(3); // input freqs
    // row 3 an 4 are norm_d21
    // 5 and 6 and const slope in s21
    
    // d21 related
    bool use_savgol_deriv = config.get_typed<bool>("finder_use_savgol_deriv");
    auto smooth_size = static_cast<Index>(config.get_typed<int>("finder_smooth_size"));
    if (use_savgol_deriv) {
        if (smooth_size < 5) {
            throw std::runtime_error("smooth size too small for SavGol deriv");
        }
        if (smooth_size % 2 == 0) {
            smooth_size += 1;
        }
    }
    auto smooth_i0 = (smooth_size - 1) / 2; // valid index
    bool pre_smooth = (smooth_size > 0);
    gram_sg::SavitzkyGolayFilter first_derivative_filter(
        smooth_i0, 0, 2, 1, sweepstep);

    grppi::map(ex, toneindices, toneindices, [&](auto c) {
        // tula::alg::gradient(data.iqs().col(c), data.fs().col(c), diqs.col(c));
            // smooth and deriv
            if  (use_savgol_deriv) {
                //  we ensured the smooth size >=5 previously
                // here we just write the deriv data to diqs
                tula::alg::convolve1d(
                    data.iqs().col(c),
                    [&data, &use_savgol_deriv, &first_derivative_filter](const auto &patch) {
                        // calculate the deriv
                        auto xx = tula::eigen_utils::to_stdvec(patch.real());
                        auto yy = tula::eigen_utils::to_stdvec(patch.imag());
                        double x = first_derivative_filter.filter(xx);
                        double y = first_derivative_filter.filter(yy);
                        return std::complex(x, y);
                    },
                    smooth_size,
                    diqs.col(c).segment(smooth_i0, nsweeps - smooth_size + 1));
            } else {
                // do smooth if needed
                Eigen::VectorXcd iqs(nsweeps);
                iqs.setConstant(std::complex(0., 0.));
               if (pre_smooth) {
                    tula::alg::convolve1d(
                        data.iqs().col(c),
                        [&data, &use_savgol_deriv, &first_derivative_filter](const auto &patch) {
                            if (use_savgol_deriv) {
                                // calculate the deriv
                                auto xx = tula::eigen_utils::to_stdvec(patch.real());
                                auto yy = tula::eigen_utils::to_stdvec(patch.imag());
                                double x = first_derivative_filter.filter(xx);
                                double y = first_derivative_filter.filter(yy);
                                return std::complex(x, y);
                            } else {
                                return std::complex(tula::alg::median(patch.real()),
                                                tula::alg::median(patch.imag()));
                            }
                        },
                        smooth_size,
                        iqs.segment(smooth_i0, nsweeps - smooth_size + 1));
                }
                // do deriv
                tula::alg::gradient(iqs, data.fs().col(c), diqs.col(c));
            }
        auto [converged, summary] = tula::alg::ceresfit::fit<D21>(
            data.fs().col(c), diqs.col(c), uncertainty.col(c), d21ps.col(c),
            {{"fp", ParamSetting::getfixed(d21ps.coeff(0, c))}});
        update_cost(d21costs.col(c), summary);
        SPDLOG_TRACE("tone #{}: fitted d21 param: {}", c, d21ps.col(c));
        // check fit convergence
        if (!converged) {
            update_flag(Flag::D21NotConverged, c);
        }
        // check large offset
        auto f0 = fs0(c);
        auto f = d21ps(2, c);
        if ((f < toneranges.coeff(0, c)) || (f > toneranges.coeff(1, c))) {
            update_flag(Flag::D21OutOfRange, c);
        }
        auto large_offset = large_offset_at_500MHz * f0 / 5e8;
        if (std::abs(f - f0) > large_offset) {
            update_flag(Flag::D21LargeOffset, c);
        }
        // check D21 Qr range
        auto d21_Qr = d21ps(1, c);
        if ((d21_Qr < lim_Qr_min) || (d21_Qr > lim_Qr_max)) {
            update_flag(Flag::D21QrOutOfRange, c);
        }

        return c;
    });
    ModelSpec modelspec{ModelSpec::gainlintrend};
    if (auto spec = ModelSpec_meta::from_name(config.get_str("modelspec"));
        spec.has_value()) {
        modelspec = spec.value().value;
    } else {
        SPDLOG_WARN("unknown model spec {}, use {} instead",
                    config.get_str("modelspec"), modelspec);
    }
    using ModelSpecs = tula::meta::cases<ModelSpec::gain, ModelSpec::gainlintrend,
                                   ModelSpec::trans, ModelSpec::translintrend>;
    // dispatch different models
    auto fitresult = tula::meta::switch_invoke<ModelSpecs>(
        [&](auto spec_) {
            constexpr auto spec = std::decay_t<decltype(spec_)>::value;
            using Model = model_t<spec>;
            SPDLOG_DEBUG("fit with model: {}", Model::model);
            Eigen::MatrixXd ps(Model::NP, ntones);
            Eigen::MatrixXd costs(2, ntones);
            // copy over basic s21 params Qr Qc fr and A
            ps.setZero();
            ps.topRows(s21ps.rows()) = s21ps;
            grppi::map(ex, toneindices, toneindices, [&](auto c) {
                // update other d21 params if flag is good
                if (!(flag[c] & Flag::D21OutOfRange)) {
                    auto qr_d21 = d21ps.coeff(1, c);
                    auto fr_d21 = d21ps.coeff(2, c);
                    auto nr_d21 = d21ps.coeff(3, c);
                    auto ni_d21 = d21ps.coeff(4, c);
                    auto sr_d21 = d21ps.coeff(5, c);
                    auto si_d21 = d21ps.coeff(6, c);
                    // update
                    ps.coeffRef(1, c) = qr_d21;
                    // (2) Qc is set from s21ps
                    ps.coeffRef(3, c) = fr_d21;
                    // (4) A is set from s21ps
                    if constexpr ((Model::model == SweepModel::S21WithGain) ||
                                  (Model::model ==
                                   SweepModel::S21WithGainLinTrend)) {
                        // (5, 6) norm
                        ps.coeffRef(5, c) =
                            -ni_d21 * fr_d21 / 2. / qr_d21 / qr_d21;
                        ps.coeffRef(6, c) =
                            +nr_d21 * fr_d21 / 2. / qr_d21 / qr_d21;
                    }
                    if constexpr (Model::model ==
                                  SweepModel::S21WithGainLinTrend) {
                        // (7, 8) slope
                        ps.coeffRef(7, c) = sr_d21;
                        ps.coeffRef(8, c) = si_d21;
                    }
                    // update uncertainty to new fr
                    if (!(flag[c] & Flag::D21LargeOffset)) {
                        tula::meta::switch_invoke<WeightTypes>(
                            [&](auto type) {
                                update_uncertainty(type, c, fr_d21);
                            },
                            weight_option);
                    }
                } else {
                    // out of range d21, we have to proceed with some guess
                    // Qr, Qc, fr, A use s21ps value
                    if constexpr ((Model::model == SweepModel::S21WithGain) ||
                                  (Model::model ==
                                   SweepModel::S21WithGainLinTrend)) {
                        // (5, 6) norm
                        ps.coeffRef(5, c) = 1.;
                        ps.coeffRef(6, c) = 0.;
                    }
                    if constexpr (Model::model ==
                                  SweepModel::S21WithGainLinTrend) {
                        // (7, 8) slope
                        ps.coeffRef(7, c) = 0.;
                        ps.coeffRef(8, c) = 0.;
                    }
                }
                auto [converged0, summary0] = tula::alg::ceresfit::fit<Model>(
                    data.fs().col(c), data.iqs().col(c), uncertainty.col(c),
                    ps.col(c),
                    {{"Qr", ParamSetting::getfixed(ps.coeff(1, c))},
                     {"Qc", ParamSetting::getfixed(1.)},
                     // {"fr", ParamSetting::getfixed(ps.coeff(3, c))},
                     {"A", ParamSetting::getfixed(0.)},
                     {"normI", ParamSetting{ps.coeff(5, c)}},
                     {"normQ", ParamSetting{ps.coeff(6, c)}},
                     {"slopeI", ParamSetting{ps.coeff(7, c)}},
                     {"slopeQ", ParamSetting{ps.coeff(8, c)}},
                     // {"interceptI",
                     // ParamSetting{data.iqs.col(c).real().minCoeff()}},
                     // {"interceptQ",
                     // ParamSetting{data.iqs.col(c).imag().minCoeff()}},
                     {"fp", ParamSetting::getfixed(ps.coeff(0, c))}});
                // 2nd pass with all parameters free
                auto [converged, summary] = tula::alg::ceresfit::fit<Model>(
                    data.fs().col(c), data.iqs().col(c), uncertainty.col(c),
                    ps.col(c),
                    {{"Qc", tula::alg::ceresfit::ParamSetting<double>::getfixed(1.)},
                     {"A", ParamSetting::getfixed(0.)},
                     {"fp", ParamSetting::getfixed(ps.coeff(0, c))}});
                // check if s21 does not improved the cost
                // here we could use a tolerance to allow equally good fit
                // 1.1 is chosen arbituraly
                if (summary.final_cost > 1.1 * summary0.final_cost) {
                    update_flag(Flag::D21FitsBetter, c);
                    // use the first pass fit result
                    converged = converged0;
                    summary = std::move(summary0);
                }
                update_cost(costs.col(c), summary);
                // check fit convergence
                if (!converged) {
                    update_flag(Flag::NotConverged, c);
                }
                auto f0 = fs0(c);
                auto f = ps(3, c);
                if ((f < toneranges.coeff(0, c)) ||
                    (f > toneranges.coeff(1, c))) {
                    update_flag(Flag::OutOfRange, c);
                }
                auto large_offset = large_offset_at_500MHz * f0 / 5e8;
                if (std::abs(f - f0) > large_offset) {
                    update_flag(Flag::LargeOffset, c);
                }
                // check S21 fitting Qr range
                auto Qr = ps(1, c);
                if ((Qr < lim_Qr_min) || (Qr > lim_Qr_max)) {
                    update_flag(Flag::QrOutOfRange, c);
                }
                // check S21 gain
                auto gain = std::sqrt(ps(5, c) * ps(5, c) + ps(6, c) * ps(6, c));
                if (gain < lim_gain_min) {
                    update_flag(Flag::LowGain, c);
                }
                return c;
            });
            SPDLOG_TRACE("best-fit ps{}", ps);
            // fill in outputs
            std::vector<std::string> colnames = {
                "f_out",
                "f_in",
                "flag",
            };
            colnames.insert(colnames.end(), Model::param_names.begin(),
                            Model::param_names.end());
            Eigen::MatrixXd output(ntones, colnames.size());
            // 0:f_out <- fr
            output.col(0) = ps.row(3);
            // 1:f_in  <- fs0
            output.col(1) = fs0; // output_1: f_in; param_2: fr
            // 2:flag <- flag
            for (Index i = 0; i < ntones; ++i) {
                auto f = flag[TULA_SIZET(i)];
                output.coeffRef(i, 2) = f.bits();
                // replace failed fit with f_out = f_in
                if (f & (Flag::OutOfRange | Flag::LargeOffset)) {
                    output.coeffRef(i, 0) = output.coeffRef(i, 1);
                }
            }
            // model parameters
            output.rightCols(Model::NP) = ps.transpose();
            // print summary
            auto flag_summary = [&flag](auto &&flags) {
                auto summary = std::unordered_map<Flag, int>{};
                for (auto f : flags) {
                    summary[f] = 0;
                }
                for (auto &[f, n] : summary) {
                    // SPDLOG_TRACE("count flag: {}", f);
                    for (auto f_i : flag) {
                        // SPDLOG_TRACE("check flag: {}", f_i);
                        if (((f == Flag::Good) && (f_i == Flag::Good)) ||
                            (f & f_i)) {
                            ++n;
                        }
                    }
                }
                return summary;
            }(tula::enum_utils::values<Flag>());
            std::stringstream flag_report;
            flag_report << fmt::format("fit result of {} tones:", ntones);
            auto fmt_entry = []() {
                auto names = tula::enum_utils::names<Flag>();
                auto max_key_width = 0;
                std::for_each(names.begin(), names.end(),
                              [&max_key_width](auto &it) mutable {
                                  if (auto size = it.size();
                                      size > max_key_width) {
                                      max_key_width = size;
                                  }
                              });
                return fmt::format("\n {{:{}s}}: {{}}", max_key_width);

            }();
            for (auto &[f, n] : flag_summary) {
                auto fstr = fmt::format("{}", f);
                flag_report << fmt::format(fmt::runtime(fmt_entry), fstr, n);
            }
            SPDLOG_INFO(flag_report.str().c_str());
            return SweepFitResult{data,
                                  window_Qr,
                                  std::move(uncertainty),
                                  std::move(diqs),
                                  std::move(d21ps),
                                  std::move(d21costs),
                                  modelspec,
                                  std::move(ps),
                                  std::move(costs),
                                  std::move(output),
                                  std::move(colnames),
                                  std::move(flag_summary)};
        },
        modelspec);
    if (fitresult) {
        return fitresult.value();
    }
    throw std::runtime_error(fmt::format("unknown model spec {}", modelspec));
}

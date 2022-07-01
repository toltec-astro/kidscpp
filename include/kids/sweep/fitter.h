#pragma once

#include "../core/kidsdata.h"
#include "model.h"
#include <tula/container.h>
#include <tula/enum.h>
#include <tula/switch_invoke.h>

namespace kids {

struct SweepFitResult;

struct SweepFitter {

    using TargetSweepData = KidsData<KidsDataKind::TargetSweep>;
    using Config = KidsData<>::meta_t;

    // clang-format off
    TULA_ENUM_DECL(WeightOption, int,
              boxcar,
              gauss,
              lorentz
              );
    TULA_ENUM_DECL(FittingOption, int,
              default_,
              fixfr
              );
    TULA_ENUM_DECL(ModelSpec, int,
              gain,
              gainlintrend,
              trans,
              translintrend);
    // clang-format on

    SweepFitter(const Config &config_);

    //     /// @brief fit with all tones assuming tone center
    //     SweepFitResult operator()(const TargetSweepData &,
    //                               const Config &config = {});

    /// @brief fit with given input position
    auto operator()(const TargetSweepData &,
                    const std::vector<double> &inputs = {},
                    const Config &config = {}) -> SweepFitResult;

    Config config;

    template <ModelSpec spec>
    using model_t = tula::meta::switch_t<
        spec,
        tula::meta::case_t<ModelSpec::gain,
                           kids::Model<SweepModel::S21WithGain>>,
        tula::meta::case_t<ModelSpec::gainlintrend,
                           kids::Model<SweepModel::S21WithGainLinTrend>>,
        tula::meta::case_t<ModelSpec::trans,
                           kids::Model<SweepModel::S21WithTrans>>,
        tula::meta::case_t<ModelSpec::translintrend,
                           kids::Model<SweepModel::S21WithTransLinTrend>>>;
};

struct SweepFitResult {
    // clang-format off
    TULA_ENUM_DECL(Flag, int,
              Good            = 0     ,
              D21FitsBetter   = 1 << 0,
              D21LargeOffset  = 1 << 1,
              D21NotConverged = 1 << 2,
              D21OutOfRange   = 1 << 3,
              D21QrOutOfRange = 1 << 4,
              LargeOffset     = 1 << 5,
              NotConverged    = 1 << 6,
              OutOfRange      = 1 << 7,
              QrOutOfRange    = 1 << 8,
              LowGain         = 1 << 9
            );
    // clang-format on
    SweepFitter::TargetSweepData data{};
    // double window_width{0};
    double window_Qr{0};
    Eigen::MatrixXcd uncertainty{};
    Eigen::MatrixXcd diqs{};
    Eigen::MatrixXd d21ps{};
    Eigen::MatrixXd d21costs{};
    SweepFitter::ModelSpec modelspec{};
    Eigen::MatrixXd ps{};
    Eigen::MatrixXd costs{};
    Eigen::MatrixXd output{};
    std::vector<std::string> colnames{};
    std::unordered_map<Flag, int> flag_summary{};

    void save(const std::string &filepath) const;
    void save_processed(const std::string &filepath) const;
    void plot() const;
    void save_plot(const std::string &filepath) const;
    auto operator()(const std::string &name) const {
        if (auto i = tula::container_utils::indexof(colnames, name);
            i.has_value()) {
            return output.col(i.value());
        }
        throw std::runtime_error(
            fmt::format("invalid colname {}. {}", name, colnames));
    }
};

TULA_ENUM_REGISTER(SweepFitter::WeightOption);
TULA_ENUM_REGISTER(SweepFitter::FittingOption);
TULA_ENUM_REGISTER(SweepFitter::ModelSpec);
TULA_ENUM_REGISTER(SweepFitResult::Flag);
TULA_BITFLAG_VALUE_MASK(SweepFitResult::Flag, 0xFFFF);

} // namespace kids

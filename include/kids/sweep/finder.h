#pragma once

#include "../core/kidsdata.h"
#include "fitter.h"
#include <tula/enum.h>

namespace kids {

struct SweepKidsFinderResult;

struct SweepKidsFinder {

    using VnaSweepData = KidsData<KidsDataKind::VnaSweep>;
    using Config = KidsData<>::meta_t;
    using Segment = std::pair<Eigen::Index, Eigen::Index>;
    using Candidate = std::tuple<Eigen::Index, Eigen::Index, Eigen::Index>;
    struct IterStep {
        double thresh;
        std::vector<Segment> segments;
        std::vector<Eigen::Index> cands;
        std::vector<Eigen::Index> icands_good;
        Eigen::VectorXd iterdata;
        Eigen::VectorXd residual;
        Eigen::VectorXd adiqsmdl;
        SweepFitResult candsfitresult;
    };

    // clang-format off
    TULA_ENUM_DECL(FittingOption, int,
              default_,
              fixfr
              );
    // clang-format on

    SweepKidsFinder(Config config_);

    /// @brief find kids in sweep file
    auto operator()(const VnaSweepData &, const Config &config = {})
        -> SweepKidsFinderResult;

    Config config;
};

struct SweepKidsFinderResult {
    // clang-format on
    SweepKidsFinder::VnaSweepData data{};
    Eigen::VectorXd rfs{};
    Eigen::VectorXd adiqs{};
    Eigen::VectorXd adiqscov{};
    Eigen::VectorXd stats_rfs{};
    Eigen::VectorXd adiqs_fcor{};
    Eigen::VectorXd adiqsmean{};
    Eigen::VectorXd adiqsstd{};
    double threshold;
    std::vector<SweepKidsFinder::IterStep> itersteps;
    std::map<Eigen::Index, std::vector<std::pair<Eigen::Index, Eigen::Index>>>
        candidates;
    std::vector<SweepKidsFinder::Candidate> unique_candidates;
    Eigen::MatrixXd output{};

    void save(const std::string &filepath) const;
    void save_d21(const std::string &filepath) const;
    void save_processed(const std::string &filepath) const;
    void plot() const;
    void save_plot(const std::string &filepath) const;
};

TULA_ENUM_REGISTER(SweepKidsFinder::FittingOption);

} // namespace kids

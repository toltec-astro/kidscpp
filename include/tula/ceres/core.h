#pragma once

#include <Eigen/Core>
#include <ceres/ceres.h>
#include <initializer_list>
#include <ranges>
#include <stdexcept>
#include <tula/eigen.h>
#include <tula/formatter/container.h>
#include <tula/formatter/matrix.h>
#include <tula/logging.h>
#include <tula/meta.h>

namespace tula::ceres_utils {

using Eigen::Dynamic;
using Eigen::Index;

using ceres::AutoDiffCostFunction;
using ceres::DynamicAutoDiffCostFunction;
using ceres::CauchyLoss;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::SubsetManifold;

/// @brief The mode of the evaluation
/// @var Residual
///     Evaluate residual.
/// @var Model
///     Evaluate model.
enum class EvalMode { Residual, Model };

/// @brief Defines properties for model parameter
template <typename Scalar = double>
struct Parameter {

    using OptScalar = std::optional<Scalar>;
    std::string name{"unnamed"};
    OptScalar value{std::nullopt};
    OptScalar lower_bound{std::nullopt};
    OptScalar upper_bound{std::nullopt};
    bool vary{true};

    /// @brief Create new parameter with bounds
    static auto Bounded(const OptScalar &lower, const OptScalar &upper,
                        const OptScalar &value) {
        return Parameter{
            .lower_bound = lower, .upper_bound = upper, .value = value};
    }

    static auto Bounded(const OptScalar &lower, const OptScalar &upper) {
        return Bounded(lower, upper, std::nullopt);
    }

    /// @brief Create new parameters that are fixed
    static auto Fixed(Scalar value) {
        return Parameter{.value = value, .vary = false};
    }
};

template <typename Scalar = double>
struct Parameters {

    using parameter_t = Parameter<Scalar>;

    Parameters() = default;

    Parameters(std::vector<parameter_t> params) : m_params{std::move(params)} {
        make_param_name_index();
    };

    auto operator[](const std::string &name) -> const parameter_t & {
        return this->operator[](m_name_index.at(name));
    }
    auto operator[](Index idx) -> const parameter_t & {
        return m_params.at(idx);
    }

    auto names() const noexcept -> decltype(auto) {
        return std::ranges::transform_view(
            m_params, [](const auto &param) { return param.name; });
    }
    auto size() -> decltype(auto) { return m_params.size(); }

    auto params() const -> const auto & { return m_params; }
    ///
    /// @brief Initialize parameter data.
    template <typename Derived>
    auto init_data(const Eigen::DenseBase<Derived> &param_data_, bool init_value = true){
        auto &param_data =
            const_cast<Eigen::DenseBase<Derived> &>(param_data_).derived();
        if constexpr (tula::eigen_utils::is_plain_v<Derived>) {
            if (auto size = param_data.size(); size != this->size()) {
                // make sure we resize the param data if not already
                param_data.resize(size);
                SPDLOG_TRACE("resize params size {} to {}", size,
                             param_data.size());
            }
        }
        // makesure params.data() is continugous up to n_params
        if (param_data.innerStride() != 1 || param_data.innerSize() < this->size()) {
            throw std::runtime_error(
                fmt::format("fitter requires params data of size {}"
                            " in continugous memory",
                            this->size()));
        }
        // set the intial values if specified
        if (init_value ) {
            for (std::size_t i = 0; i < this->m_params.size(); ++i) {
                const auto& param = m_params[i];
                if (param.value.has_value()) {
                    param_data.coeffRef(i) = param.value.value();
                    SPDLOG_TRACE("init value of param {}: {}", param, param_data.coeff(i));
                }
            }
        }
        return param_data.data();
    }
    
    /// @brief Create ceres problem
    template <typename Derived>
    auto create_problem(const Eigen::DenseBase<Derived> &param_data_,
                        bool init_data=false, bool init_value=true) {
        auto &param_data =
            const_cast<Eigen::DenseBase<Derived> &>(param_data_).derived();
        SPDLOG_TRACE("create problem for param {} init_data={} init_value={}",
                param_data,
                init_data,
                init_value
        );
        if (init_data) {
            this->init_data(param_data, init_value);
        }
        auto n_params = this->size();
        if (n_params != param_data.size()) {
            throw std::runtime_error(fmt::format(
                "param setting size {} mismatch params size {}", n_params ,param_data.size()));
        }
        auto param_value = param_data.data();
        auto problem = std::make_shared<Problem>();
        problem->AddParameterBlock(param_value, param_data.size());
        std::vector<int> fixed_params;
        // setup params
        for (std::size_t i = 0; i < n_params; ++i) {
            auto param = m_params[i];
            if (param.lower_bound.has_value()) {
                problem->SetParameterLowerBound(param_value, i,
                                                param.lower_bound.value());
            }
            if (param.upper_bound.has_value()) {
                problem->SetParameterUpperBound(param_value, i,
                                                param.upper_bound.value());
            }
            if (! param.vary) {
                fixed_params.push_back(i);
            }
        }
        if (fixed_params.size() > 0) {
            SubsetManifold *sm =
                new SubsetManifold(n_params, fixed_params);
            problem->SetManifold(param_value, sm);
        }
        // if (fixed_params.size() > 0) {
        //     SubsetParameterization *sp =
        //         new SubsetParameterization(n_params, fixed_params);
        //     problem->SetParameterization(param_value, sp);
        // }

        SPDLOG_TRACE("params init values: {}", param_data);
        SPDLOG_TRACE("fixed params: {}", fixed_params);
        return std::make_pair(problem, param_value);
    }

private:
    std::vector<parameter_t> m_params{};
    std::unordered_map<std::string, Index> m_name_index{};
    void make_param_name_index() {
        for (std::size_t i = 0; i < m_params.size(); ++i) {
            const auto &param = m_params[i];
            if (m_name_index.contains(param.name)) {
                throw std::runtime_error(
                    fmt::format("duplicated parameter name: {}", param.name));
            }
            m_name_index[param.name] = i;
        }
    }
};

template <Index NP_ = Dynamic, Index NDIM_IN_ = Dynamic,
          Index NDIM_OUT_ = Dynamic, typename Scalar_ = double>
struct ModelFitter {
    constexpr static Index NP = NP_;
    constexpr static Index NDIM_IN = NDIM_IN_;
    constexpr static Index NDIM_OUT = NDIM_OUT_;
    using Scalar = Scalar_;

    using ParamDataType = Eigen::Map<const Eigen::Matrix<Scalar, NP, 1>>;
    ModelFitter() = default;

    /// @brief Create autodiff cost function by providing number of residuals
    template <typename ResidualEvaluator>
    static auto set_residual(Problem *problem, Scalar *paramblock,
                                      ResidualEvaluator *residual_evaluator) {
        // setup residual block
        if constexpr (ResidualEvaluator::use_num_diff) {
            if constexpr (NP == Dynamic) {
                ceres::DynamicNumericDiffCostFunction <ResidualEvaluator>* cost_function  =
                new ceres::DynamicNumericDiffCostFunction<ResidualEvaluator>(residual_evaluator, ceres::Ownership::DO_NOT_TAKE_OWNERSHIP);
                cost_function->SetNumResiduals(residual_evaluator->n_y_data);
                cost_function->AddParameterBlock(problem->NumParameters());
                problem->AddResidualBlock(cost_function,
                                          residual_evaluator->loss_func(), paramblock);
            } else {
                CostFunction* cost_function
                    = new ceres::NumericDiffCostFunction<ResidualEvaluator, ceres::NumericDiffMethodType::CENTRAL, Dynamic, NP>(
                            residual_evaluator,
                        ceres::Ownership::DO_NOT_TAKE_OWNERSHIP, 
                        residual_evaluator->n_y_data);
                problem->AddResidualBlock(cost_function,
                                          residual_evaluator->loss_func(), paramblock);
            }
        }
        else if constexpr (NP == Dynamic ) {
             DynamicAutoDiffCostFunction<ResidualEvaluator, 4> *cost_function =
                new DynamicAutoDiffCostFunction<ResidualEvaluator, 4>(
                    residual_evaluator);
            cost_function->SetNumResiduals(residual_evaluator->n_y_data);
            problem->AddResidualBlock(cost_function,
                                      residual_evaluator->loss_func(), paramblock);
        
        } else {
            CostFunction *cost_function =
                new AutoDiffCostFunction<ResidualEvaluator, Dynamic, NP>(
                    residual_evaluator, residual_evaluator->n_y_data);
            problem->AddResidualBlock(cost_function,
                                      residual_evaluator->loss_func(), paramblock);
        }
    }
    
   
    template <typename ResidualEvaluator, typename ProblemParamBlock>
    static auto fit(
             ResidualEvaluator * residual_evaluator,
             const ProblemParamBlock & problem_paramblock
             ) {
        auto [problem, paramblock] = problem_paramblock;
        set_residual(problem.get(), paramblock, residual_evaluator);
        SPDLOG_TRACE("initial params {}",
                     tula::fmt_utils::pprint(paramblock, problem->NumParameters()));
        SPDLOG_TRACE("x_fit_data{}", tula::fmt_utils::pprint(residual_evaluator->x_fit_data, residual_evaluator->n_x_data));
        SPDLOG_TRACE("y_fit_data{}", tula::fmt_utils::pprint(residual_evaluator->y_fit_data, residual_evaluator->n_y_data));
        SPDLOG_TRACE("y_err_data{}", tula::fmt_utils::pprint(residual_evaluator->y_err_data, residual_evaluator->n_y_data));
        // do the fit
        Solver::Options options;
        options.linear_solver_type = ceres::DENSE_QR;
        // options.use_inner_iterations = true;
        // options.minimizer_progress_to_stdout = true;
        options.logging_type = ceres::SILENT;
        Solver::Summary summary;
        ceres::Solve(options, problem.get(), &summary);

        SPDLOG_TRACE("{}", summary.BriefReport());
        SPDLOG_TRACE("fitted params {}", tula::fmt_utils::pprint(paramblock, problem->NumParameters()));
        return std::make_tuple(summary.termination_type == ceres::CONVERGENCE,
                               std::move(summary));
        }
    };

} // namespace tula::ceres_utils

namespace fmt {

template <typename T>
struct formatter<tula::ceres_utils::Parameter<T>>
    : tula::fmt_utils::nullspec_formatter_base {
    template <typename FormatContext>
    auto format(const tula::ceres_utils::Parameter<T> &param,
                FormatContext &ctx) {
        auto it = ctx.out();
        return format_to(it, "P({},value={},lower={},upper={},vary={})",
                         param.name, param.value, param.lower_bound,
                         param.upper_bound, param.vary);
    }
};

template <typename T>
struct formatter<tula::ceres_utils::Parameters<T>>
    : tula::fmt_utils::nullspec_formatter_base {
    template <typename FormatContext>
    auto format(const tula::ceres_utils::Parameters<T> &params,
                FormatContext &ctx) {
        auto it = ctx.out();
        return format_to(it, "Ps({})", params.params());
    }
};

} // namespace fmt

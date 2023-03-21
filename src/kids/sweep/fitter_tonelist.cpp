#include "kids/sweep/fitter.h"
#include "kids/sweep/model.h"
#include <csv_parser/parser.hpp>
#include <fstream>
#include <gram_savitzky_golay/gram_savitzky_golay.h>
#include <ceres/ceres.h>
#include <tula/algorithm/ei_ceresfitter.h>
#include <tula/algorithm/ei_convolve.h>
#include <tula/algorithm/ei_numdiff.h>
#include <tula/algorithm/ei_stats.h>
#include <tula/ceres/core.h>
#include <tula/container.h>
#include <tula/ecsv/table.h>
#include <tula/eigen.h>
#include <tula/formatter/enum.h>
#include <tula/grppi.h>
#include <tula/switch_invoke.h>
#include <utility>

using Eigen::Index;
using kids::SweepFitResult;
using kids::SweepFitter;
using std::ifstream;

template <typename Config, typename meta_t>
auto loadtonelist(const Config &config, const meta_t &meta)
    -> std::optional<tula::ecsv::ECSVTable> {
    auto meta_get_int = [&](const auto &key) {
        return meta.template get_typed<int>(key);
    };
    auto pattern = fmt::format("{}{:d}_{:06d}_{:03d}_{:04d}.+_tonelist\\.ecsv",
                               meta.get_str("instru"), meta_get_int("roachid"),
                               meta_get_int("obsid"), meta_get_int("subobsid"),
                               meta_get_int("scanid"));

    std::string filepath{};
    if (config.has("tonelisttablefile")) {
        filepath = config.get_str("tonelisttablefile");
    } else if (config.has("tonelisttabledir")) {
        auto dir = config.get_str("tonelisttabledir");
        SPDLOG_TRACE("look for tonelist table dir {} with pattern {}", dir,
                     pattern);
        auto candidates = tula::filename_utils::find_regex(dir, pattern);
        if (!candidates.empty()) {
            filepath = candidates[0];
        } else {
            throw std::runtime_error(
                fmt::format("no tone list table found in {} that matches {}",
                            dir, pattern));
        }
    } else {
        throw std::runtime_error(
            fmt::format("no tonelist table location specified."));
    }
    SPDLOG_INFO("use tonelist table file {}", filepath);

    // load file

    try {
        ifstream fo(filepath);
        auto tbl = tula::ecsv::ECSVTable(tula::ecsv::ECSVHeader::read(fo));
        auto parser =
            aria::csv::CsvParser(fo).delimiter(tbl.header().delimiter());
        tbl.load_rows(parser);

        SPDLOG_DEBUG("tonelist tbl: {}", tbl);
        SPDLOG_DEBUG("tonelist tbl meta:\n{}", YAML::Dump(tbl.header().meta()));
        SPDLOG_INFO("tonelist tbl_info:\n{}", tbl.info());
        return tbl;
    } catch (std::runtime_error &e) {
        throw std::runtime_error(
            fmt::format("unable to load tonelist file: {}", e.what()));
    }
    return std::nullopt;
}

/// @brief A set of Kids Models.
template <kids::SweepModel model_>
struct KidsModelSet {
    static constexpr auto model = model_;
    using model_t = kids::internal::sweepmodel_t<model>;
    constexpr static auto n_params_per_model = model_t::param_names.size();
    using params_t = tula::ceres_utils::Parameters<double>;
    using param_t = typename params_t::parameter_t;
    using SweepModel = kids::SweepModel;
    constexpr static auto use_num_diff = true;

    // KidsModelSet(params_t params_, Index n_models_)
    //     : params{std::move(params_)}, n_models{n_models_} {}

    KidsModelSet(tula::ecsv::ECSVTable tonelist,
                 const std::vector<Index> &row_indices)
        : n_models{tula::meta::size_cast<Index>(row_indices.size())} {
        // build the model params info
        std::vector<param_t> params_init;
        constexpr auto param_names = model_t::param_names;
        Index pi_first_fr = -1;
        Index pi_first_Qr = -1;

        for (const auto &row : row_indices) {
            auto entry_get_int = [&](const auto &colname) {
                return tula::meta::size_cast<Index>(
                    tonelist.col<int64_t>(colname)(row));
            };
            auto chan_id = entry_get_int("chan_id");
            auto model_chan_id = entry_get_int("model_chan_id");
            // auto model_n_chans = tonelist.col<int>("model_n_chans")(row);

            auto tone_id = entry_get_int("tone_id");
            auto model_tone_id = entry_get_int("model_tone_id");
            // auto model_n_tones = tonelist.col<int>("model_n_tones")(row);

            auto f_fit = tonelist.col<double>("f_fit")(row);
            auto Qr_init = tonelist.col<double>("Qr_init")(row);
            constexpr double Qr_init_min = 100.;
            auto fr_init = tonelist.col<double>("fr_init")(row);
            auto fr_init_min = tonelist.col<double>("fr_init_min")(row);
            auto fr_init_max = tonelist.col<double>("fr_init_max")(row);
            auto m0_init = tonelist.col<double>("m0_init")(row);
            auto m1_init = tonelist.col<double>("m1_init")(row);

            SPDLOG_TRACE("process tonelist entry row={} chan_id={} "
                         "model_chan_id={} tone_id={} model_tone_id={}",
                         row, chan_id, model_chan_id, tone_id, model_tone_id);
            SPDLOG_TRACE("model init params: f_fit={} Qr_init={} "
                         "Qr_init_min={} fr_init={} fr_init_min={} "
                         "fr_init_max={} m0_init={} m1_init={}",
                         f_fit, Qr_init, Qr_init_min, fr_init, fr_init_min,
                         fr_init_max, m0_init, m1_init);

            // this hold the param index for each model to the param block
            std::array<Index, n_params_per_model> param_indices{};
            tula::meta::static_for<Index, 0, n_params_per_model>(
                [&, this](auto param_index_) mutable {
                    constexpr auto param_index =
                        std::decay_t<decltype(param_index_)>::value;
                    constexpr auto name = model_t::param_names[param_index];
                    auto pname = fmt::format("m{}{}_{}", model_chan_id,
                                             model_tone_id, name);
                    SPDLOG_TRACE("add model parameter {}", pname);
                    // this is the current parameter index because we are about
                    // to add the item.
                    auto pindex =
                        tula::meta::size_cast<Index>(params_init.size());
                    if constexpr (name == "fp") {
                        params_init.emplace_back(param_t{
                            .name = pname, .value = f_fit, .vary = false});
                    } else if constexpr (name == "Qr") {
                        if (pi_first_Qr < 0) {
                            pi_first_Qr = pindex;
                        } else {
                            // override pindex to be the first Qr
                            pindex = pi_first_Qr;
                        }
                        params_init.emplace_back(
                            param_t{.name = pname,
                                    .value = Qr_init,
                                    .lower_bound = Qr_init_min,
                                    .vary = true});
                    } else if constexpr (name == "Qc") {
                        params_init.emplace_back(
                            param_t{.name = pname, .value = 1., .vary = false});
                    } else if constexpr (name == "fr") {
                        if (pi_first_fr < 0) {
                            pi_first_fr = pindex;
                        } else {
                            pindex = pi_first_fr;
                        }
                        params_init.emplace_back(
                            param_t{.name = pname,
                                    .value = fr_init,
                                    .lower_bound = fr_init_min,
                                    .upper_bound = fr_init_max,
                                    .vary = true});
                    } else if constexpr (name == "A") {
                        params_init.emplace_back(
                            param_t{.name = pname, .value = 0., .vary = false});
                    } else if constexpr (name == "normI") {
                        params_init.emplace_back(
                            param_t{.name = pname, .value = 1., .vary = true});
                    } else if constexpr (name == "normQ") {
                        params_init.emplace_back(
                            param_t{.name = pname, .value = 0., .vary = true});
                    } else if constexpr (name == "slopeI") {
                        params_init.emplace_back(
                            param_t{.name = pname, .value = 0., .vary = true});
                    } else if constexpr (name == "slopeQ") {
                        params_init.emplace_back(
                            param_t{.name = pname, .value = 0., .vary = true});
                    } else if constexpr (name == "interceptI") {
                        params_init.emplace_back(
                            param_t{.name = pname, .value = m0_init, .vary = true});
                    } else if constexpr (name == "interceptQ") {
                        params_init.emplace_back(
                            param_t{.name = pname, .value = m1_init, .vary = true});
                    } else {
                        params_init.emplace_back(
                            param_t{.name = pname, .vary = true});
                    }
                    // populate param indices
                    param_indices[param_index] = pindex;
                    SPDLOG_TRACE("parameter: {}", params_init.back());
                });
            SPDLOG_TRACE("add model index {} {} {}", chan_id, tone_id,
                         param_indices);
            // populate model index vector
            this->m_model_index.push_back(
                std::tuple{chan_id, tone_id, param_indices});
        }
        this->params = params_t{std::move(params_init)};
        SPDLOG_DEBUG("initialized model parameters: {}", this->params);
    }
    
    // here residual value has layout as [real, imag, real, imag, ...]
    // template <typename T>
    // inline void eval(T const  *const * param_value, T *residual_value) const {
    template<typename T>
    inline void eval(double const * const * param_value, T* residual_value) const {
        auto n_blocks = this->x_data_blocks->size();
        // T* params_array = new T[this->n_models * this->n_params_per_model];
        Eigen::VectorXd params_array(this->n_models * this->n_params_per_model);
        params_array.setConstant(0.);
        _populate_params_array(param_value[0], params_array.data());
        // link parameter value
        // _link_param_value(param_value);
        for (Index i = 0; i < this->n_models; ++i) {
            const auto &[chan_idx, tone_idx, param_indices] =
                this->m_model_index[i];
            // // copy the params values of this model to tmp params for
            // evaluation for (Index j = 0; j < n_params_per_model; ++j) {
            //     tmp_param[j] = param_value[param_indices[j]];
            // }
            auto [x_data_offset, x_data_size] = this->x_data_blocks->at(i);
            // const auto &[x_data, residual_value_offset] =
            //     x_data_blocks[chan_idx];
            // auto nx = x_data.size();
            model_t::_eval(params_array.data() +i * n_params_per_model, 
                          residual_value + 2 * x_data_offset, // compelx
                       x_data_size, this->x_fit_data + x_data_offset);
        }
    }

    // model
    // template <typename T>
    // inline auto operator()( T const *const * param_value, T *residual_value) const
    template<typename T>
    inline auto operator()(double const *   const* param_value, T* residual_value) const 
        -> bool {
        // compute model
        eval(param_value, residual_value);
        // compute residule
        for (Index i = 0; i < this->n_x_data; ++i) {
            // compute residual
            auto y_r = this->y_fit_data[2 * i];     // S21_r
            auto y_i = this->y_fit_data[2 * i + 1]; // S21_i
            auto e_r = this->y_err_data[2 * i];     // e_S21_r
            auto e_i = this->y_err_data[2 * i + 1]; // e_S21_i
            residual_value[2 * i] = (y_r - residual_value[2 * i]) / e_r;
            residual_value[2 * i + 1] = (y_i - residual_value[2 * i + 1]) / e_i;
        }
        return true;
    }
    
    static inline auto loss_func() {
        return new ceres::CauchyLoss(2.0);
    }

    template <typename T>
    void _populate_params_array(const T *const param_values, T *params_array) const {
        // copy values in param values in array form using param_indices.
        // params_array is of shape [n_models, n_params_per_model], col
        // oriented
        for (Index i = 0; i < n_models; ++i) {
            const auto &[chan_idx, tone_idx, param_indices] = m_model_index[i];
            for (Index j = 0; j < n_params_per_model; ++j) {
                auto array_idx = i * n_params_per_model + j;
                // SPDLOG_TRACE("overwrite model param {} ({}) -> {} ({})", params_array[array_idx], array_idx, param_values[param_indices[j]], param_indices[j]);
                params_array[array_idx] = param_values[param_indices[j]];
            }
        }
    }

    // template <typename T>
    // void _link_param_value(T const  *const* param_value) const {
    //     // copy values in param values using param_indices.
    //     // params_value is of shape [n_models, n_params_per_model], col
    //     // oriented
    //     for (Index i = 0; i < n_models; ++i) {
    //         const auto &[chan_idx, tone_idx, param_indices] = m_model_index[i];
    //         for (Index j = 0; j < n_params_per_model; ++j) {
    //             auto p_idx = i * n_params_per_model + j;
    //             param_value[p_idx] = param_value[param_indices[j]];
    //         }
    //     }
    // }

    auto get_params_array(auto &&params) {
        Eigen::MatrixXd params_array(n_params_per_model, n_models);
        _populate_params_array(params.data(), params_array.data());
        return params_array;
    }
    
template <typename DerivedA, typename DerivedB,
          typename DerivedC>
auto set_data(const Eigen::DenseBase<DerivedA> &x_fit_data_,
         const Eigen::DenseBase<DerivedB> &y_fit_data_,
         const Eigen::DenseBase<DerivedC> &y_err_data_,
         std::vector<std::tuple<Index, Index>> const * x_data_blocks
         ) {
    auto &x_fit_data = x_fit_data_.derived();
    auto &y_fit_data = y_fit_data_.derived();
    auto &y_err_data = y_err_data_.derived();
    this->n_x_data = x_fit_data.size();
    this->n_y_data = y_fit_data.size() * 2;  // flattened complex
    this-> x_fit_data=reinterpret_cast<const double*>(x_fit_data.data());
    this-> y_fit_data=reinterpret_cast<const double*>(y_fit_data.data());
    this-> y_err_data=reinterpret_cast<const double*>(y_err_data.data());
    this->x_data_blocks = x_data_blocks;
     }

    params_t params;
    Index n_models;
    std::vector<std::tuple<Index, Index, std::array<Index, n_params_per_model>>>
        m_model_index;
    std::vector<std::tuple<Index, Index>> const *x_data_blocks = nullptr;
    // required by fitter
    double const *x_fit_data = nullptr;
    double const *y_fit_data = nullptr;
    double const *y_err_data = nullptr;
    Index n_x_data{0};
    Index n_y_data{0};
};

auto kids::SweepFitter::fit_tonelist(const TargetSweepData &data,
                                     const Config &config_) -> SweepFitResult {
    // merge the config object
    Config config{this->config};

    config.update(config_);

    SPDLOG_DEBUG("run kids fit_tonelist with config: {}", config.pformat());

    auto opt_tonelist = loadtonelist(config, data.meta);
    auto tonelist = opt_tonelist.value();
    SPDLOG_TRACE("tonelist int array{}",
                 tonelist.array_data<int64_t>().array());
    // split the tone list into groups by group id
    std::map<std::size_t, std::vector<Index>> tonelist_groups;
    auto col_group_id = tonelist.col<int64_t>("group_id");
    for (Index i = 0; i < tonelist.rows(); ++i) {
        auto group_id = tula::meta::size_cast<std::size_t>(col_group_id(i));
        // SPDLOG_TRACE("collect tonelist item {} group_id={}", i,
        // group_id);
        if (tonelist_groups.contains(group_id)) {
            // push to exist
            tonelist_groups[group_id].push_back(i);
            // SPDLOG_TRACE("update existing {}: {}", group_id, i);
        } else {
            // add new
            tonelist_groups[group_id] = std::vector<Index>{i};
            // SPDLOG_TRACE("add new {}: {}", group_id, i);
        }
    }
    auto n_models = tonelist.rows();
    auto n_groups = tonelist_groups.size();
    SPDLOG_DEBUG("run fit for {} tone groups", n_groups);
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
    auto groupindices = tula::container_utils::index(n_groups);
    Eigen::VectorXd fs0 = data.tones();
    // intial flag
    using Flag = SweepFitResult::Flag;
    const auto flag0 = bitmask::bitmask{Flag::Good};
    // flags container
    auto flag = std::vector(TULA_SIZET(n_models), flag0);
    // helper to update flag
    auto update_flag = [&flag](auto value, auto idx)mutable  {
        flag[idx] |= value;
        SPDLOG_TRACE("add flag {} to model #{}", value, idx);
    };

    // handle weight
    WeightOption weight_option{WeightOption::none};
    if (auto type =
            WeightOption_meta::from_name(config.get_str("weight_window_type"));
        type.has_value()) {
        SPDLOG_WARN("weight window is forced to be none in tonelist mode.");
        // weight_option = type.value().value;
    } else {
        SPDLOG_WARN("unknown weight window type {}, use {} instead",
                    config.get_str("weight_window_type"), weight_option);
    }
    using WeightTypes =
        tula::meta::cases<WeightOption::boxcar, WeightOption::gauss,
                          WeightOption::lorentz, WeightOption::none>;
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

    ModelSpec modelspec{ModelSpec::gainlintrend};
    if (auto spec = ModelSpec_meta::from_name(config.get_str("modelspec"));
        spec.has_value()) {
        modelspec = spec.value().value;
    } else {
        SPDLOG_WARN("unknown model spec {}, use {} instead",
                    config.get_str("modelspec"), modelspec);
    }
    using ModelSpecs =
        tula::meta::cases<ModelSpec::gain, ModelSpec::gainlintrend,
                          ModelSpec::trans, ModelSpec::translintrend>;
    // dispatch different models
    auto fitresult = tula::meta::switch_invoke<ModelSpecs>(
        [&](auto spec_) {
            constexpr auto spec = std::decay_t<decltype(spec_)>::value;
            using Model = model_t<spec>;
            SPDLOG_DEBUG("fit tone list with model: {}", spec);
            // loop through the tonelist to run the model fitting on a per
            // group basis make output container
            Eigen::MatrixXd output(tonelist.rows(), Model::NP);
            Index ig = 0;
            grppi::map(ex, groupindices, groupindices, [&] (auto gi) {
                    // if (gi > 1) {
                    //     return gi;
                    // }
                auto tonelist_group_id = gi;
                const auto & tonelist_row_indices = tonelist_groups[gi];

                SPDLOG_TRACE("make model for group_id={} row_indices={}",
                             tonelist_group_id, tonelist_row_indices);
                auto modelset =
                    KidsModelSet<Model::model>{tonelist, tonelist_row_indices};
                // prepare input data
                Index n_fit_data{0};
                //                    chan_id istart istop  offset length
                std::vector<std::tuple<Index, Index, Index, Index, Index>>
                    sweep_data_index;
                for (Index i = 0; i < tonelist_row_indices.size(); ++i) {
                    auto row = tonelist_row_indices[i];
                    // get data start and stop indices
                    auto entry_get_int = [&](const auto &colname) {
                        return tula::meta::size_cast<Index>(
                            tonelist.col<int64_t>(colname)(row));
                    };

                    auto chan_id = entry_get_int("chan_id");
                    auto model_chan_id = entry_get_int("model_chan_id");
                    auto tone_id = entry_get_int("tone_id");
                    auto model_tone_id = entry_get_int("model_tone_id");
                    auto index_start = entry_get_int("index_start");
                    auto index_stop = entry_get_int("index_stop");

                    SPDLOG_TRACE("process tonelist entry row={} chan_id={} "
                                 "model_chan_id={} tone_id={} model_tone_id={} "
                                 "index_start={} index_stop={}",
                                 row, chan_id, model_chan_id, tone_id,
                                 model_tone_id, index_start, index_stop);
                    Index fit_data_offset = n_fit_data;
                    Index n_data = index_stop - index_start + 1;
                    sweep_data_index.emplace_back(chan_id, index_start,
                                                  index_stop, fit_data_offset,
                                                  n_data);
                    n_fit_data += n_data;
                }
                Eigen::VectorXcd residual(n_fit_data);
                Eigen::VectorXcd y_err_data(n_fit_data);
                Eigen::VectorXcd y_fit_data(n_fit_data);
                Eigen::VectorXd x_fit_data(n_fit_data);
                std::vector<std::tuple<Index, Index>> x_data_blocks;
                for (auto [chan_id, index_start, index_stop, fit_data_offset,
                           block_size] : sweep_data_index) {
                    x_fit_data.segment(fit_data_offset, block_size) =
                        data.fs().col(chan_id).segment(index_start, block_size);
                    y_fit_data.segment(fit_data_offset, block_size) =
                        data.iqs().col(chan_id).segment(index_start,
                                                        block_size);
                    y_err_data.segment(fit_data_offset, block_size) =
                        uncertainty.col(chan_id).segment(index_start,
                                                         block_size);
                    SPDLOG_TRACE("added data segment chan_id={} index_start={} "
                                 "index_stop={} offset={} block_size={} "
                                 "total_size={}",
                                 chan_id, index_start, index_stop,
                                 fit_data_offset, block_size, n_fit_data);
                    x_data_blocks.emplace_back(fit_data_offset, block_size);
                }
                SPDLOG_TRACE("x_fit_data{}", x_fit_data);
                SPDLOG_TRACE("y_fit_data{}", y_fit_data);
                SPDLOG_TRACE("y_err_data{}", y_err_data);
                SPDLOG_TRACE("x_data_blocks{}", x_data_blocks);

                // do fit
                tula::ceres_utils::ModelFitter model_fitter{};
                modelset.set_data(x_fit_data, y_fit_data, y_err_data, &x_data_blocks);

                Eigen::VectorXd param_data(modelset.params.size());
                param_data.setConstant(0.);
                auto problem_paramblock = modelset.params.create_problem(param_data, true);
                SPDLOG_TRACE("init param_data{}", param_data);
                auto [converged, summary]  = model_fitter.fit(&modelset, problem_paramblock);
                SPDLOG_TRACE("best-fit param_data{}", param_data);
                Eigen::MatrixXd ps_out = modelset.get_params_array(param_data);
                SPDLOG_TRACE("best-fit params array{}", ps_out);
                // populate output container
                for (Index i = 0; i < tonelist_row_indices.size(); ++i) {
                    auto row = tonelist_row_indices[i];
                    output.row(row) = ps_out.col(i);
                    auto entry_get_double = [&](const auto &colname) {
                        return tonelist.col<double>(colname)(row);
                    };
                    // check fit convergence
                    if (!converged) {
                        update_flag(Flag::NotConverged, row);
                    }
                    auto f0 = ps_out(0, i);
                    auto f = ps_out(3, i);
                    if ((f < entry_get_double("f_fit_min")) ||
                        (f > entry_get_double("f_fit_max"))) {
                        update_flag(Flag::OutOfRange, row);
                    }
                    auto large_offset = large_offset_at_500MHz * f0 / 5e8;
                    if (std::abs(f - f0) > large_offset) {
                        update_flag(Flag::LargeOffset, row);
                    }
                    // check S21 fitting Qr range
                    auto Qr = ps_out(1, i);
                    if ((Qr < lim_Qr_min) || (Qr > lim_Qr_max)) {
                        update_flag(Flag::QrOutOfRange, row);
                    }
                    // check S21 gain
                    auto gain =
                        std::sqrt(ps_out(5, i) * ps_out(5, i) + ps_out(6, i) * ps_out(6, i));
                    if (gain < lim_gain_min) {
                        update_flag(Flag::LowGain, row);
                    }
                }
                return gi;
            });
            // print summary
            auto flag_summary = [&flag](auto &&flags) {
                auto summary = std::unordered_map<Flag, int>{};
                for (auto flag : flags) {
                    summary[flag] = 0;
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
            flag_report << fmt::format("fit result of n_groups={} n_models={}", n_groups, n_models);
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
            std::vector<std::string> colnames;
            for (const auto & p: Model::param_names) {
                colnames.push_back(std::string(p));
            }
            return SweepFitResult{.data=data, .window_Qr=window_Qr, .uncertainty=uncertainty, .modelspec=modelspec,
             .output = output, .colnames=colnames, .flag_summary=flag_summary};
        },
        modelspec);
    if (fitresult) {
        return fitresult.value();
    }

    throw std::runtime_error(fmt::format("unknown model spec {}", modelspec));
}
#include "kids/sweep/finder.h"
#include "kids/core/utils.h"
#include "kids/sweep/model.h"
#include <mutex>
#include <tula/algorithm/ei_ceresfitter.h>
#include <tula/algorithm/ei_convolve.h>
#include <tula/algorithm/ei_interp.h>
#include <tula/algorithm/ei_iterclip.h>
#include <tula/algorithm/ei_linspaced.h>
#include <tula/algorithm/ei_numdiff.h>
#include <tula/algorithm/ei_stats.h>
#include <tula/container.h>
#include <tula/formatter/container.h>
#include <tula/formatter/matrix.h>
#include <tula/grppi.h>
#include <tula/logging.h>
#include <gram_savitzky_golay/gram_savitzky_golay.h>
#include <tula/eigen.h>

using Eigen::Index;
using kids::SweepKidsFinder;
using kids::SweepKidsFinderResult;

SweepKidsFinder::SweepKidsFinder(Config config_) : config(std::move(config_)) {}

auto SweepKidsFinder::operator()(const VnaSweepData &data,
                                 const Config &config_)
    -> SweepKidsFinderResult {
    // merge the config object
    Config config{this->config};

    // some internal configs
    // config.set("data_smooth_size", 30);
    config.set("resample_exclude_edge", 8e3); // 8e4

    // stats
    // config.set("stats_clip_sigma", 1.5);
    config.set("stats_clip_maxiter", 5);
    config.set("stats_window_width", 5e6);
    config.set("stats_stride_width", 4e6);
    // config.set("stats_window_overlap", 0.9);

    // detect
    config.set("detect_max_niter", 30);
    // config.set("detect_min_thresh", 10.);
    config.set("detect_min_segment_size", 20);
    // s21 fit
    // this is the size of the s21 data for fitting per detector
    config.set("fit_width", 150e3);
    config.set("eval_n_fwhms", 6);
    constexpr auto fit_modelspec = SweepFitter::ModelSpec::gainlintrend;
    config.set("fitter_modelspec",
               std::string(tula::enum_utils::name(fit_modelspec)));
    // deblend
    config.set("deblend_tone_dist",
               8000.); // tones that are closer than this is merged

    config.update(config_);

    SPDLOG_DEBUG("kids finder config {}", config.pformat());
    auto ex = tula::grppi_utils::dyn_ex(config.get_str("exmode"));

    // some useful data geometry
    Index nsweeps{0};
    Index ntones{0};
    double sweepstep{0};
    Eigen::MatrixXd toneranges;
    Eigen::VectorXd rfs; // setup resampled frequency grid
    double rfstep{0.};
    Index nrfs{0};
    {
        std::tie(nsweeps, ntones) = tula::alg::shape(data.iqs());
        assert(nsweeps == data.sweeps().size());
        assert(ntones == data.tones().size());

        SPDLOG_TRACE("fs{}", data.fs());
        auto tonestep = tula::alg::step(data.tones());
        auto sweepspan = tula::alg::span(data.sweeps());
        auto coverage = int(round(sweepspan / tonestep));
        SPDLOG_TRACE("tonestep={:g}Hz sweepspan={:g}Hz coverage={}", tonestep, sweepspan, coverage);
        sweepstep = tula::alg::step(data.sweeps());
        toneranges.resize(2, ntones);
        toneranges.row(0) = data.fs().colwise().minCoeff();
        toneranges.row(1) = data.fs().colwise().maxCoeff();
        SPDLOG_TRACE(
                "toneranges{} tonespan={:g}Hz",
                toneranges,
                toneranges.row(1).maxCoeff() - toneranges.row(0).minCoeff()
                );
        // SPDLOG_TRACE(
        //   "fs{:r0} iqs{:r0} tonespan={:g}Hz tonestep={:g}Hz sweepspan={:g}Hz sweepstep={:g}Hz coverage={}",
        //   data.fs(), data.iqs(),
        //   toneranges.row(1).maxCoeff() - toneranges.row(0).minCoeff(),
        //   tonestep, sweepspan, sweepstep, coverage);
        assert(coverage == 2);
        auto [vmin, vmax, step] = tula::container_utils::parse_slice<double>(
            config.get_str("resample"));
        double rfmin = vmin.value_or(toneranges.row(0).minCoeff());
        double rfmax = vmax.value_or(toneranges.row(1).maxCoeff());
        // compute sweep coverage
        rfstep = step.value_or(sweepstep / coverage);
        SPDLOG_TRACE("resample grid [{}:{}:{}]", rfmin, rfmax, rfstep);
        rfs = tula::alg::arange(rfmin, rfmax, rfstep);
        nrfs = rfs.size();
        SPDLOG_TRACE("resampled fs{}", rfs);
    }
    // make unified adiqs
    Eigen::MatrixXcd iqs(data.iqs()); // make a copy of the iqs data in case
                                      // we need to smooth
    Eigen::VectorXd adiqs(nrfs);      // resampled abs(d21)
    Eigen::VectorXd adiqscov(nrfs);   // resampled abs(d21) coverage
    {
       bool use_savgol_deriv = config.get_typed<bool>("use_savgol_deriv");
        auto smooth_size =
            static_cast<Index>(config.get_typed<int>("data_smooth_size"));
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
        auto exclude_edge = config.get_typed<double>("resample_exclude_edge");
        auto exclude_edge_size = exclude_edge / sweepstep;

        if (pre_smooth) {
            auto smooth_i1 = smooth_size - smooth_i0;
            SPDLOG_TRACE("preprocess data with smooth_size={}", smooth_size);
            if (exclude_edge_size < smooth_i1) {
                exclude_edge_size = smooth_i1; //
                SPDLOG_DEBUG("set exclude_edge_size={} to match with "
                             "smooth_size {} ({})",
                             exclude_edge, smooth_size, smooth_i1);
            }
        }
        SPDLOG_TRACE("resample exclude_edge_size={}", exclude_edge_size);

        adiqs.setConstant(0.);
        adiqscov.setConstant(0.);

        std::mutex m_mutex; // avoid writing to same item in adiqs
        auto toneindices = tula::container_utils::index(ntones);
        grppi::map(ex, toneindices, toneindices, [&](auto c) {
            Eigen::VectorXcd diqs(nsweeps);
            diqs.setConstant(std::complex(0., 0.));
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
                    diqs.segment(smooth_i0, nsweeps - smooth_size + 1));
            } else {
                // do smooth if needed
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
                        iqs.col(c).segment(smooth_i0, nsweeps - smooth_size + 1));
                }
                // do deriv
                tula::alg::gradient(iqs.col(c), data.fs().col(c), diqs);
            }
            // interpolate onto rfs
            // get the slice of max range in rfs that is enclosed by tonerange.
            auto [il, ir1] = tula::alg::argwithin(rfs, toneranges.coeff(0, c),
                                            toneranges.coeff(1, c));
            // apply exclude_edge
            il += exclude_edge_size;
            ir1 -= exclude_edge_size;
            auto ni = ir1 - il;
            // interpolate
            SPDLOG_TRACE("interpolate tone {} on [{}:{}] ({})", c, il, ir1, ni);
            Eigen::VectorXd c_adiqs =
                tula::alg::interp(data.fs().col(c), diqs.array().abs().eval(),
                            rfs.segment(il, ni));
            {
                // need this to avoid data racing
                // since the tones overlap with each other
                std::scoped_lock lock(m_mutex);
                adiqs.segment(il, ni).array() += c_adiqs.array();
                adiqscov.segment(il, ni).array() += 1.;
            }
            return c;
        });
        // finally normalize the values by cov
        adiqs = adiqs.cwiseQuotient(adiqscov);
    }
    SPDLOG_TRACE("interpolated adiqs{} adiqscov{}", adiqs, adiqscov);

    // despike the d21
    /*
    {
        const auto sigclip = 4.5;
        const auto sigfrac = 0.3;
        const auto objlim = 5.;
        const auto maxiter = 4;
        const auto fill_value = tula::alg::median(adiqs);

        Eigen::VectorXd adiqs_unc(nrfs);
        adiqs_unc.setOnes();
        Eigen::VectorXb adiqs_mask(nrfs);
        adiqs_mask.setZero();
        auto [cleaned_data, cosmics] =
            tula::alg::lacosmic1d(adiqs, adiqs_unc, adiqs_mask, sigclip, sigfrac,
                            objlim, maxiter, fill_value, ex);
    }
    */
    // stats
    Eigen::VectorXd stats_rfs; // fs at which stats are computed
    Eigen::VectorXd adiqsmean;
    Eigen::VectorXd adiqsstd;
    Eigen::VectorXd adiqs_fcor =
        rfs.coeff(0) / rfs.array(); // account for the 1/f trend in D21 height
    std::atomic<double> adiqsmax{0.};
    double adiqsmean0{0.};
    double adiqsstd0{0.};
    {
        auto stats_clip_maxiter = config.get_typed<int>("stats_clip_maxiter");
        auto stats_clip_sigma = config.get_typed<double>("stats_clip_sigma");
        auto iterclip = tula::alg::iterclip(
            // use median-MAD for robust statistics
            [](const auto &v) { return tula::alg::nanmedmad(v); },
            // the snr cut, we only cut the positive side
            [n = stats_clip_sigma](const auto &v, const auto &c,
                                   const auto &s) { return v <= c + n * s; },
            stats_clip_maxiter);
        // infer number of stats needed
        auto stats_window_width =
            config.get_typed<double>("stats_window_width");
        auto stats_window_size = Index(stats_window_width / rfstep);
        auto stats_stride_width =
            config.get_typed<double>("stats_stride_width");
        assert(stats_window_size <= nrfs);
        auto stats_stride = Index(stats_stride_width / rfstep);
        assert(stats_stride <= nrfs);
        auto nstats = Index(nrfs / stats_stride);
        SPDLOG_TRACE("stats window width={} size={} stride={} n={}",
                     stats_window_width, stats_window_size, stats_stride,
                     nstats);
        stats_rfs.resize(nstats);
        Eigen::VectorXd stats_adiqsmean(nstats);
        Eigen::VectorXd stats_adiqsstd(nstats);
        auto statsindices = tula::container_utils::index(nstats);
        {
            tula::logging::progressbar pb0(
                [](const auto &msg) { SPDLOG_INFO("{}", msg); }, 60, "stats ");
            grppi::map(ex, statsindices, statsindices, [&](auto s) {
                pb0.count(nstats, nstats / 10);
                Index i = s * stats_stride; // index to resampled fs grid
                stats_rfs.coeffRef(s) = rfs.coeff(i);
                // compute sigma clip
                // window range
                auto [il, ir1] =
                    tula::alg::windowindex_fixedsize(i, nrfs, stats_window_size);
                SPDLOG_TRACE("compute stats at i={} window=[{}:{}]", i, il,
                             ir1);
                std::tie(std::ignore, std::ignore, stats_adiqsmean.coeffRef(s),
                         stats_adiqsstd.coeffRef(s)) =
                    iterclip(adiqs.segment(il, ir1 - il));
                if (auto m = adiqs.segment(il, ir1 - il).maxCoeff();
                    (!std::isnan(m)) && (m > adiqsmax)) {
                    adiqsmax = m;
                }
                return s;
            });
        }
        adiqsmean0 = stats_adiqsmean.mean();
        adiqsstd0 = stats_adiqsstd.mean();
        SPDLOG_TRACE("stats_fs{} stats_adiqsmean{} stats_adiqsstd{} "
                     "adiqsmax={} adiqsmean0={} adiqsstd0={}",
                     stats_rfs, stats_adiqsmean, stats_adiqsstd, adiqsmax,
                     adiqsmean0, adiqsstd0);
        // interpolate to full fs
        adiqsmean = tula::alg::interp(stats_rfs, stats_adiqsmean, rfs);
        adiqsstd = tula::alg::interp(stats_rfs, stats_adiqsstd, rfs);
    }

    // iterative find, fit and subtract
    auto initial_thresh = (adiqsmax - adiqsmean0) / adiqsstd0 / 2.;
    SPDLOG_DEBUG("detect with initial thresh={}", initial_thresh);
    auto detect_min_segment_size =
        static_cast<Index>(config.get_typed<int>("detect_min_segment_size"));
    auto deblend_tone_dist = config.get_typed<double>("deblend_tone_dist");

    auto fit_size = Index(config.get_typed<double>("fit_width") / rfstep);
    SPDLOG_TRACE("s21 fit_size={}", fit_size);
    using Flag = SweepFitResult::Flag;
    auto fitter = SweepFitter(Config{
        {"exmode", config.get_str("exmode")},
        {"weight_window_type", config.get_str("fitter_weight_window_type")},
        // {"weight_window_fwhm",
        // we make the weight wide to
        //  config.get_typed<double>("fitter_weight_window_fwhm")},
        {"weight_window_Qr",
         config.get_typed<double>("fitter_weight_window_Qr")},
        {"modelspec", config.get_str("fitter_modelspec")}});

    // this function takes the fitresult and return a list of robust entries
    auto select_tones = [](const auto &fitresult, double deblend_dist,
                           auto accept_flag =
                               Flag::Good) -> std::vector<Index> {
        const auto &output = fitresult.output;
        SPDLOG_TRACE("select from {} tones with accept_flag={} deblend_dist={}",
                     output.rows(), accept_flag, deblend_dist);
        auto fit_f = fitresult("f_out");
        auto fit_flag = fitresult("flag");
        auto fit_cost = fitresult.costs.row(1);
        auto check_f_close = [&fit_f, &deblend_dist](Index lhs, Index rhs) {
            auto f_lhs = fit_f.coeff(lhs);
            auto f_rhs = fit_f.coeff(rhs);
            return std::abs(f_lhs - f_rhs) > deblend_dist;
        };
        std::vector<Index> cands;
        // go over all rows
        for (Index c = 0; c < output.rows(); ++c) {
            auto flag = static_cast<Flag>(fit_flag.coeff(c));
            auto cost = fit_cost.coeff(c);
            [[maybe_unused]] auto f = fit_f.coeff(c);
            if (flag == accept_flag) {
                // same, keep
            } else if ((flag | accept_flag) == accept_flag) {
                // flag is contained by accept flag
                // keep
            } else {
                // flag contains flag that are not in accept_flag
                SPDLOG_TRACE("reject c={} fit_flag={} accept_flag={}", c, flag,
                             accept_flag);
                // skip and move on to next output
                continue;
            }
            // check if exists equivalent entry already
            if (auto cit = std::find_if(
                    cands.begin(), cands.end(),
                    [&](const auto &c0) { return !check_f_close(c0, c); });
                cit != cands.end()) {
                auto c0 = *cit; // the existing entry
                // compare to select the better entry
                // we compare the flags first, then the cost.
                auto flag0 = static_cast<Flag>(fit_flag.coeff(c0));
                auto cost0 = fit_cost.coeff(c0);
                [[maybe_unused]] auto f0 = fit_f.coeff(c0);
                if (((flag == Flag::Good) && (flag0 != Flag::Good)) ||
                    ((flag == flag0) && (cost < cost0)) ||
                    (static_cast<int>(flag) < static_cast<int>(flag0))) {
                    SPDLOG_TRACE("replace c={} f={} flag={} cost={} "
                                 "with c={} f={} flag={} cost={}",
                                 c0, f0, flag0, cost0, c, f, flag, cost);
                    *cit = c;
                }
            } else {
                SPDLOG_TRACE("add c={}", c);
                cands.push_back(c);
            }
        }
        std::sort(cands.begin(), cands.end());
        SPDLOG_INFO("selected {}/{} cands", cands.size(), output.rows());
        SPDLOG_TRACE("selected cands{}", cands);
        return cands;
    };
    auto iterdetect = [&](auto &iterdata, auto accept_flag,
                          auto thresh) -> std::optional<IterStep> {
        SPDLOG_INFO("iteratively detect with accept_flag={} thresh={}",
                    accept_flag, thresh);
        // detect
        std::vector<Segment> segments;
        std::vector<Index> cands;
        {
            // we need to weigth the thresh by the freqs in order
            // to account for the 1/f factor in the height of d21
            auto selectmask =
                (iterdata.array() >
                 (adiqsmean.array() + adiqsstd.array() * thresh * adiqs_fcor.array()))
                    .eval();
            SPDLOG_TRACE("selectmask{}", selectmask);
            std::vector<Index> selected; // index that is above the thresh
            selected.reserve(TULA_SIZET(nrfs));
            for (Index i = 0; i < nrfs; ++i) {
                if (selectmask.coeff(i)) {
                    selected.push_back(i);
                }
            }
            // merge the consecutive index to one segments
            for (auto i : selected) {
                if (segments.empty() || segments.back().second < i) {
                    segments.emplace_back(i, i + 1); // [begin, past-last] range
                } else {
                    ++(segments.back().second);
                }
            }
            SPDLOG_TRACE("segments{}", segments);
            // filter the segments and compute the centers using argmax
            cands.reserve(segments.size());
            segments.erase(std::remove_if(segments.begin(), segments.end(),
                                          [&](const auto &v) {
                                              const auto &[il, ir1] = v;
                                              auto n = ir1 - il;
                                              // reject small segment
                                              if (n < detect_min_segment_size) {
                                                  return true;
                                              }
                                              // i is the peak index in to
                                              // resample fs
                                              auto i =
                                                  tula::alg::argmax(iterdata.segment(
                                                      il, ir1 - il)) +
                                                  il;
                                              // reject if max is at edge
                                              if ((i == il) || (i == ir1 - 1)) {
                                                  return true;
                                              }
                                              cands.push_back(i);
                                              return false;
                                          }),
                           segments.end());
            assert(segments.size() == cands.size());
            SPDLOG_TRACE("found {} segments", segments.size());
        }
        // terminate the iteration if no candidates found
        if (cands.empty()) {
            return std::nullopt;
        }
        Index ncands = static_cast<Index>(cands.size());
        // the found candidates have a resolution of rfstep
        // do fit with found candidates, first we need
        // create data object for fitter, which is effectively
        // an targetsweep composed from the vnasweep data
        SweepFitter::TargetSweepData candsdata;
        {
            candsdata.meta = {
                {"source", fmt::format("iterdetect ncands={}", ncands)}};
            Eigen::MatrixXd toneparams(1, ncands);
            std::vector<std::string> toneparamsheader{"f_tone"};
            for (Index c = 0; c < ncands; ++c) {
                toneparams.coeffRef(0, c) = rfs.coeff(cands[TULA_SIZET(c)]);
            }
            candsdata.wcs.tone_axis = kids::ToneAxis{
                std::move(toneparams), std::move(toneparamsheader)};
            // put tone to center
            candsdata.wcs.sweep_axis = kids::SweepAxis(
                rfs.head(fit_size).array() - rfs.coeff((fit_size - 1) / 2));
            candsdata.iqs().resize(fit_size, ncands);
        }
        auto candindices = tula::container_utils::index(ncands);
        // this is to go over all the candidates, and locate the
        // original s21 data, and put them into candsdata.
        // note that the frequency grid of candsdata is aligned with rfs
        grppi::map(ex, candindices, candindices, [&](auto c) {
            // find itone that enclose the candidate tone, sorting by the
            // distance to the tone center
            // so this will be the tone with larger coverage
            auto icand = cands[c]; // index to resampled fs grid
            auto ftone = candsdata.tones().coeff(c);
            auto itone{-1};
            {
                auto sort_itone = [&](auto lhs, auto rhs) {
                    // sort itone according to distance of f_tone to its center
                    SPDLOG_TRACE("lhs={} pos={} dist={} rhs={} pos={} dist={}",
                                 lhs, tula::alg::argeq(data.fs().col(lhs), ftone),
                                 std::abs(data.tones().coeff(lhs) - ftone), rhs,
                                 tula::alg::argeq(data.fs().col(rhs), ftone),
                                 std::abs(data.tones().coeff(rhs) - ftone));
                    return std::abs(data.tones().coeff(lhs) - ftone) <
                           std::abs(data.tones().coeff(rhs) - ftone);
                };
                std::set<Index, decltype(sort_itone)> itones(sort_itone);
                for (Index i = 0; i < ntones; ++i) {
                    if ((ftone >= toneranges.coeff(0, i)) &&
                        (ftone < toneranges.coeff(1, i))) {
                        itones.insert(i);
                    }
                }
                if (itones.empty()) {
                    throw std::runtime_error(fmt::format(
                        "unable to get tone index for candidate {} f={}", c,
                        ftone));
                }
                itone = *itones.begin();
                assert(itone >= 0);
                SPDLOG_TRACE("locate candidate[{}] icand={} f_tone={} itone={} "
                             "out of tones={}",
                             c, icand, ftone, itone, itones);
            }
            // populate iqs, here we need to interpolate the iqs to
            // cdata fs
            candsdata.iqs().col(c).real() =
                tula::alg::interp(data.fs().col(itone), iqs.col(itone).real().eval(),
                            candsdata.fs().col(c));
            candsdata.iqs().col(c).imag() =
                tula::alg::interp(data.fs().col(itone), iqs.col(itone).imag().eval(),
                            candsdata.fs().col(c));
            return c;
        });
        // do the fit
        auto candsfitresult = fitter(candsdata);
        // candsfitresult.plot();

        // subtract the model and get the residual for next iteration
        Eigen::VectorXd residual(iterdata);
        // we force the slope term to zero in order to perserve the d21
        // continuum level so that the original local mean stats could stay
        // valid
        Eigen::MatrixXd s21ps(candsfitresult.ps);
        s21ps.row(7).setZero();
        s21ps.row(8).setZero();
        Eigen::MatrixXd d21ps(candsfitresult.d21ps);
        d21ps.row(5).setZero();
        d21ps.row(6).setZero();

        // compute and accumulate adiqsmdl from all the good fits, then
        // subtract from iterdata
        //
        // we need to only subtract the good, well separated fits
        // to avoid over subtracting overlapped kids when bad fit
        // occurs
        auto icands_good =
            select_tones(candsfitresult, deblend_tone_dist, accept_flag
                         // Flag::Good
            );
        Eigen::VectorXd adiqsmdl = Eigen::VectorXd::Zero(nrfs);
        for (auto c : icands_good) {
            auto i = cands[TULA_SIZET(c)]; // index to resampled fs grid
            // here we need to evaluate at a size that is proportional
            // to the FWHM
            auto fwhm =
                rfs[i] / config.get_typed<double>("fitter_weight_window_Qr");
            auto eval_size =
                Index(config.get_typed<int>("eval_n_fwhms") * fwhm / rfstep);
            auto [il, ir1] = tula::alg::windowindex_fixedsize(i, nrfs, eval_size);
            SPDLOG_TRACE("compute adiqsmdl at i={} window=[{}:{}]", i, il, ir1);
            auto fsmdl = rfs.segment(il, eval_size);
            Eigen::VectorXcd iqsmdl(eval_size);
            tula::alg::ceresfit::eval<SweepFitter::model_t<fit_modelspec>>(
                fsmdl, s21ps.col(c), iqsmdl);
            Eigen::VectorXcd diqsmdl(eval_size);
            tula::alg::gradient(iqsmdl, fsmdl, diqsmdl);
            adiqsmdl.segment(il, eval_size).array() += diqsmdl.array().abs();
        }
        // subtract from original
        residual -= adiqsmdl;
        // fix the negative value
        // these negative values is beneficial to supress the exesss of
        // of noise due to high peak and thte high noise associated
        // adiqsres = (adiqsres.array() < 0).select(0., adiqsres);
        SPDLOG_TRACE("residual{}", residual);
        Eigen::VectorXd iterdata_in(iterdata);
        iterdata = std::move(residual);
        // adiqs.noalias() = adiqsres.array().isNaN().select(adiqs, adiqsres);
        // update adiqs
        // assign residual to adiqs for next iteration
        // store the iterstep data
        return IterStep{thresh,
                        std::move(segments),
                        std::move(cands),
                        std::move(icands_good),
                        std::move(iterdata_in),
                        iterdata,
                        std::move(candsfitresult)};
    };

    // run the iterative detection
    auto detect_max_niter =
        static_cast<std::size_t>(config.get_typed<int>("detect_max_niter"));
    auto detect_min_thresh = config.get_typed<double>("detect_min_thresh");
    const double iter_thresh_factor = 0.8;
    auto next_thresh = [&](auto thresh) {
        // next thresh
        const auto eps = 0.1; // use this to signal the iter loop to end
        thresh *= iter_thresh_factor;
        if (thresh < detect_min_thresh) {
            thresh = detect_min_thresh - eps;
        }
        return thresh;
    };
    const auto iter_accept_flag =
        Flag::D21LargeOffset | Flag::D21NotConverged; // | Flag::LargeOffset
    const auto iter_accept_flag_last =
        Flag::D21LargeOffset | Flag::D21NotConverged; //  | Flag::LargeOffset;

    auto ntones_design = data.meta.get_typed<int>("ntones_design");
    auto ntones_max = data.meta.get_typed<int>("ntones_max");

    // stores per iteration data
    std::vector<IterStep> itersteps;
    // make a copy of the original adiqs for subtraction
    Eigen::VectorXd iter_adiqs(adiqs);
    // group the found candidates by its rfs index
    // iteration will stop when this list stops grow
    //       rfs index                   iter   cands index
    //        i                          t        c
    std::map<Index, std::vector<std::pair<Index, Index>>> cands;
    // logic to select what frequency to use as the output
    auto locate_cand = [&rfs](const auto &iterstep, auto c) {
        const auto &output = iterstep.candsfitresult.output;
        auto freq = output.coeff(c, 0);
        if (auto flag = static_cast<Flag>(int(output.coeff(c, 2)));
            flag & (Flag::LargeOffset)) {
            SPDLOG_TRACE("fit failed c={} flag={}", c, flag);
            freq = output.coeff(c, 1);
            SPDLOG_TRACE("failed output{}", output.row(c));
        }
        // i is the closest index in rfs to fitted freq if the fitting
        // does not have large offset. otherwise use f_in
        auto [i, eps] = tula::alg::argeq(rfs, freq);
        return std::make_pair(i, freq);
    };
    auto iter_thresh{initial_thresh};
    bool is_last_iter = false;
    while (itersteps.size() < detect_max_niter) {
        if (is_last_iter) {
            SPDLOG_INFO("this is the last iteration");
            break;
        }
        auto t = itersteps.size(); // the current iter index
        if (auto opt_it = iterdetect(iter_adiqs,
                                     is_last_iter ? iter_accept_flag_last
                                                  : iter_accept_flag,
                                     iter_thresh);
            opt_it.has_value()) {

            itersteps.push_back(opt_it.value());
            const auto &it = itersteps.back();
            auto size0 = cands.size(); // store the old size
            // go through the good cands and add to the detected list
            for (auto c : it.icands_good) {
                auto [i, f] = locate_cand(it, c);
                if (auto cit = cands.find(i); cit != cands.end()) {
                    // append to the list
                    cit->second.emplace_back(t, c);
                } else {
                    // create the list
                    cands[i] = {{t, c}};
                }
            }
            //             if (is_last_iter) {
            //                 SPDLOG_INFO("iter {} ncands {} -> {} (last
            //                 iter)", t, size0,
            //                             cands.size());
            //             }
            // check if the cands size increased
            if (auto size = cands.size(); size > size0) {
                SPDLOG_INFO("iter {} ncands {} -> {}", t, size0, size);
                if (size > TULA_SIZET(ntones_design)) {
                    SPDLOG_INFO("all candidates found, stop iter at "
                                "iter={} ncands={}",
                                t, size);
                    is_last_iter = true;
                }
            } else {
                // no new candidate
                if (iter_thresh < detect_min_thresh) {
                    SPDLOG_INFO("no more candidates could be found with "
                                "min_thresh={}, stop iter at "
                                "iter={} ncands={}",
                                detect_min_thresh, t, size);
                    is_last_iter = true;
                }
            }
        } else {
            if (iter_thresh < detect_min_thresh) {
                SPDLOG_INFO("no more candidates could be found with "
                            "min_thresh={}, stop iter at "
                            "iter={} ncands={}",
                            detect_min_thresh, t, cands.size());
                break;
            }
        }
        iter_thresh = next_thresh(iter_thresh);
    }
    if (auto niter = itersteps.size(); niter < detect_max_niter) {
        SPDLOG_INFO("detect done after {} iterations", niter);
    } else {
        SPDLOG_WARN("detect stop at maxiter={}", niter);
    }
    if (cands.size() > TULA_SIZET(ntones_max)) {
        // SPDLOG_ERROR("too many detections {} found", cands.size());
        throw std::runtime_error(
            fmt::format("too many detections {} found", cands.size()));
    }
    // itersteps.back().candsfitresult.plot();
    // group the found candidates by their frequencies
    auto sort_f_close = [&](const auto &lhs, const auto &rhs) {
        auto [il, tl, cl] = lhs;
        auto [ir, tr, cr] = rhs;
        //         auto [il_, fl] = locate(itersteps[SIZET(tl)], cl);
        //         auto [ir_, fr] = locate(itersteps[SIZET(tr)], cr);
        //         assert(il_ == il);
        //         assert(ir_ == ir);
        auto fl = itersteps[TULA_SIZET(tl)].candsfitresult.output.coeff(cl, 0);
        auto fr = itersteps[TULA_SIZET(tr)].candsfitresult.output.coeff(cr, 0);
        return std::abs(fl - fr) > deblend_tone_dist;
    };
    std::set<Candidate, decltype(sort_f_close)> unique_candidates(sort_f_close);
    auto diags = [&itersteps](const auto &e) {
        // return some useful diagnostics
        auto [i, t, c] = e;
        const auto &fitresult = itersteps[TULA_SIZET(t)].candsfitresult;
        // f_out, flag, cost
        return std::make_tuple(fitresult.output.coeff(c, 0),
                               int(fitresult.output.coeff(c, 2)),
                               fitresult.costs.coeff(1, c));
    };
    for (const auto &candidate : cands) {
        auto i = candidate.first;
        for (const auto &[t, c] : candidate.second) {
            auto e = std::make_tuple(i, t, c); // the pending new entry
            // check if exists equivalent entry already
            if (auto u = std::find_if(
                    unique_candidates.begin(), unique_candidates.end(),
                    [&](const auto &u) { return !sort_f_close(u, e); });
                u != unique_candidates.end()) {
                [[maybe_unused]] const auto &[ui, ut, uc] = *u;
                auto [ufreq, uflag, ucost] = diags(*u);
                auto [freq, flag, cost] = diags(e);
                if ((flag < uflag) || ((flag == uflag) && (cost < ucost))) {
                    SPDLOG_TRACE(
                        "replace {} {} {} freq={} flag={} cost={} with {} {} "
                        "{} freq={} flag={} cost={}",
                        ui, ut, uc, ufreq, uflag, ucost, i, t, c, freq, flag,
                        cost);
                    unique_candidates.erase(u);
                    unique_candidates.insert(e);
                }
            } else {
                SPDLOG_TRACE("add {} {} {} to set", i, t, c);
                unique_candidates.insert(e);
            }
        }
    }
    auto unique_cands = tula::container_utils::create<std::vector<Candidate>>(
        std::move(unique_candidates));
    std::sort(
        unique_cands.begin(), unique_cands.end(),
        [&itersteps](const auto lhs, const auto rhs) {
            auto [il, tl, cl] = lhs;
            auto [ir, tr, cr] = rhs;
            auto fl =
                itersteps[TULA_SIZET(tl)].candsfitresult.output.coeff(cl, 0);
            auto fr =
                itersteps[TULA_SIZET(tr)].candsfitresult.output.coeff(cr, 0);
            return fl < fr;
        });
    // make an empty iterstep in case no entry is found
    if (itersteps.empty()) {
        itersteps.emplace_back();
    }
    // generate output
    auto ncands = unique_cands.size();
    SPDLOG_INFO("unique candidates: {}/{}", ncands, ntones_design);

    Eigen::MatrixXd output(ncands,
                           itersteps.front().candsfitresult.output.cols());
    for (std::size_t j = 0; j < ncands; ++j) {
        auto [i, t, c] = unique_cands[j];
        output.row(static_cast<Index>(j)) =
            itersteps[TULA_SIZET(t)].candsfitresult.output.row(c);
    }
    SweepKidsFinderResult result;
    result.data = std::move(data);
    result.rfs = std::move(rfs);
    result.adiqs = std::move(adiqs);
    result.adiqscov = std::move(adiqscov);
    result.stats_rfs = stats_rfs;
    result.adiqs_fcor = adiqs_fcor;
    result.adiqsmean = adiqsmean;
    result.adiqsstd = adiqsstd;
    result.threshold = detect_min_thresh;
    result.itersteps = itersteps;
    result.candidates = cands;
    result.unique_candidates = unique_cands;
    result.output = output;
    return result;
}

#include "kids/sweep/fitter.h"
#include "kids/sweep/model.h"
#include <tula/algorithm/ei_ceresfitter.h>
#include <tula/algorithm/ei_stats.h>
#include <tula/datatable.h>
#include <tula/filename.h>
#include <tula/nc.h>
#include <yaml-cpp/yaml.h>
#ifdef WITH_PLOTTING
#include <matplotlibcpp.h>
#endif


void kids::SweepFitResult::save(const std::string &filepath_) const {
    auto _0 = tula::logging::scoped_timeit("save result");
    // parse path
    auto filepath = tula::filename_utils::parse_pattern(
        filepath_, data.meta.get_str("source"), fmt::arg("ext", "txt"),
        fmt::arg("suffix", ""));
    SPDLOG_TRACE("save to result to {}", filepath);
    try {
        constexpr auto format = datatable::Format::ecsv;
        std::vector<std::pair<std::string, std::string>> _{
            {"roachid", "Header.Toltec.RoachIndex"},
            {"obsid", "Header.Toltec.ObsNum"},
            {"subobsid", "Header.Toltec.SubObsNum"},
            {"scanid", "Header.Toltec.ScanNum"},
            {"accumlen", "Header.Toltec.AccumLen"}};
        // build a map with the keys replaced
        using meta_t = std::decay_t<decltype(data.meta)>;
        meta_t::storage_t meta_out;
        for (const auto &p : _) {
            meta_out[p.second] = data.meta.at(p.first);
        }
        datatable::write<format>(filepath, output, colnames, std::vector<int>{},
                                 tula::ecsv::map_to_meta(std::move(meta_out)));
        SPDLOG_INFO("finished writing file {}", filepath);
    } catch (const datatable::DumpError &e) {
        SPDLOG_ERROR("unable to write to file {}: {}", filepath, e.what());
    }
}

void kids::SweepFitResult::save_processed(const std::string &filepath_) const {
    using Eigen::Index;
    auto _0 = tula::logging::scoped_timeit("save processed data");
    // parse path
    auto filepath = tula::filename_utils::parse_pattern(
        filepath_, data.meta.get_str("source"), fmt::arg("ext", "nc"),
        fmt::arg("suffix", "_processed"));
    SPDLOG_TRACE("save processed data to {}", filepath);

    // underscore shorthand for translations
    std::unordered_map<std::string_view, std::string> _ = {
        {"is", "Data.Kids.Is"},
        {"qs", "Data.Kids.Qs"},
        {"fs", "Data.Kids.fs"},
        {"tones_orig", "Header.Toltec.ToneFreq"},
        {"tones", "Header.Kids.tones"},
        {"sweeps", "Header.Kids.sweeps"},
        // dimentions
        {"nsweeps", "nsweeps"},
        {"ntones", "ntones"},
        {"strlen", "strlen"},
        // meta
        {"kind", "Header.Kids.kind"},
        {"source", "Header.Kids.source"},
        // retained original headers
        // settings
        {"flo_center", "Header.Toltec.LoCenterFreq"},
        {"fsmp", "Header.Toltec.SampleFreq"},
        {"atten_sense", "Header.Toltec.SenseAtten"},
        {"atten_drive", "Header.Toltec.DriveAtten"},
        // info
        {"kindvar", "Header.Toltec.ObsType"},
        {"roachid", "Header.Toltec.RoachIndex"},
        {"obsid", "Header.Toltec.ObsNum"},
        {"subobsid", "Header.Toltec.SubObsNum"},
        {"scanid", "Header.Toltec.ScanNum"},
        // association
        {"cal_roachid", "Header.Toltec.TargSweepRoachIndex"},
        {"cal_obsid", "Header.Toltec.TargSweepObsNum"},
        {"cal_subobsid", "Header.Toltec.TargSweepSubObsNum"},
        {"cal_scanid", "Header.Toltec.TargSweepScanNum"},
    };

    // attach some more meta data
    typename KidsData<>::meta_t meta{data.meta};
    meta.set("kind", "processed_sweep");

    // data geometry
    Index nsweeps, ntones;
    std::tie(nsweeps, ntones) = tula::alg::shape(data.iqs());

    using netCDF::NcDim;
    using netCDF::NcFile;
    using netCDF::NcType;
    using netCDF::NcVar;
    using namespace netCDF::exceptions;

    auto make_row_major = [](auto &&m) {
        using XprType = std::decay_t<decltype(m)>;
        using RMat =
            Eigen::Matrix<typename XprType::Scalar, XprType::RowsAtCompileTime,
                          XprType::ColsAtCompileTime,
                          (XprType::ColsAtCompileTime == 1) ? Eigen::ColMajor
                                                            : Eigen::RowMajor>;
        RMat result;
        result.resize(m.rows(), m.cols());
        Eigen::Map<RMat>(result.data(), m.rows(), m.cols()) = m;
        return result;
    };
    try {
        // use replace to do clobber
        NcFile fo(filepath, NcFile::FileMode::replace,
                  NcFile::FileFormat::classic);
        NcDim d_nsweeps = fo.addDim(_["nsweeps"], TULA_SIZET(nsweeps));
        NcDim d_ntones = fo.addDim(_["ntones"], TULA_SIZET(ntones));
        //
        // put data
        SPDLOG_TRACE("write fs -> {}", _["fs"]);
        NcVar v_fs =
            fo.addVar(_["fs"], netCDF::ncDouble, {d_nsweeps, d_ntones});
        v_fs.putVar(make_row_major(data.fs()).data());
        SPDLOG_TRACE("write is -> {}", _["is"]);
        NcVar v_is =
            fo.addVar(_["is"], netCDF::ncDouble, {d_nsweeps, d_ntones});
        v_is.putVar(make_row_major(data.iqs().real()).data());

        SPDLOG_TRACE("write qs -> {}", _["qs"]);
        NcVar v_qs =
            fo.addVar(_["qs"], netCDF::ncDouble, {d_nsweeps, d_ntones});
        v_qs.putVar(make_row_major(data.iqs().imag()).data());

        {
            SPDLOG_TRACE("write sweeps -> {}", _["sweeps"]);
            NcVar v_sweeps = fo.addVar(_["sweeps"], netCDF::ncDouble, {d_nsweeps});
            Eigen::VectorXd sweeps = data.sweeps();
            v_sweeps.putVar(sweeps.data());
        }

        {
            SPDLOG_TRACE("write tones -> {}", _["tones"]);
            NcVar v_tones = fo.addVar(_["tones"], netCDF::ncDouble, {d_ntones});
            Eigen::VectorXd tones = data.tones();
            v_tones.putVar(tones.data());
            SPDLOG_INFO("tones: {}", tones);
        }

        {
            SPDLOG_TRACE("write tones_orig -> {}", _["tones_orig"]);
            NcVar v_tones_orig = fo.addVar(_["tones_orig"], netCDF::ncDouble, {d_ntones});
            Eigen::VectorXd tones_orig = data.tones().array() - meta.get_typed<double>("flo_center");
            v_tones_orig.putVar(tones_orig.data());
            SPDLOG_INFO("tones_orig: {}", tones_orig);
        }

        for (const auto &key : {"flo_center", "fsmp", "atten_sense", "atten_drive"}) {
            NcVar v = fo.addVar(_[key], netCDF::ncDouble, std::vector<NcDim>{});
            v.putVar(&std::get<double>(meta.at(key)));
        }
        NcDim d_strlen = fo.addDim(_["strlen"], 1024);
        for (const auto &key : {"kind", "source"}) {
            NcVar v = fo.addVar(_[key], netCDF::ncChar, d_strlen);
            v.putVar(std::get<std::string>(meta.at(key)).c_str());
        }
        SPDLOG_TRACE("meta: {}", meta.pformat());
        // info
        for (const auto &key :
             {"kindvar", "roachid", "obsid", "subobsid", "scanid",
              "cal_roachid", "cal_obsid", "cal_subobsid", "cal_scanid"}) {
            SPDLOG_TRACE("add meta {} -> {}", key, _[key]);
            NcVar v = fo.addVar(_[key], netCDF::ncInt, std::vector<NcDim>{});
            v.putVar(&std::get<int>(meta.at(key)));
        }

        fo.sync();
        SPDLOG_INFO("finished writing file {}", filepath);
    } catch (NcException &e) {
        SPDLOG_ERROR("{}", e.what());
        throw std::runtime_error(
            fmt::format("failed to write data to NetCDF file {}", filepath));
    }
}

void kids::SweepFitResult::plot() const {
    #ifdef WITH_PLOTTING
    namespace eiu = tula::eigen_utils;
    namespace plt = matplotlibcpp;
    using Eigen::Index;

    // data geometry
    Index nsweeps, ntones;
    std::tie(nsweeps, ntones) = tula::alg::shape(data.iqs());
    Eigen::VectorXd sfs = data.fs().col(0).array() - data.tones().coeff(0);
    [[maybe_unused]] auto sweepspan = tula::alg::span(sfs);
    [[maybe_unused]] auto sweepstep = tula::alg::step(sfs);
    [[maybe_unused]] auto tonespan = tula::alg::span(data.tones());

    // setup canvas
    auto source = data.meta.get_str("source");
    SPDLOG_DEBUG("generate plot for {}", source);
    // auto lofreq = sweepdata.emta.get<double>("lofreq");
    std::string funit = "Hz";
    std::string iqunit = "DN";
    plt::close(); // flush

    Index ntonesperpage = 6;
    Index npanels = 8;
    auto npages = ntones / ntonesperpage + 1;

    for (Index ipage = 0; ipage < npages; ++ipage) {
        plt::figure_size(200 * TULA_SIZET(npanels), 200 * TULA_SIZET(ntonesperpage));
        plt::suptitle(fmt::format("{} [{}:{}]", source, ipage * ntonesperpage,
                                  (ipage + 1) * ntonesperpage));
        for (Index itone_ = 0; itone_ < ntonesperpage; ++itone_) {
            Index itone = ipage * ntonesperpage + itone_;
            if (itone >= ntones) {
                break;
            }
            // sweep data
            auto xvec = eiu::to_stdvec(data.fs().col(itone));
            auto ivec = eiu::to_stdvec(data.iqs().col(itone).real());
            auto qvec = eiu::to_stdvec(data.iqs().col(itone).imag());
            auto wivec =
                eiu::to_stdvec(1. / uncertainty.col(itone).array().real().square());
            auto wqvec =
                eiu::to_stdvec(1. / uncertainty.col(itone).array().imag().square());
            auto avec = eiu::to_stdvec(data.iqs().col(itone).array().abs().log10());
            auto pvec = eiu::to_stdvec(data.iqs().col(itone).array().arg());
            auto advec = eiu::to_stdvec(diqs.col(itone).array().abs());
            auto fix_phase = [](auto &&m) {
                auto pi = 4. * std::atan(1);
                if ((m.maxCoeff() - m.minCoeff()) > 1.5 * pi) {
                    // fix phase angle jump
                    for (Index i = 1; i < m.size(); ++i) {
                        if (abs(m(i) - m(i - 1)) > pi) {
                            m(i) += (m(i) > m(i - 1) ? -1. : 1.) * 2. * pi;
                        }
                    }
                }
            };
            fix_phase(eiu::as_eigen(pvec));
            // other outputs
            auto fout = output(itone, 0);
            auto fin = output(itone, 1);
            [[maybe_unused]] auto flag = output(itone, 2);

            // model
            std::vector<double> imdl, qmdl, amdl, pmdl, admdl;
            auto dispatch_model = [&](auto model_) {
                constexpr auto model_v = TULA_DECAY(model_)::value;
                using Model = kids::Model<model_v>;
                // SPDLOG_DEBUG("plot with model: {}", model_v);
                Eigen::VectorXcd iqmdl(nsweeps);
                //                 Eigen::VectorXd params(Model::NP);
                //                 params << data.tfs.coeff(itone), 3e4, 5e4,
                //                     data.tfs.coeff(itone), 0., 0., 0.,
                //                     EIGEN_PI;
                //                 tula::alg::ceresfit::eval<Model>(eiu::asvec(xvec),
                //                 params, iqmdl);

                tula::alg::ceresfit::eval<Model>(eiu::as_eigen(xvec), ps.col(itone),
                                           iqmdl);
                SPDLOG_TRACE("[{}] s21 fit: {} -> {}", itone, costs(0, itone),
                             costs(1, itone));
                imdl = eiu::to_stdvec(iqmdl.real());
                qmdl = eiu::to_stdvec(iqmdl.imag());
                amdl = eiu::to_stdvec(iqmdl.array().abs().log10());
                pmdl = eiu::to_stdvec(iqmdl.array().arg());
                fix_phase(eiu::as_eigen(pmdl));
                // d21 mdl
                Eigen::VectorXcd diqmdl(nsweeps);
                tula::alg::ceresfit::eval<kids::Model<kids::SweepModel::D21>>(
                    eiu::as_eigen(xvec), d21ps.col(itone), diqmdl);
                SPDLOG_TRACE("[{}] d21 fit: {} -> {}", itone,
                             d21costs(0, itone), d21costs(1, itone));
                admdl = eiu::to_stdvec(diqmdl.array().abs());
                SPDLOG_TRACE("[{}] fitresult flag: {:l}", itone,
                             static_cast<SweepFitResult::Flag>(int(flag)));

                // plt kwargs
                std::map<std::string, std::string> vecstyle = {
                    {"markersize", "3"}, {"color", "C0"}, {"marker", "o"}};
                std::map<std::string, std::string> vec0style = {
                    {"markersize", "6"}, {"color", "C0"}, {"marker", "o"}};
                std::map<std::string, std::string> mdlstyle = {
                    {"linestyle", "-"}, {"color", "C1"}, {"marker", ""}};
                std::map<std::string, std::string> mdl0style = {
                    {"markersize", "6"}, {"color", "C1"}, {"marker", "o"}};
                std::map<std::string, std::string> resstyle = {
                    {"linestyle", "-"}, {"color", "C4"}, {"marker", ""}};
                std::map<std::string, std::string> res0style = {
                    {"markersize", "6"}, {"color", "C4"}, {"marker", "o"}};
                std::map<std::string, std::string> lastyle = {
                    {"markersize", "6"}, {"color", "C2"}, {"marker", "x"}};
                auto plot0 = [](const auto &x, const auto &y, const auto &style,
                                const auto &style0) {
                    // SPDLOG_TRACE("plot x{} y{}", x, y);
                    // plt::plot(std::vector<double>{0, 1}, std::vector<double>{0, 1}, style);
                    plt::plot(x, y, style);
                    plt::plot({x[0]}, {y[0]}, style0);
                    SPDLOG_TRACE("plot x0={} y0={}", x[0], y[0]);
                };
                auto plotf = [&]() {
                    plt::axvline(fin, 0., 1.,
                                 {{"color", "#c3c3c3"}, {"linestyle", "--"}});
                    plt::axvline(fout, 0., 1., {{"color", "r"}});
                    plt::axvline(data.tones().coeff(itone) -
                                     fin / window_Qr / 2., 0., .1,
                                 {{"color", "#e3e3e3"}, {"linestyle", ":"}});
                    plt::axvline(data.tones().coeff(itone) +
                                     fin / window_Qr / 2., 0., 1.,
                                 {{"color", "#e3e3e3"}, {"linestyle", ":"}});
                };
                // panel 1: Ampl
                plt::subplot(ntonesperpage, npanels, npanels * itone_ + 1);
                plt::ylabel(fmt::format("Log |S21| ({})", iqunit));
                plot0(xvec, avec, vecstyle, vec0style);
                plot0(xvec, amdl, mdlstyle, mdl0style);
                plotf();
                // panel 1: Ampl d21
                plt::subplot(ntonesperpage, npanels, npanels * itone_ + 2);
                plt::ylabel(fmt::format("|D21| ({})", iqunit));
                plot0(xvec, advec, vecstyle, vec0style);
                plot0(xvec, admdl, mdlstyle, mdl0style);
                plot0(xvec,
                      eiu::to_stdvec(
                          (diqs.col(itone).array() - diqmdl.array()).abs()),
                      resstyle, res0style);
                plotf();

                // panel 3: I
                plt::subplot(ntonesperpage, npanels, npanels * itone_ + 3);
                plt::ylabel(fmt::format("I"));
                plot0(xvec, ivec, vecstyle, vec0style);
                plot0(xvec, imdl, mdlstyle, mdl0style);
                plotf();
                // panel 4: Q
                plt::subplot(ntonesperpage, npanels, npanels * itone_ + 4);
                plt::ylabel(fmt::format("Q"));
                plot0(xvec, qvec, vecstyle, vec0style);
                plot0(xvec, qmdl, mdlstyle, mdl0style);
                plotf();
                // panel 5: weight I
                plt::subplot(ntonesperpage, npanels, npanels * itone_ + 5);
                plt::ylabel(fmt::format("I"));
                plot0(xvec, wivec, vecstyle, vec0style);
                plotf();
                // panel 6: weight Q
                plt::subplot(ntonesperpage, npanels, npanels * itone_ + 6);
                plt::ylabel(fmt::format("Q"));
                plot0(xvec, wqvec, vecstyle, vec0style);
                plotf();
                // panel 7: complex plane
                plt::subplot(ntonesperpage, npanels, npanels * itone_ + 7);
                plt::ylabel(fmt::format("I-Q complex plane"));
                plot0(ivec, qvec, vecstyle, vec0style);
                plot0(imdl, qmdl, mdlstyle, mdl0style);
                // get ifout and ifin
                auto ifin = tula::alg::argmin((eiu::as_eigen(xvec).array() - fin).abs());
                auto ifout =
                    tula::alg::argmin((eiu::as_eigen(xvec).array() - fout).abs());
                plt::plot({ivec[TULA_SIZET(ifin)]}, {qvec[TULA_SIZET(ifin)]},
                          {{"color", "#d3d3d3"},
                           {"marker", "o"},
                           {"markersize", "6"}});
                plt::plot(
                    {ivec[TULA_SIZET(ifout)]}, {qvec[TULA_SIZET(ifout)]},
                    {{"color", "r"}, {"marker", "o"}, {"markersize", "6"}});
            };
            using ModelSpec = SweepFitter::ModelSpec;
            switch (this->modelspec) {
            case ModelSpec::gain: {
                dispatch_model(tula::meta::scalar_t<kids::SweepModel::S21WithGain>{});
                break;
            }
            case ModelSpec::gainlintrend: {
                dispatch_model(
                    tula::meta::scalar_t<kids::SweepModel::S21WithGainLinTrend>{});
                break;
            }
            case ModelSpec::trans: {
                dispatch_model(
                    tula::meta::scalar_t<kids::SweepModel::S21WithTrans>{});
                break;
            }
            case ModelSpec::translintrend: {
                dispatch_model(
                    tula::meta::scalar_t<kids::SweepModel::S21WithTransLinTrend>{});
                break;
            }
            }
        }
        plt::xlabel(fmt::format("Freqency ({})", funit));
        plt::show();
    }
    // plt::detail::_interpreter::kill();
    #endif
}

//
// void save_plot(const std::string &filepath) const {
//     SPDLOG_TRACE("save plot to {}", filepath);
// }

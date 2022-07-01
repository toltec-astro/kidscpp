#include "kids/sweep/finder.h"
#include <tula/algorithm/ei_stats.h>
#include <tula/eigen.h>
#include <tula/datatable.h>
#include <tula/filename.h>
#include <tula/meta.h>
#ifdef WITH_PLOTTING
#include <matplotlibcpp.h>
#endif
#include <tula/nc.h>

using kids::SweepKidsFinderResult;

void SweepKidsFinderResult::save(const std::string &filepath_) const {
    using Eigen::Index;
    // tone file path
    auto filepath = tula::filename_utils::parse_pattern(
        filepath_, data.meta.get_str("source"), fmt::arg("ext", "txt"),
        fmt::arg("suffix", ""));
    SPDLOG_TRACE("save tones to {}", filepath);

    try {
        constexpr auto format = datatable::Format::ecsv;
        YAML::Node meta;
        std::vector<std::pair<std::string, std::string>> _{
            {"roachid", "Header.Toltec.RoachIndex"},
            {"obsid", "Header.Toltec.ObsNum"},
            {"subobsid", "Header.Toltec.SubObsNum"},
            {"scanid", "Header.Toltec.ScanNum"},
            {"accumlen", "Header.Toltec.AccumLen"}};
        using meta_t = std::decay_t<decltype(data.meta)>;
        meta_t::storage_t meta_out;
        for (const auto &p : _) {
            meta_out[p.second] = data.meta.at(p.first);
        }
        auto colnames = itersteps.front().candsfitresult.colnames;
        datatable::write<format>(filepath, output, colnames, std::vector<int>{},
                                 tula::ecsv::map_to_meta(std::move(meta_out)));
        SPDLOG_INFO("finished writing file {}", filepath);
    } catch (const datatable::DumpError &e) {
        SPDLOG_ERROR("unable to write to file {}: {}", filepath, e.what());
    }
}

void SweepKidsFinderResult::save_d21(const std::string &filepath_) const {
    using Eigen::Index;
    auto _0 = tula::logging::scoped_timeit("save finder d21 result");
    // parse path
    auto filepath = tula::filename_utils::parse_pattern(
        filepath_, data.meta.get_str("source"), fmt::arg("ext", "nc"),
        fmt::arg("suffix", "_d21"));
    SPDLOG_TRACE("save d21 to {}", filepath);

    // underscore shorthand for translations
    std::unordered_map<std::string_view, std::string> _ = {
        {"fs", "fs"},
        {"adiqs", "adiqs"},
        {"adiqscov", "adiqscov"},
        {"adiqsmean", "adiqsmean"},
        {"adiqsstd", "adiqsstd"},
        {"candidates", "candidates"},
        // dementions
        {"nfs", "xlen"},
        {"ncands", "clen"},
        {"strlen", "strlen"},
        // meta
        {"kind", "kind"},
        {"source", "source"},
        // retained original headers
        // settings
        {"flo_center", "Header.Toltec.LoCenterFreq"},
        {"fsmp", "Header.Toltec.SampleFreq"},
        {"atten_sense", "Header.Toltec.SenseAtten"},
        {"atten_drive", "Header.Toltec.DriveAtten"},
    };

    // attach some more meta data
    typename KidsData<>::meta_t meta{data.meta};
    meta.set("kind", "d21");

    // data geometry
    Index nfs;
    std::tie(nfs, std::ignore) = tula::alg::shape(adiqs);

    using netCDF::NcDim;
    using netCDF::NcFile;
    using netCDF::NcType;
    using netCDF::NcVar;
    using namespace netCDF::exceptions;

    try {
        // use replace to do clobber
        NcFile fo(filepath, NcFile::replace);
        NcDim d_nfs = fo.addDim(_["nfs"], TULA_SIZET(nfs));
        //
        // put data
        NcVar v_fs = fo.addVar(_["fs"], netCDF::ncDouble, {d_nfs});
        v_fs.putVar(rfs.data());
        NcVar v_adiqs = fo.addVar(_["adiqs"], netCDF::ncDouble, {d_nfs});
        v_adiqs.putVar(adiqs.data());
        NcVar v_adiqscov = fo.addVar(_["adiqscov"], netCDF::ncDouble, {d_nfs});
        v_adiqscov.putVar(adiqscov.data());
        NcVar v_adiqsmean =
            fo.addVar(_["adiqsmean"], netCDF::ncDouble, {d_nfs});
        v_adiqsmean.putVar(adiqsmean.data());
        NcVar v_adiqsstd = fo.addVar(_["adiqsstd"], netCDF::ncDouble, {d_nfs});
        v_adiqsstd.putVar(adiqsstd.data());

        NcDim d_ncands = fo.addDim(_["ncands"], TULA_SIZET(output.rows()));
        NcVar v_candidates =
            fo.addVar(_["candidates"], netCDF::ncDouble, {d_ncands});
        if (output.rows() > 0) {
            Eigen::VectorXd candidates = output.col(0);
            v_candidates.putVar(candidates.data());
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
        SPDLOG_INFO("finished writing file {}", filepath);
    } catch (NcException &e) {
        SPDLOG_ERROR("{}", e.what());
        throw std::runtime_error(
            fmt::format("failed to write data to NetCDF file {}", filepath));
    }
}

void SweepKidsFinderResult::save_processed(const std::string &filepath_) const {
    using Eigen::Index;
    auto _0 = tula::logging::scoped_timeit("save processed data");
    // parse path
    auto filepath = tula::filename_utils::parse_pattern(
        filepath_, data.meta.get_str("source"), fmt::arg("ext", "nc"),
        fmt::arg("suffix", "_processed"));
    SPDLOG_TRACE("save processed data to {}", filepath);

    // underscore shorthand for translations
    std::unordered_map<std::string_view, std::string> _ = {
        // d21
        {"d21_fs", "Data.Kids.d21_fs"},
        {"d21_adiqs", "Data.Kids.d21_adiqs"},
        {"d21_adiqscov", "Data.Kids.d21_adiqscov"},
        {"d21_adiqsmean", "Data.Kids.d21_adiqsmean"},
        {"d21_adiqsstd", "Data.Kids.d21_adiqsstd"},
        {"candidates", "Header.Kids.candidates"},
        // d21 dimensions
        {"nd21fs", "nd21fs"},
        {"ncands", "ncands"},
        // sweep
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

    Index n_d21_fs;
    std::tie(n_d21_fs, std::ignore) = tula::alg::shape(adiqs);

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
        NcFile fo(filepath, NcFile::replace,
                  NcFile::FileFormat::classic);
        NcDim d_nsweeps = fo.addDim(_["nsweeps"], TULA_SIZET(nsweeps));
        NcDim d_ntones = fo.addDim(_["ntones"], TULA_SIZET(ntones));

        NcDim d_nd21fs = fo.addDim(_["nd21fs"], TULA_SIZET(n_d21_fs));
        //
        // put data
        SPDLOG_TRACE("write d21_fs -> {}", _["d21_fs"]);
        NcVar v_d21_fs = fo.addVar(_["d21_fs"], netCDF::ncDouble, {d_nd21fs});
        v_d21_fs.putVar(rfs.data());
        SPDLOG_TRACE("write d21_adiqs -> {}", _["d21_adiqs"]);
        NcVar v_d21_adiqs =
            fo.addVar(_["d21_adiqs"], netCDF::ncDouble, {d_nd21fs});
        v_d21_adiqs.putVar(adiqs.data());
        SPDLOG_TRACE("write d21_adiqscov -> {}", _["d21_adiqscov"]);
        NcVar v_d21_adiqscov =
            fo.addVar(_["d21_adiqscov"], netCDF::ncDouble, {d_nd21fs});
        v_d21_adiqscov.putVar(adiqscov.data());
        SPDLOG_TRACE("write d21_adiqsmean -> {}", _["d21_adiqsmean"]);
        NcVar v_d21_adiqsmean =
            fo.addVar(_["d21_adiqsmean"], netCDF::ncDouble, {d_nd21fs});
        v_d21_adiqsmean.putVar(adiqsmean.data());
        SPDLOG_TRACE("write d21_adiqsstd -> {}", _["d21_adiqsstd"]);
        NcVar v_d21_adiqsstd =
            fo.addVar(_["d21_adiqsstd"], netCDF::ncDouble, {d_nd21fs});
        v_d21_adiqsstd.putVar(adiqsstd.data());
        SPDLOG_TRACE("write candicates -> {}", _["candidates"]);
        NcDim d_ncands = fo.addDim(_["ncands"], TULA_SIZET(output.rows()));
        NcVar v_candidates =
            fo.addVar(_["candidates"], netCDF::ncDouble, {d_ncands});
        if (output.rows() > 0) {
            Eigen::VectorXd candidates = output.col(0);
            v_candidates.putVar(candidates.data());
        }
        // put sweep data
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

void SweepKidsFinderResult::plot() const {
    #ifdef WITH_PLOTTING
    namespace eiu = tula::eigen_utils;
    namespace plt = matplotlibcpp;
    using Eigen::Index;

    // data geometry
    Index nfs;
    std::tie(nfs, std::ignore) = tula::alg::shape(adiqs);
    auto rfsvec = eiu::to_stdvec(rfs);
    auto adiqsvec = eiu::to_stdvec(adiqs);
    auto adiqsmeanvec = eiu::to_stdvec(adiqsmean);
    auto adiqs1sigvec = eiu::to_stdvec(adiqsmean + adiqsstd);
    // double flo = data.meta.get<double>("lofreq");

    // setup canvas
    auto source = data.meta.get_str("source");
    SPDLOG_TRACE("generate plot for {}", source);

    std::string iqunit = "$DN$";
    std::string diqunit = "$DN/Hz$";
    std::string funit = "$Hz$";
    const auto s21label = fmt::format("$|S_{{21}}|$ ({})", iqunit);
    const auto d21label = fmt::format("$|dS_{{21}}/df|$ ({})", diqunit);

    plt::close(); // flush
    plt::figure_size(2000, 1000);
    plt::suptitle(fmt::format("{}", source));

    // plt kwargs
    std::map<std::string, std::string> vecstyle = {{"linewidth", "0.1"},
                                                   {"color", "C0"}};
    std::map<std::string, std::string> vecstyle0 = {{"linewidth", "0.1"},
                                                    {"color", "C1"}};
    std::map<std::string, std::string> resstyle = {{"linewidth", "0.5"},
                                                   {"color", "#cc6600"}};
    std::map<std::string, std::string> seg0style = {{"linewidth", "1"}};
    std::map<std::string, std::string> candstyle = {
        {"linewidth", "1"}, {"color", "#cccccc"}, {"linestyle", ":"}};
    std::map<std::string, std::string> fitstyle = {
        {"linewidth", "1"}, {"color", "#cc0000"}, {"linestyle", ":"}};

    std::map<std::string, std::string> baselinestyle = {
        {"linewidth", "1"}, {"color", "#cccccc"}, {"linestyle", "-"}};
    std::map<std::string, std::string> statsstyle = {{"linestyle", "none"},
                                                     {"color", "#aaaaaa"},
                                                     {"marker", "o"},
                                                     {"markersize", "2"}};
    std::map<std::string, std::string> meanstyle = {
        {"linewidth", "1."}, {"color", "#66cc00"}, {"linestyle", "--"}};
    std::map<std::string, std::string> stdstyle = {
        {"linewidth", "1."}, {"color", "#66ccff"}, {"linestyle", "--"}};
    std::map<std::string, std::string> cutstyle = {
        {"linewidth", "1."}, {"color", "#aaaa00"}, {"linestyle", "-"}
        // {"marker", "."}
    };
    // plot each iter as one panel
    // upper two panel: plot iqs and adiqs
    // auto niter = static_cast<Index>(itersteps.size());
    // const auto iter0 = niter - 3;
    // const auto iter0 = 0;
    const auto iter0 = 0;
    const auto niter = 3;
    auto npanels = 2 + niter - iter0;
    {
        plt::subplot(npanels, 1, 1);
        plt::ylabel(s21label);
        // plt::plot(fsvec, eiu::tostd(adiqs), vecstyle);
        for (Index i = 0; i < data.tones().size(); ++i) {
            auto fsvec = eiu::to_stdvec(data.fs().col(i));
            auto iqsvec = eiu::to_stdvec(data.iqs().col(i).array().abs());
            plt::plot(fsvec, iqsvec, (i % 2 == 0) ? vecstyle : vecstyle0);
        }
        if (output.rows() > 0) {
            for (auto u : unique_candidates) {
                auto [i, t, c] = u;
                plt::axvline(
                    itersteps[TULA_SIZET(t)].candsfitresult.output.coeff(c, 0),
                    0., 1.,
                    candstyle);
                // plt::axvline(
                //    itersteps[SIZET(t)].candsfitresult.output.coeff(c, 1),
                //    fitstyle);
            }
        }
        plt::subplot(npanels, 1, 2, true);
        plt::ylabel(d21label);
        // plt::plot(fsvec, eiu::tostd(adiqs), vecstyle);
        plt::plot(rfsvec, adiqsvec, vecstyle);
        plt::axhline(0., 0., 1., baselinestyle);
        plt::plot(rfsvec, adiqsmeanvec, meanstyle);
        plt::plot(rfsvec, adiqs1sigvec, stdstyle);
        plt::plot(rfsvec,
                  eiu::to_stdvec(adiqsmean.array() + adiqsstd.array()  * threshold * adiqs_fcor.array()),
                  cutstyle);
    }
    if (output.rows() > 0) {
        // niter = 0; // disable the iterstep plots
        // plot residual
        const auto &residual = itersteps.back().residual;
        plt::plot(rfsvec, eiu::to_stdvec(residual), resstyle);
        plt::ylim(-0.1, adiqsstd.maxCoeff() * 10 + adiqsmean.maxCoeff());
        for (auto u : unique_candidates) {
            auto [i, t, c] = u;
            plt::axvline(itersteps[TULA_SIZET(t)].candsfitresult.output.coeff(c, 0),
                    0., 1.,
                         candstyle);
            plt::axvline(itersteps[TULA_SIZET(t)].candsfitresult.output.coeff(c, 1),
                    0., 1.,
                         fitstyle);
        }
        for (Index i = iter0; i < niter; ++i) {
            auto ip = i - iter0;
            plt::subplot(npanels, 1, ip + 3, true);
            const auto &it = itersteps[TULA_SIZET(i)];
            // plt::plot(fsvec, eiu::tostd(it.adiqs), vecstyle);
            //             plt::plot(eiu::tostd(stats_fs),
            //                       eiu::tostd(Eigen::VectorXd::Zero(stats_fs.size())),
            //                       statsstyle);
            plt::plot(rfsvec, adiqsmeanvec, meanstyle);
            plt::plot(rfsvec, adiqs1sigvec, stdstyle);
            plt::plot(rfsvec,
                      eiu::to_stdvec(adiqsmean.array() + adiqsstd.array() * it.thresh * adiqs_fcor.array()),
                      cutstyle);
            plt::plot(rfsvec, adiqsvec, vecstyle);
            plt::plot(rfsvec, eiu::to_stdvec(it.residual), resstyle);
            plt::plot(rfsvec, eiu::to_stdvec(it.adiqsmdl), resstyle);
            std::map<std::string, std::string> segstyle{seg0style};
            // per segment coloring
            for (std::size_t j = 0; j < it.segments.size(); ++j) {
                const auto &[si, ei] = it.segments[j];
                auto color = fmt::format("#{:02x}{:02x}{:02x}", (ip + 1) * 0x10,
                                         (j % 2 == 0) ? 0xbb : 0x55,
                                         (j % 2 == 1) ? 0xbb : 0x55);
                // SPDLOG_TRACE("color={}", color);
                segstyle["color"] = color;
                plt::plot(eiu::to_stdvec(rfs.segment(si, ei - si)),
                          eiu::to_stdvec(it.iterdata.segment(si, ei - si)),
                          segstyle);
                plt::axvline(rfs.coeff(it.cands[j]), 0., 1., candstyle);
                // plot fitresult
                plt::axvline(it.candsfitresult.output.coeff(j, 1), 0., 1., fitstyle);
                plt::text(
                    it.candsfitresult.output.coeff(j, 1), -10.,
                    fmt::format("{}",
                                static_cast<SweepFitResult::Flag>(
                                    it.candsfitresult.output.coeff(j, 2))));
            }
            plt::text(4.5e8, 20., fmt::format("ncand={}", it.segments.size()));
        }
    }
    plt::show();
    // itersteps.back().fitresult.plot();
    #endif
}

void SweepKidsFinderResult::save_plot(const std::string &filepath) const {}

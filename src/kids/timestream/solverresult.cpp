#include "kids/timestream/solver.h"
#include <tula/algorithm/ei_stats.h>
#include <tula/filename.h>
#include <tula/formatter/matrix.h>
#include <tula/logging.h>
// #include <tula/matplotlibcpp.h>

namespace {
// underscore shorthand for translations
std::unordered_map<std::string_view, std::string> _ = {
    // processed ts data
    {"xs", "Data.Kids.xs"},
    {"rs", "Data.Kids.rs"},
    {"its", "Data.Kids.its"},
    {"qts", "Data.Kids.qts"},
    {"phs", "Data.Kids.phs"},
    {"tones", "Header.Toltec.ToneFreq"},
    {"fs", "Header.Kids.tones"},
    {"ts", "Data.Toltec.Ts"},
    {"psd_fbins", "Header.Kids.PsdFreq"},
    {"is_psd", "Data.Kids.ispsd"},
    {"qs_psd", "Data.Kids.qspsd"},
    {"phs_psd", "Data.Kids.phspsd"},
    {"xs_psd", "Data.Kids.xspsd"},
    {"rs_psd", "Data.Kids.rspsd"},

    // calibration data
    // {"fitreport", "fitreport"},
    {"modelspec", "modelspec"},
    // {"fitreportheader", "fitreportheader"},
    // dementions
    {"ntimes", "ntimes"},
    {"ntones", "ntones"},
    {"npsdfs", "npsdfs"},
    {"strlen", "strlen"},
    {"ntimecols", "ntimecols"},
    // meta
    {"kind", "kind"},
    {"source", "source"},
    {"exit_0", "exit_0"},
    // retained ts headers
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
    // assoc
    {"cal_file", "cal_file"},
    {"cal_roachid", "Header.Toltec.TargSweepRoachIndex"},
    {"cal_obsid", "Header.Toltec.TargSweepObsNum"},
    {"cal_scanid", "Header.Toltec.TargSweepScanNum"},
    {"cal_subobsid", "Header.Toltec.TargSweepSubObsNum"}};

using netCDF::NcDim;
using netCDF::NcFile;
using netCDF::NcType;
using netCDF::NcVar;

auto add_var = [](auto &fo, const std::string &key, auto &&type_,
                  auto &&...args) {
    if (const auto &v = fo.getVar(_[key]); !v.isNull()) {
        SPDLOG_INFO("get var {} {}", key, _[key]);
        return v;
    }
    SPDLOG_INFO("add var {} {}", key, _[key]);
    return std::forward<decltype(fo)>(fo).addVar(
        _[key], std::forward<decltype(type_)>(type_),
        std::vector<netCDF::NcDim>{std::forward<decltype(args)>(args)...});
};

} // namespace

void kids::TimeStreamSolverResult::save(const std::string &filepath_) const {
    // parse path
    auto filepath = tula::filename_utils::parse_pattern(
        filepath_, data.meta.get_str("source"), fmt::arg("ext", "nc"),
        fmt::arg("suffix", "_processed"));
    SPDLOG_TRACE("save result to {}", filepath);
    NcFileIO io{filepath};
    save_to_nc(io, false);
}

kids::TimeStreamSolverResult::NcFileIO &
kids::TimeStreamSolverResult::save_to_nc(NcFileIO &io, bool header_only) const {

    // attach some more meta data
    typename KidsData<>::meta_t meta{data_out.meta};
    // meta.set("cal",
    //         fmt::format("toltec{:d}_{:06d}_{:02d}_{:04d}", data.cal[0],
    //                     data.cal[1], data.cal[2], data.cal[3]));
    auto &fo = io.file_obj();
    const auto &filepath = io.filepath();
    meta.set("source", filepath);

    using Eigen::Index;

    // data geometry
    Index ntimes, ntones;
    std::tie(ntimes, ntones) = tula::alg::shape(data_out.xs());
    Index ntimecols = data.wcs.time_axis.data.cols();

    using namespace netCDF::exceptions;

    auto check_data_range = [](auto name, const auto &m) {
        auto r = m.maxCoeff() - m.minCoeff();
        SPDLOG_INFO("data range of {}: {}", name, r);
    };
    try {
        // define the dimensions.
        // NcDim d_ntimes = fo.addDim(_["ntimes"], SIZET(ntimes));
        NcDim d_ntimes = fo.addDim(_["ntimes"]); // unlimited
        NcDim d_ntones = fo.addDim(_["ntones"], TULA_SIZET(ntones));
        NcDim d_ntimecols = fo.addDim(_["ntimecols"], TULA_SIZET(ntimecols));
        NcDim d_strlen = fo.addDim(_["strlen"], 1024);

        NcVar v_fs = add_var(fo, "fs", netCDF::ncDouble, d_ntones);
        check_data_range("fs", fs);
        v_fs.putVar(fs.data());
        fo.sync();

        {
            NcVar v_tones = add_var(fo, "tones", netCDF::ncDouble, d_ntones);
            Eigen::VectorXd tones =
                fs.array() - meta.get_typed<double>("flo_center");
            check_data_range("tones", tones);
            v_tones.putVar(tones.data());
            fo.sync();
        }

        if (extra_output) {
            Index npsdfs = psd_fbins.size();
            NcDim d_npsdfs = fo.addDim(_["npsdfs"], TULA_SIZET(npsdfs));

            NcVar v_fs_psd =
                add_var(fo, "psd_fbins", netCDF::ncDouble, d_npsdfs);
            check_data_range("psd_fbins", psd_fbins);
            v_fs_psd.putVar(psd_fbins.data());
            fo.sync();

            NcVar v_is_psd =
                add_var(fo, "is_psd", netCDF::ncDouble, d_npsdfs, d_ntones);
            check_data_range("is_psd", is_psd);
            v_is_psd.putVar(is_psd.data());
            fo.sync();
            NcVar v_qs_psd =
                add_var(fo, "qs_psd", netCDF::ncDouble, d_npsdfs, d_ntones);
            check_data_range("qs_psd", qs_psd);
            v_qs_psd.putVar(qs_psd.data());
            fo.sync();
            NcVar v_phs_psd =
                add_var(fo, "phs_psd", netCDF::ncDouble, d_npsdfs, d_ntones);
            check_data_range("phs_psd", phs_psd);
            v_phs_psd.putVar(phs_psd.data());
            fo.sync();
            NcVar v_xs_psd =
                add_var(fo, "xs_psd", netCDF::ncDouble, d_npsdfs, d_ntones);
            check_data_range("xs_psd", xs_psd);
            v_xs_psd.putVar(xs_psd.data());
            fo.sync();
            NcVar v_rs_psd =
                add_var(fo, "rs_psd", netCDF::ncDouble, d_npsdfs, d_ntones);
            check_data_range("rs_psd", rs_psd);
            v_rs_psd.putVar(rs_psd.data());
            fo.sync();
        }

        for (const auto &key :
             {"flo_center", "fsmp", "atten_sense", "atten_drive"}) {
            NcVar v = add_var(fo, key, netCDF::ncDouble);
            v.putVar(&std::get<double>(meta.at(key)));
            fo.sync();
        }

        for (const auto &key : {"kind", "modelspec", "source", "cal_file"}) {
            NcVar v = add_var(fo, key, netCDF::ncChar, d_strlen);
            v.putVar(std::get<std::string>(meta.at(key)).data());
            fo.sync();
        }
        for (const auto &key :
             {"kindvar", "roachid", "obsid", "subobsid", "scanid",
              "cal_roachid", "cal_obsid", "cal_subobsid", "cal_scanid"}) {
            NcVar v = add_var(fo, key, netCDF::ncInt);
            v.putVar(&std::get<int>(meta.at(key)));
            fo.sync();
        }
        if (header_only) {
            SPDLOG_INFO("finished writing file header {}", filepath);
            fo.sync();
            return io;
        }
        append_to_nc(io);

        NcVar v = add_var(fo, "exit_0", netCDF::ncInt);
        int exit_0 = 0;
        v.putVar(&exit_0);
        SPDLOG_INFO("finished writing file {}", filepath);
    } catch (NcException &e) {
        SPDLOG_ERROR("{}", e.what());
        throw std::runtime_error(
            fmt::format("failed to write data to NetCDF file {}", filepath));
    }
    fo.sync();
    return io;
}

kids::TimeStreamSolverResult::NcFileIO &
kids::TimeStreamSolverResult::append_to_nc(NcFileIO &io) const {

    using Eigen::Index;

    using netCDF::NcDim;
    using netCDF::NcFile;
    using netCDF::NcType;
    using netCDF::NcVar;
    using namespace netCDF::exceptions;

    auto &fo = io.file_obj();

    try {

        // define the dimensions.
        // NcDim d_ntimes = fo.addDim(_["ntimes"], SIZET(ntimes));
        NcDim d_ntimes = fo.getDim(_["ntimes"]); // unlimited
        NcDim d_ntones = fo.getDim(_["ntones"]);
        NcDim d_ntimecols = fo.getDim(_["ntimecols"]);
        if (d_ntimes.isNull() || d_ntones.isNull() || d_ntimecols.isNull()) {
            save_to_nc(io, true);
            d_ntimes = fo.getDim(_["ntimes"]);
            d_ntones = fo.getDim(_["ntones"]);
            d_ntimecols = fo.getDim(_["ntimecols"]);
        }
        auto ntimes_exists = d_ntimes.getSize();
        SPDLOG_INFO("append to existing ntimes={}", ntimes_exists);

        // std::vector<std::size_t> s{SIZET(ntimes), SIZET(ntones)};
        NcVar v_xs = add_var(fo, "xs", netCDF::ncDouble, d_ntimes, d_ntones);
        NcVar v_rs = add_var(fo, "rs", netCDF::ncDouble, d_ntimes, d_ntones);
        NcVar v_ts = add_var(fo, "ts", netCDF::ncInt, d_ntimes, d_ntimecols);

        std::vector<std::size_t> i0{ntimes_exists, 0};
        std::vector<std::size_t> s_d{1, d_ntones.getSize()};
        std::vector<std::size_t> s_t{1, d_ntimecols.getSize()};

        const auto &tdata = data.wcs.time_axis.data;
        assert(d_ntimecols.getSize() == SIZET(tdata.cols()));
        assert(d_ntones.getSize() == SIZET(data_out.xs.data.cols()));
        auto ntimes = data_out.xs.data.rows();
        {
            tula::logging::progressbar pb0(
                [](const auto &msg) { SPDLOG_INFO("{}", msg); }, 60,
                "write data ");

            auto put_data = [&](auto ii) mutable {
                i0[0] = ntimes_exists + ii;
                v_xs.putVar(i0, s_d, data_out.xs().row(ii).data());
                v_rs.putVar(i0, s_d, data_out.rs().row(ii).data());
                v_ts.putVar(i0, s_t, tdata.row(ii).data());
            };

            if (extra_output) {
                NcVar v_its =
                    add_var(fo, "its", netCDF::ncDouble, d_ntimes, d_ntones);
                NcVar v_qts =
                    add_var(fo, "qts", netCDF::ncDouble, d_ntimes, d_ntones);
                NcVar v_phs =
                    add_var(fo, "phs", netCDF::ncDouble, d_ntimes, d_ntones);
                // v_its.putVar(i, s, its.data());
                // v_qts.putVar(i, s, qts.data());
                // v_phs.putVar(i, s, phs.data());

                for (std::size_t ii = 0; ii < TULA_SIZET(ntimes); ++ii) {
                    put_data(ii);
                    v_its.putVar(i0, s_d, its.row(ii).data());
                    v_qts.putVar(i0, s_d, qts.row(ii).data());
                    v_phs.putVar(i0, s_d, phs.row(ii).data());
                    pb0.count(ntimes, ntimes / 10);
                }
                /*
                 */
            } else {
                // v_xs.putVar(i, s, data_out.xs().data());
                // v_rs.putVar(i, s, data_out.rs().data());

                for (std::size_t ii = 0; ii < TULA_SIZET(ntimes); ++ii) {
                    put_data(ii);
                    pb0.count(ntimes, ntimes / 10);
                }
                /*
                 */
            }
        }
    } catch (NcException &e) {
        SPDLOG_ERROR("{}", e.what());
        throw std::runtime_error(fmt::format("failed to append data to file {}",
                                             tula::nc_utils::pprint(fo)));
    }
    fo.sync();
    return io;
}

void kids::TimeStreamSolverResult::plot() const {
    /*
    namespace eiu = eigen_utils;
    namespace plt = matplotlibcpp;
    using Eigen::Index;

    // data geometry
    Index ntimes, ntones;
    std::tie(ntimes, ntones) = tula::alg::shape(data.is());
    double fsmp = data.meta.get_typed<double>("fsmp");
    Eigen::VectorXd times =
        Eigen::VectorXd::LinSpaced(ntimes, 0, (ntimes - 1) / fsmp);

    // setup canvas
    auto source = data.meta.get_str("source");
    SPDLOG_TRACE("generate plot for {}", source);

    std::string tunit = "s";
    std::string iqunit = "DN";
    std::string xunit = "";
    std::string funit = "Hz";
    std::string psdunit = "DN^2/Hz";
    std::string psdxunit = "Hz^-1";
    // limit the points in the time domain to this number
    Index npts_plot = 2000;
    Index npts_plot_ratio = ntimes / npts_plot;
    std::vector<double> tvec(SIZET(npts_plot));
    eiu::asvec(tvec).setLinSpaced(npts_plot, 0, (ntimes - 1) / fsmp);
    auto fxxvec = eiu::tostd(psd_fbins);

    plt::close(); // flush

    Index ntonesperpage = 6;
    Index npanels = 8;
    auto npages = ntones / ntonesperpage + 1;
    for (Index ipage = 0; ipage < npages; ++ipage) {
        plt::figure_size(800 * SIZET(npanels), 200 * SIZET(ntonesperpage));
        plt::suptitle(fmt::format("{} [{}:{}]", source, ipage * ntonesperpage,
                                  (ipage + 1) * ntonesperpage));
        for (Index irow = 0; irow < ntonesperpage; ++irow) {
            Index itone = ipage * ntonesperpage + irow;
            // iqs data, down sampled
            std::vector<double> ivec(SIZET(npts_plot));
            eiu::asvec(ivec) =
                Eigen::Map<const Eigen::VectorXd, 0, Eigen::InnerStride<>>(
                    &data.is().coeff(0, itone), npts_plot,
                    Eigen::InnerStride<>(ntones * npts_plot_ratio));
            std::vector<double> qvec(SIZET(npts_plot));
            eiu::asvec(qvec) =
                Eigen::Map<const Eigen::VectorXd, 0, Eigen::InnerStride<>>(
                    &data.qs().coeff(0, itone), npts_plot,
                    Eigen::InnerStride<>(ntones * npts_plot_ratio));
            // xs data, down sampled
            std::vector<double> xvec(SIZET(npts_plot));
            eiu::asvec(xvec) =
                Eigen::Map<const Eigen::VectorXd, 0, Eigen::InnerStride<>>(
                    &data_out.xs().coeff(0, itone), npts_plot,
                    Eigen::InnerStride<>(ntones * npts_plot_ratio));
            std::vector<double> rvec(SIZET(npts_plot));
            eiu::asvec(rvec) =
                Eigen::Map<const Eigen::VectorXd, 0, Eigen::InnerStride<>>(
                    &data_out.rs().coeff(0, itone), npts_plot,
                    Eigen::InnerStride<>(ntones * npts_plot_ratio));
            // psd
            auto xpxxvec = eiu::tostd(xs_psd.col(itone));
            auto rpxxvec = eiu::tostd(rs_psd.col(itone));
            auto ipxxvec = eiu::tostd(is_psd.col(itone));
            auto qpxxvec = eiu::tostd(qs_psd.col(itone));

            // plt kwargs
            std::map<std::string, std::string> vecstyle = {
                {"linewidth", "0.1"},
                // {"markersize", "0.1"},
                // {"marker", "."},
                {"color", "C0"}};
            std::map<std::string, std::string> psdstyle = {{"linewidth", "0.5"},
                                                           {"color", "C0"}};
            std::map<std::string, std::string> mdlstyle = {
                {"linestyle", "-"}, {"color", "C1"}, {"marker", ""}};
            std::map<std::string, std::string> mdl0style = {
                {"markersize", "6"}, {"color", "C1"}, {"marker", "o"}};
            std::map<std::string, std::string> lastyle = {
                {"markersize", "6"}, {"color", "C2"}, {"marker", "x"}};
            // plot
            // I - t
            plt::subplot(ntonesperpage, npanels, npanels * irow + 1);
            plt::xlabel(fmt::format("Time ({})", tunit));
            plt::ylabel(fmt::format("I ({})", iqunit));
            plt::plot(tvec, ivec, vecstyle);
            // PSD (I)
            plt::subplot(ntonesperpage, npanels, npanels * irow + 2);
            plt::yscale("log");
            plt::xlabel(fmt::format("Frequency ({})", funit));
            plt::ylabel(fmt::format("I_PSD ({})", psdunit));
            plt::plot(fxxvec, ipxxvec, psdstyle);
            // Q - t
            plt::subplot(ntonesperpage, npanels, npanels * irow + 3);
            plt::ylabel(fmt::format("Q ({})", iqunit));
            plt::plot(tvec, qvec, vecstyle);
            // PSD (Q)
            plt::subplot(ntonesperpage, npanels, npanels * irow + 4);
            plt::yscale("log");
            plt::xlabel(fmt::format("Frequency ({})", funit));
            plt::ylabel(fmt::format("Q_PSD ({})", psdunit));
            plt::plot(fxxvec, qpxxvec, psdstyle);
            // x - t
            plt::subplot(ntonesperpage, npanels, npanels * irow + 5);
            plt::xlabel(fmt::format("Time ({})", tunit));
            plt::ylabel(fmt::format("x ({})", xunit));
            plt::plot(tvec, xvec, vecstyle);
            // PSD (x)
            plt::subplot(ntonesperpage, npanels, npanels * irow + 6);
            plt::yscale("log");
            plt::xlabel(fmt::format("Frequency ({})", funit));
            plt::ylabel(fmt::format("x_PSD ({})", psdxunit));
            plt::plot(fxxvec, xpxxvec, psdstyle);
            // r - t
            plt::subplot(ntonesperpage, npanels, npanels * irow + 7);
            plt::xlabel(fmt::format("Time ({})", tunit));
            plt::ylabel(fmt::format("r ({})", xunit));
            plt::plot(tvec, rvec, vecstyle);
            // PSD (r)
            plt::subplot(ntonesperpage, npanels, npanels * irow + 8);
            plt::yscale("log");
            plt::xlabel(fmt::format("Frequency ({})", funit));
            plt::ylabel(fmt::format("r_PSD ({})", psdxunit));
            plt::plot(fxxvec, rpxxvec, psdstyle);
        }
        plt::show();
    }
    */
}

void kids::TimeStreamSolverResult::save_plot(
    const std::string &filepath) const {
    SPDLOG_ERROR("same plot to {} not implemented", filepath);
}

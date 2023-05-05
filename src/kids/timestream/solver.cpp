#include "kids/timestream/solver.h"
#include "kids/sweep/model.h"
#include "kids/timestream/solver_psd.h"
#include <tula/algorithm/ei_stats.h>
#include <tula/algorithm/index.h>
#include <tula/container.h>
#include <tula/datatable.h>
#include <tula/ecsv/core.h>
#include <tula/filename.h>
#include <tula/grppi.h>
#include <tula/logging.h>

using Eigen::Index;

kids::TimeStreamSolver::TimeStreamSolver(
    const kids::TimeStreamSolver::Config &config_)
    : config(config_) {}

namespace {

// key mapper to retrieve the meta from table header
std::vector<std::pair<std::string, std::string>> _{
    {"roachid", "Header.Toltec.RoachIndex"},
    {"obsid", "Header.Toltec.ObsNum"},
    {"subobsid", "Header.Toltec.SubObsNum"},
    {"scanid", "Header.Toltec.ScanNum"},
    {"accumlen", "Header.Toltec.AccumLen"}};

template <typename Config, typename meta_t>
auto loadfitreport(const Config &config, const meta_t &meta) {
    namespace fs = std::filesystem;
    // this is the pattern expected for the same obsnum
    auto m = [&](const auto &key) { return meta.template get_typed<int>(key); };
    auto pattern0 =
        fmt::format("{}{:d}_{:06d}_{:03d}_0001.+\\.txt",
                    meta.get_str("instru"), m("roachid"),
                    m("obsid"), m("subobsid"));
    // this is the pattern expected for the cal obsnum
    auto pattern1 = meta.get_str("cal_file");
    std::string filepath{};
    if (config.has("fitreportfile")) {
        filepath = config.get_str("fitreportfile");
    } else if (config.has("fitreportdir")) {
        auto dir = config.get_str("fitreportdir");
        SPDLOG_TRACE("look for fitreport dir {} with pattern {}", dir, pattern0);
        auto candidates0 = tula::filename_utils::find_regex(dir, pattern0);
        if (!candidates0.empty()) {
            filepath = candidates0[0];
        } else {
            SPDLOG_TRACE("look for fitreport dir {} with cal pattern {}", dir, pattern1);
            auto candidates1 = tula::filename_utils::find_regex(dir, pattern1);
            if (!candidates1.empty()) {
                filepath = candidates1[0];
            } else {
                throw std::runtime_error(fmt::format(
                    "no fit report found in {} that matches {} or {}", dir, pattern0, pattern1));
            }
        }
    } else {
        throw std::runtime_error(
            fmt::format("no fit report location specified."));
    }
    SPDLOG_INFO("use fitreport file {}", filepath);
    std::vector<std::string> header;
    Eigen::MatrixXd table;
    meta_t meta_cal{};

    try {
        YAML::Node meta_;
        table = datatable::read<double, datatable::Format::ecsv>(
            filepath, &header, &meta_);
        auto meta_map =
            tula::ecsv::meta_to_map<typename meta_t::storage_t::key_type,
                                    typename meta_t::storage_t::mapped_type>(
                meta_, &meta_);
        meta_cal = meta_t{std::move(meta_map)};
        if (!meta_.IsNull()) {
            SPDLOG_WARN("un recongnized meta:\n{}", YAML::Dump(meta_));
        }
    } catch (datatable::ParseError &e) {
        SPDLOG_WARN("unable to read fitreport file as ECSV {}: {}", filepath,
                    e.what());
        try {
            table = datatable::read<double, datatable::Format::ascii>(filepath,
                                                                      &header);
        } catch (datatable::ParseError &e) {
            SPDLOG_WARN("unable to read fitreport file as ASCII {}: {}",
                        filepath, e.what());
            throw e;
        }
    }
    SPDLOG_INFO("meta_cal: {}", meta_cal.pformat());
    return std::tuple{
        kids::ToneAxis(std::move(table).transpose(), std::move(header)),
        std::move(meta_cal)};
}

} // namespace

kids::TimeStreamSolverResult kids::TimeStreamSolver::operator()(
    const kids::TimeStreamSolver::RawTimeStreamData &data,
    const kids::TimeStreamSolver::Config &config_) {
    // merge the config object
    Config config{this->config};
    config.update(config_);
    SPDLOG_DEBUG("run kids solver with config: {}", config.pformat());
    using meta_t = KidsData<>::meta_t;
    // try load fitreport
    kids::ToneAxis tone_axis;
    meta_t meta_cal;
    try {
        // tula::logging::scoped_loglevel<spdlog::level::trace> l0;
        std::tie(tone_axis, meta_cal) = loadfitreport(config, data.meta);
    } catch (std::runtime_error &e) {
        SPDLOG_WARN("unable to load fitreport file: {}", e.what());
        SPDLOG_INFO("use built-in model params");
        tone_axis = data.wcs.tone_axis;
    }
    auto &tone_data = tone_axis.data;
    auto &tone_data_header = tone_axis.row_labels.labels();
    auto fs = tone_axis("f_in");

    SPDLOG_TRACE("tone data header: {}", tone_data_header);
    SPDLOG_TRACE("tone data{}", tone_data);
    SPDLOG_TRACE("tone fs{}", fs);
    SPDLOG_TRACE("meta_cal{}", meta_cal.pformat());
    // data geometry
    Index ntimes, ntones;
    std::tie(ntimes, ntones) = tula::alg::shape(data.is());
    double fsmp = data.meta.get_typed<double>("fsmp");
    const Index nfft = 2048;
    auto [npts, nfs, df] = ::internal::fftstat(nfft, fsmp);

    // ex policy
    auto toneindices = tula::container_utils::index(ntones);
    auto ex = tula::grppi_utils::dyn_ex(config.get_str("exmode"));
    SPDLOG_TRACE("use grppi ex mode {}", config.get_str("exmode"));

    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    using internal::RMatrixXd;
    // data objects
    RMatrixXd xs(ntimes, ntones);
    RMatrixXd rs(ntimes, ntones);
    RMatrixXd its(ntimes, ntones);
    RMatrixXd qts(ntimes, ntones);

    // optional
    auto extra_output = config.get_typed<bool>("extra_output");
    SPDLOG_INFO("compute extra output: {}", extra_output ? "yes" : "no");
    RMatrixXd phs;

    VectorXd fxx;
    MatrixXd ipxx;
    MatrixXd qpxx;
    MatrixXd phpxx;
    MatrixXd xpxx;
    MatrixXd rpxx;

    constexpr auto modelspec = kids::SweepModel::S21WithGainLinTrend;
    using Model = kids::Model<modelspec>;
    // make mdl param on axis 1
    auto mdl = tone_data.array().bottomRows(Model::NP);
    //  0  1   2   3  4  5      6       7       8       9          10
    // fp Qr, Qc, fr, A, normI, normQ, slopeI, slopeQ, interceptI, interceptQ
    SPDLOG_TRACE("mdl{}", mdl);

    {
        tula::logging::scoped_timeit t0{"calc x"};

        if (extra_output) {
            // raw phase angle
            phs.resize(ntimes, ntones);
        }

        // make chunks
        std::vector<std::pair<Index, Index>> chunks;
        const Index chunksize = std::min(4880l, ntimes);
        for (Index i = 0; i + chunksize <= ntimes; i += chunksize) {
            chunks.emplace_back(i, i + chunksize);
        }
        if (auto i = chunks.back().second; i < ntimes) {
            chunks.back().second = ntimes;
        }
        auto nchunks = chunks.size();
        SPDLOG_DEBUG("chunks{}", chunks);

        auto atan2 = [](double a, double b) { return std::atan2(a, b); };
        auto norm2 = (mdl.row(5) * mdl.row(5) + mdl.row(6) * mdl.row(6)).eval();
        auto mi = (mdl.row(9) + (fs.array() - mdl.row(0)) * mdl.row(7)).eval();
        auto mq = (mdl.row(10) + (fs.array() - mdl.row(0)) * mdl.row(8)).eval();
        // check accumu len
        double iq_norm{};
        if (meta_cal.is_set("accumlen")) {
            iq_norm = data.meta.get_typed<int>("accumlen") /
                      meta_cal.get_typed<int>("accumlen");
        } else {
            auto accumlen_def = 1 << 19; // 20200326: This is to address the
            // repeated data on taco. This 2^19 is for sample freq of 488 Hz.
            SPDLOG_WARN("found no accumlen in fitreport, use default {}",
                        accumlen_def);
            iq_norm = data.meta.get_typed<int>("accumlen") / accumlen_def;
        }
        // 20200326: we need to check the obsnum, and make sure the old, scaled
        // data does not get scaled again. This applies to obsnum <=10899 and
        // master != repeat (3)
        if ((data.meta.get_typed<int>("obsid") <= 10899) &&
            (data.meta.get_typed<int>("master") != 3)) {
            iq_norm = 1.;
            SPDLOG_WARN("force iq_norm={} for legacy data", iq_norm);
        }
        SPDLOG_INFO("iq_norm={}", iq_norm);

        tula::logging::progressbar pb0(
            [](const auto &msg) { SPDLOG_INFO("{}", msg); }, 60, "calc x ");
        grppi::map(ex, chunks, chunks, [&](auto chunk) {
            auto [i0, i1] = chunk;
            auto block = [i0 = i0, i1 = i1, &ntones](auto &&m) {
                return std::forward<decltype(m)>(m).array().block(
                    i0, 0, i1 - i0, ntones);
            };
            auto rblock = [&block](auto &&m) {
                return block(std::forward<decltype(m)>(m)).rowwise();
            };

            // detrend s21
            auto is1 = (rblock(data.is() / iq_norm) - mi).eval();
            auto qs1 = (rblock(data.qs() / iq_norm) - mq).eval();
            // deroated iqs
            auto is0 = block(its);
            auto qs0 = block(qts);
            is0 = (is1.rowwise() * mdl.row(5) + qs1.rowwise() * mdl.row(6))
                      .rowwise() /
                  norm2;
            qs0 = (qs1.rowwise() * mdl.row(5) - is1.rowwise() * mdl.row(6))
                      .rowwise() /
                  norm2;
            // compute x
            auto denom =
                ((is0 * is0 + qs0 * qs0).rowwise() * mdl.row(2) * 2.).eval();
            block(xs) = (is0.rowwise() * mdl.row(4) - qs0).cwiseQuotient(denom);
            block(rs) = (qs0.rowwise() * mdl.row(4) + is0).cwiseQuotient(denom);
            if (extra_output) {
                // raw phase angle
                block(phs) =
                    block(data.qs()).binaryExpr(block(data.is()), atan2);
            }
            pb0.count(nchunks, nchunks / 10);
            return chunk;
        });
    }
    // calculate psd

    if (extra_output) {
        tula::logging::scoped_timeit t0{"calc psd"};
        fxx = ::internal::fftfs(npts, nfs, df);
        ipxx.resize(nfs, ntones);
        qpxx.resize(nfs, ntones);
        phpxx.resize(nfs, ntones);
        xpxx.resize(nfs, ntones);
        rpxx.resize(nfs, ntones);

        // make chunks
        auto chunks = tula::alg::indexchunks(0l, ntones, 32l, 0l);
        auto nchunks = chunks.size();
        SPDLOG_DEBUG("chunks{}", chunks);

        tula::logging::progressbar pb0(
            [](const auto &msg) { SPDLOG_INFO("{}", msg); }, 60, "calc psd ");
        grppi::map(ex, chunks, chunks, [&](auto chunk) {
            auto [i0, i1] = chunk;
            for (auto c = i0; c < i1; ++c) {
                ::internal::welch(data.is().col(c), ipxx.col(c), fsmp, nfft);
                ::internal::welch(data.qs().col(c), qpxx.col(c), fsmp, nfft);
                ::internal::welch(phs.col(c), phpxx.col(c), fsmp, nfft);
                ::internal::welch(xs.col(c), xpxx.col(c), fsmp, nfft);
                ::internal::welch(rs.col(c), rpxx.col(c), fsmp, nfft);
            }
            pb0.count(nchunks, nchunks / 10);
            return chunk;
        });
    } else {
        // release its and iqs
        its.resize(0, 0);
        qts.resize(0, 0);
    }
    // SPDLOG_INFO("its{:r10c10} qts{:r10c10} phs{:r10c10}", its, qts, phs);
    // SPDLOG_INFO("xs{:r10c10} rs{:r10c10}", xs, rs);

    /*
    {
        logging::progressbar pb0(
            [](const auto &msg) { SPDLOG_INFO("{}", msg); }, 60, "solving ");
        grppi::map(ex, toneindices, toneindices, [&](auto c) {
            pb0.count(ntones, ntones / 10);
            // raw phase angle
            phs.col(c) = data.qs().col(c).binaryExpr(data.is().col(c), atan2);
            // detrend s21
            assert(fs(c) == mdl(0, c));
            auto is1 = (data.is().col(c).array() - mdl(9, c) -
                        mdl(7, c) * (fs(c) - mdl(0, c)));
            auto qs1 = (data.qs().col(c).array() - mdl(10, c) -
                        mdl(8, c) * (fs(c) - mdl(0, c)));
            // divide norm out
            auto norm2 = mdl(5, c) * mdl(5, c) + mdl(6, c) * mdl(6, c);
            // deroated iqs
            its.col(c) = (is1 * mdl(5, c) + qs1 * mdl(6, c)) / norm2;
            qts.col(c) = (qs1 * mdl(5, c) - is1 * mdl(6, c)) / norm2;
            auto is0 = its.col(c).array();
            auto qs0 = qts.col(c).array();
            // compute x
            auto denom = (2. * mdl(2, c) * (is0 * is0 + qs0 * qs0)).eval();
            xs.col(c) = (mdl(4, c) * is0 - qs0).cwiseQuotient(denom);
            rs.col(c) = (mdl(4, c) * qs0 + is0).cwiseQuotient(denom);
            ::internal::welch(data.is().col(c), ipxx.col(c), fsmp, nfft);
            ::internal::welch(data.qs().col(c), qpxx.col(c), fsmp, nfft);
            ::internal::welch(phs.col(c), phpxx.col(c), fsmp, nfft);
            ::internal::welch(xs.col(c), xpxx.col(c), fsmp, nfft);
            ::internal::welch(rs.col(c), rpxx.col(c), fsmp, nfft);
            return c;
        });
    }
    */
    // SPDLOG_INFO("its{:r10c10} qts{:r10c10} phs{:r10c10}", its, qts, phs);
    // SPDLOG_INFO("xs{:r10c10} rs{:r10c10}", xs, rs);

    KidsData<>::meta_t meta{data.meta};
    meta.set("modelspec", std::string(tula::enum_utils::name(modelspec)));
    meta.set("kind", std::string(tula::enum_utils::name(
                         KidsDataKind::SolvedTimeStream)));
    // KidsData<KidsDataKind::SolvedTimeStream> result;
    TimeStreamSolverResult result(data);
    {
        tula::logging::scoped_timeit t0{"populate result"};
        result.data_out.meta = std::move(meta);
        result.data_out.xs.data = std::move(xs);
        result.data_out.rs.data = std::move(rs);
        result.extra_output = extra_output;
        result.its = std::move(its);
        result.qts = std::move(qts);
        result.phs = std::move(phs);
        result.fs = std::move(fs);
        result.psd_fbins = std::move(fxx);
        result.is_psd = std::move(ipxx);
        result.qs_psd = std::move(qpxx);
        result.phs_psd = std::move(phpxx);
        result.xs_psd = std::move(xpxx);
        result.rs_psd = std::move(rpxx);
    }
    return result;
    /*
    return TimeStreamSolverResult{std::move(data),
                                  std::move(result),
                                  extra_output,
                                  std::move(its),
                                  std::move(qts),
                                  std::move(phs),
                                  std::move(fs),
                                  std::move(fxx),
                                  std::move(ipxx),
                                  std::move(qpxx),
                                  std::move(phpxx),
                                  std::move(xpxx),
                                  std::move(rpxx)};
    */
}

#pragma once

#include "../core/kidsdata.h"
#include <tula/config/flatconfig.h>
#include <tula/nc.h>

namespace kids {

struct TimeStreamSolverResult;

struct TimeStreamSolver {
    using RawTimeStreamData = KidsData<KidsDataKind::RawTimeStream>;
    using SolvedTimeStreamData = KidsData<KidsDataKind::SolvedTimeStream>;
    using Config = KidsData<>::meta_t;

    TimeStreamSolver(const Config &config_);

    /// @brief solve timestream
    TimeStreamSolverResult operator()(const RawTimeStreamData &,
                                      const Config &config = {});

    Config config;
};

struct TimeStreamSolverResult {

	TimeStreamSolverResult(const TimeStreamSolver::RawTimeStreamData& data_): data(data_) {}

    const TimeStreamSolver::RawTimeStreamData& data;
    TimeStreamSolver::SolvedTimeStreamData data_out{};
    using RMatrixXd =
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    bool extra_output{false};
    RMatrixXd its{};
    RMatrixXd qts{};
    RMatrixXd phs{};
    Eigen::VectorXd fs{};
    RMatrixXd psd_fbins{};
    RMatrixXd is_psd{};
    RMatrixXd qs_psd{};
    RMatrixXd phs_psd{};
    RMatrixXd xs_psd{};
    RMatrixXd rs_psd{};

    class NcFileIO {
    public:
        using NcFile = netCDF::NcFile;
        NcFileIO() : m_filepath{}, m_fo{NcFile()} {}
        NcFileIO(const std::string &filepath)
            : NcFileIO(filepath, NcFile::FileMode::replace,
                       NcFile::FileFormat::classic) {}
        template <typename... Args>
        NcFileIO(const std::string &filepath, Args &&... args)
            : m_filepath{filepath}, m_fo{NcFile(m_filepath, std::forward<decltype(args)>(args)...)} {}
        auto &file_obj() { return m_fo; }
        const auto &filepath() const { return m_filepath; }

    private:
        std::string m_filepath;
        netCDF::NcFile m_fo;
    };

    NcFileIO &append_to_nc(NcFileIO &io) const;
    NcFileIO &save_to_nc(NcFileIO &io, bool header_only = false) const;
    void save(const std::string &filepath) const;
    void plot() const;
    void save_plot(const std::string &filepath) const;
};

}  // namespace kids

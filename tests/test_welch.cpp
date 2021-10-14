#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <kids/timestream/solver_psd.h>

namespace {

using namespace ::testing;

TEST(timestream, welch) {
    using namespace Eigen;
    VectorXd xdata;
    xdata.setLinSpaced(16, 0, 15);

    int psdsize = 4;
    VectorXd psddata;
    SPDLOG_INFO("before xdata{} psddata{}", xdata, psddata);
    spdlog::set_level(spdlog::level::trace);
    auto stat = ::internal::welch<::internal::Window::Hann>(xdata, psddata,
                                                            500., psdsize);
    SPDLOG_INFO("xdata{} psddata{}", xdata, psddata);
    auto [npts, nfs, df] = stat;
    auto fs = ::internal::fftfs(npts, nfs, df);
    SPDLOG_INFO("psd fs{}", fs);
}

} // namespace

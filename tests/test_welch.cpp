#include <gtest/gtest.h>
#include <kids/timestream/solver_psd.h>

TEST(TimestreamPsd, FftShapeUsesEvenSampleCount) {
    const auto [npts, nfs, df] = ::internal::fftstat(17, 80.0);
    EXPECT_EQ(npts, 16);
    EXPECT_EQ(nfs, 9);
    EXPECT_DOUBLE_EQ(df, 5.0);

    const auto frequencies = ::internal::fftfs(npts, nfs, df);
    ASSERT_EQ(frequencies.size(), nfs);
    EXPECT_DOUBLE_EQ(frequencies(0), 0.0);
    EXPECT_DOUBLE_EQ(frequencies(nfs - 1), 40.0);
}

TEST(TimestreamPsd, WelchProducesFiniteOneSidedSpectrum) {
    Eigen::VectorXd samples(16);
    for (Eigen::Index i = 0; i < samples.size(); ++i) {
        samples(i) = static_cast<double>(i);
    }

    Eigen::VectorXd spectrum;
    const auto [npts, nfs, df] =
        ::internal::welch<::internal::Hann>(samples, spectrum, 8.0, 8);

    EXPECT_EQ(npts, 8);
    EXPECT_EQ(nfs, 5);
    EXPECT_DOUBLE_EQ(df, 1.0);
    ASSERT_EQ(spectrum.size(), nfs);
    EXPECT_TRUE(spectrum.array().isFinite().all());
    EXPECT_TRUE((spectrum.array() >= 0.0).all());
}

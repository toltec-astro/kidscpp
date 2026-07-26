#include <gtest/gtest.h>
#include <kids/timestream/solver.h>

#include <Eigen/Core>

#include <string>
#include <vector>

namespace {

auto make_raw_timestream() -> kids::TimeStreamSolver::RawTimeStreamData
{
    using kids::internal::RMatrixXd;

    kids::TimeStreamSolver::RawTimeStreamData data;
    Eigen::MatrixXd tone_data(12, 1);
    tone_data.col(0) << 100.0, // f_in
        100.0,                // fp
        2.0,                  // Qr
        4.0,                  // Qc
        100.0,                // fr
        2.0,                  // A
        1.0,                  // normI
        0.0,                  // normQ
        0.0,                  // slopeI
        0.0,                  // slopeQ
        0.0,                  // interceptI
        0.0;                  // interceptQ
    data.wcs.tone_axis = kids::ToneAxis{
        std::move(tone_data),
        std::vector<std::string>{
            "f_in",
            "fp",
            "Qr",
            "Qc",
            "fr",
            "A",
            "normI",
            "normQ",
            "slopeI",
            "slopeQ",
            "interceptI",
            "interceptQ",
        }};
    data.is.data = RMatrixXd::Constant(4, 1, 1.0);
    data.qs.data = RMatrixXd::Zero(4, 1);
    data.meta.set("fsmp", 488.0);
    data.meta.set("accumlen", 1 << 19);
    data.meta.set("obsid", 20'000);
    data.meta.set("master", 0);
    return data;
}

} // namespace

TEST(TimestreamSolver, ConvertsIqToDetuningAndDissipation)
{
    const auto data = make_raw_timestream();
    const kids::TimeStreamSolver solver{
        {{"exmode", "seq"}, {"extra_output", false}}};

    const auto result = solver(data);

    ASSERT_EQ(result.data_out.xs.data.rows(), 4);
    ASSERT_EQ(result.data_out.xs.data.cols(), 1);
    EXPECT_TRUE(result.data_out.xs.data.isConstant(0.25));
    EXPECT_TRUE(result.data_out.rs.data.isConstant(0.125));
    EXPECT_FALSE(result.extra_output);
    EXPECT_EQ(
        result.data_out.meta.get_str("modelspec"),
        "S21WithGainLinTrend");
    EXPECT_EQ(result.data_out.meta.get_str("kind"), "SolvedTimeStream");
}

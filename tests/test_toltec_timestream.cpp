#include <kids/toltec/timestream.h>

#include <gtest/gtest.h>

#include <cstdlib>
#include <filesystem>

namespace {

constexpr auto fixture_relative_path =
    "data_lmt/toltec/ics/toltec0/"
    "toltec0_018230_111_0000_2024_05_20_13_40_08_timestream.nc";

[[nodiscard]] auto fixture_path() -> std::filesystem::path
{
    return std::filesystem::path{std::getenv("TOLTECA_TEST_DATA_ROOT")} /
           fixture_relative_path;
}

TEST(ToltecTimeStream, ReadsRealMetadataAndSlice)
{
    if (std::getenv("TOLTECA_TEST_DATA_ROOT") == nullptr) {
        GTEST_SKIP()
            << "set TOLTECA_TEST_DATA_ROOT to run the real-file test";
    }
    const auto path = fixture_path();
    ASSERT_TRUE(std::filesystem::is_regular_file(path)) << path;

    const auto meta = kids::toltec::get_raw_timestream_meta(path);
    EXPECT_EQ(meta.get_typed<int>("roachid"), 0);
    EXPECT_EQ(meta.get_typed<int>("obsid"), 18230);
    EXPECT_EQ(meta.get_typed<int>("subobsid"), 111);
    EXPECT_EQ(meta.get_typed<int>("scanid"), 0);
    EXPECT_EQ(meta.get_typed<int>("ntimes_all"), 610);
    EXPECT_EQ(meta.get_typed<int>("ntones"), 649);
    EXPECT_DOUBLE_EQ(meta.get_typed<double>("fsmp"), 122.0703125);
    EXPECT_EQ(
        meta.get_str("cal_file"),
        R"(toltec0_018228_000_0000.+\.txt)");

    const auto data = kids::toltec::read_raw_timestream_slice(
        path, kids::toltec::SampleSlice{0, 2, 1});
    EXPECT_EQ(data.is.data.rows(), 2);
    EXPECT_EQ(data.is.data.cols(), 649);
    EXPECT_EQ(data.qs.data.rows(), 2);
    EXPECT_EQ(data.wcs.time_axis.data.rows(), 2);
    EXPECT_EQ(data.wcs.time_axis.data.cols(), 6);
    EXPECT_EQ(data.wcs.tone_axis.data.rows(), 15);
    EXPECT_EQ(data.wcs.tone_axis.data.cols(), 649);
    EXPECT_EQ(data.wcs.tone_axis.row_labels.label(0), "f_tone");
    EXPECT_EQ(data.wcs.tone_axis.row_labels.label(3), "f_in");

    EXPECT_DOUBLE_EQ(data.is.data(0, 0), -88454.0);
    EXPECT_DOUBLE_EQ(data.qs.data(0, 0), -217279.0);
    EXPECT_EQ(data.wcs.time_axis.data(0, 0), 1715644457);
    EXPECT_EQ(data.meta.get_typed<int>("sample_slice_start"), 0);
    EXPECT_EQ(data.meta.get_typed<int>("sample_slice_stop"), 2);
    EXPECT_EQ(data.meta.get_typed<int>("sample_slice_step"), 1);
    EXPECT_EQ(data.meta.get_typed<int>("sample_slice_size"), 2);
}

TEST(ToltecTimeStream, RejectsNonPositiveStride)
{
    if (std::getenv("TOLTECA_TEST_DATA_ROOT") == nullptr) {
        GTEST_SKIP()
            << "set TOLTECA_TEST_DATA_ROOT to run the real-file test";
    }
    EXPECT_THROW(
        static_cast<void>(kids::toltec::read_raw_timestream_slice(
            fixture_path(), kids::toltec::SampleSlice{0, 2, 0})),
        kids::toltec::RawTimeStreamIOError);
}

} // namespace

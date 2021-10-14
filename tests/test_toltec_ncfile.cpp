#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <kids/core/kidsdata.h>
#include <kids/toltec/toltec.h>

namespace {

using namespace ::testing;

class KidsData : public Test {
public:
    auto read(const std::string &filepath) {
        using kids::KidsData;
        using kids::KidsDataKind;
        namespace spec = kids::toltec;
        SPDLOG_INFO("use data spec: {}", spec::name);
        // read data
        auto kidsdata = spec::read_data<>(filepath);
        std::visit(
            [](const auto &data) {
                using Data = DECAY(data);
                constexpr auto kind = Data::kind();
                SPDLOG_INFO("{}", data);
                if constexpr (kind & KidsDataKind::Sweep) {
                    SPDLOG_INFO("{}", data.wcs.sweep_axis.data);
                }
            },
            kidsdata);
        return kidsdata;
    }
};

TEST_F(KidsData, ncfile_pre191206_targsweep) {
    const auto filepath =
        "data/pre191206/"
        "toltec0_007565_00_0000_2019_12_03_23_27_15_targsweep.nc";
    auto kidsdata = read(filepath);
}

TEST_F(KidsData, ncfile_lofreq_per_sample_targsweep) {
    const auto filepath =
        "data/lofreq_per_sample/"
        "toltec3_000444_00_0000_2019_12_04_21_42_32_targsweep.nc";
    auto kidsdata = read(filepath);
}

TEST_F(KidsData, ncfile_lofreq_per_sample_timestream) {
    const auto filepath =
        "data/lofreq_per_sample/"
        "toltec3_000445_00_0000_2019_12_04_21_44_02_timestream.nc";
    auto kidsdata = read(filepath);
}

TEST_F(KidsData, ncfile_atten_drive_sense) {
    const auto filepath =
        "data/atten_drive_sense/"
        "toltec12_014659_000_0000_2021_04_23_22_10_43_timestream.nc";
    auto kidsdata = read(filepath);
}

} // namespace

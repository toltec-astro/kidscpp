#pragma once

#include <tula/enum.h>

#include <array>
#include <string_view>

namespace kids {

TULA_ENUM(
    CalibrationModel,
    int,
    D21,
    S21Basic,
    S21WithGain,
    S21WithGainLinTrend,
    S21WithTrans,
    S21WithTransLinTrend);

template <CalibrationModel>
struct CalibrationModelTraits;

template <>
struct CalibrationModelTraits<CalibrationModel::S21WithGainLinTrend> {
    static constexpr std::array<std::string_view, 11> parameter_names{
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
    };
    static constexpr auto parameter_count = parameter_names.size();
};

} // namespace kids

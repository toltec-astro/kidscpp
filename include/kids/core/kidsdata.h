#pragma once

#include "wcs.h"

#include <Eigen/Core>
#include <tula/config/flatconfig.h>
#include <tula/enum.h>
#include <tula/nddata/core.h>
#include <tula/nddata/labelmapper.h>

namespace kids {

TULA_ENUM(
    KidsDataKind,
    int,
    RawTimeStream = 1 << 4,
    SolvedTimeStream = 1 << 5);

namespace internal {

using RMatrixXd =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using RMatrixXi =
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

} // namespace internal

struct ToneAxis : wcs::Axis<ToneAxis, wcs::CoordsKind::Column>,
                  wcs::LabeledData<ToneAxis>,
                  tula::nddata::EigenData<Eigen::MatrixXd> {
    using wcs::LabeledData<ToneAxis>::operator();
    using tula::nddata::EigenData<Eigen::MatrixXd>::operator();

    ToneAxis() = default;
    ToneAxis(
        Eigen::MatrixXd data_,
        tula::nddata::LabelMapper<ToneAxis> row_labels_)
        : tula::nddata::EigenData<Eigen::MatrixXd>{std::move(data_)},
          row_labels{std::move(row_labels_)}
    {
    }

    std::string_view name{"tone"};
    tula::nddata::LabelMapper<ToneAxis> row_labels;
};

struct TimeAxis : wcs::Axis<TimeAxis, wcs::CoordsKind::Row>,
                  wcs::LabeledData<TimeAxis>,
                  tula::nddata::EigenData<internal::RMatrixXi> {
    using wcs::LabeledData<TimeAxis>::operator();
    using tula::nddata::EigenData<internal::RMatrixXi>::operator();

    TimeAxis() = default;
    TimeAxis(
        internal::RMatrixXi data_,
        tula::nddata::LabelMapper<TimeAxis> col_labels_)
        : tula::nddata::EigenData<internal::RMatrixXi>{std::move(data_)},
          col_labels{std::move(col_labels_)}
    {
    }

    std::string_view name{"time"};
    tula::nddata::LabelMapper<TimeAxis> col_labels;
};

struct TimeStreamFrame
    : wcs::Frame2D<TimeStreamFrame, TimeAxis, ToneAxis> {
    TimeAxis time_axis;
    ToneAxis tone_axis;

    [[nodiscard]] auto row_axis() const noexcept -> const TimeAxis &
    {
        return time_axis;
    }

    [[nodiscard]] auto col_axis() const noexcept -> const ToneAxis &
    {
        return tone_axis;
    }
};

template <KidsDataKind kind>
struct KidsDataBase {
    static constexpr auto kind_value = kind;
    using meta_t = tula::config::FlatConfig;

    meta_t meta;
    TimeStreamFrame wcs;
};

template <KidsDataKind>
struct KidsData;

template <>
struct KidsData<KidsDataKind::RawTimeStream>
    : KidsDataBase<KidsDataKind::RawTimeStream> {
    tula::nddata::EigenData<internal::RMatrixXd> is;
    tula::nddata::EigenData<internal::RMatrixXd> qs;
};

template <>
struct KidsData<KidsDataKind::SolvedTimeStream>
    : KidsDataBase<KidsDataKind::SolvedTimeStream> {
    tula::nddata::EigenData<internal::RMatrixXd> rs;
    tula::nddata::EigenData<internal::RMatrixXd> xs;
};

} // namespace kids

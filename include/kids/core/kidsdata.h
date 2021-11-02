#pragma once

#include "wcs.h"
#include <Eigen/Core>
#include <tula/config/flatconfig.h>
#include <tula/enum.h>
#include <tula/formatter/container.h>
#include <tula/formatter/enum.h>
#include <tula/formatter/matrix.h>
#include <tula/formatter/utils.h>
#include <tula/logging.h>
#include <tula/nddata/cacheddata.h>
#include <tula/nddata/core.h>
#include <tula/nddata/labelmapper.h>

namespace kids {

// clang-format off
// NOLINTNEXTLINE
TULA_BITFLAG(KidsDataKind, int,        0xFFFF,
         VnaSweep                = 1 << 0,
         TargetSweep             = 1 << 1,
         Sweep                   = VnaSweep | TargetSweep,
         RawTimeStream           = 1 << 4,
         SolvedTimeStream        = 1 << 5,
         TimeStream              = RawTimeStream | SolvedTimeStream,
         Any                     = Sweep | TimeStream
         );
// clang-format on

/// @brief Kids data class.
template <KidsDataKind kind_ = KidsDataKind::Any>
struct KidsData;
} // namespace kids

namespace std {

// Register KidsData as a variant type
// Below are mandatory to inherit from variant on gcc
template <kids::KidsDataKind kind>
struct variant_size<kids::KidsData<kind>>
    : variant_size<typename kids::KidsData<kind>::variant_t> {};
template <size_t N, auto kind>
struct variant_alternative<N, kids::KidsData<kind>>
    : variant_alternative<N, typename kids::KidsData<kind>::variant_t> {};

#if defined(__GNUC__) && !defined(__clang__)
#if (__GNUC__ >= 9)
// this is need to allow inherit from std::variant on GCC
namespace __detail {
namespace __variant {

template <typename _Ret, typename _Visitor, auto kind, size_t __first>
struct _Multi_array<_Ret (*)(_Visitor, kids::KidsData<kind>), __first>
    : _Multi_array<_Ret (*)(_Visitor, typename kids::KidsData<kind>::variant_t),
                   __first> {
    static constexpr int __do_cookie = 0;
};
template <typename _Maybe_variant_cookie, auto kind>
struct _Extra_visit_slot_needed<_Maybe_variant_cookie, kids::KidsData<kind>>
    : _Extra_visit_slot_needed<_Maybe_variant_cookie,
                               typename kids::KidsData<kind>::variant_t> {};

template <typename _Maybe_variant_cookie, auto kind>
struct _Extra_visit_slot_needed<_Maybe_variant_cookie, kids::KidsData<kind> &>
    : _Extra_visit_slot_needed<_Maybe_variant_cookie,
                               typename kids::KidsData<kind>::variant_t &> {};
} // namespace __variant
} // namespace __detail
#else
#endif
#endif

} // namespace std

namespace kids {

namespace internal {

using RMatrixXd =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using RMatrixXi =
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

template <typename Derived>
struct impl_traits {
    using meta_t = tula::config::FlatConfig;
    constexpr static auto has_tones = requires { Derived::tones; };
    constexpr static auto has_sweeps = requires { Derived::sweeps; };
    constexpr static auto has_fs = requires { Derived::fs; };
    constexpr static auto has_iqs = requires { Derived::iqs; };
    constexpr static auto has_eiqs = requires { Derived::eiqs; };
    constexpr static auto has_cal = requires { Derived::cal; };
    constexpr static auto has_is = requires { Derived::is; };
    constexpr static auto has_qs = requires { Derived::qs; };
    constexpr static auto has_rs = requires { Derived::rs; };
    constexpr static auto has_xs = requires { Derived::xs; };
};

/// @brief CRTP base class for KidsData type
template <typename Derived>
struct KidsDataBase;

// Primitive kind
template <KidsDataKind kind_>
requires(!tula::enum_utils::is_compound_v<kind_>) struct KidsDataBase<
    KidsData<kind_>> {
    static constexpr auto kind() { return kind_; }
    using meta_t = typename internal::impl_traits<KidsData<kind_>>::meta_t;
    meta_t meta;
    static constexpr bool is_variant{false};
};

// Compound kind
template <KidsDataKind kind_>
requires(tula::enum_utils::is_compound_v<kind_>) struct KidsDataBase<
    KidsData<kind_>> {
    static constexpr auto kind() { return kind_; }
    using meta_t = typename internal::impl_traits<KidsData<kind_>>::meta_t;
    static constexpr bool is_variant{true};
};

} // namespace internal

// WCS objects

struct ToneAxis : wcs::Axis<ToneAxis, wcs::CoordsKind::Column>,
                  wcs::LabeledData<ToneAxis>,
                  tula::nddata::EigenData<Eigen::MatrixXd> {
    using wcs::LabeledData<ToneAxis>::operator();
    using tula::nddata::EigenData<Eigen::MatrixXd> ::operator();
    ToneAxis() = default;
    ToneAxis(Eigen::MatrixXd data_,
             tula::nddata::LabelMapper<ToneAxis> row_labels_)
        : tula::nddata::EigenData<Eigen::MatrixXd>{std::move(data_)},
          row_labels{std::move(row_labels_)} {}
    std::string_view name{"tone"};
    tula::nddata::LabelMapper<ToneAxis> row_labels;
};
struct SweepAxis : wcs::Axis<SweepAxis, wcs::CoordsKind::Scalar>,
                   tula::nddata::EigenData<Eigen::VectorXd> {
    SweepAxis() = default;
    SweepAxis(Eigen::VectorXd data_)
        : tula::nddata::EigenData<Eigen::VectorXd>{std::move(data_)} {}
    std::string_view name{"sweep"};
};

/// @brief The wcs frame for sweep data.
struct SweepFrame;

struct SweepFrame : wcs::Frame2D<SweepFrame, SweepAxis, ToneAxis> {
    SweepAxis sweep_axis;
    ToneAxis tone_axis;

    // Frame2D impl
    auto row_axis() const noexcept -> const SweepAxis & { return sweep_axis; }
    auto col_axis() const noexcept -> const ToneAxis & { return tone_axis; }
    auto fs() const -> const auto & { return m_fs(*this); }

private:
    constexpr static auto fs_evaluator = [](const auto &frame) {
        auto tfs = frame.tone_axis("f_tone");
        auto sfs = frame.sweep_axis();
        return sfs.replicate(1, tfs.size()) + tfs.replicate(sfs.size(), 1);
    };
    struct tula::nddata::CachedData < Eigen::MatrixXd, fs_evaluator> m_fs;
};

struct TimeAxis : wcs::Axis<TimeAxis, wcs::CoordsKind::Row>,
                  wcs::LabeledData<TimeAxis>,
                  tula::nddata::EigenData<internal::RMatrixXi> {
    TimeAxis() = default;
    TimeAxis(internal::RMatrixXi data_,
             tula::nddata::LabelMapper<TimeAxis> col_labels_)
        : tula::nddata::EigenData<internal::RMatrixXi>{std::move(data_)},
          col_labels{std::move(col_labels_)} {}
    std::string_view name{"tone"};
    tula::nddata::LabelMapper<TimeAxis> col_labels;
};

struct TimeStreamFrame : wcs::Frame2D<TimeStreamFrame, TimeAxis, ToneAxis> {
    TimeAxis time_axis;
    ToneAxis tone_axis;
    // Frame2D impl
    [[nodiscard]] auto row_axis() const noexcept -> const TimeAxis & {
        return time_axis;
    }
    [[nodiscard]] auto col_axis() const noexcept -> const ToneAxis & {
        return tone_axis;
    }
};

// Data objects

/// @brief Base class for frequency sweep data
template <typename Derived>
struct Sweep : internal::KidsDataBase<Derived>,
               tula::nddata::NDData<Sweep<Derived>> {
    SweepFrame wcs;

    auto derived() const noexcept -> const auto & {
        return static_cast<const Derived &>(*this);
    }
    auto tones() const noexcept { return derived().wcs.tone_axis("f_tone"); }
    auto sweeps() const noexcept -> const auto & {
        return derived().wcs.sweep_axis();
    }
    auto fs() const noexcept -> const auto & { return derived().wcs.fs(); }
    // data objects
    tula::nddata::EigenData<Eigen::MatrixXcd> iqs;
    tula::nddata::EigenData<Eigen::MatrixXcd> eiqs; // uncertainty
};

template <>
struct KidsData<KidsDataKind::VnaSweep>
    : Sweep<KidsData<KidsDataKind::VnaSweep>> {};
template <>
struct KidsData<KidsDataKind::TargetSweep>
    : Sweep<KidsData<KidsDataKind::TargetSweep>> {};

/// @brief Base class for time stream data
template <typename Derived>
struct TimeStream : internal::KidsDataBase<Derived>,
                    tula::nddata::NDData<TimeStream<Derived>> {
    TimeStreamFrame wcs;
};

template <>
struct KidsData<KidsDataKind::RawTimeStream>
    : TimeStream<KidsData<KidsDataKind::RawTimeStream>> {
    using Base = TimeStream<KidsData<KidsDataKind::RawTimeStream>>;
    tula::nddata::EigenData<internal::RMatrixXd> is;
    tula::nddata::EigenData<internal::RMatrixXd> qs;
};

template <>
struct KidsData<KidsDataKind::SolvedTimeStream>
    : TimeStream<KidsData<KidsDataKind::SolvedTimeStream>> {
    using Base = TimeStream<KidsData<KidsDataKind::SolvedTimeStream>>;
    tula::nddata::EigenData<internal::RMatrixXd> rs;
    tula::nddata::EigenData<internal::RMatrixXd> xs;
};

/// @brief Kids data class of runtime variant kind.
template <KidsDataKind kind_>
requires(tula::enum_utils::is_compound_v<kind_>) struct KidsData<kind_>
    : tula::enum_utils::enum_to_variant_t<kind_, KidsData>,
      internal::KidsDataBase<KidsData<kind_>> {
    // using Base = tula::enum_utils::enum_to_variant_t<kind_, KidsData>;
    using variant_t = tula::enum_utils::enum_to_variant_t<kind_, KidsData>;
    // note that this class does not contain any data on itself.
    // construct from primitive type
    // template <KidsDataKind kind1,
    // REQUIRES_V(!enum_utils::is_compound_v<kind1>)> KidsData(KidsData<kind1>
    // other) : Base(std::move(other)) {}
    auto variant() const -> const variant_t & { return *this; }
};

} // namespace kids

namespace fmt {

template <>
struct formatter<kids::ToneAxis>
    : formatter<kids::wcs::FrameBase<kids::ToneAxis>> {
    using Base = formatter<kids::wcs::FrameBase<kids::ToneAxis>>;
    template <typename FormatContext>
    auto format(const kids::ToneAxis &data, FormatContext &ctx) {
        auto it = Base::format(data, ctx);
        return it = format_to(it, " labels={}", data.row_labels());
    }
};
template <>
struct formatter<kids::SweepAxis>
    : formatter<kids::wcs::FrameBase<kids::SweepAxis>> {
    using Base = formatter<kids::wcs::FrameBase<kids::SweepAxis>>;
    template <typename FormatContext>
    auto format(const kids::SweepAxis &data, FormatContext &ctx) {
        return Base::format(data, ctx);
    }
};
template <>
struct formatter<kids::TimeAxis>
    : formatter<kids::wcs::FrameBase<kids::TimeAxis>> {
    using Base = formatter<kids::wcs::FrameBase<kids::TimeAxis>>;
    template <typename FormatContext>
    auto format(const kids::TimeAxis &data, FormatContext &ctx) {
        return Base::format(data, ctx);
    }
};

template <kids::KidsDataKind kind_>
struct formatter<kids::KidsData<kind_>>
    : tula::fmt_utils::charspec_formatter_base<'l', 's'> {
    // s: the short form
    // l: the long form
    using Data = kids::KidsData<kind_>;

    template <typename FormatContext>
    auto format(const Data &data, FormatContext &ctx) {
        using data_traits = kids::internal::impl_traits<Data>;
        auto it = ctx.out();
        constexpr auto kind = Data::kind();
        if constexpr (tula::enum_utils::is_compound_v<kind>) {
            // format compound kind as variant
            return format_to(it, "({}) {:0}", kind, data.variant());
        } else {
            /// format simple kind type
            auto spec = spec_handler();
            // meta
            switch (spec) {
            case 's': {
                it = format_to(it, "({})", kind);
                if constexpr (data_traits::has_iqs) {
                    return format_member(it, "iqs", data.iqs());
                } else if constexpr (data_traits::has_xs) {
                    return format_member(it, "xs", data.xs());
                }
            }
            case 'l': {
                it =
                    format_to(it, "kind={} meta={}", kind, data.meta.pformat());
                bool sep = false;
                if constexpr (data_traits::has_tones) {
                    it = format_member(it, "tones", data.tones(), &sep);
                }
                if constexpr (data_traits::has_sweeps) {
                    it = format_member(it, "sweeps", data.sweeps(), &sep);
                }
                if constexpr (data_traits::has_fs) {
                    it = format_member(it, "fs", data.fs(), &sep);
                }
                if constexpr (data_traits::has_iqs) {
                    it = format_member(it, "iqs", data.iqs(), &sep);
                }
                if constexpr (data_traits::has_eiqs) {
                    it = format_member(it, "eiqs", data.eiqs(), &sep);
                }
                if constexpr (data_traits::has_is) {
                    it = format_member(it, "is", data.is(), &sep);
                }
                if constexpr (data_traits::has_qs) {
                    it = format_member(it, "qs", data.qs(), &sep);
                }
                if constexpr (data_traits::has_xs) {
                    it = format_member(it, "xs", data.xs(), &sep);
                }
                if constexpr (data_traits::has_rs) {
                    it = format_member(it, "rs", data.rs(), &sep);
                }
                return it;
            }
            }
            return it;
        }
    }
    template <typename T, typename FormatContextOut>
    auto format_member(FormatContextOut &it, std::string_view name, const T &m,
                       bool *sep = nullptr) {
        auto spec = spec_handler();
        switch (spec) {
        case 's': {
            if ((m.rows() == 1) || (m.cols() == 1)) {
                return format_to(it, "[{}]", m.size());
            }
            return format_to(it, "[{}, {}]", m.rows(), m.cols());
        }
        case 'l': {
            it = format_to(it, fmt::runtime("{}{}={:r0}"), sep ? " " : "", name,
                           m);
            // as soon as this it called, sep is set
            *sep = true;
            return it;
        }
        }
        return it;
    }
};
} // namespace fmt

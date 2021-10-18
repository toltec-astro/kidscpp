#pragma once

#include <Eigen/Core>
#include <tula/formatter/matrix.h>
#include <tula/formatter/utils.h>
#include <tula/meta.h>
#include <tula/nddata/core.h>
#include <tula/nddata/eigen.h>

// This unfortunately requires c++20
// template <typename DataType, auto evaluator>
// struct CachedAttr : NDDataAttr<CachedAttr<DataType, evaluator>> {
//     mutable DataType data;
//     mutable bool initialized{false};
//     template <typename T> const auto &get(const T &derived) const {
//         if (!initialized) {
//             data = evaluator(derived);
//             initialized = true;
//         }
//         return this->data;
//     }
// };
//

namespace kids::wcs {

template <typename Derived>
struct LabeledData : tula::nddata::NDData<Derived> {

    using Base = tula::nddata::NDData<Derived>;
    using Base::derived;
    using index_t = typename tula::nddata::type_traits<Derived>::index_t;
    using label_t = typename tula::nddata::type_traits<Derived>::label_t;

    template <typename U = Derived>
    requires requires { U::row_labels; }
    auto operator()(const label_t &name) const {
        return derived().data.row(derived().row_labels.index(name));
    }
    template <typename U = Derived>
    requires requires { U::row_labels; }
    auto operator()(const label_t &name) {
        return derived().data.row(derived().row_labels.index(name));
    }

    template <typename U = Derived>
    requires requires { U::col_labels; }
    auto operator()(const label_t &name) const {
        return derived().data.col(derived().col_labels.index(name));
    }
    template <typename U = Derived>
    requires requires { U::col_labels; }
    auto operator()(const label_t &name) {
        return derived().data.col(derived().col_labels.index(name));
    }
};

namespace internal {

struct frame_impl_traits : tula::nddata::type_traits<frame_impl_traits> {
    using Base = tula::nddata::type_traits<frame_impl_traits>;
    using index_t = Base::index_t;
    using physical_type_t = Base::physical_type_t;
    using unit_t = Base::unit_t;
    template <typename data_t>
    static constexpr auto is_valid_data_type =
        std::is_base_of_v<tula::nddata::NDData<data_t>, data_t>;
};

} // namespace internal

/// @brief Frame base interface.
template <typename Derived>
struct FrameBase {
    using index_t = internal::frame_impl_traits::index_t;
    [[nodiscard]] constexpr auto array_n_dim() const noexcept -> index_t {
        return Derived::frame_dimension;
    }
    auto derived() const noexcept -> const auto & {
        return static_cast<const Derived &>(*this);
    }
    auto derived() noexcept -> auto & { return static_cast<Derived &>(*this); }
};

enum class CoordsKind { Scalar, Row, Column };

/// @brief Axis base interface that represent one single dimension.
template <typename Derived, CoordsKind coords_kind>
struct Axis : FrameBase<Derived> {
    using Base = FrameBase<Derived>;
    using index_t = typename Base::index_t;
    using Base::derived;

    // FrameBase impl
    static constexpr index_t frame_dimension = 1;

    // WCS impl
    auto array_index_to_world_values(index_t i) const -> const auto & {
        if constexpr (coords_kind == CoordsKind::Scalar) {
            return derived().data.coeff(i);
        } else if constexpr (coords_kind == CoordsKind::Row) {
            return derived().data.row(i);
        } else if constexpr (coords_kind == CoordsKind::Column) {
            return derived().data.col(i);
        }
    }
    void resize(index_t n) {
        if constexpr (coords_kind == CoordsKind::Scalar) {
            return derived().data.resize(n);
        } else if constexpr (coords_kind == CoordsKind::Row) {
            return derived().data.resize(n, world_n_dim());
        } else if constexpr (coords_kind == CoordsKind::Column) {
            return derived().data.resize(world_n_dim(), n);
        }
    }
    auto array_shape() const -> index_t {
        if constexpr (coords_kind == CoordsKind::Scalar) {
            return derived().data.size();
        } else if constexpr (coords_kind == CoordsKind::Row) {
            return derived().data.rows();
        } else if constexpr (coords_kind == CoordsKind::Column) {
            return derived().data.cols();
        }
    }
    auto world_n_dim() const -> index_t {
        constexpr index_t one = 1;
        if constexpr (coords_kind == CoordsKind::Scalar) {
            return one;
        } else if constexpr (coords_kind == CoordsKind::Row) {
            return derived().data.cols();
        } else if constexpr (coords_kind == CoordsKind::Column) {
            return derived().data.rows();
        }
    }
};

/// @brief Two dimensional frame that consists of two axes.
template <typename Derived, typename RowAxis, typename ColAxis>
struct Frame2D : FrameBase<Derived> {
    using Base = FrameBase<Derived>;
    using index_t = typename Base::index_t;
    using Base::derived;

    // FrameBase impl
    static constexpr index_t dimension = 2;

    auto array_index_to_world_values(index_t i, index_t j) const {
        return std::make_tuple(
            derived().row_axis().array_index_to_world_values(i),
            derived().col_axis().array_index_to_world_values(j));
    }
    auto array_shape() -> std::pair<index_t, index_t> {
        return {derived().row_axis().array_shape(),
                derived().col_axis().array_shape()};
    }
    auto world_n_dim() -> index_t {
        return derived().row_axis().world_n_dim() +
               derived().col_axis().world_n_dim();
    }
};

} // namespace kids::wcs

namespace fmt {

template <typename Derived, typename Char>
struct formatter<kids::wcs::FrameBase<Derived>, Char>
    : tula::fmt_utils::charspec_formatter_base<'l', 's'> {
    // s: the short form
    // l: the long form

    template <typename FormatContext>
    auto format(const kids::wcs::FrameBase<Derived> &data_,
                FormatContext &ctx) {
        const auto &data = data_.derived();
        auto it = ctx.out();
        /// format simple kind type
        auto spec = spec_handler();
        std::string name = "unnamed";
        if constexpr (requires { Derived::name; }) {
            name = data.name;
        }
        // meta
        switch (spec) {
        case 's': {
            it = format_to(it, "{}({}d{})[{}]", name, data.array_n_dim(),
                           data.world_n_dim(), data.array_shape());
        }
        case 'l': {
            it = format_to(it, "wcs{}d name={} shape=[{}] coords_dim={}",
                           data.array_n_dim(), name, data.array_shape(),
                           data.world_n_dim());
        }
        }
        return it;
    }
};

} // namespace fmt

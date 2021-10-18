#pragma once

#include <unsupported/Eigen/FFT>
#include <tula/algorithm/ei_stats.h>
#include <tula/container.h>
#include <tula/formatter/matrix.h>
#include <tula/formatter/utils.h>

namespace internal {

using Eigen::Index;
using Eigen::VectorXd;

//                         npts   nf,    df
using FFTStat = std::tuple<Index, Index, double>;

FFTStat fftstat(Index npts_, double fsmp) {
    // make an even number of data points by rounding-down
    Index npts = npts_;
    if (npts % 2 == 1) {
        --npts;
    }
    // prepare containers in frequency domain
    Index nfs = npts / 2 + 1; // number of one sided freq bins
    double df = fsmp / npts;
    SPDLOG_TRACE("using npts={} (out of n={})", npts, npts_);
    SPDLOG_TRACE("fsmp={} df={} nfs={}", fsmp, df, nfs);
    return {npts, nfs, df};
}

// fft frequency bins
VectorXd fftfs(Index npts, Index nfs, double df) {
    return df * VectorXd::LinSpaced(nfs, 0, npts / 2);
}

enum Window { NoWindow = 0, Hann = 1, HannScaled = 2 };

inline auto hann(Index npts) {
    // hann window
    // N numbers starting from 0 to (include) 2pi/N * (N-1)
    // 0.5 - 0.5 cos([0, 2pi/N,  2pi/N * 2, ... 2pi/N * (N - 1)])
    return (0.5 -
            0.5 *
                Eigen::ArrayXd::LinSpaced(npts, 0,
                                          2 * static_cast<double>(EIGEN_PI) /
                                              npts * (npts - 1))
                    .cos())
        .matrix();
}

template <typename... Args> inline auto hann_scaled(Args &&... args) {
    // NENBW = 1.5 * df therefore we devide by 1.5 here to get the
    // equivalent scale as if no window is used
    return hann(FWD(args)...) / 1.5;
}

// one sided psd [0, fsmp/2];
template <Window win = Hann, typename DerivedA, typename DerivedB>
FFTStat psd(const Eigen::DenseBase<DerivedA> &data_,
            Eigen::DenseBase<DerivedB> const &psddata_, double fsmp) {
    // decltype(auto) scan = _scan.derived();
    // decltype(auto) forward the return type of derived() so it declares a
    // refernce as expected if scan has storage, this is equivalent to:
    decltype(auto) data = data_.derived();
    decltype(auto) psddata =
        const_cast<Eigen::DenseBase<DerivedB> &>(psddata_).derived();

    auto stat = fftstat(data.size(), fsmp);
    auto [npts, nfs, df] = stat;

    // prepare fft
    Eigen::FFT<double> fft;
    fft.SetFlag(Eigen::FFT<double>::HalfSpectrum);
    fft.SetFlag(Eigen::FFT<double>::Unscaled);
    Eigen::VectorXcd buf;

    double scale = 1. / fsmp; // unit of Hz^{-1}
    auto dispatch_window = [](auto npts) {
        if constexpr (win == Hann) {
            return hann(npts);
        } else if constexpr (win == HannScaled) {
            return hann_scaled(npts);
        }
    };
    // make copy of the data
    auto _data = data.head(npts).array().eval();
    auto trend_const = _data.mean();
    // subtrace mean
    _data -= trend_const;

    // branch according to whether applying window
    if constexpr (win == NoWindow) {
        SPDLOG_TRACE("fft with no window");
        fft.fwd(buf, _data.matrix());
        scale /= npts;
    } else if constexpr ((win == Hann) || (win == HannScaled)) {
        // we need the eval() here per the requirement of fft.fwd()
        auto _win = dispatch_window(npts);
        SPDLOG_TRACE("fft with win={} {}", win, _win);
        _data *= _win.array();
        fft.fwd(buf, _data.matrix());
        scale /= _win.squaredNorm();
    } else {
        throw std::runtime_error("unknwon window type");
        // static_assert(meta::always_false<Window>{}, "UNKNOWN WINDOW TYPE");
    } // note: at this point the freqdata is not normalized to NEBW yet

    // calcualte psd
    // normalize to frequency resolution
    SPDLOG_TRACE("fft data{} result buf{} scale{}", data, buf, scale);
    psddata = buf.cwiseAbs2() * scale;
    // account for the negative frequencies by an extra factor of 2. note: first
    // and last are 0 and nquist freq, so they only appear once
    psddata.segment(1, nfs - 2) *= 2.;
    SPDLOG_TRACE("psds{}", psddata);
    return stat;
}

template <Window win = Hann, typename A, typename B, typename DerivedC>
FFTStat psd(A &&data, B &&psddata, Eigen::DenseBase<DerivedC> const &fs_,
            double fsmp) {

    auto stat = psd<win>(FWD(data), FWD(psddata), fsmp);
    decltype(auto) fs = const_cast<Eigen::DenseBase<DerivedC> &>(fs_).derived();
    auto [npts, nfs, df] = stat;
    fs = fftfs(npts, nfs, df);
    SPDLOG_TRACE("psd fs{}", fs);
    return stat;
}

/*
 * This mimics scipy.welch.
 */
template <Window win = Hann, typename A, typename DerivedB>
FFTStat welch(A &&data, Eigen::DenseBase<DerivedB> const &psddata_, double fsmp,
              int nperseg, double overlap_fraction = 0.5) {

    std::vector<std::pair<Index, Index>> segments;
    auto noverlap = int(overlap_fraction * nperseg);
    for (Index i = 0; i + nperseg <= data.size(); i += nperseg - noverlap) {
        segments.emplace_back(i, i + nperseg);
    }
    SPDLOG_TRACE("segments{}", segments);
    auto stat = fftstat(nperseg, fsmp);
    auto [npts, nfs, df] = stat;

    auto nsegs = segments.size();
    auto calc_psds = [&, npts = npts](auto &out) {
        Eigen::VectorXd buf;
        out.setConstant(0.);
        for (std::size_t i = 0; i < nsegs; ++i) {
            psd<win>(data.segment(segments[i].first, npts), buf, fsmp);
            out += buf;
        }
        out /= double(nsegs);
    };

    decltype(auto) psddata =
        const_cast<Eigen::DenseBase<DerivedB> &>(psddata_).derived();
    if (psddata.size() > 0) {
        // modify psddata in place
        calc_psds(psddata);
    } else {
        Eigen::VectorXd _psddata(nfs);
        calc_psds(_psddata);
        psddata = std::move(_psddata);
    }
    SPDLOG_TRACE("psd data{}", psddata);
    return stat;
}

} // namespace internal

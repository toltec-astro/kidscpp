#include <kids/sweep/model.h>
#include <utils/logging.h>
#include <Eigen/Dense>


int main(int argc, char* argv[]) {
    logging::init<>(true);
    {
        logging::scoped_timeit timer("calc x demo");

        using Model = kids::Model<kids::SweepModel::S21WithGainLinTrend>;
        // load model params
        // the model params has to be an Eigen::Matrix of shape
        // [Model::NP, n_detectors]
        const auto n_detectors = 2;
        Eigen::MatrixXd mdl(Model::NP, n_detectors);
        // let's create some example model params
        // the data is from the first two row of tone file
        // toltec3_012758_000_0000_2020_09_24_17_04_06_tune.txt
        mdl <<
            677912089.84375,     678221171.875,
            29357.46891504806,   23666.72598789425,
            1,                   1,
            677911368.0522399,   678220061.5168798,
            0,                   0,
            15.33526832825879,   107.9779834980814,
            66.49090362276078,   -56.63154753085788,
            -1.009214019834231,  -0.4645663804972067,
            1.222792635775451,   -0.1812642778482604,
            -4029207.455497678,  -5894015.522510329,
            -6036638.110360871,  3422826.505198508;
        //  0  1   2   3  4  5      6       7       8       9          10
        // fp Qr, Qc, fr, A, normI, normQ, slopeI, slopeQ, interceptI, interceptQ
        SPDLOG_INFO("mdl{}", mdl);

        // now we are ready to evaluate the model
        // note that this should be properly encapsulated in the
        // next version of kidscpp
        {
            // this is the ratio of accum len between timestream and the tune.
            // for example, if tiemstream is downsampled by 4x, i.e., 122Hz,
            // iqnorm should be 4 because we always take tune at 488Hz.
            const double iq_norm = 1;

            // the function to compute r,x from I, Q
            auto iq2rx = [&mdl, &iq_norm] (const auto& is, const auto& qs, const auto& fs) {
                auto m = mdl.array();
                static const auto norm2 = (m.row(5) * m.row(5) + m.row(6) * m.row(6)).eval();
                auto mi = (m.row(9) + (fs.array() - m.row(0)) * m.row(7)).eval();
                auto mq = (m.row(10) + (fs.array() - m.row(0)) * m.row(8)).eval();
                // detrend s21
                auto is1 = ((is.array() / iq_norm) - mi);
                auto qs1 = ((qs.array() / iq_norm) - mq);
                // deroated iqs
                auto is0 = (is1.rowwise() * m.row(5) + qs1.rowwise() * m.row(6))
                          .rowwise() / norm2;
                auto qs0 = (qs1.rowwise() * m.row(5) - is1.rowwise() * m.row(6))
                          .rowwise() / norm2;
                // compute x
                auto denom =
                    ((is0 * is0 + qs0 * qs0).rowwise() * m.row(2) * 2.).eval();
                auto xs = (is0.rowwise() * m.row(4) - qs0).cwiseQuotient(denom).eval();
                auto rs = (qs0.rowwise() * m.row(4) + is0).cwiseQuotient(denom).eval();
                return std::tuple{rs, xs};
            };

            auto rx2iq = [&mdl, &iq_norm] (const auto& rs, const auto& xs, const auto& fs) {
                // the function to get I, Q from r, x
                // make a copy
                Eigen::MatrixXd m = mdl;
                //  0  1   2   3  4  5      6       7       8       9          10
                // fp Qr, Qc, fr, A, normI, normQ, slopeI, slopeQ, interceptI, interceptQ
                // update fr and Qr for each sample and evalute at fp
                // the evaluation returns complex number so we use MatrixXcd
                // to store is and qs.
                auto n_samples = rs.rows();
                auto n_detectors = rs.cols();
                Eigen::MatrixXcd iqs(n_samples, n_detectors);
                for (int i = 0; i < n_samples; ++i)  {
                    // Qr = 1 / (2r)
                    m.row(1) = 1. / (2. * rs.row(i).array());
                    // fr = fp / (x + 1)
                    m.row(3) = fs.array() / (xs.array() + 1);
                    for (int j = 0; j < n_detectors; ++j) {
                    tula::alg::ceresfit::eval<Model>(
                        fs.segment(j, 1), m.col(j), iqs.block(i, j, 1, 1));
                    }
                }
                // now copy out the is and qs
                auto is = iqs.real().eval();
                auto qs = iqs.imag().eval();
                return std::tuple(is, qs);
            };

            // the is and qs has to be MatrixXd of shape [n_samples, n_detectors]
            // the fs has to be *row* vector of shape [n_detectors] which is the
            // probe tone frequency in Hz
            Eigen::MatrixXd is(1, n_detectors);
            Eigen::MatrixXd qs(1, n_detectors);
            Eigen::RowVectorXd fs(n_detectors);
            is << 12345., 22222;
            qs << 54321., 55555;
            fs << 5e8, 6e8; // Hz
            auto [rs, xs] = iq2rx(is, qs, fs);
            SPDLOG_INFO("I{}, Q{} ->\nr{}, x{}", is, qs, rs, xs);
            // round trip
            auto [is2, qs2] = rx2iq(rs, xs, fs);
            SPDLOG_INFO("r{}, x{} ->\nI_roundtrip{}, Q_roundtrip{}", rs, xs, is2, qs2);
       }
    }
}

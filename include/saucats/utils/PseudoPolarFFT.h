#ifndef SAUCATS_UTILS_PSEUDO_POLAR_FFT_H
#define SAUCATS_UTILS_PSEUDO_POLAR_FFT_H

#include <saucats/utils/FFT.h>
#include <iostream>



namespace saucats {
	class PseudoPolarFFT {
	public:
		static Eigen::VectorXcd
		frft(const Eigen::VectorXcd& A, double alpha) {
			const std::complex<double> j(0., 1.);
    
			Eigen::Index N = A.size();
		
			Eigen::VectorXcd K(2 * N);
			for(Eigen::Index i = 0; i < N; ++i) {
				K[i] = i;
				K[2 * N - 1 - i] = -(i + 1); 
			}
			K = ((-M_PI * alpha * j) * K.array().square()).exp();

			Eigen::VectorXcd X(2 * N);
			X.head(N) = A;
			X.tail(N) = Eigen::VectorXcd::Zero(N);
			X.array() *= K.array();

			Eigen::VectorXcd XX = FFT::fft_1d(X);
			XX.array() *= FFT::fft_1d(K.conjugate()).array();
			Eigen::VectorXcd Y = FFT::ifft_1d(XX);
			Y.array() *= K.array();

			return Y.head(N);
		}



		static Eigen::VectorXcd
		centered_frft(const Eigen::VectorXcd& A, double alpha) {
			const std::complex<double> j(0., 1.);

			Eigen::Index N = A.size();

			Eigen::VectorXcd K(N);
			for(Eigen::Index i = 0; i < N; ++i)
				K[i] = i;
			K *= M_PI * alpha * N * j;
			
			Eigen::VectorXcd Ap = A.array() * K.array().exp();
			return frft(Ap, alpha);
		}

		

		static void
		fftshift_rows(Eigen::MatrixXcd& X) {
			for(Eigen::Index i = 0; i < X.rows() / 2; ++i) {
				Eigen::Index ip = i + (X.rows() / 2);
				for(Eigen::Index j = 0; j < X.cols(); ++j)
					std::swap(X(i, j), X(ip, j));
			}
		}

		static void
		fftshift_cols(Eigen::MatrixXcd& X) {
			for(Eigen::Index i = 0; i < X.rows(); ++i) {
				for(Eigen::Index j = 0; j < X.cols() / 2; ++j)
					std::swap(X(i, j), X(i, j + X.cols() / 2));
			}
		}


		static Eigen::MatrixXcd
		pseudo_polar_fft_2d(const Eigen::MatrixXd& A) {
			const std::complex<double> j(0., 1.);

			Eigen::Index N = A.rows();
			double N2 = N * N;

			Eigen::MatrixXcd ret(2 * N, 2 * N);
			ret.array() = 0.;

			// Construction of quadrant 1 and 3
      {
				Eigen::MatrixXcd Ap(2 * N, N);
				Ap.topRows(N) = A;
				Ap.bottomRows(N).array() = 0.;

				Eigen::MatrixXcd U(2 * N, N);
				for(Eigen::Index i = 0; i < N; ++i)
					U.col(i) = FFT::fft_1d(Ap.col(i));
				fftshift_rows(U);

				for(Eigen::Index i = -N; i < N; ++i)
					ret.row(i+N).head(N) = centered_frft(U.row(i+N), i / N2).reverse();
			}

			// Construction of quadrant 2 and 4
			{
				Eigen::MatrixXcd Ap(N, 2 * N);
				Ap.leftCols(N) = A;
				Ap.rightCols(N).array() = 0.;

				Eigen::MatrixXcd U(N, 2 * N);
				for(Eigen::Index i = 0; i < N; ++i)
					U.row(i) = FFT::fft_1d(Ap.row(i));
				fftshift_cols(U);
				U.transposeInPlace();

				for(Eigen::Index i = -N; i < N; ++i) {
					Eigen::VectorXcd K = (((N / 2 - 1) * (i / (N2)) * 2 * M_PI * j) * Eigen::VectorXcd::LinSpaced(N, 0., N - 1)).array().exp();
					K.array() *= U.row(i+N).array();
					ret.row(i+N).tail(N) = frft(K, i / N2);			
				}
			}

			// Job done
			return ret;
		}
	}; // class PseudoPolarFFT
} // namespace saucats



#endif // SAUCATS_UTILS_PSEUDO_POLAR_FFT_H

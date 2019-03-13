#ifndef SAUCATS_UTILS_PSEUDO_POLAR_FFT_H
#define SAUCATS_UTILS_PSEUDO_POLAR_FFT_H

#include <saucats/utils/FFT.h>
#include <iostream>



namespace saucats {
	class PseudoPolarFFT2d {
	public:
		PseudoPolarFFT2d(Eigen::Index size) :
			m_fft_1d(2 * size),
			m_N(size) { } 



		Eigen::MatrixXcd
		operator () (const Eigen::MatrixXd& A) {
			const std::complex<double> j(0., 1.);

			double N2 = m_N * m_N;

			Eigen::MatrixXcd ret(2 * m_N, 2 * m_N);
			ret.array() = 0.;

			// Construction of quadrant 1 and 3
      {
				Eigen::MatrixXcd Ap(2 * m_N, m_N);
				Ap.topRows(m_N) = A;
				Ap.bottomRows(m_N).array() = 0.;

				Eigen::MatrixXcd U(2 * m_N, m_N);
				for(Eigen::Index i = 0; i < m_N; ++i)
					U.col(i) = m_fft_1d.fft(Ap.col(i));
				
				fftshift_rows(U);

				for(Eigen::Index i = -m_N; i < m_N; ++i)
					ret.row(i+m_N).head(m_N) = centered_frft(U.row(i+m_N), i / N2).reverse();
			}

			// Construction of quadrant 2 and 4
			{
				Eigen::MatrixXcd Ap(m_N, 2 * m_N);
				Ap.leftCols(m_N) = A;
				Ap.rightCols(m_N).array() = 0.;

				Eigen::MatrixXcd U(m_N, 2 * m_N);
				for(Eigen::Index i = 0; i < m_N; ++i)
					U.row(i) = m_fft_1d.fft(Ap.row(i));

				fftshift_cols(U);
				U.transposeInPlace();

				for(Eigen::Index i = -m_N; i < m_N; ++i) {
					Eigen::VectorXcd K = (((m_N / 2 - 1) * (i / (N2)) * 2 * M_PI * j) * Eigen::VectorXcd::LinSpaced(m_N, 0., m_N - 1)).array().exp();
					K.array() *= U.row(i+m_N).array();
					ret.row(i+m_N).tail(m_N) = frft(K, i / N2);			
				}
			}

			// Job done
			return ret;
		}

	private:
		Eigen::VectorXcd
		frft(const Eigen::VectorXcd& A, double alpha) {
			const std::complex<double> j(0., 1.);
		
			Eigen::VectorXcd K(2 * m_N);
			for(Eigen::Index i = 0; i < m_N; ++i) {
				K[i] = i;
				K[2 * m_N - 1 - i] = -(i + 1); 
			}
			K = ((-M_PI * alpha * j) * K.array().square()).exp();

			Eigen::VectorXcd X(2 * m_N);
			X.head(m_N) = A;
			X.tail(m_N) = Eigen::VectorXcd::Zero(m_N);
			X.array() *= K.array();

			Eigen::VectorXcd XX = m_fft_1d.fft(X);
			XX.array() *= m_fft_1d.fft(K.conjugate()).array();
			Eigen::VectorXcd Y = m_fft_1d.ifft(XX);
			Y.array() *= K.array();

			return Y.head(m_N);
		}



		Eigen::VectorXcd
		centered_frft(const Eigen::VectorXcd& A, double alpha) {
			const std::complex<double> j(0., 1.);

			Eigen::VectorXcd K(m_N);
			for(Eigen::Index i = 0; i < m_N; ++i)
				K[i] = i;
			K *= M_PI * alpha * m_N * j;
			
			Eigen::VectorXcd Ap = A.array() * K.array().exp();
			return frft(Ap, alpha);
		}



		FFT1d m_fft_1d;
		Eigen::Index m_N;		
	}; // class PseudoPolarFFT
} // namespace saucats



#endif // SAUCATS_UTILS_PSEUDO_POLAR_FFT_H

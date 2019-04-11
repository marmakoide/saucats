#ifndef SAUCATS_UTILS_PSEUDO_POLAR_FFT_H
#define SAUCATS_UTILS_PSEUDO_POLAR_FFT_H

#include <saucats/utils/FFT.h>
#include <iostream>



namespace saucats {
	void
	fftshift_rows(Eigen::MatrixXcd& X) {
		for(Eigen::Index i = 0; i < X.rows() / 2; ++i) {
			Eigen::Index ip = i + (X.rows() / 2);
			for(Eigen::Index j = 0; j < X.cols(); ++j)
				std::swap(X(i, j), X(ip, j));
		}
	}

	void
	fftshift_cols(Eigen::MatrixXcd& X) {
		for(Eigen::Index i = 0; i < X.rows(); ++i)
			for(Eigen::Index j = 0; j < X.cols() / 2; ++j)
				std::swap(X(i, j), X(i, j + X.cols() / 2));
	}



	class FractionalFFT1d {
	public:
		FractionalFFT1d(Eigen::Index size) :
			m_fft_1d(2 * size),
			m_N(size) { } 

		Eigen::VectorXcd
		frft(const Eigen::VectorXcd& A, double alpha) {
			const std::complex<double> imag_j(0., 1.);
		
			Eigen::VectorXcd K(2 * m_N);
			K.head(m_N) = Eigen::VectorXcd::LinSpaced(m_N, 0, m_N - 1);
			K.tail(m_N) = Eigen::VectorXcd::LinSpaced(m_N, -m_N, -1);
			K = ((-M_PI * alpha * imag_j) * K.array().square()).exp();

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
			const std::complex<double> imag_j(0., 1.);

			Eigen::VectorXcd K = (M_PI * alpha * m_N * imag_j) * Eigen::VectorXcd::LinSpaced(m_N, 0, m_N - 1);
			Eigen::VectorXcd Ap = A.array() * K.array().exp();

			return frft(Ap, alpha);
		}
		
	private:
		FFT1d m_fft_1d;
		Eigen::Index m_N;	
	}; // class FractionalFFT1d



	/*
	 * 2d pseudo-polar FFT, see "Pseudo-polar based estimation of large 
	 * translations rotations and scalings in images" by Y. Keller, A. Averbuch, 
	 * and M. Isreali for reference
	 */

	class PseudoPolarFFT2d {
	public:
		PseudoPolarFFT2d(Eigen::Index size) :
			m_fft_1d(2 * size),
			m_frac_fft_1d(size),
			m_N(size) { } 



		Eigen::MatrixXcd
		operator () (const Eigen::MatrixXd& A) {
			const std::complex<double> j(0., 1.);

			double N2 = m_N * m_N;

			Eigen::MatrixXcd U(2 * m_N, m_N);
			Eigen::MatrixXcd ret = Eigen::MatrixXcd::Zero(2 * m_N, 2 * m_N);

			// Construction of quadrant 1 and 3
			U.topRows(m_N) = A;
			U.bottomRows(m_N).array() = 0.;

			for(Eigen::Index i = 0; i < m_N; ++i)
				U.col(i) = m_fft_1d.fft(U.col(i));
				
			fftshift_rows(U);

			for(Eigen::Index i = -m_N; i < m_N; ++i)
				ret.row(i + m_N).head(m_N) = m_frac_fft_1d.centered_frft(U.row(i + m_N), i / N2).reverse();

			// Construction of quadrant 2 and 4
			U.topRows(m_N) = A.transpose();
			U.bottomRows(m_N).array() = 0.;	

			for(Eigen::Index i = 0; i < m_N; ++i)
				U.col(i) = m_fft_1d.fft(U.col(i));

			fftshift_rows(U);

			for(Eigen::Index i = -m_N; i < m_N; ++i) {
				Eigen::VectorXcd K = (((m_N / 2 - 1) * (i / (N2)) * 2 * M_PI * j) * Eigen::VectorXcd::LinSpaced(m_N, 0., m_N - 1)).array().exp();
				K.array() *= U.row(i + m_N).array();
				ret.row(i + m_N).tail(m_N) = m_frac_fft_1d.frft(K, i / N2);
			}

			// Job done
			return ret;
		}

	private:
		FFT1d m_fft_1d;
		FractionalFFT1d m_frac_fft_1d;
		Eigen::Index m_N;		
	}; // class PseudoPolarFFT2d



	class PseudoPolarCylindricalFFT3d {
	public:
		typedef RealFFT3d::real_data_type real_data_type;
		typedef RealFFT3d::complex_data_type complex_data_type;



		PseudoPolarCylindricalFFT3d(Eigen::Index size) :
			m_xy_fft_1d(2 * size),
			m_z_fft_1d(size),
			m_frac_fft_1d(size),
			m_N(size) { }

		complex_data_type
		operator () (const real_data_type& A) {
			const std::complex<double> imag_j(0., 1.); 

			double N2 = m_N * m_N;

			complex_data_type ret(2 * m_N, 2 * m_N, m_N);
			ret.setZero();			

			// 1d FFT along Z axis
			complex_data_type U(m_N, m_N, m_N);

			{
				Eigen::VectorXcd A_slice(m_N);
				Eigen::VectorXcd U_slice(m_N);

				for(Eigen::Index i = 0; i < m_N; ++i)
					for(Eigen::Index j = 0; j < m_N; ++j) {
						for(Eigen::Index k = 0; k < m_N; ++k)
							A_slice(k) = A(i, j, k);
	
						U_slice = m_z_fft_1d.fft(A_slice);

						for(Eigen::Index k = 0; k < m_N; ++k)
							U(i, j, k) = U_slice(k);
					}
			}

			// 2d pseudo-polar FFT on each XY slice
			Eigen::MatrixXcd V(2 * m_N, m_N);

			for(Eigen::Index i = 0; i < m_N; ++i) {
				// Construction of quadrant 1 and 3
				{
					V.bottomRows(m_N).array() = 0.;
					for(Eigen::Index j = 0; j < m_N; ++j)
						for(Eigen::Index k = 0; k < m_N; ++k)
							V(k, j) = U(j, k, i); 

					for(Eigen::Index j = 0; j < m_N; ++j)
						V.col(j) = m_xy_fft_1d.fft(V.col(j));
				
					fftshift_rows(V);

					for(Eigen::Index j = -m_N; j < m_N; ++j) {
						Eigen::VectorXcd W = m_frac_fft_1d.centered_frft(V.row(j + m_N), j / N2).reverse();
						for(Eigen::Index k = 0; k < m_N; ++k)
							ret(j + m_N, k, i) = W(k);
					}
				}

				// Construction of quadrant 2 and 4
				{
					V.bottomRows(m_N).array() = 0.;
					for(Eigen::Index j = 0; j < m_N; ++j)
						for(Eigen::Index k = 0; k < m_N; ++k)
							V(k, j) = U(k, j, i); 

					for(Eigen::Index j = 0; j < m_N; ++j)
						V.col(j) = m_xy_fft_1d.fft(V.col(j));

					fftshift_rows(V);

					for(Eigen::Index j = -m_N; j < m_N; ++j) {
						Eigen::VectorXcd K = (((m_N / 2 - 1) * (j / (N2)) * 2 * M_PI * imag_j) * Eigen::VectorXcd::LinSpaced(m_N, 0., m_N - 1)).array().exp();
						K.array() *= V.row(j + m_N).array();
						Eigen::VectorXcd W = m_frac_fft_1d.frft(K, j / N2);
						for(Eigen::Index k = 0; k < m_N; ++k)
							ret(j + m_N, k + m_N, i) = W(k);			
					}
				}
			}

			// Job done
			return ret;
		}

	private:
		FFT1d m_z_fft_1d;
		FFT1d m_xy_fft_1d;
		FractionalFFT1d m_frac_fft_1d;
		Eigen::Index m_N;		
	}; //  class PseudoPolarCylindricalFFT3d
} // namespace saucats



#endif // SAUCATS_UTILS_PSEUDO_POLAR_FFT_H

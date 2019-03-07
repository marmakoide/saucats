#ifndef SAUCATS_UTILS_FFT_H
#define SAUCATS_UTILS_FFT_H

#include <fftw3.h>
#include <Eigen/Dense>



namespace saucats {
	/*
	 * FFT routines, interfacing with the fftw3 library
	 */

	class FFT {
	public:
		static Eigen::VectorXcd
		real_fft_1d(const Eigen::VectorXd& A) {
			// Allocate ressources
			double* A_data = (double*)fftw_malloc(sizeof(double) * A.size());
			fftw_complex* U_data = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (1 + A.size() / 2));

			// Mapping of raw data for Eigen library
			Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1> > A_data_map(A_data, A.size());
			Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> > U_data_map((std::complex<double>*)U_data, 1 + A.size() / 2);
	
			// Setup FFT computations
			fftw_plan fftw_plan = fftw_plan_dft_r2c_1d(A.size(), A_data, U_data, FFTW_ESTIMATE);

			// Compute FFT
			A_data_map = A;
			fftw_execute(fftw_plan);
			Eigen::VectorXcd U = U_data_map;

			// Free ressources	
			fftw_destroy_plan(fftw_plan);
			fftw_free(U_data);
			fftw_free(A_data);

			// Job done
			return U;
		}

		static Eigen::VectorXd
		real_ifft_1d(const Eigen::VectorXcd& U) {
			// Allocate ressources
			double* A_data = (double*)fftw_malloc(sizeof(double) * (2 * (U.size() - 1)));
			fftw_complex* U_data = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * U.size());

			// Mapping of raw data for Eigen library
			Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1> > A_data_map(A_data, 2 * (U.size() - 1));
			Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> > U_data_map((std::complex<double>*)U_data, U.size());
	
			// Setup FFT computations
			fftw_plan fftw_plan = fftw_plan_dft_c2r_1d(2 * (U.size() - 1), U_data, A_data, FFTW_ESTIMATE);

			// Compute inverse FFT
			U_data_map = U;
			fftw_execute(fftw_plan);
			Eigen::VectorXd A = A_data_map;
			A /= A.size();

			// Free ressources	
			fftw_destroy_plan(fftw_plan);
			fftw_free(U_data);
			fftw_free(A_data);

			// Job done
			return A;
		}

		static Eigen::MatrixXcd
		real_fft_2d(const Eigen::MatrixXd& A) {
			// Allocate ressources
			double* A_data = (double*)fftw_malloc(sizeof(double) * A.rows() * A.cols());
			fftw_complex* U_data = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * A.rows() * (1 + A.cols() / 2));

			// Mapping of raw data for Eigen library
			Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > A_data_map(A_data, A.rows(), A.cols());
			Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > U_data_map((std::complex<double>*)U_data, A.rows(), 1 + A.cols() / 2);
	
			// Setup FFT computations
			fftw_plan fftw_plan = fftw_plan_dft_r2c_2d(A.cols(), A.rows(), A_data, U_data, FFTW_ESTIMATE);

			// Compute FFT
			A_data_map = A;
			fftw_execute(fftw_plan);
			Eigen::MatrixXcd U = U_data_map;

			// Free ressources	
			fftw_destroy_plan(fftw_plan);
			fftw_free(U_data);
			fftw_free(A_data);

			// Job done
			return U;
		}



		static Eigen::MatrixXd
		real_ifft_2d(const Eigen::MatrixXcd& U) {
			// Allocate ressources
			double* A_data = (double*)fftw_malloc(sizeof(double) * U.rows() * (2 * (U.cols() - 1)));
			fftw_complex* U_data = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * U.rows() * U.cols());

			// Mapping of raw data for Eigen library
			Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > A_data_map(A_data, U.rows(), 2 * (U.cols() - 1));
			Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > U_data_map((std::complex<double>*)U_data, U.rows(), U.cols());
	
			// Setup FFT computations
			fftw_plan fftw_plan = fftw_plan_dft_c2r_2d(2 * (U.cols() - 1), U.rows(), U_data, A_data, FFTW_ESTIMATE);

			// Compute inverse FFT
			U_data_map = U;
			fftw_execute(fftw_plan);
			Eigen::MatrixXd A = A_data_map;
			A /= A.rows() * A.cols();

			// Free ressources	
			fftw_destroy_plan(fftw_plan);
			fftw_free(U_data);
			fftw_free(A_data);

			// Job done
			return A;
		}
	}; // class FFT
} // namespace saucats



#endif // SAUCATS_UTILS_FFT_H

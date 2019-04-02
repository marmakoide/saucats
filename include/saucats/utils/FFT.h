#ifndef SAUCATS_UTILS_FFT_H
#define SAUCATS_UTILS_FFT_H

#include <fftw3.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>



namespace saucats {
	/*
	 * FFT routines, interfacing with the fftw3 library
	 */

	class FFT1d {
	public:
		typedef Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> > mapping_type;



		FFT1d(Eigen::Index size) :
			m_A_data(allocate_array(size)),
			m_U_data(allocate_array(size)),
			m_A_map((std::complex<double>*)m_A_data, size),
			m_U_map((std::complex<double>*)m_U_data, size)
		{
			m_fwd_plan = fftw_plan_dft_1d(size, m_A_data, m_U_data, FFTW_FORWARD, FFTW_ESTIMATE);		
			m_bwd_plan = fftw_plan_dft_1d(size, m_U_data, m_A_data, FFTW_BACKWARD, FFTW_ESTIMATE);
		}

		~FFT1d() {
			fftw_destroy_plan(m_fwd_plan);
			fftw_destroy_plan(m_bwd_plan);
			fftw_free(m_U_data);
			fftw_free(m_A_data);
		}

		mapping_type&
		fft(const Eigen::VectorXcd& A) {
			m_A_map = A;
			fftw_execute(m_fwd_plan);
			return m_U_map;
		}

		mapping_type&
		ifft(const Eigen::VectorXcd& U) {
			m_U_map = U;
			fftw_execute(m_bwd_plan);
			return m_A_map;
		}

	private:
		static fftw_complex*
		allocate_array(Eigen::Index size) {
			return (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size);
		}



		fftw_complex* m_A_data;
		fftw_complex* m_U_data;
		fftw_plan m_fwd_plan;
		fftw_plan m_bwd_plan;
		mapping_type m_A_map;
		mapping_type m_U_map;
	}; // class FFT1d



	class FFT2d {
	public:
		typedef Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > mapping_type;



		FFT2d(Eigen::Index size) :
			m_A_data(allocate_array(size)),
			m_U_data(allocate_array(size)),
			m_A_map((std::complex<double>*)m_A_data, size, size),
			m_U_map((std::complex<double>*)m_U_data, size, size)
		{
			m_fwd_plan = fftw_plan_dft_2d(size, size, m_A_data, m_U_data, FFTW_FORWARD, FFTW_ESTIMATE);		
			m_bwd_plan = fftw_plan_dft_2d(size, size, m_U_data, m_A_data, FFTW_BACKWARD, FFTW_ESTIMATE);
		}

		~FFT2d() {
			fftw_destroy_plan(m_fwd_plan);
			fftw_destroy_plan(m_bwd_plan);
			fftw_free(m_U_data);
			fftw_free(m_A_data);
		}

		mapping_type&
		fft(const Eigen::MatrixXcd& A) {
			m_A_map = A;
			fftw_execute(m_fwd_plan);
			return m_U_map;
		}

		mapping_type&
		ifft(const Eigen::MatrixXcd& U) {
			m_U_map = U;
			fftw_execute(m_bwd_plan);
			return m_A_map;
		}

	private:
		static fftw_complex*
		allocate_array(Eigen::Index size) {
			return (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size * size);
		}



		fftw_complex* m_A_data;
		fftw_complex* m_U_data;
		fftw_plan m_fwd_plan;
		fftw_plan m_bwd_plan;
		mapping_type m_A_map;
		mapping_type m_U_map;
	}; // class FFT2d




	class RealFFT1d {
	public:
		typedef Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1> > real_mapping_type;
		typedef Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> > complex_mapping_type;



		RealFFT1d(Eigen::Index size) :
			m_A_data(allocate_real_array(size)),
			m_U_data(allocate_complex_array(size)),
			m_A_map(m_A_data, size),
			m_U_map((std::complex<double>*)m_U_data, 1 + size / 2)
		{
			m_fwd_plan = fftw_plan_dft_r2c_1d(size, m_A_data, m_U_data, FFTW_ESTIMATE);
			m_bwd_plan = fftw_plan_dft_c2r_1d(size, m_U_data, m_A_data, FFTW_ESTIMATE);
		}

		~RealFFT1d() {
			fftw_destroy_plan(m_fwd_plan);
			fftw_destroy_plan(m_bwd_plan);
			fftw_free(m_U_data);
			fftw_free(m_A_data);
		}

		complex_mapping_type&
		fft(const Eigen::VectorXd& A) {
			m_A_map = A;
			fftw_execute(m_fwd_plan);
			return m_U_map;
		}

		real_mapping_type&
		ifft(const Eigen::VectorXcd& U) {
			m_U_map = U;
			fftw_execute(m_bwd_plan);
			return m_A_map;
		}

	private:
		static double*
		allocate_real_array(Eigen::Index size) {
			return (double*)fftw_malloc(sizeof(double) * size);
		}

		static fftw_complex*
		allocate_complex_array(Eigen::Index size) {
			return (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (1 + size / 2));
		}



		double* m_A_data;
		fftw_complex* m_U_data;
		fftw_plan m_fwd_plan;
		fftw_plan m_bwd_plan;
		real_mapping_type m_A_map;
		complex_mapping_type m_U_map;
	}; // class RealFFT1d



	class RealFFT2d {
	public:
		typedef Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > real_mapping_type;
		typedef Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > complex_mapping_type;



		RealFFT2d(Eigen::Index size) :
			m_A_data(allocate_real_array(size)),
			m_U_data(allocate_complex_array(size)),
			m_A_map(m_A_data, size, size),
			m_U_map((std::complex<double>*)m_U_data, size, 1 + size / 2)
		{
			m_fwd_plan = fftw_plan_dft_r2c_2d(size, size, m_A_data, m_U_data, FFTW_ESTIMATE);
			m_bwd_plan = fftw_plan_dft_c2r_2d(size, size, m_U_data, m_A_data, FFTW_ESTIMATE);
		}

		~RealFFT2d() {
			fftw_destroy_plan(m_fwd_plan);
			fftw_destroy_plan(m_bwd_plan);
			fftw_free(m_U_data);
			fftw_free(m_A_data);
		}

		complex_mapping_type&
		fft(const Eigen::MatrixXd& A) {
			m_A_map = A;
			fftw_execute(m_fwd_plan);
			return m_U_map;
		}

		real_mapping_type&
		ifft(const Eigen::MatrixXcd& U) {
			m_U_map = U;
			fftw_execute(m_bwd_plan);
			return m_A_map;
		}

	private:
		static double*
		allocate_real_array(Eigen::Index size) {
			return (double*)fftw_malloc(sizeof(double) * size * size);
		}

		static fftw_complex*
		allocate_complex_array(Eigen::Index size) {
			return (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size * (1 + size / 2));
		}



		double* m_A_data;
		fftw_complex* m_U_data;
		fftw_plan m_fwd_plan;
		fftw_plan m_bwd_plan;
		real_mapping_type m_A_map;
		complex_mapping_type m_U_map;
	}; // class RealFFT2d



	class RealFFT3d {
	public:
		typedef Eigen::Tensor<double, 3, Eigen::RowMajor> real_data_type;
		typedef Eigen::Tensor<std::complex<double>, 3, Eigen::RowMajor> complex_data_type;
		typedef Eigen::TensorMap<real_data_type> real_mapping_type;
		typedef Eigen::TensorMap<complex_data_type> complex_mapping_type;


		RealFFT3d(Eigen::Index size) :
			m_A_data(allocate_real_array(size)),
			m_U_data(allocate_complex_array(size)),
			m_A_map(m_A_data, size, size, size),
			m_U_map((std::complex<double>*)m_U_data, size, size, 1 + size / 2)
		{
			m_fwd_plan = fftw_plan_dft_r2c_3d(size, size, size, m_A_data, m_U_data, FFTW_ESTIMATE);
			m_bwd_plan = fftw_plan_dft_c2r_3d(size, size, size, m_U_data, m_A_data, FFTW_ESTIMATE);
		}

		~RealFFT3d() {
			fftw_destroy_plan(m_fwd_plan);
			fftw_destroy_plan(m_bwd_plan);
		}

		complex_mapping_type&
		fft(const real_data_type& A) {
			m_A_map = A;
			fftw_execute(m_fwd_plan);
			return m_U_map;
		}

		real_mapping_type&
		ifft(const complex_data_type& U) {
			m_U_map = U;
			fftw_execute(m_bwd_plan);
			return m_A_map;
		}

	private:
		static double*
		allocate_real_array(Eigen::Index size) {
			return (double*)fftw_malloc(sizeof(double) * size * size * size);
		}

		static fftw_complex*
		allocate_complex_array(Eigen::Index size) {
			return (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size * size * (1 + size / 2));
		}



		double* m_A_data;
		fftw_complex* m_U_data;
		fftw_plan m_fwd_plan;
		fftw_plan m_bwd_plan;
		real_mapping_type m_A_map;
		complex_mapping_type m_U_map;
	}; // class RealFFT3d
} // namespace saucats



#endif // SAUCATS_UTILS_FFT_H

#ifndef SAUCATS_SDF_UTILS_REGISTRATION_H
#define SAUCATS_SDF_UTILS_REGISTRATION_H

#include <saucats/geometry/Sampling.h>
#include <saucats/utils/TensorMath.h>
#include <saucats/utils/PseudoPolarFFT.h>



namespace saucats {
	/*
	 * Find the optimal pose so that a shape B overlap/dock into a shaoe A
	 */

	namespace registration {
		/*
		 * From a signed distance field, define a 2d scalar field called 
		 * "registration field".
		 */

		template <class func_type>
		class RegistrationFieldFunctor {
		public:
			typedef typename func_type::scalar_type scalar_type;
			typedef typename func_type::vector_type vector_type;
			typedef typename func_type::sphere_type sphere_type;
			typedef Eigen::Matrix<scalar_type, Eigen::Dynamic, vector_type::RowsAtCompileTime> sample_array_type;



			RegistrationFieldFunctor(const func_type& func,
	    		                     const sample_array_type& sample_array) :
				m_func(func),
				m_sample_array(sample_array) { }

			scalar_type 
			operator () (const vector_type& P) const {
				const scalar_type epsilon = .125;

				scalar_type P_dist = m_func.dist(P);
				scalar_type P_radius = std::fabs(P_dist);
	
				scalar_type sum = 0.;
				for(Eigen::Index i = 0; i < m_sample_array.rows(); ++i) {	
					sphere_type search_sphere(P, P_radius + epsilon);
					LineT<vector_type> search_line(m_sample_array.row(i), P);
					auto it = get_intersection_iterator(m_func, search_line, search_sphere);		
					if (!it.completed())
						sum += 1;
				}

			return std::copysign(sum / m_sample_array.rows(), P_dist);
		}

		private:
			const func_type& m_func;
			const sample_array_type& m_sample_array;
		}; // class RegistrationField2dFunctor



		template <class func_type>
		RegistrationFieldFunctor<func_type>
		get_registration_field_functor(const func_type& func,
		                               const Eigen::Matrix<typename func_type::scalar_type, Eigen::Dynamic, func_type::vector_type::RowsAtCompileTime>& sample_array) {
			return RegistrationFieldFunctor<func_type>(func, sample_array);
		}



		/*
		 * Augment a scalar field with a mask
		 */


		template <class func_type, class shape_type>
		class MaskedFunctor {
		public:
			typedef typename func_type::scalar_type scalar_type;
			typedef typename func_type::vector_type vector_type;



			MaskedFunctor(const func_type& func,
			              const shape_type& shape) :
				m_func(func),
				m_shape(shape) { }

			inline scalar_type
			operator() (const vector_type& P) const {
				if (m_shape.contains(P))
					return m_func(P);
				return 0.;
			} 
	
			private:
				const func_type& m_func;
				const shape_type& m_shape;
		}; // class MaskedFunctor

 

		template <class func_type, class shape_type>
		MaskedFunctor<func_type, shape_type>
		get_masked_functor(const func_type& func,
		                   const shape_type& shape) {
			return MaskedFunctor<func_type, shape_type>(func, shape);
		}



		/*
		 * Find the translation so that B cross-correlation with A is maximal
		 * A and B should have the same dimensions.
		 */

		Eigen::Vector2d
		solve_translation_2d(const Eigen::MatrixXd& A,
		                     const Eigen::MatrixXd& B,
		                     double& score) {
			RealFFT2d real_fft_2d(A.rows());

			// Compute C as cross-correlation of A and B
			Eigen::MatrixXcd U = real_fft_2d.fft(A);
			Eigen::MatrixXcd V = real_fft_2d.fft(B).conjugate();
			U.array() *= V.array();
			Eigen::MatrixXd C = real_fft_2d.ifft(U);	

			// Search for the maximum element of C
			Eigen::Index row_id = 0;
			Eigen::Index col_id = 0;
			score = C.maxCoeff(&row_id, &col_id); 

			// Adjust coordinates
			if (row_id > (A.rows() / 2))
				row_id -= A.rows();

			if (col_id > (A.cols() / 2))
				col_id -= A.cols();

			// Job done
			return Eigen::Vector2d(col_id, row_id);
		}

		/*
		 * Find the translation so that B cross-correlation with A is maximal
		 * A and B should have the same dimensions.
		 */

		Eigen::Vector3d
		solve_translation_3d(const RealFFT3d::real_data_type& A,
		                     const RealFFT3d::real_data_type& B,
		                     double& score) {
			// Compute C as cross-correlation of A and B
			RealFFT3d real_fft_3d(A.dimension(0));

			RealFFT3d::complex_data_type U = real_fft_3d.fft(A);
			RealFFT3d::complex_data_type V = real_fft_3d.fft(B).conjugate();
			U *= V;
			RealFFT3d::real_data_type C = real_fft_3d.ifft(U);	
	
			// Search for the maximum element of C
			Eigen::Array3i arg_max;
			get_3d_real_tensor_arg_max(C, arg_max);

			// Adjust coordinates
			for(int i = 0; i < 3; i++)
				if (arg_max.coeff(i) > (A.dimension(i) / 2))
					arg_max.coeffRef(i) -= A.dimension(i);

			// Job done
			return arg_max.matrix().cast<double>();
		}



		template <class func_type>
		void
		fill_density_matrix(const func_type& func,
                        Eigen::MatrixXd& M,
                        const typename func_type::sphere_type& domain,
		                    const Eigen::Matrix<typename func_type::scalar_type, Eigen::Dynamic, func_type::vector_type::RowsAtCompileTime>& sample_array,
                        double resolution) {
			// Compute the box domain
			auto box_domain = BoxT<typename func_type::vector_type>(Eigen::Vector2d(resolution * M.cols(), resolution * M.rows()), domain.center());

			// Fills the matrix
			sample_function_2d(box_domain, M, registration::get_masked_functor(registration::get_registration_field_functor(func, sample_array), domain));
		}



		template <class func_type>
		void
		fill_density_tensor(const func_type& func,
                        RealFFT3d::real_data_type& M,
                        const typename func_type::sphere_type& domain,
		                    const Eigen::Matrix<typename func_type::scalar_type, Eigen::Dynamic, func_type::vector_type::RowsAtCompileTime>& sample_array,
                        double resolution) {
			// Compute the box domain
			auto box_domain = BoxT<typename func_type::vector_type>(Eigen::Vector3d(resolution * M.dimension(0), resolution * M.dimension(1), resolution * M.dimension(2)), domain.center());

			// Fills the matrix
			sample_function_3d(box_domain, M, registration::get_masked_functor(registration::get_registration_field_functor(func, sample_array), domain));
		}
	} // namespace registration



	template <class left_func_type, class right_func_type>
	Eigen::Transform<double, 2, Eigen::Affine>
	registrate_2d(const left_func_type& left_func,
  	            const right_func_type& right_func,
  	            double resolution,
	              bool complement) {
		auto unit_circle_samples = get_circle_points<double>(256);

		// Compute left and right sdf bounding sphere
		auto left_bounding_sphere = left_func.get_bounding_sphere();
		auto right_bounding_sphere = right_func.get_bounding_sphere();

		// Compute left and right domains
		Sphere2d left_domain;
		if (complement)
			left_domain = Sphere2d(left_bounding_sphere.center(), left_bounding_sphere.radius() + 2 * right_bounding_sphere.radius());
		else
			left_domain = Sphere2d(left_bounding_sphere.center(), std::max(left_bounding_sphere.radius(), right_bounding_sphere.radius()));

		Sphere2d right_domain = Sphere2d(right_bounding_sphere.center(), right_bounding_sphere.radius());

		// Compute matrix size
		double domain_radius = std::max(left_domain.radius(), right_domain.radius());
		Eigen::Index matrix_size = get_smallest_greater_power_of_2(Eigen::Index(std::ceil(2 * domain_radius / resolution)));

		// Discretize the left and right function
		Eigen::MatrixXd left_M(matrix_size, matrix_size);
		Eigen::MatrixXd right_M(matrix_size, matrix_size);

		registration::fill_density_matrix(left_func, left_M, left_domain, unit_circle_samples, resolution);
		registration::fill_density_matrix(right_func, right_M, right_domain, unit_circle_samples, resolution);

		// Switch sign of left function
		if (complement)
			left_M = -left_M;

		// Solve the rotation component
		PseudoPolarFFT2d pseudo_polar_fft(matrix_size);
		FFT2d fft_2d(2 * matrix_size);

		// Compute U, pseudo-polar FFT for left function (required only once)
		Eigen::MatrixXcd U = pseudo_polar_fft(left_M).cwiseAbs();
		Eigen::MatrixXcd A = fft_2d.fft(U);
		A.array() /= A.cwiseAbs().array();
	
		std::vector<Eigen::Index> angle_list;
		angle_list.push_back(0);

		auto rotated_right_func = get_sdf_rotation2(right_func, 0.);
		bool converged = false;
		int max_epoch = 100;
		for(int epoch = 0; (epoch < max_epoch) and (!converged); ++epoch) {
			// Compute C, the pseudo FFT of pseudo-polar FFT for right function
			Eigen::MatrixXcd V = pseudo_polar_fft(right_M).cwiseAbs();
			Eigen::MatrixXcd B = fft_2d.fft(V);

			// Compute correlation of A and B
			Eigen::MatrixXcd C = B.conjugate();
			C.array() *= A.array();
			C.array() /= B.cwiseAbs().array();
			Eigen::MatrixXcd W = fft_2d.ifft(C);
			Eigen::MatrixXd W_abs = W.cwiseAbs();

			// Look for the maximum correlation coefficient
			Eigen::Index row_id = 0;
			Eigen::Index col_id = 0;
			W_abs.maxCoeff(&row_id, &col_id);	

			// Compute the angle increment
			Eigen::Index theta_delta = -col_id;
			if (col_id >= matrix_size)
				theta_delta += 2 * matrix_size;
			
			// Compute the angle
			Eigen::Index angle = (angle_list.back() + theta_delta) % (4 * matrix_size);

			// Check that we didn't find this angle before
			if (std::find(angle_list.begin(), angle_list.end(), angle) != angle_list.end())
				converged = true;

			angle_list.push_back(angle);

			// Resample the right function
			rotated_right_func = get_sdf_rotation2(right_func, -(.5 * M_PI / matrix_size) * angle);
			right_bounding_sphere = rotated_right_func.get_bounding_sphere();
			right_domain = Sphere2d(right_bounding_sphere.center(), right_bounding_sphere.radius());
			registration::fill_density_matrix(rotated_right_func, right_M, right_domain, unit_circle_samples, resolution);
 		}
		
		// Convert the angle to radian 
		double angle = (.5 * M_PI / matrix_size) * angle_list.back();
		
		// Solve the translation component
		double hi_score = 0.;
		Eigen::Vector2d hi_T = registration::solve_translation_2d(left_M, right_M, hi_score);

		matrix_rotate_180(right_M);

		double lo_score = 0.;
		Eigen::Vector2d lo_T = registration::solve_translation_2d(left_M, right_M, lo_score);

		Eigen::Vector2d T;
		if (hi_score > lo_score)
			T = hi_T;
		else {
			T = lo_T;
			angle += M_PI;
			rotated_right_func = get_sdf_rotation2(right_func, -angle);
			right_bounding_sphere = rotated_right_func.get_bounding_sphere();
		}

		T *= -resolution;
		T += right_bounding_sphere.center() - left_bounding_sphere.center(); 
	
		// Job done
		return
			Eigen::Translation<double, 2>(-T) *
			Eigen::Rotation2D<double>(-angle);
	}



	template <class left_func_type, class right_func_type>
	Eigen::Transform<double, 3, Eigen::Affine>
	registrate_3d(const left_func_type& left_func,
  	            const right_func_type& right_func,
  	            double resolution,
	              bool complement) {
		auto unit_sphere_samples = get_spiral_sphere<double>(256);

		// Compute left and right sdf bounding sphere
		auto left_bounding_sphere = left_func.get_bounding_sphere();
		auto right_bounding_sphere = right_func.get_bounding_sphere();

		// Compute left and right domains
		Sphere3d left_domain;
		if (complement)
			left_domain = Sphere3d(left_bounding_sphere.center(), left_bounding_sphere.radius() + 2 * right_bounding_sphere.radius());
		else
			left_domain = Sphere3d(left_bounding_sphere.center(), std::max(left_bounding_sphere.radius(), right_bounding_sphere.radius()));

		Sphere3d right_domain = Sphere3d(right_bounding_sphere.center(), right_bounding_sphere.radius());

		// Compute matric size
		double domain_radius = std::sqrt(std::max(left_domain.squared_radius(), right_domain.squared_radius()));
		Eigen::Index matrix_size = get_smallest_greater_power_of_2(Eigen::Index(std::ceil(2 * domain_radius / resolution)));

		// Discretize the left and right function
		RealFFT3d::real_data_type left_M(matrix_size, matrix_size, matrix_size);
		RealFFT3d::real_data_type right_M(matrix_size, matrix_size, matrix_size);

		registration::fill_density_tensor(left_func, left_M, left_domain, unit_sphere_samples, resolution);
		registration::fill_density_tensor(right_func, right_M, right_domain, unit_sphere_samples, resolution);

		// Switch sign of right function
		if (complement)
			right_M = -right_M;
	
		// Solve the Z rotation angle
		PseudoPolarCylindricalFFT3d pseudo_polar_cylindrical_fft(matrix_size);
		FFT3d fft_3d(2 * matrix_size);

		RealFFT3d::complex_data_type U = pseudo_polar_cylindrical_fft(left_M).abs();
		RealFFT3d::complex_data_type A = fft_3d.fft(U);
		RealFFT3d::real_data_type A_abs = A.abs();

		std::vector<Eigen::Index> angle_list;
		angle_list.push_back(0);

		Eigen::Transform<double, 3, Eigen::Affine> transform;
		auto rotated_right_func = get_sdf_affine_transform(right_func, transform);

		int max_epoch = 10;
		for(int epoch = 0; epoch < max_epoch; ++epoch) {
			RealFFT3d::complex_data_type V = pseudo_polar_cylindrical_fft(right_M).abs();
			RealFFT3d::complex_data_type B = fft_3d.fft(V);
			RealFFT3d::real_data_type B_abs = B.abs();

			RealFFT3d::complex_data_type C = B.conjugate();
			C *= A;
			C /= A_abs;
			C /= B_abs;
			RealFFT3d::complex_data_type W = fft_3d.ifft(C);
			RealFFT3d::real_data_type W_abs = W.abs();

			Eigen::Array3i arg_max;
			get_3d_real_tensor_arg_max(W_abs, arg_max);

			Eigen::Index theta;
			if (arg_max.coeff(0) < matrix_size)
				theta = -arg_max.coeff(0);
			else
				theta = 2 * matrix_size - arg_max.coeff(0);
		
			Eigen::Index angle = (angle_list.back() + theta) % (4 * matrix_size);

			if (std::find(angle_list.begin(), angle_list.end(), angle) != angle_list.end()) {
				angle_list.push_back(angle);
				break;
			}

			angle_list.push_back(angle);

			transform = Eigen::AngleAxisd(Eigen::AngleAxisd(-(.5 * M_PI / matrix_size) * angle, Eigen::Vector3d::UnitZ()));
			rotated_right_func = get_sdf_affine_transform(right_func, transform);
			std::cout << "angle = " << (.5 * M_PI / matrix_size) * angle << std::endl;

			right_bounding_sphere = rotated_right_func.get_bounding_sphere();
			right_domain = Sphere3d(right_bounding_sphere.center(), right_bounding_sphere.radius());
			registration::fill_density_tensor(rotated_right_func, right_M, right_domain, unit_sphere_samples, resolution);

			if (complement)
				right_M = -right_M;
 		}

		double angle = (.5 * M_PI / matrix_size) * angle_list.back();
		
		// Job done
		Eigen::Transform<double, 3, Eigen::Affine> ret;
		ret.setIdentity();
		ret = Eigen::AngleAxisd(Eigen::AngleAxisd(-angle, Eigen::Vector3d::UnitZ()));
		return ret;
	}

	/*
	template <class left_func_type, class right_func_type>
	Eigen::Transform<double, 3, Eigen::Affine>
	registrate_3d(const left_func_type& left_func,
  	            const right_func_type& right_func,
  	            double resolution,
	              bool complement) {
		auto unit_sphere_samples = get_spiral_sphere<double>(256);

		// Compute left and right sdf bounding sphere
		auto left_bounding_sphere = left_func.get_bounding_sphere();
		auto right_bounding_sphere = right_func.get_bounding_sphere();

		// Compute left and right domains
		auto left_domain = Sphere3d(left_bounding_sphere.center(), left_bounding_sphere.radius() + 2 * right_bounding_sphere.radius());
		auto right_domain = Sphere3d(right_bounding_sphere.center(), right_bounding_sphere.radius());

		// Compute matric size
		double domain_radius = std::sqrt(std::max(left_domain.squared_radius(), right_domain.squared_radius()));
		Eigen::Index matrix_size = get_smallest_greater_power_of_2(Eigen::Index(std::ceil(2 * domain_radius / resolution)));

		// Discretize the left and right function
		RealFFT3d::real_data_type left_M(matrix_size, matrix_size, matrix_size);
		RealFFT3d::real_data_type right_M(matrix_size, matrix_size, matrix_size);

		registration::fill_density_tensor(left_func, left_M, left_domain, unit_sphere_samples, resolution);
		registration::fill_density_tensor(right_func, right_M, right_domain, unit_sphere_samples, resolution);

		// Switch sign of right function
		if (complement)
			right_M = -right_M;
	
		// Solve the translation component
		double score = 0.;
		Eigen::Vector3d T = registration::solve_translation_3d(left_M, right_M, score);

		T *= -2 * domain_radius / matrix_size;
		T += right_bounding_sphere.center() - left_bounding_sphere.center(); 
		
		// Job done
		Eigen::Transform<double, 3, Eigen::Affine> ret;
		ret = Eigen::Translation<double, 3>(-T);
		return ret;
	}
	*/
} // namespace saucats



#endif // SAUCATS_SDF_UTILS_REGISTRATION_H

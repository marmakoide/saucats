#ifndef SAUCATS_SDF_UTILS_REGISTRATION_H
#define SAUCATS_SDF_UTILS_REGISTRATION_H

#include <saucats/utils/PseudoPolarFFT.h>



namespace saucats {
	/*
	 * Find the optimal pose so that a shape B overlap/dock into a shaoe A
	 */

	namespace registration {
		/*
		 * From a signed distance field, define a 2dscalar field called "partial
		 * occlusion field".
		 */

		template <class func_type>
		class PartialOcclusionFunctor {
		public:
			typedef typename func_type::scalar_type scalar_type;
			typedef typename func_type::vector_type vector_type;
			typedef typename func_type::sphere_type sphere_type;



			PartialOcclusionFunctor(const func_type& func,
	    		                    Eigen::Index sample_count) :
				m_U(sample_count, 2), 
				m_func(func) {
				for(Eigen::Index i = 0; i < sample_count; ++i) {
					scalar_type theta = (2 * M_PI / sample_count) * i;
					m_U.row(i) << std::cos(theta), std::sin(theta); 
				}
			}

			scalar_type 
			operator () (const vector_type& P) const {
				const scalar_type epsilon = .125;

				scalar_type P_dist = m_func.dist(P);
				scalar_type P_radius = std::fabs(P_dist);
	
				scalar_type sum = 0.;
				for(Eigen::Index i = 0; i < m_U.rows(); ++i) {	
					sphere_type search_sphere(P, P_radius + epsilon);
					LineT<vector_type> search_line(m_U.row(i), P);
					auto it = get_intersection_iterator(m_func, search_line, search_sphere);		
					if (!it.completed())
						sum += 1;
				}

			return std::copysign(sum / m_U.rows(), P_dist);
		}

		private:
			Eigen::Matrix<scalar_type, Eigen::Dynamic, 2> m_U;
			const func_type& m_func;
		}; // class PartialOcclusionFunctor



		/*
		 * Augment a scalar field with a mask
		 */

		template <class func_type>
		PartialOcclusionFunctor<func_type>
		get_partial_occlusion_functor(const func_type& func,
		                              Eigen::Index sample_count) {
			return PartialOcclusionFunctor<func_type>(func, sample_count);
		}



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
		solve_translation(const Eigen::MatrixXd& A,
		                  const Eigen::MatrixXd& B,
		                  double& score) {
			// Compute C as cross-correlation of A and B
			RealFFT2d real_fft_2d(A.rows());

			Eigen::MatrixXcd U = real_fft_2d.fft(A);
			Eigen::MatrixXcd V = real_fft_2d.fft(B).conjugate();
			U.array() *= V.array();
			Eigen::MatrixXd C = real_fft_2d.ifft(U);	
	
			// Search for the maximum element of C
			Eigen::Index i = 0;
			Eigen::Index j = 0;
			score = C.maxCoeff(&i, &j); 

			// Adjust coordinates
			if (i > (A.rows() / 2))
				i -= A.rows();

			if (j > (A.cols() / 2))
				j -= A.cols();

			// Job done
			return Eigen::Vector2d(i, j);
		}



		template <class func_type>
		void
		fill_density_matrix(const func_type& func,
                        Eigen::MatrixXd& M,
                        const typename func_type::sphere_type& domain, 
                        double resolution) {
			// Compute the box domain
			auto box_domain = BoxT<typename func_type::vector_type>(Eigen::Vector2d(resolution * M.rows(), resolution * M.cols()), domain.center());

			// Fills the matrix
			sample_function_2d(box_domain, M, registration::get_masked_functor(registration::get_partial_occlusion_functor(func, 256), domain));
		}
	} // namespace registration



	template <class left_func_type, class right_func_type>
	Eigen::Transform<double, 2, Eigen::Affine>
	registrate_2d(const left_func_type& left_func,
  	            const right_func_type& right_func,
  	            double resolution,
	              bool complement) {
		// Compute left and right sdf bounding sphere
		auto left_bounding_sphere = left_func.get_bounding_sphere();
		auto right_bounding_sphere = right_func.get_bounding_sphere();

		// Compute left and right domains
		auto left_domain = Sphere2d(left_bounding_sphere.center(), left_bounding_sphere.radius() + 2 * right_bounding_sphere.radius());
		auto right_domain = Sphere2d(right_bounding_sphere.center(), right_bounding_sphere.radius());

		// Compute matric size
		double domain_radius = std::sqrt(std::max(left_domain.squared_radius(), right_domain.squared_radius()));
		Eigen::Index matrix_size = get_smallest_greater_power_of_2(Eigen::Index(std::ceil(2 * domain_radius / resolution)));

		// Discretize the left and right function
		Eigen::MatrixXd left_M(matrix_size, matrix_size);
		Eigen::MatrixXd right_M(matrix_size, matrix_size);

		registration::fill_density_matrix(left_func, left_M, left_domain, resolution);
		registration::fill_density_matrix(right_func, right_M, right_domain, resolution);

		// Switch sign of right function
		if (complement)
			right_M = -right_M;

		// Solve the rotation component
		PseudoPolarFFT2d pseudo_polar_fft(matrix_size);
		FFT2d fft_2d(2 * matrix_size);

		Eigen::MatrixXcd U = pseudo_polar_fft(left_M).cwiseAbs();
		Eigen::MatrixXcd A = fft_2d.fft(U);
		Eigen::MatrixXd A_abs = A.cwiseAbs();
	
		std::vector<Eigen::Index> angle_list;
		angle_list.push_back(0);

		auto rotated_right_func = get_sdf_rotation2(right_func, 0.);
		int max_epoch = 100;
		for(int epoch = 0; epoch < max_epoch; ++epoch) {
			Eigen::MatrixXcd V = pseudo_polar_fft(right_M).cwiseAbs();
			Eigen::MatrixXcd B = fft_2d.fft(V);
			Eigen::MatrixXd B_abs = B.cwiseAbs();

			Eigen::MatrixXcd C = B.conjugate();
			C.array() *= A.array();
			C.array() /= A_abs.array();
			C.array() /= B_abs.array();
			Eigen::MatrixXcd W = fft_2d.ifft(C);
			Eigen::MatrixXd W_abs = W.cwiseAbs();

			Eigen::Index i = 0;
			Eigen::Index j = 0;
			W_abs.maxCoeff(&i, &j);	

			Eigen::Index theta = j;
			if (theta >= matrix_size)
				theta = j - 2 * matrix_size;
		
			Eigen::Index angle = angle_list.back() + theta;
			if (angle < 0)
				angle += 4 * matrix_size;
			if (angle >= 4 * matrix_size)
				angle -= 4 * matrix_size;
		
			if (std::find(angle_list.begin(), angle_list.end(), angle) != angle_list.end())
				break;

			angle_list.push_back(angle);

			rotated_right_func = get_sdf_rotation2(right_func, -(.5 * M_PI / matrix_size) * angle);
			right_bounding_sphere = rotated_right_func.get_bounding_sphere();
			right_domain = Sphere2d(right_bounding_sphere.center(), right_bounding_sphere.radius());
			registration::fill_density_matrix(rotated_right_func, right_M, right_domain, resolution);

			if (complement)
				right_M = -right_M;
 		}

		double angle = (.5 * M_PI / matrix_size) * angle_list.back();
	
		// Solve the translation component
		double hi_score = 0.;
		Eigen::Vector2d hi_T = registration::solve_translation(left_M, right_M, hi_score);

		matrix_rotate_180(right_M);

		double lo_score = 0.;
		Eigen::Vector2d lo_T = registration::solve_translation(left_M, right_M, lo_score);

		Eigen::Vector2d T;
		if (hi_score > lo_score)
			T = hi_T;
		else {
			T = lo_T;
			angle += M_PI;
			rotated_right_func = get_sdf_rotation2(right_func, -angle);
			right_bounding_sphere = rotated_right_func.get_bounding_sphere();
		}

		T *= -2 * domain_radius / matrix_size;
		T += right_bounding_sphere.center() - left_bounding_sphere.center(); 
	
		// Job done
		return
			Eigen::Translation<double, 2>(-T) *
			Eigen::Rotation2D<double>(-angle);
	}
} // namespace saucats



#endif // SAUCATS_SDF_UTILS_REGISTRATION_H

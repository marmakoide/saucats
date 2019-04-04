#ifndef SAUCATS_UTILS_ALGORITHM_H
#define SAUCATS_UTILS_ALGORITHM_H

#include <algorithm>
#include <random>



namespace saucats {
	/*
	 * Return the smallest power of 2 that is greater than x
	 */

	Eigen::Index
	get_smallest_greater_power_of_2(Eigen::Index x) {
		Eigen::Index ret;
		for(ret = 1; ret < x; ret *= 2);
		return ret;
	}



	/*
	 * Return an iterator on an element of a collection
	 */

	template<typename iterator_type, typename rng_type>
	iterator_type
	select_randomly(iterator_type start, iterator_type end, rng_type& rng) {
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(rng));
    return start;
	}



	/*
	 * 180 degrees rotation of a matrix
	 */

	template <class matrix_type>
	void
	matrix_rotate_180(matrix_type& U) {
		U.rowwise().reverseInPlace();
		U.colwise().reverseInPlace();
	}



	/*
	 * Sample a 2d function
	 */

	template <class func_type, class matrix_type>
	void
	sample_function_2d(const BoxT<typename func_type::vector_type>& domain,
                     matrix_type& matrix,
                     const func_type& func) {
		for(Eigen::Index i = 0; i < matrix.rows(); ++i) {
			for(Eigen::Index j = 0; j < matrix.cols(); ++j) {
				Eigen::Vector2d P;

				P.x() = (j + .5) / matrix.cols();			
				P.y() = (i + .5) / matrix.rows();
				P.array() -= .5;
				P.array() *= 2 * domain.half_extent().array();
				P += domain.center();

				matrix.coeffRef(i, j) = func(P);
			}
		}
	}


	
	/*
	 * Sample a 3d function
	 */

	template <class func_type, class tensor_type>
	void
	sample_function_3d(const BoxT<typename func_type::vector_type>& domain,
                     tensor_type& tensor,
                     const func_type& func) {
		for(Eigen::Index i = 0; i < tensor.dimension(0); ++i) {
			for(Eigen::Index j = 0; j < tensor.dimension(1); ++j) {
				for(Eigen::Index k = 0; k < tensor.dimension(2); ++k) {
					Eigen::Vector3d P;
					P << 
					(i + .5) / tensor.dimension(0),
					(j + .5) / tensor.dimension(1),
					(k + .5) / tensor.dimension(2);
				
					P.array() -= .5;
					P.array() *= 2 * domain.half_extent().array();
					P += domain.center();

					tensor.coeffRef(i, j, k) = func(P);
				}
			}
		}
	}
} // namespace saucats



#endif // SAUCATS_UTILS_ALGORITHM_H

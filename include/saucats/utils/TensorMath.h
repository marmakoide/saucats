#ifndef SAUCATS_UTILS_TENSOR_MATH_H
#define SAUCATS_UTILS_TENSOR_MATH_H

#include <unsupported/Eigen/CXX11/Tensor>



namespace saucats {
	/*
	 * Return index of maximum element in a 3d tensor
	 */

	template<class T>
	T
	get_3d_real_tensor_arg_max(const Eigen::Tensor<T, 3, Eigen::RowMajor>& U,
	                           Eigen::Array3i& arg_max) {
		arg_max = Eigen::Array3i::Zero();
		T coeff_max = U(arg_max.coeff(0), arg_max.coeff(1), arg_max.coeff(2));

		for(Eigen::Index i = 0; i < U.dimension(0); ++i)
			for(Eigen::Index j = 0; j < U.dimension(1); ++j)
				for(Eigen::Index k = 0; k < U.dimension(2); ++k)
					if (U(i, j, k) > coeff_max) {
						coeff_max = U(i, j, k);
						arg_max << i, j, k;
					}

		return coeff_max;
	}
} // namespace saucats



#endif // SAUCATS_UTILS_TENSOR_MATH_H

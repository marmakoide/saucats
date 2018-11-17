#ifndef SAUCATS_SDF_UTILS_GRADIENT_H
#define SAUCATS_SDF_UTILS_GRADIENT_H



namespace saucats {
	/*
	 * Estimate the gradient of a signed distance function at given point
	 * Uses central finite difference, evaluating the sdf 2 * N times
	 */


	template <class FuncT, typename VectorT>
	VectorT
	get_sdf_gradient(const FuncT& func,
	                 const VectorT& X,
	                 typename VectorT::Scalar eps) {
		VectorT G;
		for(Eigen::Index i = 0; i < VectorT::RowsAtCompileTime; ++i) {
			VectorT U = X;
			VectorT V = X;
			U[i] += eps;
			V[i] -= eps;
			G[i] = (func.dist(U) - func.dist(V)) / (2 * eps);
		}

		return G;
	}
}; // namespace saucats



#endif // SAUCATS_SDF_UTILS_GRADIENT_H

#ifndef SAUCATS_SDF_AFFINE_TRANSFORM_H
#define SAUCATS_SDF_AFFINE_TRANSFORM_H

#include <Eigen/Geometry>
#include <saucats/geometry/Sphere.h>



namespace saucats {
	/*
	 * Define a signed distance function as the affine transform of a given 
	 * signed distance function
	 */

	template <class func_type>
	class SDFAffineTransform {
	public:
		typedef typename func_type::scalar_type scalar_type;
		typedef typename func_type::vector_type vector_type;
		//typedef Eigen::Transform<scalar_type, 2, Eigen::Affine> transform_type;
		typedef Eigen::Transform<scalar_type, vector_type::RowsAtCompileTime, Eigen::Affine> transform_type;
		typedef SphereT<vector_type> sphere_type;
		


		SDFAffineTransform(const func_type& func,
		                   const transform_type transform) :
			m_func(func),
			m_transform(transform) {
			m_inv_transform = transform.inverse();
		}

		inline transform_type transform() const {
			return m_transform;
		}

		template <typename VectorT>
		inline typename VectorT::Scalar
		dist(const VectorT& X) const {
			VectorT Xp = m_inv_transform * X;
			return m_func.dist(Xp);
		}
	
		inline sphere_type get_bounding_sphere() const {
			sphere_type S = m_func.get_bounding_sphere();
			return sphere_type(m_transform * S.center(), S.radius());
		}

	private:
		func_type m_func;
		transform_type m_transform;
		transform_type m_inv_transform;
	}; // class SDFAffineTransform


 
	template <class func_type>
	SDFAffineTransform<func_type>
	get_sdf_affine_transform(const func_type& func,
	                         const Eigen::Transform<typename func_type::scalar_type, func_type::vector_type::RowsAtCompileTime, Eigen::Affine>& T) {
		return SDFAffineTransform<func_type>(func, T);
	}
} // namespace saucats



#endif // SAUCATS_SDF_AFFINE_TRANSFORM_H

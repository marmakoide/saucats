#ifndef SAUCATS_SDF_ROTATION_H
#define SAUCATS_SDF_ROTATION_H

#include <Eigen/Geometry>
#include <saucats/geometry/Sphere.h>



namespace saucats {
	/*
	 * Define a signed distance function as the rotation of a given signed distance 
	 * function
	 */

	template <class func_type>
	class SDFRotation2 {
	public:
		typedef typename func_type::scalar_type scalar_type;
		typedef typename func_type::vector_type vector_type;
		typedef Eigen::Matrix<scalar_type, 2, 2> matrix_type;
		typedef SphereT<vector_type> sphere_type;



		SDFRotation2(const func_type& func,
		             scalar_type angle) :
			m_func(func),
			m_angle(angle) {
			m_inv_rot_matrix = Eigen::Rotation2D<scalar_type>(-m_angle);
		}

		inline scalar_type angle() const {
			return m_angle;
		}

		template <typename VectorT>
		inline typename VectorT::Scalar
		dist(const VectorT& X) const {
			VectorT Xp = m_inv_rot_matrix * X;
			return m_func.dist(Xp);
		}
	
		inline sphere_type get_bounding_sphere() const {
			sphere_type S = m_func.get_bounding_sphere();
			return sphere_type(Eigen::Rotation2D<scalar_type>(m_angle) * S.center(), S.radius());
		}

	private:
		func_type m_func;
		scalar_type m_angle;
		matrix_type m_inv_rot_matrix;
	}; // class SDFRotation2


 
	template <class func_type>
	SDFRotation2<func_type>
	get_sdf_rotation2(const func_type& func,
	                  typename func_type::scalar_type angle) {
		return SDFRotation2<func_type>(func, angle);
	}



	/*
	 * Define a signed distance function as the rotation of a given signed distance 
	 * function
	 */

	template <class func_type>
	class SDFRotation3 {
	public:
		typedef typename func_type::scalar_type scalar_type;
		typedef typename func_type::vector_type vector_type;
		typedef Eigen::Matrix<scalar_type, 3, 3> matrix_type;
		typedef Eigen::Quaternion<scalar_type> quaternion_type;
		typedef SphereT<vector_type> sphere_type;



		SDFRotation3(const func_type& func,
		             const quaternion_type& Q) :
			m_func(func),
			m_Q(Q) {
			m_inv_rot_matrix = Q.inverse();
		}

		inline quaternion_type Q() const {
			return m_Q;
		}

		template <typename VectorT>
		inline typename VectorT::Scalar
		dist(const VectorT& X) const {
			VectorT Xp = m_inv_rot_matrix * X;
			return m_func.dist(Xp);
		}
	
		inline sphere_type get_bounding_sphere() const {
			sphere_type S = m_func.get_bounding_sphere();
			return sphere_type(m_Q * S.center(), S.radius());
		}

	private:
		func_type m_func;
		quaternion_type m_Q;
		matrix_type m_inv_rot_matrix;
	}; // class SDFRotation3


 
	template <class func_type>
	SDFRotation3<func_type>
	get_sdf_rotation3(const func_type& func,
	                  const Eigen::Quaternion<typename func_type::scalar_type>& Q) {
		return SDFRotation3<func_type>(func, Q);
	}

	template <class func_type>
	SDFRotation3<func_type>
	get_sdf_rotation3(const func_type& func,
	                  const Eigen::AngleAxis<typename func_type::scalar_type>& AA) {
		Eigen::Quaternion<typename func_type::scalar_type> Q;
		Q = AA;
		return SDFRotation3<func_type>(func, Q);
	}
} // namespace saucats



#endif // SAUCATS_SDF_ROTATION_H

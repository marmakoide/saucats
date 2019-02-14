#ifndef SAUCATS_SDF_OFFSET_H
#define SAUCATS_SDF_OFFSET_H

#include <saucats/geometry/Sphere.h>



namespace saucats {
	/*
	 * Define a signed distance function as the offset of a given signed distance 
	 * function
	 */

	template <class func_type>
	class SDFOffset {
	public:
		typedef typename func_type::scalar_type scalar_type;
		typedef typename func_type::vector_type vector_type;
		typedef SphereT<vector_type> sphere_type;



		SDFOffset(const func_type& func,
		          scalar_type offset) :
			m_func(func),
			m_offset(offset) { }

		inline scalar_type& offset() {
			return m_offset;
		}

		inline scalar_type offset() const {
			return m_offset;
		}

		template <typename VectorT>
		inline typename VectorT::Scalar
		dist(const VectorT& X) const {
			return m_func.dist(X) + m_offset;
		}
	
		inline sphere_type get_bounding_sphere() const {
			sphere_type bsphere = m_func.get_bounding_sphere();
			return sphere_type(bsphere.center(), bsphere.radius() - m_offset);
		}

	private:
		func_type m_func;
		scalar_type m_offset;
	}; // class SDFOffset


 
	template <class func_type>
	SDFOffset<func_type>
	get_sdf_offset(const func_type& func,
	               typename func_type::scalar_type offset) {
		return SDFOffset<func_type>(func, offset);
	}
} // namespace saucats



#endif // SAUCATS_SDF_OFFSET_H

#ifndef SAUCATS_SDF_TRANSLATION_H
#define SAUCATS_SDF_TRANSLATION_H

#include <Eigen/Geometry>
#include <saucats/geometry/Sphere.h>



namespace saucats {
	/*
	 * Define a signed distance function as the translation of a given signed distance 
	 * function
	 */

	template <class func_type>
	class SDFTranslation {
	public:
		typedef typename func_type::scalar_type scalar_type;
		typedef typename func_type::vector_type vector_type;
		typedef SphereT<vector_type> sphere_type;



		SDFTranslation(const func_type& func,
		               const vector_type& offset) :
			m_func(func),
			m_offset(offset) { }

		inline vector_type offset() const {
			return m_offset;
		}

		template <typename VectorT>
		inline typename VectorT::Scalar
		dist(const VectorT& X) const {
			return m_func.dist((X - m_offset).eval());
		}
	
		inline sphere_type get_bounding_sphere() const {
			sphere_type S = m_func.get_bounding_sphere();
			return sphere_type(S.center() + m_offset, S.radius());
		}

	private:
		func_type m_func;
		vector_type m_offset;
	}; // class SDFTranslation


 
	template <class func_type>
	SDFTranslation<func_type>
	get_sdf_translation(const func_type& func,
	                    typename func_type::vector_type T) {
		return SDFTranslation<func_type>(func, T);
	}
} // namespace saucats



#endif // SAUCATS_SDF_TRANSLATION_H

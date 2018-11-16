#ifndef SAUCATS_SDF_SPHERE_SDF_H
#define SAUCATS_SDF_SPHERE_SDF_H

#include <saucats/geometry/bounds/SphereBounds.h>



namespace saucats {
	/*
	 Euclidean signed distance to a sphere
	 */

	template <class VectorT>
	class SphereSDF {
	public:
		typedef typename VectorT::Scalar scalar_type;
		typedef VectorT vector_type;
		typedef SphereT<VectorT> sphere_type;



		inline SphereSDF(const sphere_type& sphere)
			: m_radius(sphere.radius()),
			  m_sphere(sphere) { }

		inline const sphere_type& sphere() const {
			return m_sphere;
		}

		template <typename InVectorT>
		inline typename InVectorT::Scalar
		dist(const InVectorT& X) const {
			return (m_sphere.center() - X).norm() - m_radius;
		}

		inline sphere_type bounding_sphere() const {
			return get_bounding_sphere(m_sphere);
		}

	private:
		scalar_type m_radius;
		sphere_type m_sphere;
	}; // class SphereSDF



	template <class VectorT>
	SphereSDF<VectorT>
	get_sphere_sdf(const SphereT<VectorT>& sphere) {
		return SphereSDF<VectorT>(sphere);
	}
} // namespace saucats



#endif // SAUCATS_SDF_SPHERE_SDF_H

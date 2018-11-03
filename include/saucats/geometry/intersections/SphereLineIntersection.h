#ifndef SAUCATS_GEOMETRY_INTERSECTION_SPHERE_LINE_H
#define SAUCATS_GEOMETRY_INTERSECTION_SPHERE_LINE_H

#include <saucats/geometry/Line.h>
#include <saucats/geometry/Sphere.h>



namespace saucats {
	/*
	 * The SphereLineIntersectionT represents the intersection between a
	 * N-dimensional sphere and an N-dimensional line.
	 */

	template <class VectorT>
	class SphereLineIntersectionT {
	public:
		typedef typename VectorT::Scalar scalar_type;



		inline SphereLineIntersectionT() :
			m_valid(false),
			m_tmin(0.),
			m_tmax(0.) { }

		inline SphereLineIntersectionT(scalar_type tmin, scalar_type tmax) :
			m_valid(true),
			m_tmin(tmin),
			m_tmax(tmax) { }

		inline SphereLineIntersectionT(const SphereLineIntersectionT& other) :
			m_valid(other.m_valid),
			m_tmin(other.m_tmin),
			m_tmax(other.m_tmax) { }

		inline bool is_valid() const {
			return m_valid;
		}

		inline scalar_type tmin() const {
			return m_tmin;
		}

		inline scalar_type tmax() const {
			return m_tmax;
		}

	private:
		bool m_valid;
		scalar_type m_tmin, m_tmax;
	}; // class SphereLineIntersectionT



	template <class VectorT>
	SphereLineIntersectionT<VectorT> 
	get_sphere_line_intersection(const SphereT<VectorT>& sphere,
	                             const LineT<VectorT>& line) {
		typedef typename VectorT::Scalar scalar_type;

		VectorT diff = sphere.center() - line.origin();
		scalar_type t = line.direction().dot(diff);
		scalar_type d_sqr = sphere.radius_sqr() - diff.squaredNorm() + t * t;
		if (d_sqr < 0.)
			return SphereLineIntersectionT<VectorT>();

		scalar_type d = std::sqrt(d_sqr);
		return SphereLineIntersectionT<VectorT>(t - d, t + d);
	}
} // namespace saucats



#endif // SAUCATS_GEOMETRY_INTERSECTION_SPHERE_LINE_H

#ifndef SAUCATS_GEOMETRY_INTERSECTION_SPHERE_SPHERE_H
#define SAUCATS_GEOMETRY_INTERSECTION_SPHERE_SPHERE_H

#include <saucats/geometry/Circle.h>
#include <saucats/geometry/Sphere.h>



namespace saucats {
/*
   Sphere/Sphere intersection

   The intersection is a circle of center M, embedded on a plane centered on M
   and normal N. Even if there's no intersection, the separator plane is computed.
   This gives a plane guaranted to separate two spheres.
 */

	template <class VectorT>
	class SphereSphereIntersectionT {
	public:
		typedef typename VectorT::Scalar scalar_type;
		typedef CircleT<VectorT> circle_type;



		inline SphereSphereIntersectionT() :
			m_valid(false),
			m_MA_dist(0),
			m_MB_dist(0) { }

		inline SphereSphereIntersectionT(const circle_type& circle,
		                                 scalar_type MA_dist,
		                                 scalar_type MB_dist,
			                               bool valid) :
			m_valid(valid),
			m_circle(circle),
			m_MA_dist(MA_dist),
			m_MB_dist(MB_dist) { }

		inline SphereSphereIntersectionT(const SphereSphereIntersectionT& other) :
			m_valid(other.m_valid),
			m_circle(other.m_circle),
			m_MA_dist(other.m_MA_dist),
			m_MB_dist(other.m_MB_dist) { }

		inline bool is_valid() const {
			return m_valid;
		}

		inline const circle_type& circle() const {
			return m_circle;
		}

		inline scalar_type MA_dist() const {
			return m_MA_dist;
		}

		inline scalar_type MB_dist() const {
			return m_MB_dist;
		}

	private:
		bool m_valid;          // Tell if we have an intersection or not
		circle_type m_circle;  // The intersection itself. Even if it's invalid, center and normal will have consistent values
		scalar_type m_MA_dist; // Distance from A center to separator plane origin
		scalar_type m_MB_dist; // Distance from A center to separator plane origin
	}; // class SphereSphereIntersectionT



	template <class VectorT>
	SphereSphereIntersectionT<VectorT> 
	get_sphere_sphere_intersection(const SphereT<VectorT>& A,
	                               const SphereT<VectorT>& B) {
		typedef typename VectorT::Scalar scalar_type;
		typedef CircleT<VectorT> circle_type;

		VectorT AB = B.center() - A.center();
		scalar_type AB_norm_sqr = AB.squaredNorm();
		scalar_type AB_norm = std::sqrt(AB_norm_sqr);
		Eigen::Vector3d U = AB / AB_norm;

		scalar_type k = (A.radius_sqr() - B.radius_sqr() + AB_norm_sqr) / (2 * AB_norm);
		VectorT M = A.center() + k * U;
		scalar_type h2 = A.radius_sqr() - k * k;

		SphereSphereIntersectionT<VectorT> ret;
		if (h2 >= 0.)
			ret = SphereSphereIntersectionT<VectorT>(circle_type(M, U, std::sqrt(h2)), k, AB_norm - k, true);
		else
			ret = SphereSphereIntersectionT<VectorT>(circle_type(M, U, 0), k, AB_norm - k, false);

		return ret;
	}
} // namespace saucats



#endif // SAUCATS_GEOMETRY_INTERSECTION_SPHERE_SPHERE_H

#ifndef SAUCATS_GEOMETRY_SPHERE_H
#define SAUCATS_GEOMETRY_SPHERE_H

#include <Eigen/Dense>



namespace saucats {
	/*
	 * The SphereT class represents a N-dimensional sphere
	 */

	template <class VectorT>
	class SphereT {
	public:
		typedef typename VectorT::Scalar scalar_type;
		typedef VectorT vector_type;



		inline SphereT() :
			m_radius_sqr(0),
			m_center(vector_type::Zero()) { }

		inline SphereT(scalar_type radius) :
			m_radius_sqr(radius * radius),
			m_center(vector_type::Zero()) { }

		inline SphereT(const vector_type& center,
		               scalar_type radius) :
			m_radius_sqr(radius * radius),
			m_center(center) { }

		inline SphereT(const vector_type& center) :
			m_radius_sqr(0),
			m_center(center) { }

		inline SphereT(const SphereT& sphere) :
			m_radius_sqr(sphere.m_radius_sqr),
			m_center(sphere.m_center) { }

		inline SphereT& operator = (const SphereT& sphere) {
			m_radius_sqr = sphere.m_radius_sqr;
			m_center = sphere.m_center;
			return *this;
		}

		// Returns the radius of the sphere
		inline scalar_type radius() const {
			return std::sqrt(m_radius_sqr);
		}

		// Returns the square of the radius of the sphere
		inline scalar_type squared_radius() const {
			return m_radius_sqr;
		}

		// Returns the center of the sphere
		inline vector_type& center() {
			return m_center;
		}

		// Returns the center of the sphere
		inline const vector_type& center() const {
			return m_center;
		}
		
		// Returns true if a point is inside the sphere
		inline bool contains(const vector_type& X) const {
			return (X - m_center).squaredNorm() <= m_radius_sqr;
		}

	private:
		scalar_type m_radius_sqr;
		vector_type m_center;
	}; // class SphereT



	typedef SphereT<Eigen::Vector2f> Sphere2f;
	typedef SphereT<Eigen::Vector2d> Sphere2d;
	typedef SphereT<Eigen::Vector3f> Sphere3f;
	typedef SphereT<Eigen::Vector3d> Sphere3d;
} // namespace saucats



#endif // SAUCATS_GEOMETRY_SPHERE_H

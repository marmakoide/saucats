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

		// Surface area of the sphere (in 2d, it equals to the perimeter)
		inline scalar_type surface_area() const {
			return UNIT_BALL_SURFACE_AREA * std::pow(std::sqrt(m_radius_sqr), vector_type::RowsAtCompileTime - 1);
		}

		// Volume of the sphere (in 2d, it equals to the surface area)
		inline scalar_type volume() const {
			return UNIT_BALL_VOLUME * std::pow(std::sqrt(m_radius_sqr), vector_type::RowsAtCompileTime);
		}

	private:
		static const double UNIT_BALL_SURFACE_AREA;
		static const double UNIT_BALL_VOLUME;

		scalar_type m_radius_sqr;
		vector_type m_center;
	}; // class SphereT



	typedef SphereT<Eigen::Vector2f> Sphere2f;
	typedef SphereT<Eigen::Vector2d> Sphere2d;
	typedef SphereT<Eigen::Vector3f> Sphere3f;
	typedef SphereT<Eigen::Vector3d> Sphere3d;



	/*
	 * Compilation time computation of the volume and surface areas constants
	 */

	constexpr double get_unit_ball_surface_area(int n);
	constexpr double get_unit_ball_volume(int n);

	constexpr double
	get_unit_ball_surface_area(int n) {
		return n == 0 ? 2. : (2 * M_PI) * get_unit_ball_volume(n - 1);
	}

	constexpr double
	get_unit_ball_volume(int n) {
		return n == 0 ? 1. : get_unit_ball_surface_area(n - 1) / n;
	}



	template <class VectorT> const double
	SphereT<VectorT>::UNIT_BALL_SURFACE_AREA = get_unit_ball_surface_area(VectorT::RowsAtCompileTime - 1);

	template <class VectorT> const double
	SphereT<VectorT>::UNIT_BALL_VOLUME = get_unit_ball_volume(VectorT::RowsAtCompileTime);
} // namespace saucats



#endif // SAUCATS_GEOMETRY_SPHERE_H

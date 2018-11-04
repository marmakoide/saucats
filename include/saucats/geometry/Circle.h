#ifndef SAUCATS_GEOMETRY_CIRCLE_H
#define SAUCATS_GEOMETRY_CIRCLE_H

#include <Eigen/Dense>



// Represents a circle embedded on a straight plane
namespace saucats {
	/*
	 * The CircleT class represents a circle embedded on a N-dimensional plane
	 */

	template <class VectorT>
	class CircleT {
	public:
		typedef typename VectorT::Scalar scalar_type;
		typedef VectorT vector_type;



		inline CircleT() :
			m_radius(0),
			m_center(vector_type::Zero()),
			m_normal(vector_type::Unit(vector_type::RowsAtCompileTime, vector_type::RowsAtCompileTime - 1)) { }

		inline CircleT(scalar_type radius) :
			m_radius(radius),
			m_center(vector_type::Zero()),
			m_normal(vector_type::Unit(vector_type::RowsAtCompileTime, vector_type::RowsAtCompileTime - 1)) { }

		inline CircleT(const vector_type& center, scalar_type radius) :
			m_radius(radius),
			m_center(center),
			m_normal(vector_type::Unit(vector_type::RowsAtCompileTime, vector_type::RowsAtCompileTime - 1)) { }

		inline CircleT(const vector_type& center) :
			m_radius(0),
			m_center(center),
			m_normal(vector_type::Unit(vector_type::RowsAtCompileTime, vector_type::RowsAtCompileTime - 1)) { }

		inline CircleT(const vector_type& center, const vector_type& normal, scalar_type radius) :
			m_radius(radius),
			m_center(center),
			m_normal(normal) { }

		inline CircleT(const vector_type& center, const vector_type& normal) :
			m_radius(0),
			m_center(center),
			m_normal(normal) { }

		inline CircleT(const CircleT& circle) :
			m_radius(circle.m_radius),
			m_center(circle.m_center),
			m_normal(circle.m_normal) { }

		inline CircleT& operator = (const CircleT& circle) {
			m_radius = circle.m_radius;
			m_center = circle.m_center;
			m_normal = circle.m_normal;
			return *this;
		}

		// Returns the radius of the circle
		inline scalar_type radius() const {
			return m_radius;
		}

		// Returns the center of the circle
		inline vector_type& center() {
			return m_center;
		}

		// Returns the center of the circle
		inline const vector_type& center() const {
			return m_center;
		}

		// Returns the normal of the plane the circle is embedded on
		inline vector_type& normal() {
			return m_normal;
		}

		// Returns the normal of the plane the circle is embedded on
		inline const vector_type& normal() const {
			return m_normal;
		}

		// Returns the square of the Euclidean distance from a point to the circle
		scalar_type
		squared_dist(const vector_type& X) const {
			vector_type Y = X - m_center;
			scalar_type d = m_normal.dot(Y);
			vector_type H = Y - d * m_normal;
			Eigen::Matrix<scalar_type, 2, 1> I(H.norm() - m_radius, d);
			return I.squaredNorm();
		}

		// Returns the Euclidean distance from a point to the circle
		inline scalar_type
		dist(const vector_type& X) const {
			return std::sqrt(squared_dist(X));
		}

	private:
		scalar_type m_radius;
		vector_type m_center;
		vector_type m_normal;
	}; // class CircleT



	typedef CircleT<Eigen::Vector3f> Circle3f;
	typedef CircleT<Eigen::Vector3d> Circle3d;
} // namespace saucats



#endif // SAUCATS_GEOMETRY_CIRCLE_H

#ifndef SAUCATS_GEOMETRY_PLANE_H
#define SAUCATS_GEOMETRY_PLANE_H

#include <Eigen/Dense>



// Represents an infinite plane in N dimension, ie. an hyperplane
namespace saucats {
	template <class VectorT>
	class PlaneT {
	public:
		typedef typename VectorT::Scalar scalar_type;
		typedef VectorT vector_type;



		inline PlaneT() :
			m_normal(vector_type::Unit(vector_type::RowsAtCompileTime, vector_type::RowsAtCompileTime - 1)),
			m_constant(0) { }

		inline PlaneT(const vector_type& normal) :
			m_normal(normal),
			m_constant(0) { }

		inline PlaneT(const vector_type& normal,
		              const vector_type& origin) :
			m_normal(normal),
			m_constant(-origin.dot(normal)) { }

		inline PlaneT(const PlaneT& plane) :
			m_normal(plane.m_normal),
			m_constant(plane.m_constant) { }

		inline PlaneT& operator = (const PlaneT& plane) {
			m_normal = plane.m_normal;
			m_constant = plane.m_constant;
			return *this;
		}

		inline const vector_type& normal() const {
			return m_normal;
		}

		inline scalar_type signed_dist(const vector_type& X) const {
			return X.dot(m_normal) + m_constant;
		}

		/*
		// Build a plane centered on A, and with unit normal AB ^ AC 
		static PlaneT from_3_points(const Eigen::Vector3d& A,
		                           const Eigen::Vector3d& B,
		                           const Eigen::Vector3d& C) {
			return Plane(((B - A).cross(C - A)).normalized(), A);
		}
		*/

	private:
		vector_type m_normal;
		scalar_type m_constant;
	}; // class PlaneT



	typedef PlaneT<Eigen::Vector3d> Plane3d;
} // namespace saucats



#endif // SAUCATS_GEOMETRY_PLANE_H

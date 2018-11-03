#ifndef SAUCATS_GEOMETRY_LINE_H
#define SAUCATS_GEOMETRY_LINE_H

#include <Eigen/Dense>



namespace saucats {
	/*
	 * The LineT class represents a N-dimensional infinite line
	 */

	template <class VectorT>
	class LineT {
	public:
		typedef typename VectorT::Scalar scalar_type;
		typedef VectorT vector_type;



		inline LineT() :
			m_direction(vector_type::Unit(vector_type::RowsAtCompileTime, vector_type::RowsAtCompileTime - 1)),
			m_origin(vector_type::Zero()) { }

		inline LineT(const vector_type& direction) :
			m_direction(direction),
			m_origin(vector_type::Zero()) { }

		inline LineT(const vector_type& direction,
		             const vector_type& origin) :
			m_direction(direction),
			m_origin(origin) { }

		inline LineT(const LineT& line) :
			m_direction(line.m_direction),
			m_origin(line.m_origin) { }

		inline LineT& operator = (const LineT& line) {
			m_direction = line.m_direction;
			m_origin = line.m_origin;
			return *this;
		}

		// Return the normalized direction of the line
		inline const vector_type& direction() const {
			return m_direction;
		}

		// Return the origin of the line
		inline vector_type& origin() {
			return m_origin;
		}

		// Return the origin of the line
		inline const vector_type& origin() const {
			return m_origin;
		}

		// Return the point on the line with parameter t
		inline vector_type point(scalar_type t) const {
			return m_origin + t * m_direction;
		}

		// Build a line with origin A and direction AB 
		static LineT from_2_points(const vector_type& A,
		                           const vector_type& B) {
			return LineT((B - A).normalized(), A);
		}

	private:
		vector_type m_direction;
		vector_type m_origin;
	}; // class LineT



	typedef LineT<Eigen::Vector2f> Line2f;
	typedef LineT<Eigen::Vector2d> Line2d;
	typedef LineT<Eigen::Vector3f> Line3f;
	typedef LineT<Eigen::Vector3d> Line3d;
} // namespace saucats



#endif // SAUCATS_GEOMETRY_LINE_H

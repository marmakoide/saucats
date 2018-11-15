#ifndef SAUCATS_GEOMETRY_SEGMENT_H
#define SAUCATS_GEOMETRY_SEGMENT_H

#include <saucats/geometry/Line.h>



namespace saucats {
	/*
	 * The SegmentT class represent a segment from an infinite straight line in N dimension
	 */

	template <class VectorT>
	class SegmentT {
	public:
		typedef typename VectorT::Scalar scalar_type;
		typedef VectorT vector_type;
		typedef LineT<VectorT> line_type;



		inline SegmentT() :
			m_line(line_type()),
			m_half_length(0) { }

		inline SegmentT(scalar_type half_length) :
			m_line(line_type()),
			m_half_length(half_length) { }

		inline SegmentT(const line_type& line,
	  	              scalar_type half_length) :
			m_line(line),
			m_half_length(half_length) { }

		inline SegmentT(const SegmentT& segment) :
			m_line(segment.m_line),
			m_half_length(segment.m_half_length) { }

		inline SegmentT& operator = (const SegmentT& segment) {
			m_line = segment.m_line;
			m_half_length = segment.m_half_length;
			return *this;
		}

		inline const line_type& line() const {
			return m_line;
		}

		inline scalar_type half_length() const {
			return m_half_length;
		}

		// Build a segment with end point A (left) and B (right) 
		static SegmentT from_2_points(const vector_type& A,
		                              const vector_type& B) {
			vector_type U = B - A;
			scalar_type U_norm = U.norm();
			return SegmentT(line_type(U / U_norm, .5 * (A + B)), .5 * U_norm);
		}

	private:
		line_type m_line;
		scalar_type m_half_length;
	}; // class SegmentT



	typedef SegmentT<Eigen::Vector2f> Segment2f;
	typedef SegmentT<Eigen::Vector2d> Segment2d;
	typedef SegmentT<Eigen::Vector3f> Segment3f;
	typedef SegmentT<Eigen::Vector3d> Segment3d;
} // namespace saucats



#endif // SAUCATS_GEOMETRY_SEGMENT_H

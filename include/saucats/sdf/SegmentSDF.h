#ifndef SAUCATS_SDF_SEGMENT_SDF_H
#define SAUCATS_SDF_SEGMENT_SDF_H

#include <saucats/geometry/Segment.h>



namespace saucats {
	/*
	 Euclidean signed distance to a segment
	 */

	template <class VectorT>
	class SegmentSDF {
	public:
		typedef typename VectorT::Scalar scalar_type;
		typedef VectorT vector_type;
		typedef SegmentT<VectorT> segment_type;
		typedef SphereT<VectorT> sphere_type;



		inline SegmentSDF(const segment_type& segment)
			: m_segment(segment) { }

		template <typename InVectorT>
		inline typename InVectorT::Scalar
		dist(const InVectorT& X) const {
			InVectorT U = X - m_segment.line().origin();
			typename InVectorT::Scalar x = m_segment.line().direction().dot(U);
			typename InVectorT::Scalar k = std::copysign(std::min(std::fabs(x), m_segment.half_length()), x);
			return std::sqrt(U.squaredNorm() + k * (k - 2 * x));				
		}

		inline sphere_type get_bounding_sphere() const {
			return sphere_type(m_segment.line().origin(), m_segment.half_length());
		}

	private:
		segment_type m_segment;
	}; // class SegmentSDF



	template <class VectorT>
	SegmentSDF<VectorT>
	get_segment_sdf(const SegmentT<VectorT>& segment) {
		return SegmentSDF<VectorT>(segment);
	}
} // namespace saucats



#endif // SAUCATS_SDF_SEGMENT_SDF_H

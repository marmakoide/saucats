#ifndef SAUCATS_SDF_SEGMENT_SDF_H
#define SAUCATS_SDF_SEGMENT_SDF_H

#include <saucats/geometry/bounds/SegmentBounds.h>



namespace saucats {
	/*
	 Euclidean signed distance to a segment
	 */

	template <class VectorT>
	class SegmentSDF {
	public:
		typedef typename VectorT::Scalar scalar_type;
		typedef VectorT vector_type;
		typedef SegmentT<vector_type> segment_type;
		typedef SphereT<vector_type> sphere_type;


		inline SegmentSDF() { }

		inline SegmentSDF(const segment_type& segment)
			: m_segment(segment) { }

		inline const segment_type& segment() const {
			return m_segment;
		}

		template <typename InVectorT>
		inline typename InVectorT::Scalar
		dist(const InVectorT& X) const {
			return std::sqrt(squared_dist(X));				
		}

		template <typename InVectorT>
		inline typename InVectorT::Scalar
		squared_dist(const InVectorT& X) const {
			InVectorT U = X - m_segment.line().origin();
			typename InVectorT::Scalar x = m_segment.line().direction().dot(U);
			typename InVectorT::Scalar k = std::copysign(std::min(std::fabs(x), m_segment.half_length()), x);
			return U.squaredNorm() + k * (k - 2 * x);				
		}

		inline sphere_type get_bounding_sphere() const {
			return get_bounding_sphere(m_segment);
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

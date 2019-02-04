#ifndef SAUCATS_SDF_POLYGON_SDF_H
#define SAUCATS_SDF_POLYGON_SDF_H

#include <saucats/utils/ArrayT.h>
#include <saucats/sdf/SegmentSDF.h>
#include <saucats/geometry/bounds/PointListBounds.h>



namespace saucats {
	/*
	 Euclidean signed distance to a 2d polygon
	 */
	template <class VectorT>
	class PolygonSDF {
	public:
		typedef typename VectorT::Scalar scalar_type;
		typedef VectorT vector_type;
		typedef SegmentT<VectorT> segment_type;
		typedef SegmentSDF<VectorT> segment_sdf_type;
		typedef SphereT<VectorT> sphere_type;

		typedef Eigen::Matrix<scalar_type, Eigen::Dynamic, 2> vertex_array_type;
		typedef Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 2> edge_array_type;


		
		inline PolygonSDF(const vertex_array_type& vertex_array,
		                  const edge_array_type& edge_array) :
			m_segment_sdf_array(edge_array.rows()) {
			// Compute the segments
			for(Eigen::Index i = 0; i < edge_array.rows(); ++i)
				m_segment_sdf_array[i] = segment_sdf_type(segment_type::from_2_points(vertex_array.row(edge_array(i, 0)), vertex_array.row(edge_array(i, 1))));

			// Compute the bounding disc
			m_bounding_sphere = point_collection::get_bounding_sphere(vertex_array.rowwise().begin(), vertex_array.rowwise().end());
		}
		


		template <typename InVectorT>
		inline typename InVectorT::Scalar
		dist(const InVectorT& X) const {
			double dist = 
				std::accumulate(std::next(m_segment_sdf_array.begin()), m_segment_sdf_array.end(),
				                std::fabs(m_segment_sdf_array.begin()->dist(X)),
				                [&X] (const segment_sdf_type& sdf, scalar_type min_dist) {
													return std::min(min_dist, std::fabs(sdf.dist(X)));
												});
			
			if (is_inside(X))
				dist = -dist;

			return dist;
		}

		inline sphere_type get_bounding_sphere() const {
			return m_bounding_sphere;
		}

	private:
		// Check whether a point is inside or outside a polygon
		template <typename InVectorT>
		inline bool
		is_inside(const InVectorT& P) const {	
			// Compute winding number
			Eigen::Index winding_number = 0;
			for(const segment_sdf_type& seg_sdf : m_segment_sdf_array) {
				const segment_type& seg = seg_sdf.segment();

				vector_type A = seg.line().point(-seg.half_length());
				vector_type B = seg.line().point( seg.half_length());
				vector_type U = (2 * seg.half_length()) * seg.line().direction();

				if (A.y() <= P.y()) {
					if (B.y() > P.y())
						if (U.cross(P - A) > 0)
							winding_number += 1;
				}
				else {
				if (B.y() <= P.y())
						if (U.cross(P - A) < 0)
							winding_number -= 1;
				}
			}
			
			// Job done
			return winding_number != 0;	
		}



		sphere_type m_bounding_sphere;
		ArrayT<segment_sdf_type> m_segment_sdf_array;
	}; // class PolygonSDF



	template <class VectorT>
	PolygonSDF<VectorT>
	get_polygon_sdf(const BoxT<VectorT>& box) {
		return PolygonSDF<VectorT>(box);
	}
} // namespace saucats



#endif // SAUCATS_SDF_POLYGON_SDF_H

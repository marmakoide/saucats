#ifndef SAUCATS_SDF_POLYGON_SDF_H
#define SAUCATS_SDF_POLYGON_SDF_H

#include <saucats/utils/ArrayT.h>
#include <saucats/sdf/SegmentSDF.h>
#include <saucats/geometry/bounds/PointListBounds.h>



namespace saucats {
	/*
	 Euclidean signed distance to a 2d polygon
	 */
	template <class ScalarT>
	class PolygonSDF {
	public:
		typedef ScalarT scalar_type;
		typedef Eigen::Matrix<scalar_type, 2, 1> vector_type;
		typedef SegmentT<vector_type> segment_type;
		typedef SegmentSDF<vector_type> segment_sdf_type;
		typedef SphereT<vector_type> sphere_type;

		typedef Eigen::Matrix<scalar_type, Eigen::Dynamic, 2> vertex_array_type;
		typedef Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 2> edge_array_type;


		
		inline PolygonSDF(const vertex_array_type& vertex_array,
		                  const edge_array_type& edge_array) :
			m_segment_sdf_array(edge_array.rows()) {
			// Compute the segments
			for(Eigen::Index i = 0; i < edge_array.rows(); ++i)
				m_segment_sdf_array[i] = segment_sdf_type(segment_type::from_2_points(vertex_array.row(edge_array(i, 0)), vertex_array.row(edge_array(i, 1))));

			// Compute the bounding disc
			// TODO : soon-to-be-release version of Eigen will interface directly with STL routines
			//m_bounding_sphere = point_collection::get_bounding_sphere(vertex_array.rowwise().begin(), vertex_array.rowwise().end());

			std::vector<vector_type> points(vertex_array.size());
			for(Eigen::Index i = 0; i < vertex_array.rows(); ++i)
				points[i] = vertex_array.row(i);
			m_bounding_sphere = point_collection::get_bounding_sphere(points.begin(), points.end());
		}
		


		template <typename InVectorT>
		inline typename InVectorT::Scalar
		dist(const InVectorT& X) const {
			auto seg_it = m_segment_sdf_array.begin();
			scalar_type dist = std::fabs(seg_it->dist(X));
			for(++seg_it; seg_it != m_segment_sdf_array.end(); ++seg_it)
				dist = std::min(dist, std::fabs(seg_it->dist(X)));

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

				const vector_type& U = seg.line().direction();
				scalar_type ly = seg.half_length() * U.y();

				if (-ly <= P.y() - seg.line().origin().y()) {
					if (ly > P.y() - seg.line().origin().y()) {
						vector_type V = P - seg.line().point(-seg.half_length());
						if (U.x() * V.y() - U.y() * V.x() > 0)
							winding_number += 1;
					}
				}
				else {
					if (ly <= P.y() - seg.line().origin().y()) {
						vector_type V = P - seg.line().point(-seg.half_length());
						if (U.x() * V.y() - U.y() * V.x() < 0)
							winding_number -= 1;
					}
				}
			}
			
			// Job done
			return winding_number != 0;	
		}



		sphere_type m_bounding_sphere;
		ArrayT<segment_sdf_type> m_segment_sdf_array;
	}; // class PolygonSDF



	template <class ScalarT>
	PolygonSDF<ScalarT>
	get_polygon_sdf(const Eigen::Matrix<ScalarT, Eigen::Dynamic, 2>& vertex_array,
	                const Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 2>& edge_array) {
		return PolygonSDF<ScalarT>(vertex_array, edge_array);
	}


	template <class ScalarT>
	PolygonSDF<ScalarT>
	get_polygon_sdf(const Eigen::Matrix<ScalarT, Eigen::Dynamic, 2>& vertex_array) {
		Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 2> edge_array(vertex_array.rows(), 2);
		for(Eigen::Index i = 0; i < vertex_array.rows() - 1; ++i) {
			edge_array.coeffRef(i, 0) = i;
			edge_array.coeffRef(i, 1) = i + 1;
		}
		edge_array.coeffRef(vertex_array.rows() - 1, 0) = vertex_array.rows() - 1;
		edge_array.coeffRef(vertex_array.rows() - 1, 1) = 0;


		return PolygonSDF<ScalarT>(vertex_array, edge_array);
	}
} // namespace saucats



#endif // SAUCATS_SDF_POLYGON_SDF_H

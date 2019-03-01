#ifndef SAUCATS_SDF_POLYGON2_SDF_H
#define SAUCATS_SDF_POLYGON2_SDF_H

#include <saucats/utils/ArrayT.h>
#include <saucats/sdf/SegmentSDF.h>
#include <saucats/geometry/BoxTree.h>
#include <saucats/geometry/bounds/PointListBounds.h>



namespace saucats {
	/*
	 Euclidean signed distance to a 2d polygon
	 */
	template <class ScalarT>
	class Polygon2SDF {
	public:
		typedef ScalarT scalar_type;
		typedef Eigen::Matrix<scalar_type, 2, 1> vector_type;
		typedef SegmentT<vector_type> segment_type;
		typedef SegmentSDF<vector_type> segment_sdf_type;
		typedef SphereT<vector_type> sphere_type;

		typedef BoxT<vector_type> box_type;
		typedef BoxTreeT<box_type, std::size_t> boxtree_type;

		typedef Eigen::Matrix<scalar_type, Eigen::Dynamic, 2> vertex_array_type;
		typedef Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 2> edge_array_type;



		Polygon2SDF(const vertex_array_type& vertex_array,
		            const edge_array_type& edge_array) :
			m_segment_sdf_array(edge_array.rows()) {
			// Compute the segments
			for(Eigen::Index i = 0; i < edge_array.rows(); ++i)
				m_segment_sdf_array[i] = segment_sdf_type(segment_type::from_2_points(vertex_array.row(edge_array(i, 0)), vertex_array.row(edge_array(i, 1))));

			// Build the boxtree
			std::vector<std::size_t> id_list(edge_array.rows());
			std::iota(id_list.begin(), id_list.end(), 0);
			
			m_boxtree =
				boxtree_type::from_kd_construction(id_list.begin(), id_list.end(), 
					[this](std::size_t i) -> box_type { return get_bounding_box(m_segment_sdf_array[i].segment()); }
				);

			// Compute the bounding disc
			// TODO : soon-to-be-release version of Eigen will interface directly with STL routines
			//m_bounding_sphere = point_collection::get_bounding_sphere(vertex_array.rowwise().begin(), vertex_array.rowwise().end());

			std::vector<vector_type> points(vertex_array.rows());
			for(Eigen::Index i = 0; i < vertex_array.rows(); ++i)
				points[i] = vertex_array.row(i);
			m_bounding_sphere = point_collection::get_bounding_sphere(points.begin(), points.end());
		}
		
		template <typename InVectorT>
		inline typename InVectorT::Scalar
		dist(const InVectorT& X) const {
			typedef typename boxtree_type::Node const* node_ptr_type;
			// Find the distance to the closest triangle
			bool first_hit = true;
			typename InVectorT::Scalar min_dist = 0;

			std::stack<node_ptr_type> stack;
			stack.push(m_boxtree.root());
			while(!stack.empty()) {
				// Pop a node
				node_ptr_type node = stack.top();
				stack.pop();

				// Skip nodes that can't possibly closer than the closest segment found so far
				if ((!first_hit) and (node->box().dist(X) >= min_dist))
					continue;

				// If we reached a leaf node, update closest segment found so far
				if (node->is_leaf()) {
					if (first_hit) {
						min_dist = m_segment_sdf_array[node->value()].dist(X);
						first_hit = false;
					}
					else
						min_dist = std::min(min_dist, m_segment_sdf_array[node->value()].dist(X));
				}
				// If we are not in a leaf node 
				else {
					node_ptr_type left_child = node->children().first;
					node_ptr_type right_child = node->children().second;

					if ((left_child->box().center() - X).squaredNorm() > (right_child->box().center() - X).squaredNorm())
						std::swap(left_child, right_child);

					stack.push(right_child);
					stack.push(left_child);
				}
			}

			// Negative distance if inside the mesh
			if (is_inside(X))
				min_dist = -min_dist;

			// Job done
			return min_dist;
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



		ArrayT<segment_sdf_type> m_segment_sdf_array;
		sphere_type m_bounding_sphere;
		boxtree_type m_boxtree;
	}; // class Polygon2SDF



	template <class ScalarT>
	Polygon2SDF<ScalarT>
	get_polygon2_sdf(const Eigen::Matrix<ScalarT, Eigen::Dynamic, 2>& vertex_array,
	                const Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 2>& edge_array) {
		return Polygon2SDF<ScalarT>(vertex_array, edge_array);
	}


	template <class ScalarT>
	Polygon2SDF<ScalarT>
	get_polygon2_sdf(const Eigen::Matrix<ScalarT, Eigen::Dynamic, 2>& vertex_array) {
		Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 2> edge_array(vertex_array.rows(), 2);
		for(Eigen::Index i = 0; i < vertex_array.rows() - 1; ++i) {
			edge_array.coeffRef(i, 0) = i;
			edge_array.coeffRef(i, 1) = i + 1;
		}
		edge_array.coeffRef(vertex_array.rows() - 1, 0) = vertex_array.rows() - 1;
		edge_array.coeffRef(vertex_array.rows() - 1, 1) = 0;


		return Polygon2SDF<ScalarT>(vertex_array, edge_array);
	}
} // namespace saucats



#endif // SAUCATS_SDF_POLYGON2_SDF_H

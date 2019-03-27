#ifndef SAUCATS_SDF_TRIMESH_SDF_H
#define SAUCATS_SDF_TRIMESH_SDF_H

#include <saucats/geometry/bounds/PointListBounds.h>
#include <saucats/sdf/Triangle3SDF.h>
#include <saucats/geometry/BoxTree.h>
#include <saucats/utils/KahanSum.h>



namespace saucats {
	template <class ScalarT>
	class TriMeshSDF {
	public:
		typedef ScalarT scalar_type;
		typedef Eigen::Matrix<scalar_type, 3, 1> vector_type;
		typedef Eigen::Matrix<scalar_type, 3, 3> triangle_type;
		typedef Triangle3SDF<scalar_type> triangle_sdf_type;
		typedef SphereT<vector_type> sphere_type;
		typedef BoxT<vector_type> box_type;
		typedef BoxTreeT<box_type, std::size_t> boxtree_type;



		TriMeshSDF(const std::vector<triangle_type>& triangle_list) {
			// Copy the triangle list
			m_triangle_list.reserve(triangle_list.size());
			std::copy(triangle_list.begin(), triangle_list.end(), std::back_inserter(m_triangle_list));

			// Build the boxtree
			std::vector<std::size_t> id_list(triangle_list.size());
			std::iota(id_list.begin(), id_list.end(), 0);
			
			m_boxtree =
				boxtree_type::from_kd_construction(id_list.begin(), id_list.end(), 
					[&triangle_list](size_t i) -> box_type { return triangle_bounding_box(triangle_list[i]); }
				);
			
			// Compute individual triangle signed distance functions
			m_triangle_sdf_list.reserve(triangle_list.size());
			for(const auto& t : triangle_list)
				m_triangle_sdf_list.push_back(triangle_sdf_type(t.row(0), t.row(1), t.row(2)));

			// Compute bounding sphere
			std::vector<vector_type> vertex_list;
			vertex_list.reserve(3 * triangle_list.size());
			for(const auto& t : triangle_list)
				for(Eigen::Index i = 0; i < 3; ++i)
					vertex_list.push_back(t.row(i));

			m_bounding_sphere = point_collection::get_bounding_sphere(vertex_list.begin(), vertex_list.end());
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
						min_dist = m_triangle_sdf_list[node->value()].dist(X);
						first_hit = false;
					}
					else
						min_dist = std::min(min_dist, m_triangle_sdf_list[node->value()].dist(X));
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
		static box_type
		triangle_bounding_box(const triangle_type& t) {
			std::vector<vector_type> triplet(3);
			triplet[0] = t.row(0);
			triplet[1] = t.row(1);
			triplet[2] = t.row(2);
			return point_collection::get_bounding_box(triplet.begin(), triplet.end());
		}



		// Check whether a point is inside or outside a polygon
		template <typename InVectorT>
		inline bool
		is_inside(const InVectorT& P) const {
			KahanSum<scalar_type> winding_number;

			// For each triangle
			for(const triangle_type& t : m_triangle_list) {
				triangle_type M = t;

				M.row(0) -= P;
				M.row(1) -= P;
				M.row(2) -= P;

				scalar_type a = M.row(0).norm();
				scalar_type b = M.row(1).norm();
				scalar_type c = M.row(2).norm();

				winding_number += std::atan2(M.determinant(), (a * b * c) + c * M.row(0).dot(M.row(1)) + a * M.row(1).dot(M.row(2)) + b * M.row(2).dot(M.row(0)));
			}

			// Job done
			return winding_number() >= 2 * M_PI - 1e-7;
		}



		std::vector<triangle_type> m_triangle_list;
		std::vector<triangle_sdf_type> m_triangle_sdf_list;
		sphere_type m_bounding_sphere;
		boxtree_type m_boxtree;
	}; // class TriMeshSDF



	template <class scalar_type>
	TriMeshSDF<scalar_type>
	get_trimesh_sdf(const std::vector<Eigen::Matrix<scalar_type, 3, 3> >& triangle_list) {
		return TriMeshSDF<scalar_type>(triangle_list);
	}
} // namespace saucats



#endif // SAUCATS_SDF_TRIMESH_SDF_H

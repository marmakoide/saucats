#ifndef SAUCATS_SDF_TRIMESH_SDF_H
#define SAUCATS_SDF_TRIMESH_SDF_H

#include <saucats/geometry/bounds/PointListBounds.h>
#include <saucats/sdf/Triangle3SDF.h>



namespace saucats {
	template <class ScalarT>
	class TriMeshSDF {
	public:
		typedef ScalarT scalar_type;
		typedef Eigen::Matrix<scalar_type, 3, 1> vector_type;
		typedef Eigen::Matrix<scalar_type, 3, 3> triangle_type;
		typedef Triangle3SDF<scalar_type> triangle_sdf_type;
		typedef SphereT<vector_type> sphere_type;



		TriMeshSDF(const std::vector<triangle_type>& triangle_list) {
			// Copy the triangle list
			m_triangle_list.reserve(triangle_list.size());
			std::copy(triangle_list.begin(), triangle_list.end(), m_triangle_list.begin());

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
			// Find the distance to the closest triangle
			typename std::vector<triangle_sdf_type>::const_iterator it = m_triangle_sdf_list.begin();
			typename InVectorT::Scalar min_dist = (*it).dist(X);
			for(++it; it != m_triangle_sdf_list.end(); ++it)
				min_dist = std::min(min_dist, (*it).dist(X));

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
			scalar_type winding_number = 0;

			// For each triangle
			for(const triangle_type& t : m_triangle_list) {
				triangle_type M = t;
				//M.array().rowwise() -= P.array();
				M.row(0) -= P;
				M.row(1) -= P;
				M.row(2) -= P;

				scalar_type a = M.row(0).norm();
				scalar_type b = M.row(1).norm();
				scalar_type c = M.row(2).norm();
				scalar_type omega = M.determinant() /  ((a * b * c) + c * M.row(0).dot(M.row(1)) + a * M.row(1).dot(M.row(2)) + b * M.row(2).dot(M.row(0)));
				winding_number += std::atan(omega);
			}

			// Job done
			winding_number /= 2 * M_PI;
			return winding_number > .5;
		}


		std::vector<triangle_type> m_triangle_list;
		std::vector<triangle_sdf_type> m_triangle_sdf_list;
		sphere_type m_bounding_sphere;
	}; // class TriMeshSDF



	template <class scalar_type>
	TriMeshSDF<scalar_type>
	get_trimesh_sdf(const std::vector<Eigen::Matrix<scalar_type, 3, 3> >& triangle_list) {
		return TriMeshSDF<scalar_type>(triangle_list);
	}
} // namespace saucats



#endif // SAUCATS_SDF_TRIMESH_SDF_H

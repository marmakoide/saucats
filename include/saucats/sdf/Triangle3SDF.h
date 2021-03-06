#ifndef SAUCATS_SDF_TRIANGLE3_SDF_H
#define SAUCATS_SDF_TRIANGLE3_SDF_H

#include <vector>
#include <saucats/sdf/SegmentSDF.h>



namespace saucats {
	/*
	 Euclidean signed distance to a 3d triangle
	 */
	template <class ScalarT>
	class Triangle3SDF {
	public:
		typedef ScalarT scalar_type;
		typedef Eigen::Matrix<scalar_type, 3, 1> vector_type;
		typedef SegmentT<vector_type> segment_type;
		typedef SegmentSDF<vector_type> segment_sdf_type;
		typedef SphereT<vector_type> sphere_type;



		Triangle3SDF(const vector_type& A,
		             const vector_type& B,
		             const vector_type& C) {
			m_A = A;
			m_N = (B - A).cross(C - A).normalized();

			Eigen::Matrix<scalar_type, 3, 3> M;
			M.col(0) = B - A;
			M.col(1) = C - A;
			M.col(2) = m_N;
			Eigen::Matrix<scalar_type, 3, 3> M_inv = M.inverse();
			m_M.row(0) = M_inv.row(0);
			m_M.row(1) = M_inv.row(1);

			m_U = get_segment_sdf(segment_type::from_2_points(A, B));
			m_V = get_segment_sdf(segment_type::from_2_points(B, C));
			m_W = get_segment_sdf(segment_type::from_2_points(C, A));

			// Compute the circumpsphere
			std::vector<vector_type> triplet(3);
			triplet[0] = A;
			triplet[1] = B;
			triplet[2] = C;
			m_circumsphere = get_circumsphere(triplet.begin(), triplet.end());
		}

		template <typename InVectorT>
		inline typename InVectorT::Scalar
		dist(const InVectorT& X) const {
			// Check if X projected on the triangle plane is inside the triangle
			vector_type P = X - m_A;
			auto K = (m_M * P).eval();
			if ((K.x() + K.y() < 1) and (K.x() > 0) and (K.y() > 0))
				return std::fabs(m_N.dot(P));

			// Return distance to the closest edge
			return std::sqrt(std::min(m_U.squared_dist(X), std::min(m_V.squared_dist(X), m_W.squared_dist(X))));
		}

		inline sphere_type get_bounding_sphere() const {
			return m_circumsphere;
		}

	private:
		vector_type m_A;      // Vertices of the triangle
		vector_type m_N;      // Normal of the triangle	
		Eigen::Matrix<scalar_type, 2, 3> m_M;
		segment_sdf_type m_U; // Edges of the triangle
		segment_sdf_type m_V;
		segment_sdf_type m_W;

		sphere_type m_circumsphere; // Circumsphere of the triangle
	}; // class Triangle3SDF



	template <class scalar_type>
	Triangle3SDF<scalar_type>
	get_triangle3_sdf(const Eigen::Matrix<scalar_type, 3, 1>& A,
		                const Eigen::Matrix<scalar_type, 3, 1>& B,
		                const Eigen::Matrix<scalar_type, 3, 1>& C) {
		return Triangle3SDF<scalar_type>(A, B, C);
	}
} // namespace saucats



#endif // SAUCATS_SDF_TRIANGLE3_SDF_H

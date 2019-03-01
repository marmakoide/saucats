#ifndef SAUCATS_GEOMETRY_BOX_H
#define SAUCATS_GEOMETRY_BOX_H

#include <Eigen/Dense>



namespace saucats {
	/*
	 * The BoxT class represents a N-dimensional axis-aligned box
	 */

	template <class VectorT>
	class BoxT {
	public:
		typedef typename VectorT::Scalar scalar_type;
		typedef VectorT vector_type;

		static constexpr unsigned int CORNER_COUNT = 1 << vector_type::RowsAtCompileTime;



		BoxT(const vector_type& extent = vector_type::Zero(),
		     const vector_type& center = vector_type::Zero()) :
			m_half_extent(extent / 2),
			m_center(center) { } 

		inline BoxT(const BoxT& box) :
			m_half_extent(box.m_half_extent),
			m_center(box.m_center) { }

		inline BoxT& operator = (const BoxT& box) {
			m_half_extent = box.m_half_extent;
			m_center = box.m_center;
			return *this;
		}

		// Returns the center of the box
		inline vector_type& center() {
			return m_center;
		}

		// Returns the center of the box
		inline const vector_type& center() const {
			return m_center;
		}

		// Returns the size of the box on each dimension, divided by two
		inline const vector_type& half_extent() const {
			return m_half_extent;
		}

		// Returns the corner with the lowest coordinates
		inline vector_type min_corner() const {
			return m_center - m_half_extent;
		}

		// Returns the corner with the highest coordinates
		inline vector_type max_corner() const {
			return m_center + m_half_extent;
		}

		// Returns true if a point is inside the box
		inline bool contains(const vector_type& X) const {
			return ((X - m_center).cwiseAbs() - m_half_extent).maxCoeff() <= 0;
		}

		// Returns signed euclidean distance
		template <typename InVectorT>
		inline typename InVectorT::Scalar
		dist(const InVectorT& X) const {
			InVectorT M = ((X - m_center).cwiseAbs() - m_half_extent);
	
			typename InVectorT::Scalar max_coeff = M.maxCoeff();
			if (max_coeff < 0)
				return max_coeff;
	
			return M.cwiseMax(0).norm();
		}

		// Returns the corners of the box
		Eigen::Matrix<scalar_type, CORNER_COUNT, vector_type::RowsAtCompileTime> corners() const {
			Eigen::Matrix<scalar_type, CORNER_COUNT, vector_type::RowsAtCompileTime> ret;
			for(unsigned int i = 0; i < CORNER_COUNT; ++i)
				for(Eigen::Index j = 0; j < vector_type::RowsAtCompileTime; ++j)
					ret(i, j) = i & (1 << j) ? m_half_extent(j) : -m_half_extent(j);

			ret.rowwise() += m_center.transpose();
	
			return ret;
		}

		// Returns the volume (or surface area in 2d) of the box
		inline scalar_type volume() const {
			return std::pow(scalar_type(2), scalar_type(vector_type::RowsAtCompileTime)) * m_half_extent.prod();
		}

	private:
		vector_type m_half_extent;
		vector_type m_center;
	}; // class BoxT



	typedef BoxT<Eigen::Vector2f> Box2f;
	typedef BoxT<Eigen::Vector2d> Box2d;
	typedef BoxT<Eigen::Vector3f> Box3f;
	typedef BoxT<Eigen::Vector3d> Box3d;
} // namespace saucats



#endif // SAUCATS_GEOMETRY_BOX_H

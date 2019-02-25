#ifndef SAUCATS_RENDER_PERSPECTIVE_PROJECTION_H
#define SAUCATS_RENDER_PERSPECTIVE_PROJECTION_H

#include <Eigen/Geometry>
#include <saucats/geometry/Line.h>



namespace saucats {
	template <class ScalarT>
	class PerspectiveProjection {
	public:
		typedef ScalarT scalar_type;
		typedef Eigen::Matrix<scalar_type, 3, 1> vector_type;
		typedef Eigen::Matrix<scalar_type, 3, 3> matrix_type;
		typedef Eigen::Matrix<scalar_type, 2, 1> uv_coord_type;
		typedef LineT<vector_type> line_type;



		PerspectiveProjection(const vector_type& eye_pos = vector_type::Zero(),
		                      const matrix_type& head_rot = matrix_type::Identity(),
	                        double view_angle = .5 * M_PI) :
			m_eye_pos(eye_pos),
			m_head_rot(head_rot),
			m_focal_length(focal_length_from_angle(view_angle)) { }

		PerspectiveProjection(const PerspectiveProjection& other) :
			m_eye_pos(other.m_eye_pos),
			m_head_rot(other.m_head_rot),
			m_focal_length(other.m_focal_length) { }

		inline vector_type& eye_pos() {
			return m_eye_pos;
		}

		inline const vector_type& rot() const {
			return m_eye_pos;
		}

		inline matrix_type& head_rot() {
			return m_head_rot;
		}

		inline const matrix_type& head_rot() const {
			return m_head_rot;
		}

		line_type get_eye_line(const uv_coord_type& UV) const {
			vector_type eye_dir(UV.x() - scalar_type(.5), UV.y() - scalar_type(.5), m_focal_length);
			eye_dir = m_head_rot * eye_dir;
			return line_type(eye_dir.normalized(), m_eye_pos);
		}

	private:
		static scalar_type focal_length_from_angle(double angle) {
			return 1 / (2 * std::tan(angle / 2));
		}

		vector_type m_eye_pos;
		matrix_type m_head_rot;
		scalar_type m_focal_length;
	}; // class PerspectiveProjection



	typedef PerspectiveProjection<float> PerspectiveProjectionf;
	typedef PerspectiveProjection<double> PerspectiveProjectiond;
} // namespace saucats




#endif // SAUCATS_RENDER_PERSPECTIVE_PROJECTION_H

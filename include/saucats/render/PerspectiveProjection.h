#ifndef SAUCATS_RENDER_PERSPECTIVE_PROJECTION_H
#define SAUCATS_RENDER_PERSPECTIVE_PROJECTION_H

#include <Eigen/Geometry>
#include <saucats/geometry/Line.h>



namespace saucats {
	class PerspectiveProjection {
	public:
		PerspectiveProjection(const Eigen::Vector3d& eye_pos = Eigen::Vector3d::Zero(),
		                      const Eigen::Matrix3d& head_rot = Eigen::Matrix3d::Identity(),
	                        double view_angle = .5 * M_PI) :
			m_eye_pos(eye_pos),
			m_head_rot(head_rot),
			m_focal_length(focal_length_from_angle(view_angle)) { }

		PerspectiveProjection(const PerspectiveProjection& other) :
			m_eye_pos(other.m_eye_pos),
			m_head_rot(other.m_head_rot),
			m_focal_length(other.m_focal_length) { }

		inline Eigen::Vector3d& eye_pos() {
			return m_eye_pos;
		}

		inline const Eigen::Vector3d& rot() const {
			return m_eye_pos;
		}

		inline Eigen::Matrix3d& head_rot() {
			return m_head_rot;
		}

		inline const Eigen::Matrix3d& head_rot() const {
			return m_head_rot;
		}

		Line3d get_eye_line(const Eigen::Vector2d& UV) const {
			Eigen::Vector3d eye_dir(UV.x() - .5, UV.y() - .5, m_focal_length);
			eye_dir = m_head_rot * eye_dir;
			return Line3d(eye_dir.normalized(), m_eye_pos);
		}

	private:
		static double focal_length_from_angle(double angle) {
			return .5 / std::tan(.5 * angle);
		}

		Eigen::Vector3d m_eye_pos;
		Eigen::Matrix3d m_head_rot;
		double m_focal_length;
	}; // class PerspectiveProjection
} // namespace saucats




#endif // SAUCATS_RENDER_PERSPECTIVE_PROJECTION_H

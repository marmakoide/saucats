#ifndef SAUCATS_RENDER_COLOR_RAMP_SHADER_H
#define SAUCATS_RENDER_COLOR_RAMP_SHADER_H

#include <Eigen/Dense>



namespace saucats {
	/*
	 * Implements a shader that renders a vertical color ramp.
	 *
	 * Its intended purpose is to be used as a background for 3d renders.
	 */

	template <class color_ramp_type>
	class ColorRampShader {
	public:
		typedef typename color_ramp_type::color_type color_type;
		typedef typename color_ramp_type::scalar_type scalar_type;
		typedef Eigen::Matrix<typename color_ramp_type::scalar_type, 2, 1> uv_coord_type;



		ColorRampShader(const color_ramp_type& color_ramp,
		                scalar_type start = 1,
		                scalar_type end = 0) :
			m_start(end),
			m_scale(start - end),
			m_color_ramp(color_ramp) { }

		inline color_type
		get_color(uv_coord_type UV) const {
			return m_color_ramp.get_color(m_scale * UV.y() + m_start);
		}

	private:
		scalar_type m_start;
		scalar_type m_scale;
		color_ramp_type m_color_ramp;
	}; // class ColorRampShader



	template <class color_ramp_type>
	ColorRampShader<color_ramp_type>
	get_color_ramp_shader(const color_ramp_type& color_ramp,
	                      typename color_ramp_type::scalar_type start = 1,
	                      typename color_ramp_type::scalar_type end = 0) {
		return ColorRampShader<color_ramp_type>(color_ramp, start, end);
	}
} // namespace saucats



#endif // SAUCATS_RENDER_COLOR_RAMP_SHADER_H

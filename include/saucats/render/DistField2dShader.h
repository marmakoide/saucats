#ifndef SAUCATS_RENDER_DIST_FIELD_2D_SHADER_H
#define SAUCATS_RENDER_DIST_FIELD_2D_SHADER_H



namespace saucats {
	/*
	 * Implements a shader for 2d signed distance field. It renders each pixel
	 * as a function of the distance.
	 *
	 * When initialized, the view is automatically scaled and centered to show the 
	 * shape.
	 */

	template <class func_type, class dist_shader_type>
	class DistField2dShaderT {
	public:
		typedef Eigen::Matrix<typename func_type::scalar_type, 2, 1> uv_coord_type;
		typedef typename dist_shader_type::color_type color_type;



		DistField2dShaderT(const func_type& func,
		                   const dist_shader_type& shader) :
			m_func(func),
			m_shader(shader) {
			auto bounding_sphere = func.get_bounding_sphere();
			m_center = bounding_sphere.center();
			m_scale = Eigen::Vector2d(2.5 * bounding_sphere.radius(), -2.5 * bounding_sphere.radius());
		}

		inline uv_coord_type& center() {
			return m_center;
		}

		inline const uv_coord_type& center() const {
			return m_center;
		}

		inline uv_coord_type& scale() {
			return m_scale;
		}

		inline const uv_coord_type& scale() const {
			return m_scale;
		}

		color_type
		get_color(uv_coord_type UV) const {
			UV.array() -= .5;
			UV = m_scale.asDiagonal() * UV;
			UV += m_center;
			return m_shader.get_color(m_func.dist(UV));
		}

	private:
		uv_coord_type m_center;
		uv_coord_type m_scale;
		const func_type& m_func;
		const dist_shader_type& m_shader;
	}; // class DistField2dShaderT



	// Helper function to get an instance of a shader
	template <class func_type, class dist_shader_type>
	DistField2dShaderT<func_type, dist_shader_type>
	get_dist_field_2d_shader(const func_type& func,
	                         const dist_shader_type& shader) {
		return DistField2dShaderT<func_type, dist_shader_type>(func, shader);
	}



	/*
	 * Implements a simple distance shader using a color ramp
	 */

	template <class color_ramp_type>
	class ColorRampDistShader {
	public:
		typedef typename color_ramp_type::color_type color_type;

		ColorRampDistShader(const color_ramp_type& color_ramp) :
			m_color_ramp(color_ramp) { }

		color_type get_color(double dist) const {
			return m_color_ramp.get_color(dist);
		}

	private:
		color_ramp_type m_color_ramp;
	}; // class ColorRampDistShader



	template <class color_ramp_type>
	ColorRampDistShader<color_ramp_type>
	get_color_ramp_dist_shader(const color_ramp_type& color_ramp) {
		return ColorRampDistShader<color_ramp_type>(color_ramp);
	}
} // namespace saucats



#endif // SAUCATS_RENDER_DIST_FIELD_2D_SHADER_H

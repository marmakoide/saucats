#ifndef SAUCATS_RENDER_ISOSURFACE_3D_SHADER_H
#define SAUCATS_RENDER_ISOSURFACE_3D_SHADER_H

#include <saucats/sdf/utils/Intersection.h>



namespace saucats {
	/*
	 * Implements a shader for 3d signed distance field. It renders the
	 * isosurface.
	 *
	 * When initialized, the view is automatically scaled and centered to show the 
	 * isosurface.
	 */

	template <class func_type,
	          class fragment_shader_type,
	          class bg_shader_type,
	          class projection_type>
	class Isosurface3dShaderT {
	public:
		typedef Eigen::Matrix<typename func_type::scalar_type, 2, 1> uv_coord_type;
		typedef typename func_type::vector_type vector_type;
		typedef typename func_type::sphere_type sphere_type;
		typedef typename fragment_shader_type::color_type color_type;



		Isosurface3dShaderT(const func_type& func,
		                    const fragment_shader_type& frag_shader,
		                    const bg_shader_type& bg_shader,
		                    const projection_type& projection) :
			m_bounding_sphere(func.get_bounding_sphere()),
			m_func(func),
			m_frag_shader(frag_shader),
			m_bg_shader(bg_shader),
			m_projection(projection) {
		}

		color_type
		get_color(uv_coord_type UV) const {
			// Get the eye line
			auto eye_line = m_projection.get_eye_line(UV);
			
			// Search for an intersection
			auto intersection_it = get_intersection_iterator(m_func, eye_line, m_bounding_sphere);

			// No intersection found
			if (intersection_it.completed())
				return m_bg_shader.get_color(UV);

			// Compute the color of the fragment at the intersection
			return m_frag_shader.get_color(m_func, eye_line, m_bounding_sphere, intersection_it.get_dist());
		}

	private:
		sphere_type m_bounding_sphere;         // Bounding sphere of the isosurface
		projection_type m_projection;          // Screen projection (ie. a camera)
		const func_type& m_func;               // SDF of the isosurface
		const fragment_shader_type& m_frag_shader; // Shader used to draw isosurface
		const bg_shader_type& m_bg_shader;     // Shader used to draw background
	}; // class Isosurface3dShaderT



	// Helper function to get an instance of a shader
	template <class func_type,
	          class fragment_shader_type,
	          class bg_shader_type,
	          class projection_type>
	Isosurface3dShaderT<func_type, fragment_shader_type, bg_shader_type, projection_type>
	get_isosurface_3d_shader(const func_type& func,
	                         const fragment_shader_type& frag_shader,
		                       const bg_shader_type& bg_shader,
	                         const projection_type& projection) {
		return
			Isosurface3dShaderT<func_type, fragment_shader_type, bg_shader_type, projection_type>(
				func,
				frag_shader,
				bg_shader,
				projection);
	}
} // namespace saucats



#endif // SAUCATS_RENDER_ISOSURFACE_3D_SHADER_H

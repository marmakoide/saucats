#ifndef SAUCATS_RENDER_RENDERER_3D_T_H
#define SAUCATS_RENDER_RENDERER_3D_T_H

#include <saucats/utils/StopWatch.h>



namespace saucats {
	/*
	 * Implements the rendering of a mapping of 3d coordinates to colors (ie. a shader)
	 * with support for supersampling, and generic render target
	 *
	 * The [0, 1] unit cube is rendered.
	 */

	template <class shader_type, class target_type>
	class Renderer3dT {
	public:
		typedef typename shader_type::uvw_coord_type uvw_coord_type;



		Renderer3dT(const shader_type& shader,
		            target_type& target) :
			m_shader(shader),
			m_target(target),
			m_render_time(0) { }

		inline unsigned int render_time() const {
			return m_render_time;
		}

		void render() {
			// Initialize render
			m_target.start_render();
			precise_stopwatch stop_watch;

			// Get the render target size
			int w = m_target.width();
			int h = m_target.height(); 
			int d = m_target.depth(); 

			// Render each pixel one by one, in scanline order
			uvw_coord_type UVW;
			for(int i = 0; i < d; ++i) {
				UVW.z() = (i + .5) / d;
				for(int j = 0; j < h; ++j) {
					UVW.y() = (j + .5) / h;
					for(int k = 0; k < w; ++k) {
						UVW.x() = (k + .5) / w;

						// Compute voxel type at current position
						std::uint8_t voxel_type_id;
						if (m_shader.get_voxel_type(UVW, voxel_type_id))
							m_target.set_voxel(k, j, i, voxel_type_id);
					}
				}
			}

			// Finalize render
			m_render_time = stop_watch.elapsed_time<unsigned int, std::chrono::milliseconds>();
			m_target.end_render();
		}

	private:
		const shader_type& m_shader;
		target_type& m_target;
		unsigned int m_render_time;
	}; // class Renderer3d



	// Helper function to get an instance of a renderer
	template <class shader_type, class target_type>
	Renderer3dT<shader_type, target_type>
	get_renderer_3d(const shader_type& shader,
	                target_type& target) {
		return Renderer3dT<shader_type, target_type>(shader, target);
	}
} // namespace saucats



#endif // SAUCATS_RENDER_RENDERER_3D_T_H

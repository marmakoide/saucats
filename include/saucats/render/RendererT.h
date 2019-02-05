#ifndef SAUCATS_RENDER_RENDERER_T_H
#define SAUCATS_RENDER_RENDERER_T_H

#include <saucats/utils/StopWatch.h>



namespace saucats {
	/*
	 * Implements the rendering of a mapping of 2d coordinates to colors (ie. a shader)
	 * with support for supersampling, and generic render target
	 *
	 * The [0, 1] unit square is rendered.
	 */

	template <class shader_type, class target_type>
	class RendererT {
	public:
		typedef typename shader_type::uv_coord_type uv_coord_type;
		typedef typename shader_type::color_type color_type;



		RendererT(const shader_type& shader,
		          target_type& target,
			        int supersampling_level = 1) :
			m_shader(shader),
			m_target(target),
			m_render_time(0),
			m_supersampling_level(supersampling_level) { }

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

			// Render each pixel one by one, in scanline order
			uv_coord_type UV;
			for(int i = 0; i < h; ++i) {
				for(int j = 0; j < w; ++j) {
					// Compute pixel color
					color_type color = color_type::Zero();
					for(int is = 0; is < m_supersampling_level; ++is) {
						UV.y() = (((i + ((is + .5) / m_supersampling_level - .5)) + .5) / h);
						for(int js = 0; js < m_supersampling_level; ++js) {
							UV.x() = (((j + ((js + .5) / m_supersampling_level - .5)) + .5) / w);

							// Compute color at current position
							color += m_shader.get_color(UV);
						}
					}
					color /= m_supersampling_level * m_supersampling_level;

					// Write pixel
					m_target.set_pixel(j, i, color);
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
		int m_supersampling_level;
	}; // class Renderer



	// Helper function to get an instance of a renderer
	template <class shader_type, class target_type>
	RendererT<shader_type, target_type>
	get_renderer(const shader_type& shader, target_type& target) {
		return RendererT<shader_type, target_type>(shader, target);
	}
} // namespace saucats



#endif // SAUCATS_RENDER_RENDERER_T_H

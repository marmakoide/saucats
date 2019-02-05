#include <saucats/Geometry>
#include <saucats/Macros>
#include <saucats/Render>
#include <saucats/SDF>

using namespace saucats;



/*
 * A shader showing a shape using two colors, background and foreground
 */

template <class func_type>
class Shader {
public:
	typedef Eigen::Vector2d uv_coord_type;
	typedef Eigen::Vector3d color_type;



	Shader(const func_type& func) :
		m_func(func),
		m_bounding_sphere(func.get_bounding_sphere()),
		m_fg_color(1., 1., 1.),
		m_bg_color(0., 0., 0.) { }

	color_type
	get_color(uv_coord_type UV) const {
		if (m_bounding_sphere.contains(UV) and (m_func.dist(UV) < 0.))
			return m_fg_color;
		return m_bg_color;
	}

private:
	const func_type& m_func;
	typename func_type::sphere_type m_bounding_sphere;
	color_type m_fg_color;
	color_type m_bg_color;
}; // class Shader



template <class func_type>
Shader<func_type>
get_shader(const func_type& func) {
	return Shader<func_type>(func);
}


 
int
main(int UNUSED_PARAM(argc), char** UNUSED_PARAM(argv)) {
	// Setup the distance field
	auto sdf = get_sphere_sdf(Sphere2d(1.));

	// Setup the render target
	PNGRenderTarget render_target(256, 256, "out.png");

	// Setup the shader
	auto shader = get_shader(sdf);

	// Render the distance field
	get_renderer(shader, render_target).render();

	// Job done
	return EXIT_SUCCESS;
}

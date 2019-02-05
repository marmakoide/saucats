#include <saucats/Geometry>
#include <saucats/Macros>
#include <saucats/Render>
#include <saucats/SDF>

using namespace saucats;


 
int
main(int UNUSED_PARAM(argc), char** UNUSED_PARAM(argv)) {
	// Setup the distance field
	auto sdf = get_box_sdf(Box2d(Eigen::Vector2d(1., 1.)));

	// Setup the render target
	PNGRenderTarget render_target(256, 256, "out.png");

	// Setup the shader
	auto color_ramp = abs_color_ramp(get_color_ramp(ColorMapd::get_RdYlBu_map(), LinearInterpolation<Eigen::Vector3d>()));
	auto dist_shader = get_color_ramp_dist_shader(color_ramp);
	auto shader = get_dist_field_2d_shader(sdf, dist_shader);

	// Render the distance field
	get_renderer(shader, render_target).render();

	// Job done
	return EXIT_SUCCESS;
}

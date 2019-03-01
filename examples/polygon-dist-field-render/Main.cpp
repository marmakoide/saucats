#include <saucats/Geometry>
#include <saucats/Macros>
#include <saucats/Render>
#include <saucats/SDF>

using namespace saucats;


 
int
main(int UNUSED_PARAM(argc), char** UNUSED_PARAM(argv)) {
	// Setup the distance field
	Eigen::Matrix<double, Eigen::Dynamic, 2> vertex_array(12, 2);
	vertex_array <<
		 0., 2.,
		 1., 2.,
		 1., 1.,
		 2., 1.,
		 2., 2.,
		 3., 2.,
		 3., 1.,
		 4., 1.,
		 4., 2.,
		 5., 2.,
		 5., 0.,
     0., 0.;

	auto sdf = get_polygon2_sdf(vertex_array);

	// Setup the render target
	PNGRenderTarget render_target("out.png", 512, 512);

	// Setup the shader
	auto color_ramp =
		get_signed_color_ramp(get_color_ramp(ColorMapd::get_Blues_map()),
		                      get_color_ramp(ColorMapd::get_Oranges_map()) );
	auto dist_shader = get_color_ramp_dist_shader(color_ramp);
	auto shader = get_dist_field_2d_shader(sdf, dist_shader);

	// Render the distance field
	get_renderer(shader, render_target).render();

	// Job done
	return EXIT_SUCCESS;
}

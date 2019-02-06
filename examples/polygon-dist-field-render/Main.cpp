#include <saucats/Geometry>
#include <saucats/Macros>
#include <saucats/Render>
#include <saucats/SDF>

using namespace saucats;


 
int
main(int UNUSED_PARAM(argc), char** UNUSED_PARAM(argv)) {
	// Setup the distance field
	Eigen::Matrix<double, Eigen::Dynamic, 2> vertex_array(8, 2);
	vertex_array <<
		 1.,  1.,
		-1.,  1.,
		-2., -1.,
		-1., -1.,
		-1.,  0.,
		 1.,  0.,
		 1., -1.,
		 2., -1.;

	Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 2> edge_array(8, 2);
	edge_array <<
		0, 1,
		1, 2,
		2, 3,
		3, 4,
		4, 5,
		5, 6,
		6, 7,
		7, 0;

	auto sdf = get_polygon_sdf(vertex_array, edge_array);

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

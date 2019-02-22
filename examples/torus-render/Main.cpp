#include <saucats/Geometry>
#include <saucats/Macros>
#include <saucats/Render>
#include <saucats/SDF>

using namespace saucats;



int
main(int UNUSED_PARAM(argc), char** UNUSED_PARAM(argv)) {
	// Setup the distance field
	auto sdf = get_sdf_offset(
		get_sdf_rotation3(
			get_circle_sdf(Circle3d(1.)),
			Eigen::AngleAxisd(M_PI / 3., Eigen::Vector3d(0., -1., 0.))
		), 
	-.2);

	// Setup the render target
	PNGRenderTarget render_target("out.png", 256, 256);

	// Setup the shader
	PerspectiveProjection projection;
	Sphere3d bound_sphere = sdf.get_bounding_sphere();
	projection.eye_pos() = bound_sphere.center() + Eigen::Vector3d(0., 0., -1.5 * bound_sphere.radius());

	auto frag_shader = get_phong_fragment_shader();
	auto bg_shader = get_color_ramp_shader(get_color_ramp(ColorMapd::get_Blues_map()), 1., .5);
	auto shader = get_isosurface_3d_shader(sdf, frag_shader, bg_shader, projection);

	// Render the distance field
	get_renderer(shader, render_target).render();

	// Job done
	return EXIT_SUCCESS;
}

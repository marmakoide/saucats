#include <saucats/Geometry>
#include <saucats/Macros>
#include <saucats/Render>
#include <saucats/SDF>
#include <STLParser.h>


using namespace std;
using namespace saucats;



int
main(int UNUSED_PARAM(argc), char** UNUSED_PARAM(argv)) {
	// Load triangle mesh
	std::vector<TriMeshSDF<double>::triangle_type> triangle_list;
	STLParser<double>::load("./data/fox.stl", std::back_inserter(triangle_list));
	
	// Setup the distance field
	auto sdf = get_sdf_rotation3(
		get_trimesh_sdf(triangle_list),
		Eigen::AngleAxisd(M_PI / 2., Eigen::Vector3d(0., -1., 0.)) * Eigen::AngleAxisd(M_PI / 2., Eigen::Vector3d(1., 0., 0.))
	);

	// Setup the render target
	PNGRenderTarget render_target("out.png", 32, 32);

	// Setup the shader
	PerspectiveProjectiond projection;
	auto bound_sphere = sdf.get_bounding_sphere();
	projection.eye_pos() = bound_sphere.center() + Eigen::Vector3d(0., 0., -1.5 * bound_sphere.radius());

	PhongFragmentShaderd frag_shader;
	auto bg_shader = get_color_ramp_shader(get_color_ramp(ColorMapd::get_Blues_map()), 1., .5);
	auto shader = get_isosurface_3d_shader(sdf, frag_shader, bg_shader, projection);

	// Render the distance field
	auto renderer = get_renderer(shader, render_target);
	renderer.render();
	cout << "render time = " << renderer.render_time() << endl;

	// Job done
	return EXIT_SUCCESS;
}

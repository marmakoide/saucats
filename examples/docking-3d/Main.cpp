#include <saucats/Geometry>
#include <saucats/Macros>
#include <saucats/Render>
#include <saucats/SDF>
#include <STLParser.h>

#include <tclap/CmdLine.h>



using namespace std;
using namespace saucats;



/*
 * Define the shapes to dock together
 */

auto
get_left_sdf() {
	/*
	std::vector<TriMeshSDF<double>::triangle_type> triangle_list;
	STLParser<double>::load("./data/L-shape-container.stl", std::back_inserter(triangle_list));
	return get_trimesh_sdf(triangle_list);
	*/
	std::vector<TriMeshSDF<double>::triangle_type> triangle_list;
	STLParser<double>::load("./data/L-shape.stl", std::back_inserter(triangle_list));
	//return get_trimesh_sdf(triangle_list);
	return get_trimesh_sdf(triangle_list);	
}



auto
get_right_sdf(double angle) {
	std::vector<TriMeshSDF<double>::triangle_type> triangle_list;
	STLParser<double>::load("./data/L-shape.stl", std::back_inserter(triangle_list));
	//return get_trimesh_sdf(triangle_list);
	return get_sdf_rotation3(get_trimesh_sdf(triangle_list), Eigen::AngleAxisd((M_PI / 180.) * angle, Eigen::Vector3d::UnitZ()));	
}



/*
 * Rendering code
 */

template <class left_func_type, class right_func_type>
class Shader {
public:
	typedef Eigen::Vector3d uvw_coord_type;



	Shader(const left_func_type& left_func,
	       const right_func_type& right_func) :
		m_left_func(left_func),
		m_right_func(right_func),
		m_left_bounding_sphere(left_func.get_bounding_sphere()),
		m_right_bounding_sphere(right_func.get_bounding_sphere()) { }

	bool
	get_voxel_type(uvw_coord_type& UVW, std::uint8_t& voxel_id) const {
		Eigen::Vector3d P = UVW;
		P.array() -= .5;
		P *= 2 * (m_left_bounding_sphere.radius() + m_right_bounding_sphere.radius());
		P += m_left_bounding_sphere.center();

		bool is_left = m_left_bounding_sphere.contains(P) and (m_left_func.dist(P) < 0.);
		bool is_right = m_right_bounding_sphere.contains(P) and (m_right_func.dist(P) < 0.);
			
		if (is_left and is_right) {
			voxel_id = 3;
			return true;
		}

		if (is_left) {
			voxel_id = 1;
			return true;
		}

		if (is_right) {
			voxel_id = 2;
			return true;
		}

		return false;				
	}

private:
	const left_func_type& m_left_func;
	const right_func_type& m_right_func;

	typename left_func_type::sphere_type m_left_bounding_sphere;
	typename right_func_type::sphere_type m_right_bounding_sphere;
}; // class Shader



template <class left_func_type, class right_func_type>
Shader<left_func_type, right_func_type>
get_shader(const left_func_type& left_func,
           const right_func_type& right_func) {
	return Shader<left_func_type, right_func_type>(left_func, right_func);
}



/*
 * Main entry point
 */

void
process(double angle,
        double sampling_resolution,
        bool complementary_mode) {
	// Define the shapes to dock together
	auto left_sdf = get_left_sdf();
	auto right_sdf = get_right_sdf(angle);

	// Compute the transformation to dock 'right' shape into 'left' shape
	Eigen::Transform<double, 3, Eigen::Affine> transform = 
		registrate_3d(left_sdf, right_sdf, sampling_resolution, complementary_mode);

	// Setup the render target
	VOXRenderTarget render_target("out.vox", 128, 128, 128);

	// Setup the shader
	auto transformed_right_sdf = get_sdf_affine_transform(right_sdf, transform);
	auto shader = get_shader(left_sdf, transformed_right_sdf);

	// Render the distance field
	get_renderer_3d(shader, render_target).render();
}



int
main(int UNUSED_PARAM(argc), char** UNUSED_PARAM(argv)) {
	try {
		// Define command line
		TCLAP::CmdLine parser("3d docking demo", ' ', "1.0");

		
		TCLAP::SwitchArg complementary_arg("c", "complementary-mode", "Switch from alignment to docking mode", false);
		parser.add(complementary_arg);

		TCLAP::ValueArg<double> resolution_arg("r", "resolution", "Sampling resolution", false, .1, "positive number");
		parser.add(resolution_arg);

		TCLAP::ValueArg<double> angle_arg("a", "angle", "Rotation angle in degrees", false, 0., "positive number");
		parser.add(angle_arg);

		// Parse command line
		parser.parse(argc, argv);	

		// Process stuffs
		process(angle_arg.getValue(),
		        resolution_arg.getValue(),
		        complementary_arg.getValue());
	}
	catch (TCLAP::ArgException& e) {
		std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
	}

	// Job done
	return EXIT_SUCCESS;
}



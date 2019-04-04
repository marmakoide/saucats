#include <saucats/Geometry>
#include <saucats/Macros>
#include <saucats/Render>
#include <saucats/SDF>

#include <saucats/utils/Algorithm.h>

#include <tclap/CmdLine.h>

#include <fstream>
#include <iostream>

using namespace std;
using namespace saucats;



/*
 * Define the shapes to dock together
 */

auto
get_left_sdf() {
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

	return get_polygon2_sdf(vertex_array);
}



auto
get_right_sdf(double angle) {
	Eigen::Matrix<double, Eigen::Dynamic, 2> vertex_array(8, 2);
	vertex_array <<
		 0.,  0.,
		 1.,  0.,
		 1.,  1.,
		 2.,  1.,
		 2.,  0.,
		 3.,  0.,
		 3.,  2.,
		 0.,  2.;

	return get_sdf_rotation2(get_polygon2_sdf(vertex_array), (M_PI / 180.) * angle);	
}



/*
 * Rendering code
 */

template <class left_func_type, class right_func_type>
class Shader {
public:
	typedef Eigen::Vector2d uv_coord_type;
	typedef Eigen::Vector3d color_type;



	Shader(const left_func_type& left_func,
	       const right_func_type& right_func) :
		m_left_func(left_func),
		m_right_func(right_func),
		m_left_bounding_sphere(left_func.get_bounding_sphere()),
		m_right_bounding_sphere(right_func.get_bounding_sphere()),
		m_a_color(.552, .827, .780),
		m_b_color(1., 1., .701),
		m_c_color(.745, .729, .854),
		m_d_color(.984, .501, .447) { }

	color_type
	get_color(uv_coord_type UV) const {
		uv_coord_type P = UV;
		P.array() -= .5;
		P *= 2 * (m_left_bounding_sphere.radius() + m_right_bounding_sphere.radius());
		P += m_left_bounding_sphere.center();
		
		bool is_left = m_left_bounding_sphere.contains(P) and (m_left_func.dist(P) < 0.);
		bool is_right = m_right_bounding_sphere.contains(P) and (m_right_func.dist(P) < 0.);
			
		if (is_left and is_right)
			return m_a_color;

		if (is_left)
			return m_b_color;

		if (is_right)
			return m_c_color;

		return m_d_color;				
	}

private:
	const left_func_type& m_left_func;
	const right_func_type& m_right_func;

	typename left_func_type::sphere_type m_left_bounding_sphere;
	typename right_func_type::sphere_type m_right_bounding_sphere;

	color_type m_a_color;
	color_type m_b_color;
	color_type m_c_color;
	color_type m_d_color;
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
	Eigen::Transform<double, 2, Eigen::Affine> transform = 
		registrate_2d(left_sdf, right_sdf, sampling_resolution, complementary_mode);

	// Output transformation parameters
	{
		Eigen::Matrix2d M = transform.rotation();
		Eigen::Vector2d T = transform.translation();
		cout << "rotation = " << (180. / M_PI) * std::atan2(M(1, 0), M(0, 0)) << endl;
		cout << "translation = [" << T.x() << ", " << T.y() << "]" << endl;
	}

	// Setup the render target
	PNGRenderTarget render_target("out.png", 512, 512);

	// Setup the shader
	auto transformed_right_sdf = get_sdf_affine_transform(right_sdf, transform);
	auto shader = get_shader(left_sdf, transformed_right_sdf);

	// Render the distance field
	get_renderer(shader, render_target).render();
}



int
main(int argc, char* argv[]) {
	try {
		// Define command line
		TCLAP::CmdLine parser("2d docking demo", ' ', "1.0");

		
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


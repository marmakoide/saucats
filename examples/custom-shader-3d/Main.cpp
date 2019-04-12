#include <saucats/Geometry>
#include <saucats/Macros>
#include <saucats/Render>
#include <saucats/SDF>

using namespace saucats;




template <class func_type>
class Shader {
public:
	typedef Eigen::Vector3d uvw_coord_type;



	Shader(const func_type& func) :
		m_func(func),
		m_bounding_sphere(func.get_bounding_sphere()) { }

	bool
	get_voxel_type(uvw_coord_type& UVW, std::uint8_t& voxel_id) const {
		Eigen::Vector3d P = UVW;
		P.array() -= .5;
		P *= 2 * m_bounding_sphere.radius();
		P += m_bounding_sphere.center();

		if (m_bounding_sphere.contains(P) and (m_func.dist(P) < 0.)) {
			voxel_id = 1;
			return true;
		}

		return false;
	}

private:
	const func_type& m_func;
	typename func_type::sphere_type m_bounding_sphere;
}; // class Shader



template <class func_type>
Shader<func_type>
get_shader(const func_type& func) {
	return Shader<func_type>(func);
}



/*
 * Main entry point
 */

int
main(int UNUSED_PARAM(argc), char** UNUSED_PARAM(argv)) {
	// Setup the distance field
	auto sdf = get_sphere_sdf(Sphere3d(1.));

	// Setup the render target
	VOXRenderTarget render_target("out.vox", 16, 16, 16);

	// Setup the shader
	auto shader = get_shader(sdf);

	// Render the distance field
	get_renderer_3d(shader, render_target).render();

	// Job done
	return EXIT_SUCCESS;
}

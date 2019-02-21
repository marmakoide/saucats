#ifndef SAUCATS_RENDER_PHONG_FRAGMENT_SHADER_H
#define SAUCATS_RENDER_PHONG_FRAGMENT_SHADER_H

#include <saucats/sdf/utils/Gradient.h>
#include <saucats/sdf/utils/Intersection.h>



namespace saucats {
	class PhongFragmentShader {
	public:
		typedef Eigen::Vector3d vector_type;
		typedef Eigen::Vector3d color_type;



		PhongFragmentShader() :
			m_light_dir(vector_type(-.2, .2, 1.).normalized()) { }

		inline vector_type& light_dir() {
			return m_light_dir;
		}

		inline const vector_type& light_dir() const {
			return m_light_dir;
		}

		template <class func_type>
		double
		get_shadow(const func_type& sdf,
		           const Line3d& line,
		           const typename func_type::sphere_type& bounding_sphere,
		           double k = 2) const {
			const double tol = 1e-6;
	
			auto inter = get_sphere_line_intersection(bounding_sphere, line);
			double t_max = inter.tmax();

			double res = 1.0;
			double ph = 1e20;
			for(double t = 1e-4; t < t_max; ) {
				double h = sdf.dist(line.point(t));
				if (h < tol)
					return 0.;

				double y = h * h / (2. * ph);
				double d = std::sqrt(h * h - y * y);
				res = std::fmin(res, k * d / std::fmax(0., t - y));
				ph = h;
				t += h;
			}
	
			return res;
		}

		template <class func_type>
		color_type
		get_color(const func_type& sdf,
		          const Line3d& eye_line,
		          const typename func_type::sphere_type& bounding_sphere,
		          double dist) const {
			typedef typename func_type::vector_type vector_type;

			// Intersection position
			vector_type P = eye_line.point(dist);

			// Shadow
			double u = get_shadow(sdf, Line3d(-m_light_dir, P), bounding_sphere);

			// Normal vector
			vector_type N = get_sdf_gradient(sdf, P, 1e-6).normalized();

			// Light
			vector_type L = m_light_dir;
			double l = std::fmax(0, -N.dot(L));
			vector_type color(1., .7, 0.);
			return std::fmax(u * l, .1) * color; //m_texture.get_color(P);
		}

		private:
			vector_type m_light_dir;
	}; // class PhongFragmentShader



	PhongFragmentShader
	get_phong_fragment_shader() {
		return PhongFragmentShader();
	}
} // namespace saucats



#endif // SAUCATS_RENDER_PHONG_FRAGMENT_SHADER_H

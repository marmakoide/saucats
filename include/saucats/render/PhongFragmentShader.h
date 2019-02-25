#ifndef SAUCATS_RENDER_PHONG_FRAGMENT_SHADER_H
#define SAUCATS_RENDER_PHONG_FRAGMENT_SHADER_H

#include <saucats/sdf/utils/Gradient.h>
#include <saucats/sdf/utils/Intersection.h>



namespace saucats {
	template <class ScalarT>
	class PhongFragmentShaderT {
	public:
		typedef ScalarT scalar_type;
		typedef Eigen::Matrix<scalar_type, 3, 1> vector_type;
		typedef Eigen::Matrix<scalar_type, 3, 1> color_type;
		typedef LineT<vector_type> line_type;
		typedef SphereT<vector_type> sphere_type;



		PhongFragmentShaderT() :
			m_light_dir(vector_type(-.2, .2, 1.).normalized()) { }

		inline vector_type& light_dir() {
			return m_light_dir;
		}

		inline const vector_type& light_dir() const {
			return m_light_dir;
		}

		template <class func_type>
		scalar_type
		get_shadow(const func_type& sdf,
		           const line_type& line,
		           const sphere_type& bounding_sphere,
		           scalar_type k = 2) const {
			const scalar_type tol = 1e-6;
	
			auto inter = get_sphere_line_intersection(bounding_sphere, line);
			scalar_type t_max = inter.tmax();

			scalar_type res = 1;
			scalar_type ph = 1e20;
			for(scalar_type t = 1e-4; t < t_max; ) {
				scalar_type h = sdf.dist(line.point(t));
				if (h < tol)
					return 0.;

				scalar_type y = h * h / (2 * ph);
				scalar_type d = std::sqrt(h * h - y * y);
				res = std::fmin(res, k * d / std::fmax(scalar_type(0), t - y));
				ph = h;
				t += h;
			}
	
			return res;
		}

		template <class func_type>
		color_type
		get_color(const func_type& sdf,
		          const line_type& eye_line,
		          const sphere_type& bounding_sphere,
		          scalar_type dist) const {
			// Intersection position
			vector_type P = eye_line.point(dist);

			// Shadow
			scalar_type u = get_shadow(sdf, line_type(-m_light_dir, P), bounding_sphere);

			// Normal vector
			vector_type N = get_sdf_gradient(sdf, P, 1e-6).normalized();

			// Light
			vector_type L = m_light_dir;
			scalar_type l = std::fmax(0, -N.dot(L));
			color_type color;
			color << 1., .7, 0.;
			return std::fmax(u * l, .1) * color; //m_texture.get_color(P);
		}

		private:
			vector_type m_light_dir;
	}; // class PhongFragmentShaderT



	typedef PhongFragmentShaderT<float> PhongFragmentShaderf;
	typedef PhongFragmentShaderT<double> PhongFragmentShaderd;
} // namespace saucats



#endif // SAUCATS_RENDER_PHONG_FRAGMENT_SHADER_H

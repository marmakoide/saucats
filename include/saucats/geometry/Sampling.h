#ifndef SAUCATS_GEOMETRY_SAMPLING_H
#define SAUCATS_GEOMETRY_SAMPLING_H

#include <random>

#include <saucats/geometry/Box.h>
#include <saucats/geometry/Sphere.h>



namespace saucats {
	/*
	 *  Random uniform sampling of a box volume
	 */

	template <class VectorT>
	class BoxVolumeSampler {
	public:
		typedef VectorT vector_type;
		typedef BoxT<VectorT> box_type;



		BoxVolumeSampler(const box_type& box) :
			m_box(box),
			m_ud(-1., 1.) { }

		template <class rng_type>
		vector_type
		get_sample(rng_type& rng) {
			VectorT P;
			for(Eigen::Index i = 0; i < vector_type::RowsAtCompileTime; ++i)
				P.coeffRef(i) = m_ud(rng);

			P.array() *= m_box.half_extent().array();
			P += m_box.center();

			return P;
		}

	private:
		box_type m_box;
		std::uniform_real_distribution<> m_ud;
	}; // class BoxSampler

	template <class VectorT>
	BoxVolumeSampler<VectorT>
	get_box_volume_sampler(const BoxT<VectorT>& box) {
		return BoxVolumeSampler<VectorT>(box);
	}



	/*
	 * Random uniform sampling of a sphere surface
	 */

	// Generic implementation
	template <class VectorT>
	class SphereSurfaceSampler {
	public:
		typedef VectorT vector_type;
		typedef SphereT<VectorT> sphere_type;
		typedef typename VectorT::Scalar value_type;

		SphereSurfaceSampler(const sphere_type& sphere) :
			m_center(sphere.center()),
			m_radius(sphere.radius()) { }		

		template <class rng_type>
		vector_type
		get_sample(rng_type& rng) {
			vector_type P;
			for(Eigen::Index i = 0; i < vector_type::RowsAtCompileTime; ++i)
				P.coeffRef(i) = m_gd(rng);

			P *= m_radius / P.norm();
			P += m_center;
			return P;
		}

	private:
		vector_type m_center;
		value_type m_radius;
		std::normal_distribution<value_type> m_gd;
	}; // class SphereSurfaceSampler
	


	// Specialized 2d implementation
	template <class ValueT>
	class SphereSurfaceSampler<Eigen::Matrix<ValueT, 2, 1> > {
	public:
		typedef Eigen::Matrix<ValueT, 2, 1> vector_type;
		typedef SphereT<vector_type> sphere_type;
		typedef ValueT value_type;

		SphereSurfaceSampler(const sphere_type& sphere) :
			m_center(sphere.center()),
			m_radius(sphere.radius()),
			m_r_d(0, 1) { }

		template <class rng_type>
		vector_type
		get_sample(rng_type& rng) {
			value_type theta = 2 * M_PI * m_r_d(rng);
			vector_type P(std::cos(theta), std::sin(theta));
			P *= m_radius;
			P += m_center;
			return P;
		}

	private:
		vector_type m_center;
		value_type m_radius;
		std::uniform_real_distribution<value_type> m_r_d;
	}; // class SphereSurfaceSampler
	


	// Specialized 3d implementation
	template <class ValueT>
	class SphereSurfaceSampler<Eigen::Matrix<ValueT, 3, 1> > {
	public:
		typedef Eigen::Matrix<ValueT, 3, 1> vector_type;
		typedef SphereT<vector_type> sphere_type;
		typedef ValueT value_type;

		SphereSurfaceSampler(const sphere_type& sphere) :
			m_center(sphere.center()),
			m_radius(sphere.radius()),
			m_r_d(0, 1) { }

		template <class rng_type>
		vector_type
		get_sample(rng_type& rng) {
			value_type theta = 2 * M_PI * m_r_d(rng);
			value_type u = 2 * m_r_d(rng) - 1;
			value_type v = std::sqrt(1 - u * u);
			vector_type P(v * std::cos(theta), v * std::sin(theta), u);
			P *= m_radius;
			P += m_center;
			return P;
		}

	private:
		vector_type m_center;
		value_type m_radius;
		std::uniform_real_distribution<value_type> m_r_d;
	}; // class SphereSurfaceSampler



	template <class VectorT>
	SphereSurfaceSampler<VectorT>
	get_sphere_surface_sampler(const SphereT<VectorT>& sphere) {
		return SphereSurfaceSampler<VectorT>(sphere);
	}




	/*
	 * Random uniform sampling of a sphere volume
	 */

	template <class VectorT>
	class SphereVolumeSampler {
	public:
		typedef VectorT vector_type;
		typedef SphereT<VectorT> sphere_type;
		typedef typename VectorT::Scalar value_type;

		SphereVolumeSampler(const sphere_type& sphere) :
			m_center(sphere.center()),
			m_radius(sphere.radius()),
			m_r_d(0, 1) { }		

		template <class rng_type>
		vector_type
		get_sample(rng_type& rng) {
			vector_type P;
			for(Eigen::Index i = 0; i < vector_type::RowsAtCompileTime; ++i)
				P.coeffRef(i) = m_gd(rng);

			P *= (m_radius * std::sqrt(m_r_d(rng))) / P.norm();
			P += m_center;
			return P;
		}

	private:
		vector_type m_center;
		value_type m_radius;
		std::normal_distribution<value_type> m_gd;
		std::uniform_real_distribution<value_type> m_r_d;
	}; // class SphereSurfaceVolume



	// Specialized 2d implementation
	template <class ValueT>
	class SphereVolumeSampler<Eigen::Matrix<ValueT, 2, 1> > {
	public:
		typedef Eigen::Matrix<ValueT, 2, 1> vector_type;
		typedef SphereT<vector_type> sphere_type;
		typedef ValueT value_type;

		SphereVolumeSampler(const sphere_type& sphere) :
			m_center(sphere.center()),
			m_radius(sphere.radius()),
			m_r_d(0, 1) { }

		template <class rng_type>
		vector_type
		get_sample(rng_type& rng) {
			value_type theta = (2 * M_PI) * m_r_d(rng);
			vector_type P(std::cos(theta), std::sin(theta));
			P *= m_radius * std::sqrt(m_r_d(rng));
			P += m_center;
			return P;
		}

	private:
		vector_type m_center;
		value_type m_radius;
		std::uniform_real_distribution<value_type> m_r_d;
	}; // class SphereSurfaceSampler



	// Specialized 3d implementation
	template <class ValueT>
	class SphereVolumeSampler<Eigen::Matrix<ValueT, 3, 1> > {
	public:
		typedef Eigen::Matrix<ValueT, 3, 1> vector_type;
		typedef SphereT<vector_type> sphere_type;
		typedef ValueT value_type;

		SphereVolumeSampler(const sphere_type& sphere) :
			m_center(sphere.center()),
			m_radius(sphere.radius()),
			m_r_d(0, 1) { }

		template <class rng_type>
		vector_type
		get_sample(rng_type& rng) {
			value_type theta = 2 * M_PI * m_r_d(rng);
			value_type u = 2 * m_r_d(rng) - 1;
			value_type v = std::sqrt(1 - u * u);
			vector_type P(v * std::cos(theta), v * std::sin(theta), u);
			P *= m_radius * std::sqrt(m_r_d(rng));
			P += m_center;
			return P;
		}

	private:
		vector_type m_center;
		value_type m_radius;
		std::uniform_real_distribution<value_type> m_r_d;
	}; // class SphereVolumeSampler



	template <class VectorT>
	SphereVolumeSampler<VectorT>
	get_sphere_volume_sampler(const SphereT<VectorT>& sphere) {
		return SphereVolumeSampler<VectorT>(sphere);
	}



	/*
	 * Approximatively regular uniform sampling of a sphere surface
	 */

	template <class ValueT>
	Eigen::Matrix<ValueT, Eigen::Dynamic, 3>
	get_spiral_sphere(Eigen::Index point_count) {
		typedef Eigen::Array<ValueT, Eigen::Dynamic, 1> array_1d_type;

		const ValueT golden_angle = M_PI * (3. - std::sqrt(5.));

		array_1d_type theta = array_1d_type::LinSpaced(point_count, 0., golden_angle * (point_count - 1));
		array_1d_type Z = array_1d_type::LinSpaced(point_count, 1. - (1. / point_count), (1. / point_count) - 1.);
		array_1d_type R = (array_1d_type::Ones(point_count) - Z.square()).sqrt();

		auto ret = Eigen::Matrix<ValueT, Eigen::Dynamic, 3>(point_count, 3);
		ret.col(0) = R * theta.cos();
		ret.col(1) = R * theta.sin();
		ret.col(2) = Z;
		return ret;
	}
} // namespace saucats



#endif // SAUCATS_GEOMETRY_SAMPLING_H

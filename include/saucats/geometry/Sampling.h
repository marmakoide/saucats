#ifndef SAUCATS_GEOMETRY_SAMPLING_H
#define SAUCATS_GEOMETRY_SAMPLING_H

#include <random>

#include <saucats/geometry/Box.h>
#include <saucats/geometry/Sphere.h>



namespace saucats {
	// Uniform sampling of a box volume
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
	 * Uniform sampling of a sphere surface
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
			m_theta_d(0, 2 * M_PI) { }

		template <class rng_type>
		vector_type
		get_sample(rng_type& rng) {
			value_type theta = m_theta_d(rng);
			vector_type P(std::cos(theta), std::sin(theta));
			P *= m_radius;
			P += m_center;
			return P;
		}

	private:
		vector_type m_center;
		value_type m_radius;
		std::uniform_real_distribution<> m_theta_d;
	}; // class SphereSurfaceSampler
	


	// Specialized 2d implementation
	template <class ValueT>
	class SphereSurfaceSampler<Eigen::Matrix<ValueT, 3, 1> > {
	public:
		typedef Eigen::Matrix<ValueT, 3, 1> vector_type;
		typedef SphereT<vector_type> sphere_type;
		typedef ValueT value_type;

		SphereSurfaceSampler(const sphere_type& sphere) :
			m_center(sphere.center()),
			m_radius(sphere.radius()),
			m_theta_d(0, 2 * M_PI),
			m_u_d(-1, 1) { }

		template <class rng_type>
		vector_type
		get_sample(rng_type& rng) {
			value_type theta = m_theta_d(rng);
			value_type u = m_u_d(rng);
			value_type v = std::sqrt(1 - u * u);
			vector_type P(v * std::cos(theta), v * std::sin(theta), u);
			P *= m_radius;
			P += m_center;
			return P;
		}

	private:
		vector_type m_center;
		value_type m_radius;
		std::uniform_real_distribution<> m_theta_d;
		std::uniform_real_distribution<> m_u_d;
	}; // class SphereSurfaceSampler



	template <class VectorT>
	SphereSurfaceSampler<VectorT>
	get_sphere_surface_sampler(const SphereT<VectorT>& sphere) {
		return SphereSurfaceSampler<VectorT>(sphere);
	}
} // namespace saucats



#endif // SAUCATS_GEOMETRY_SAMPLING_H

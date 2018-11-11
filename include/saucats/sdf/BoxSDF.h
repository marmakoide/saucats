#ifndef SAUCATS_SDF_BOX_SDF_H
#define SAUCATS_SDF_BOX_SDF_H

#include <saucats/geometry/Box.h>
#include <saucats/geometry/Sphere.h>



namespace saucats {
	/*
	 Euclidean signed distance to a box
	 */

	template <class VectorT>
	class BoxSDF {
	public:
		typedef typename VectorT::Scalar scalar_type;
		typedef VectorT vector_type;
		typedef BoxT<VectorT> box_type;
		typedef SphereT<VectorT> sphere_type;



		inline BoxSDF(const box_type& box)
			: m_box(box) { }

		template <typename InVectorT>
		inline typename InVectorT::Scalar
		dist(const InVectorT& X) const {
			vector_type M = ((X - m_box.center()).cwiseAbs() - m_box.half_extent());
	
			scalar_type max_coeff = M.maxCoeff();
			if (max_coeff < 0)
				return max_coeff;
	
			return M.cwiseMax(0).norm();
		}

		inline sphere_type get_bounding_sphere() const {
			return sphere_type(m_box.half_extent().norm(), m_box.center());
		}

	private:
		box_type m_box;
	}; // classBoxSDF



	template <class VectorT>
	BoxSDF<VectorT>
	get_box_sdf(const BoxT<VectorT>& box) {
		return BoxSDF<VectorT>(box);
	}
} // namespace saucats



#endif // SAUCATS_SDF_BOX_SDF_H

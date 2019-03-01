#ifndef SAUCATS_SDF_BOX_SDF_H
#define SAUCATS_SDF_BOX_SDF_H

#include <saucats/geometry/bounds/BoxBounds.h>



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
			return m_box.dist(X);
		}

		inline sphere_type get_bounding_sphere() const {
			return saucats::get_bounding_sphere(m_box);
		}

	private:
		box_type m_box;
	}; // class BoxSDF



	template <class VectorT>
	BoxSDF<VectorT>
	get_box_sdf(const BoxT<VectorT>& box) {
		return BoxSDF<VectorT>(box);
	}
} // namespace saucats



#endif // SAUCATS_SDF_BOX_SDF_H

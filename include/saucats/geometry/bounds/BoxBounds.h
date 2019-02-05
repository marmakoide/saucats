#ifndef SAUCATS_GEOMETRY_BOUNDS_BOX_BOUNDS_H
#define SAUCATS_GEOMETRY_BOUNDS_BOX_BOUNDS_H

#include <saucats/geometry/Box.h>
#include <saucats/geometry/Sphere.h>



namespace saucats {
	// Bounding box of a box (yeah... for completleness)
	template <class VectorT>
	BoxT<VectorT>
	get_bounding_box(const BoxT<VectorT>& box) {
		return BoxT<VectorT>(box);
	}

	// Bounding sphere of a box
	template <class VectorT>
	SphereT<VectorT>
	get_bounding_sphere(const BoxT<VectorT>& box) {
		return SphereT<VectorT>(box.center(), box.half_extent().norm());
	}
} // namespace saucats



#endif // SAUCATS_GEOMETRY_BOUNDS_BOX_BOUNDS_H

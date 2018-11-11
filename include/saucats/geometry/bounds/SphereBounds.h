#ifndef SAUCATS_GEOMETRY_BOUNDS_SPHERE_BOUNDS_H
#define SAUCATS_GEOMETRY_BOUNDS_SPHERE_BOUNDS_H

#include <saucats/geometry/Box.h>
#include <saucats/geometry/Sphere.h>



namespace saucats {
	// Bounding box of a sphere
	template <class VectorT>
	BoxT<VectorT>
	get_bounding_box(const SphereT<VectorT>& sphere) {
		return BoxT<VectorT>(VectorT::Constant(sphere.radius()), sphere.center());
	}

	// Bounding sphere of a sphere (yeah... for completleness)
	template <class VectorT>
	SphereT<VectorT>
	get_bounding_sphere(const SphereT<VectorT>& sphere) {
		return SphereT<VectorT>(sphere);
	}
} // namespace saucats



#endif // SAUCATS_GEOMETRY_BOUNDS_SPHERE_BOUNDS_H

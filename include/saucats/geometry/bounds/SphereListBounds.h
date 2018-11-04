#ifndef SAUCATS_GEOMETRY_BOUNDS_SPHERE_LIST_BOUNDS_H
#define SAUCATS_GEOMETRY_BOUNDS_SPHERE_LIST_BOUNDS_H

#include <saucats/geometry/Sphere.h>



namespace saucats {
	/*
	 * Compute the minimum volume axis-aligned box bounding of a collection of spheres
	 */

	template <class sphere_collection_type>
	BoxT<typename sphere_collection_type::value_type::vector_type>
	get_bounding_box(const sphere_collection_type& sphere_collection) {
		typedef typename sphere_collection_type::value_type sphere_type;
		typedef typename sphere_type::vector_type vector_type;

		vector_type min_corner, max_corner;

		bool first_item = true;
		for(const sphere_type& sphere : sphere_collection) {
			vector_type lo = sphere.center().array() - sphere.radius();
			vector_type hi = sphere.center().array() + sphere.radius();

			if (first_item) {
				min_corner = lo;
				max_corner = hi;
				first_item = false;
			}
			else {
				min_corner = min_corner.cwiseMin(lo);
				max_corner = max_corner.cwiseMax(hi);
			}
		}

		return
			BoxT<vector_type>(max_corner - min_corner,
			                  (min_corner + max_corner) / 2);
	}
} // namespace saucats



#endif // SAUCATS_GEOMETRY_BOUNDS_SPHERE_LIST_BOUNDS_H

#ifndef SAUCATS_GEOMETRY_BOUNDS_BOX_LIST_BOUNDS_H
#define SAUCATS_GEOMETRY_BOUNDS_BOX_LIST_BOUNDS_H

#include <saucats/geometry/Box.h>



namespace saucats {
	namespace box_collection {
		/*
		 * Compute the minimum volume axis-aligned box bounding of a collection of boxes
		 */

		template <class iterator_type>
		BoxT<typename iterator_type::value_type::vector_type>
		get_bounding_box(iterator_type start_it,
		                 iterator_type end_it) {
			typedef typename iterator_type::value_type box_type;
			typedef typename box_type::vector_type vector_type;
	
			vector_type min_corner = start_it->min_corner();
			vector_type max_corner = start_it->max_corner();
			
			for(++start_it; start_it != end_it; ++start_it) {
				min_corner = min_corner.cwiseMin(start_it->min_corner());
				max_corner = max_corner.cwiseMax(start_it->max_corner());
			}
			
			return
				BoxT<vector_type>(max_corner - min_corner,
				                  (min_corner + max_corner) / 2);
		}
	} // namespace box_collection
} // namespace saucats



#endif // SAUCATS_GEOMETRY_BOUNDS_BOX_LIST_BOUNDS_H

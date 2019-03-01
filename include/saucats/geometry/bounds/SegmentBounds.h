#ifndef SAUCATS_GEOMETRY_BOUNDS_SEGMENT_BOUNDS_H
#define SAUCATS_GEOMETRY_BOUNDS_SEGMENT_BOUNDS_H

#include <saucats/geometry/Box.h>
#include <saucats/geometry/Segment.h>
#include <saucats/geometry/Sphere.h>



namespace saucats {
	// Bounding box of a segment
	template <class VectorT>
	BoxT<VectorT>
	get_bounding_box(const SegmentT<VectorT>& segment) {
		VectorT extent = (2 * segment.half_length()) * segment.line().direction().cwiseAbs();
		return BoxT<VectorT>(extent, segment.line().origin());
	}

	// Bounding sphere of a segment
	template <class VectorT>
	SphereT<VectorT>
	get_bounding_sphere(const SegmentT<VectorT>& segment) {
		return sphere_type(segment.line().origin(), segment.half_length());
	}
} // namespace saucats



#endif // SAUCATS_GEOMETRY_BOUNDS_SEGMENT_BOUNDS_H

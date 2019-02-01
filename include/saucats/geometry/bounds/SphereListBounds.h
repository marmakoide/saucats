#ifndef SAUCATS_GEOMETRY_BOUNDS_SPHERE_LIST_BOUNDS_H
#define SAUCATS_GEOMETRY_BOUNDS_SPHERE_LIST_BOUNDS_H

#include <vector>
#include <saucats/geometry/Sphere.h>



namespace saucats {
	namespace sphere_collection {
		/*
		 * Compute the minimum volume axis-aligned box bounding of a collection of spheres
		 */

		template <class iterator_type>
		BoxT<typename iterator_type::value_type::vector_type>
		get_bounding_box(iterator_type start_it,
		                 iterator_type end_it) {
			typedef typename iterator_type::value_type sphere_type;
			typedef typename sphere_type::vector_type vector_type;
	
			vector_type min_corner = *start_it;
			vector_type max_corner = *start_it;
			
			for(++start_it; start_it != end_it; ++start_it) {
				const sphere_type& sphere = *start_it;
				min_corner = min_corner.cwiseMin(sphere.center().array() - sphere.radius());
				max_corner = max_corner.cwiseMax(sphere.center().array() + sphere.radius());
			}

			return
				BoxT<vector_type>(max_corner - min_corner,
				                  (min_corner + max_corner) / 2);
		}



		/*
		 * Compute a minimum volume bounding sphere for a collection of spheres
		 */

		template <class sphere_collection_type>
		class SphereCollectionBoundingBallT {
		public:
			typedef typename sphere_collection_type::value_type sphere_type;
			typedef typename sphere_type::vector_type vector_type;
			typedef typename std::vector<vector_type> vector_collection_type;
			typedef typename sphere_type::scalar_type scalar_type;



			struct ball_type {
				inline ball_type(const vector_type& center, scalar_type radius) :
					m_center(center),
					m_radius(radius) { }

				vector_type m_center;
				scalar_type m_radius;
			}; // struct ball_type



			static sphere_type
			get_bounding_sphere(const sphere_collection_type& sphere_collection,
			                    scalar_type epsilon) {
				if (sphere_collection.empty())
					return sphere_type();

				if (sphere_collection.size() == 1)
					return sphere_type(sphere_collection.front());

				ball_type qp_hp = find_farthest_point(sphere_collection.front().center() , sphere_collection);
				ball_type q_h = find_farthest_point(qp_hp.m_center, sphere_collection);
				ball_type c_r = ball_type((qp_hp.m_center + q_h.m_center) / 2, q_h.m_radius / 2);
		
				vector_collection_type C;
				C.push_back(qp_hp.m_center);
				C.push_back(q_h.m_center);

				unsigned int epoch_count = (unsigned int)std::ceil(2. / epsilon);
				for(unsigned int i = 0; i < epoch_count; ++i) {
					q_h = find_farthest_point(c_r.m_center, sphere_collection);
					if (q_h.m_radius < c_r.m_radius * (1. + epsilon))
						break;

					c_r.m_radius = ((c_r.m_radius * c_r.m_radius) / q_h.m_radius + q_h.m_radius) / 2;
					c_r.m_center = q_h.m_center + (c_r.m_radius / q_h.m_radius) * (c_r.m_center - q_h.m_center);

					C.push_back(q_h.m_center);
					c_r = solve_ball(c_r, C, epsilon / 2);
				}
	
				return sphere_type(c_r.m_center, q_h.m_radius);
			}

		private:
			static ball_type
			find_farthest_point(const vector_type& c, const vector_collection_type& P) {
				bool first = true;
				scalar_type max_dist_sqr = 0;
				vector_type max_point;

				for(const vector_type& p : P) {
					scalar_type dist_sqr = (p - c).squaredNorm();
					if (first or (max_dist_sqr < dist_sqr)) {
						first = false;
						max_dist_sqr = dist_sqr;
						max_point = p;
					}
				}
				return ball_type(max_point, std::sqrt(max_dist_sqr));
			}

			static ball_type
			find_farthest_point(const vector_type& c, const sphere_collection_type& sphere_collection) {
				bool first = true;
				scalar_type max_dist = 0;
				sphere_type max_sphere;

				for(const sphere_type& sphere : sphere_collection) {
					scalar_type dist = (sphere.center() - c).norm() + sphere.radius();
					if (first or (max_dist < dist)) {
						first = false;
						max_dist = dist;
						max_sphere = sphere;
					}
				}
	
				vector_type q = max_sphere.center() + (max_sphere.radius() / (max_dist - max_sphere.radius())) * (max_sphere.center() - c);
				return ball_type(q, max_dist);
			}

			static ball_type
			solve_ball(const ball_type& B, const vector_collection_type& P, scalar_type delta) {
				ball_type ret = B;

				unsigned int epoch_count = (unsigned int)std::ceil(2 / delta);
				for(unsigned int i = 0; i < epoch_count; ++i) {
					ball_type q_h = find_farthest_point(ret.m_center, P);
					if (q_h.m_radius < ret.m_radius * (1 + delta))
						break;

					ret.m_radius = ((ret.m_radius * ret.m_radius) / q_h.m_radius + q_h.m_radius) / 2;
					ret.m_center = q_h.m_center + (ret.m_radius / q_h.m_radius) * (ret.m_center - q_h.m_center);
				}

				return ret;
			}
		}; // class SphereCollectionBoundingBallT



		template <class sphere_collection_type>
		static SphereT<typename sphere_collection_type::value_type::vector_type>
		get_bounding_sphere(const sphere_collection_type& sphere_collection,
  		                  typename sphere_collection_type::value_type::scalar_type epsilon) {
			return SphereCollectionBoundingBallT<sphere_collection_type>().get_bounding_sphere(sphere_collection, epsilon);
		}
	} // namespace sphere_collection
} // namespace saucats



#endif // SAUCATS_GEOMETRY_BOUNDS_SPHERE_LIST_BOUNDS_H

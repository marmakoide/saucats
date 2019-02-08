#ifndef SAUCATS_SDF_UTILS_INTERSECTION_H
#define SAUCATS_SDF_UTILS_INTERSECTION_H

#include <saucats/geometry/intersections/SphereLineIntersection.h>
#include <saucats/utils/KahanSum.h>



namespace saucats {
	/*
	 * Enumerate the intersections between the isosurface f = 0 of a signed 
	 * distance function and a line, assuming we have a bounding sphere for the 
	 * isosurface
	 */

	template <class FuncT>
	class IntersectionIterator {
	public:
		typedef FuncT func_type;
		typedef typename FuncT::scalar_type scalar_type;
		typedef typename FuncT::vector_type vector_type;
		typedef typename FuncT::sphere_type sphere_type;
		typedef LineT<vector_type> line_type;
		typedef KahanSum<scalar_type> sum_type;


		IntersectionIterator(const func_type& sdf,
		                     const line_type& line,
		                     const sphere_type& bounding_sphere,
		                     std::size_t max_step_count = 1024,
  	                     scalar_type tol = 1e-6) :
			m_sdf(sdf),
			m_line(line),
			m_completed(false),
			m_inside_bound(false),
			m_max_dist(0),
			m_step_count(0),
			m_max_step_count(max_step_count),
			m_tol(tol) {

			if (!setup(bounding_sphere))
				m_completed = true;
			else
				m_completed = !march();
		}

		inline void next() {
			m_dist_sum += m_tol;
			m_completed = !march();
		}

		inline bool completed() const {
			return m_completed;
		}

		inline const line_type& line() const {
			return m_line;
		}

		inline scalar_type get_dist() const {
			return m_dist_sum();
		}

	private:
		bool setup(const sphere_type& bounding_sphere) {
			// Check if the line intersects the bounding sphere
			auto clint = get_sphere_line_intersection(bounding_sphere, m_line);
			if ((!clint.is_valid()) || (clint.tmax() < 0.))
				return false;

			// Setup the initial and maximum marching distance
			m_max_dist = clint.tmax();
			m_dist_sum = sum_type(std::fmax(0., clint.tmin()));
		
			// Job done
			return true;
		}

		bool march() {
			for(; (m_step_count < m_max_step_count) and (m_dist_sum() < m_max_dist); ++m_step_count) {
				scalar_type dist = std::fabs(m_sdf.dist(m_line.point(m_dist_sum())));

				if (m_inside_bound) {
					if (dist > m_tol)
						m_inside_bound = false;
					else
						dist = 2 * m_tol;
				}
				else if ((!m_inside_bound) and (dist < m_tol)) {
					m_inside_bound = true;
					return true;
				}

				m_dist_sum += dist;
			}
			return false;
		}

		const func_type& m_sdf;
		const line_type& m_line;
		bool m_completed;
		bool m_inside_bound;
		scalar_type m_max_dist;
		sum_type m_dist_sum;
		std::size_t m_step_count;
		std::size_t m_max_step_count;
		scalar_type m_tol;
	}; // class IntersectionIterator



	template <class FuncT>
	IntersectionIterator<FuncT>
	get_intersection_iterator(const FuncT& sdf,
  	                        const LineT<typename FuncT::vector_type>& line,
  	                        const typename FuncT::sphere_type& bounding_sphere) {
		return IntersectionIterator<FuncT>(sdf, line, bounding_sphere);
	}
} // namespace saucats



#endif // SAUCATS_SDF_UTILS_INTERSECTION_H

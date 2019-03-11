#ifndef SAUCATS_UTILS_ALGORITHM_H
#define SAUCATS_UTILS_ALGORITHM_H

#include <algorithm>
#include <random>



namespace saucats {
	/*
	 * Return an iterator on an element of a collection
	 */

	template<typename iterator_type, typename rng_type>
	iterator_type
	select_randomly(iterator_type start, iterator_type end, rng_type& rng) {
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(rng));
    return start;
	}



	/*
	 * 180 degrees rotation of a matrix
	 */

	template <class matrix_type>
	void
	matrix_rotate_180(matrix_type& U) {
		U.rowwise().reverseInPlace();
		U.colwise().reverseInPlace();
	}	
} // namespace saucats



#endif // SAUCATS_UTILS_ALGORITHM_H

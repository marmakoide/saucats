#ifndef SAUCATS_RENDER_COLOR_MAP_H
#define SAUCATS_RENDER_COLOR_MAP_H

#include <Eigen/Dense>



namespace saucats {
	/*
	 * Defines color maps
	 */

	template <class ScalarT>
	struct ColorMap {
		typedef ScalarT scalar_type;
		typedef Eigen::Matrix<ScalarT, Eigen::Dynamic, 3> color_array_type;



		static color_array_type&
		get_RdYlBu_map() {
			static ColorMap::color_array_type ret = 
			(Eigen::Matrix<scalar_type, 11, 3>() <<
				0.647, 0.0,   0.149,
				0.843, 0.188, 0.153,
				0.957, 0.427, 0.263,
				0.992, 0.682, 0.380,
				0.996, 0.878, 0.565,
				1.,    1.,    0.749,
				0.878, 0.953, 0.972,
				0.670, 0.850, 0.914,
				0.455, 0.678, 0.819,
				0.270, 0.459, 0.706,
				0.192, 0.212, 0.584
			).finished();

			return ret;
		}
	}; // struct ColorMap



	typedef ColorMap<float> ColorMapf;
} // namespace saucats



#endif // SAUCATS_RENDER_COLOR_MAP_H


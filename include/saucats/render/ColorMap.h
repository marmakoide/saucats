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
			static color_array_type ret = 
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

		static color_array_type&
		get_PiYG_map() {
			static color_array_type ret = 
			(Eigen::Matrix<scalar_type, 11, 3>() <<
				0.556, 0.004, 0.321,
				0.772, 0.105, 0.490,
				0.870, 0.466, 0.682,
				0.945, 0.713, 0.854,
				0.992, 0.878, 0.937,
				0.968, 0.968, 0.968,
				0.901, 0.960, 0.815,
				0.721, 0.882, 0.525,
				0.498, 0.737, 0.254,
				0.301, 0.572, 0.129,
				0.152, 0.392, 0.098
			).finished();

			return ret;
		}

		static color_array_type&
		get_Blues_map() {
			static color_array_type ret = 
			(Eigen::Matrix<scalar_type, 9, 3>() <<
				0.968, 0.984, 1.,
				0.870, 0.921, 0.968,
				0.776, 0.858, 0.937,
				0.619, 0.792, 0.882,
				0.419, 0.682, 0.839,
				0.258, 0.572, 0.776,
				0.129, 0.443, 0.709,
				0.031, 0.317, 0.611,
				0.031, 0.188, 0.419
			).finished();

			return ret;
		}

		static color_array_type&
		get_Oranges_map() {
			static color_array_type ret = 
			(Eigen::Matrix<scalar_type, 9, 3>() <<
				1.,    0.960, 0.921,
				0.996, 0.901, 0.807,
				0.992, 0.815, 0.635,
				0.992, 0.682, 0.419,
				0.992, 0.552, 0.235,
				0.945, 0.411, 0.074,
				0.850, 0.282, 0.004,
				0.650, 0.211, 0.011,
				0.498, 0.152, 0.015
			).finished();

			return ret;
		}
	}; // struct ColorMap



	typedef ColorMap<float> ColorMapf;
	typedef ColorMap<double> ColorMapd;
} // namespace saucats



#endif // SAUCATS_RENDER_COLOR_MAP_H


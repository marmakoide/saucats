#include <saucats/Geometry>
#include <saucats/Macros>
#include <saucats/SDF>
#include <saucats/utils/Algorithm.h>

#include <fstream>
#include <iostream>

using namespace std;
using namespace saucats;



// --- Define the shapes to dock together -------------------------------------

auto
get_left_sdf() {
	Eigen::Matrix<double, Eigen::Dynamic, 2> vertex_array(12, 2);
	vertex_array <<
		 0., 2.,
		 1., 2.,
		 1., 1.,
		 2., 1.,
		 2., 2.,
		 3., 2.,
		 3., 1.,
		 4., 1.,
		 4., 2.,
		 5., 2.,
		 5., 0.,
     0., 0.;

	return get_polygon2_sdf(vertex_array);
}



auto
get_right_sdf() {
	Eigen::Matrix<double, Eigen::Dynamic, 2> vertex_array(8, 2);
	vertex_array <<
		 0.,  0.,
		 1.,  0.,
		 1.,  1.,
		 2.,  1.,
		 2.,  0.,
		 3.,  0.,
		 3.,  2.,
		 0.,  2.;

	return get_sdf_rotation2(get_polygon2_sdf(vertex_array), (M_PI / 180.) * 42.);	
}



// --- Rendering code ---------------------------------------------------------

template <class left_func_type, class right_func_type>
class DebugFunctor {
public:
	typedef typename left_func_type::scalar_type scalar_type;
	typedef typename left_func_type::vector_type vector_type;



	DebugFunctor(const left_func_type& left_func,
	             const right_func_type& right_func) :
		m_left_func(left_func),
		m_right_func(right_func) { }

	inline scalar_type
	operator() (const vector_type& P) const {
		bool a = m_left_func.dist(P) < 0;
		bool b = m_right_func.dist(P) < 0;
											
		if (a) {
			if (b)
				return 2.;
			else
				return 1.5;
		}
		else {
			if (b)
				return .5;
			else
				return 0.;
			}
	} 
	
private:
	const left_func_type& m_left_func;
	const right_func_type& m_right_func;
}; // class DebugFunctor

 

template <class left_func_type, class right_func_type>
DebugFunctor<left_func_type, right_func_type>
get_debug_functor(const left_func_type& left_func,
                  const right_func_type& right_func) {
	return DebugFunctor<left_func_type, right_func_type>(left_func, right_func);
}



template <class matrix_type>
void
output_matrix(matrix_type& matrix, std::ostream& out) {
	for(Eigen::Index i = 0; i < matrix.rows(); ++i) {
		for(Eigen::Index j = 0; j < matrix.cols(); ++j) {
			if (j > 0)
				out << " ";
			out << matrix(i, j);
		}
		out << endl;
	}
}



// --- Main entry point -------------------------------------------------------

int
main(int UNUSED_PARAM(argc), char** UNUSED_PARAM(argv)) {
	// Define the shapes to dock together
	auto left_sdf = get_left_sdf();
	auto right_sdf = get_right_sdf();

	// Compute the transformation to dock 'right' shape into 'left' shape
	Eigen::Transform<double, 2, Eigen::Affine> transform = registrate(left_sdf, right_sdf, .1, true);

	// Output the result
	auto left_bounding_sphere = left_sdf.get_bounding_sphere();
	auto right_bounding_sphere = right_sdf.get_bounding_sphere();

	Eigen::MatrixXd Z(256, 256);
	Box2d box_domain((2 * (left_bounding_sphere.radius() + right_bounding_sphere.radius())) * Eigen::Vector2d::Ones(), left_bounding_sphere.center());
	
	sample_function_2d(box_domain, Z, get_debug_functor(left_sdf, get_sdf_affine_transform(right_sdf, transform)));	
	std::ofstream out_file("out.txt");
	output_matrix(Z, out_file);	

	// Job done
	return EXIT_SUCCESS;
}


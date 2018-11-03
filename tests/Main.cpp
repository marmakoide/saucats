#include <minunit.h>

#include <saucats/Geometry>
#include <saucats/Macros>

#include <iostream>
using namespace std;
using namespace saucats;



// --- Utilities --------------------------------------------------------------

template <class Scalar>
Scalar
abs_distance(Scalar a, Scalar b) {
	return std::fabs(a - b);
}



// --- Sphere class checks ----------------------------------------------------

template <class VectorT>
int
sphere_distance_consistency_check() {
	SphereT<VectorT> S(1.);

	const unsigned int step_count = 16;
	for(unsigned int i = 0; i < step_count; ++i) {
		VectorT A = VectorT::Zero();
		A(0) = (2. * i) / step_count;
		mu_assert(abs_distance(std::pow(S.signed_dist(A), (typename VectorT::Scalar)2), S.squared_dist(A)) <= std::numeric_limits<typename VectorT::Scalar>::epsilon(), "inconsistent values returned for SphereT::squared_dist and SphereT::signed_dist");
	}

	// Job done
	return 0;
}



// --- Entry point ------------------------------------------------------------

int
all_tests() {	
	// Sphere class tests
	mu_run_test(sphere_distance_consistency_check<Eigen::Vector2f>);
	mu_run_test(sphere_distance_consistency_check<Eigen::Vector2d>);
	mu_run_test(sphere_distance_consistency_check<Eigen::Vector3f>);
	mu_run_test(sphere_distance_consistency_check<Eigen::Vector3d>);

	// Job done
	return 0;
}



int
main(int UNUSED_PARAM(argc), char** UNUSED_PARAM(argv)) {
	// Enforce a deterministic run
	srand(1953);

	// Run the tests
	int failed = all_tests();

	// Print diagnostic
	if (!failed)
		printf("ALL TESTS PASSED\n");

	printf("Tests run: %u\n", mu_tests_run);

	// Job done
	return failed;
}



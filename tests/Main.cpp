#include <minunit.h>

#include <saucats/Geometry>

using namespace saucats;



template <class Scalar>
bool
relative_distance(Scalar a, Scalar b) {
	return std::fabs((a - b) / std::min(a, b));
}



// --- Sphere class checks ----------------------------------------------------

int
sphere_distance_consistency_check() {
	Sphere3d S(1.);

	Sphere3d::vector_type A(.5, 0., 0.);
	mu_assert(relative_distance(std::pow(S.signed_dist(A), 2), S.squared_dist(A)) <= std::numeric_limits<float>::epsilon(), "inconsistent values returned for SphereT::squared_dist and SphereT::signed_dist");

	Sphere3d::vector_type B(2., 0., 0.);
	mu_assert(relative_distance(std::pow(S.signed_dist(A), 2), S.squared_dist(A)) <= std::numeric_limits<float>::epsilon(), "inconsistent values returned for SphereT::squared_dist and SphereT::signed_dist");

	// Job done
	return 0;
}



// --- Entry point ------------------------------------------------------------

int
all_tests() {	
	// Sphere class tests
	mu_run_test(sphere_distance_consistency_check);

	// Job done
	return 0;
}



int
main(int argc, char* argv[]) {
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



#include <minunit.h>

#include <saucats/Geometry>
#include <saucats/Macros>
#include <saucats/SDF>

#include <iostream>

using namespace std;
using namespace saucats;



// --- Utilities --------------------------------------------------------------

template <class Scalar>
Scalar
abs_distance(Scalar a, Scalar b) {
	return std::fabs(a - b);
}



// --- Box class checks -------------------------------------------------------

template <class VectorT>
int
box_contains_consistency_check() {
	BoxT<VectorT> B(VectorT::Constant(1.));
	auto sdf = get_box_sdf(B);

	const unsigned int step_count = 16;
	for(unsigned int i = 0; i < step_count; ++i) {
		VectorT A = VectorT::Zero();
		A(0) = (2. * i) / step_count;
		mu_assert(B.contains(A) == (sdf.dist(A) <= 0.), "inconsistent return returned for BoxT::contains and BoxSDF::dist");
	}

	// Job done
	return 0;
}



// --- Sphere class checks ----------------------------------------------------

template <class VectorT>
int
sphere_contains_consistency_check() {
	SphereT<VectorT> S(1.);
	auto sdf = get_sphere_sdf(S);

	const unsigned int step_count = 16;
	for(unsigned int i = 0; i < step_count; ++i) {
		VectorT A = VectorT::Zero();
		A(0) = (2. * i) / step_count;
		mu_assert(S.contains(A) == (sdf.dist(A) <= 0.), "inconsistent return returned for SphereT::contains and SphereSDF::dist");
	}

	// Job done
	return 0;
}



template <class VectorT>
int
sphere_volume_check() {
	// Setup a unit sphere
	SphereT<VectorT> S(1.);
	auto sdf = get_sphere_sdf(S);
	auto bound_box = get_bounding_box(sdf.bounding_sphere());

	// Estimate the volume with Monte-Carlo integration
	const unsigned int sample_count = 10000000;
	unsigned int sample_inside_count = 0;
	for(unsigned int i = 0; i < sample_count; ++i) {
		VectorT X = VectorT::Random();
		sample_inside_count += sdf.dist(X) <= 0;
	}

	typename VectorT::Scalar volume = bound_box.volume();
	volume *= sample_inside_count;
	volume /= sample_count;

	// Check that the estimate is close enough
	mu_assert(std::fabs(volume - S.volume()) <= 1e-2, "possibly incorrect sphere volume");

	// Job done
	return 0;
}



// --- Entry point ------------------------------------------------------------

int
all_tests() {	
	// Sphere class tests
	mu_run_test(box_contains_consistency_check<Eigen::Vector2f>);
	mu_run_test(box_contains_consistency_check<Eigen::Vector2d>);
	mu_run_test(box_contains_consistency_check<Eigen::Vector3f>);
	mu_run_test(box_contains_consistency_check<Eigen::Vector3d>);

	// Sphere class tests
	mu_run_test(sphere_contains_consistency_check<Eigen::Vector2f>);
	mu_run_test(sphere_contains_consistency_check<Eigen::Vector2d>);
	mu_run_test(sphere_contains_consistency_check<Eigen::Vector3f>);
	mu_run_test(sphere_contains_consistency_check<Eigen::Vector3d>);

	mu_run_test(sphere_volume_check<Eigen::Vector2d>);
	mu_run_test(sphere_volume_check<Eigen::Vector3d>);
	mu_run_test(sphere_volume_check<Eigen::Vector4d>);	

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



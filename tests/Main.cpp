#include <minunit.h>

#include <saucats/Geometry>
#include <saucats/Macros>
#include <saucats/SDF>
#include <saucats/utils/ArrayT.h>

#include <random>
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



template <class VectorT>
int
sphere_circumsphere_check() {
	typedef typename VectorT::Scalar value_type;

	std::default_random_engine rng;
	std::normal_distribution<value_type> gd;
	std::uniform_real_distribution<> ud(0., 1.);

	for(int trial_count = 0; trial_count < 16; ++trial_count) {
		// Generate the radius and the center of the disc
		value_type radius = ud(rng);
		VectorT center = VectorT::Random();

		// Generate n + 1 points on the disc
		std::vector<VectorT> points(VectorT::RowsAtCompileTime + 1);
		for(std::size_t i = 0; i < points.size(); ++i) {
			VectorT P;
			for(Eigen::Index j = 0; j < VectorT::RowsAtCompileTime; ++j)
				P.coeffRef(j) = gd(rng);

			points[i] = radius * P.normalized() + center;
		}

		// Compute circumcircle
		SphereT<VectorT> circle = get_circumsphere(points.begin(), points.end());

		// Check the circumcircle
		mu_assert(std::fabs(radius - circle.radius()) <= 1e-3, "wrong circumsphere radius");
		mu_assert((center - circle.center()).squaredNorm() <= 1e-5, "wrong circumsphere center");
	}

	// Job done
	return 0;
}



// --- Smallest bounding sphere test ------------------------------------------

template <class VectorT>
int
smallest_bounding_sphere_check() {
	typedef typename VectorT::Scalar value_type;

	std::default_random_engine rng;
	std::normal_distribution<value_type> gd;
	std::uniform_real_distribution<> ud(0., 1.);

	for(int trial_count = 0; trial_count < 16; ++trial_count) {
		for(Eigen::Index count = VectorT::RowsAtCompileTime + 2; count < VectorT::RowsAtCompileTime + 30; ++count) {
			// Generate a support sphere from n+1 points
			std::vector<VectorT> support_points(VectorT::RowsAtCompileTime + 1);
			std::generate(support_points.begin(), support_points.end(), [] () -> VectorT { return VectorT::Random(); });
			SphereT<VectorT> support_sphere = get_bounding_sphere(support_points.begin(), support_points.end());

			// Generate points inside the support sphere
			std::vector<VectorT> points(count);
			std::copy(support_points.begin(), support_points.end(), points.begin());
			std::generate(points.begin() + support_points.size(), points.end(), [&support_sphere, &rng, &gd, &ud] () -> VectorT {
				VectorT P;
				for(Eigen::Index j = 0; j < VectorT::RowsAtCompileTime; ++j)
					P.coeffRef(j) = gd(rng);

				value_type radius = ud(rng);
				P *= support_sphere.radius() * radius / P.norm();
				P += support_sphere.center();
				return P;
				//return support_sphere.center();
			});

			// Get the bounding sphere
			SphereT<VectorT> bounding_sphere = get_bounding_sphere(points.begin(), points.end());

			// Check that the bounding sphere and the support sphere are equivalent up to machine precision
			mu_assert(std::fabs(bounding_sphere.radius() - support_sphere.radius()) <= 1e-3, "wrong circumsphere radius");
			mu_assert((bounding_sphere.center() - support_sphere.center()).squaredNorm() <= 1e-5, "wrong circumsphere center");
		}
	}

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
	
	mu_run_test(sphere_circumsphere_check<Eigen::Vector2f>);
	mu_run_test(sphere_circumsphere_check<Eigen::Vector2d>);
	mu_run_test(sphere_circumsphere_check<Eigen::Vector3f>);
	mu_run_test(sphere_circumsphere_check<Eigen::Vector3d>);
	mu_run_test(sphere_circumsphere_check<Eigen::Vector4f>);
	mu_run_test(sphere_circumsphere_check<Eigen::Vector4d>);

	// Smallest bounding sphere test
	mu_run_test(smallest_bounding_sphere_check<Eigen::Vector2d>);
	mu_run_test(smallest_bounding_sphere_check<Eigen::Vector3d>);
	mu_run_test(smallest_bounding_sphere_check<Eigen::Vector4d>);

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



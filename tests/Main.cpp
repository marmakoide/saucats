#include <minunit.h>

#include <saucats/Geometry>
#include <saucats/Macros>
#include <saucats/Render>
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



// --- Sphere sampling checks -------------------------------------------------

template <class VectorT>
int
sphere_surface_sampling_consistency_check() {
	VectorT center = VectorT::Random();

	SphereT<VectorT> sphere(center, 1);
	auto sampler = get_sphere_surface_sampler(sphere);

	std::default_random_engine rng;
	for(unsigned int i = 0; i < 32; ++i) {
		VectorT X = sampler.get_sample(rng);
		mu_assert(((X - center).squaredNorm() - sphere.squared_radius()) < 1e-6, "inconsistent distance to sphere center");
	}

	// Job done
	return 0;
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
	auto bound_box = get_bounding_box(sdf.get_bounding_sphere());

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
	std::default_random_engine rng;
	std::uniform_real_distribution<> ud(0., 1.);

	for(int trial_count = 0; trial_count < 16; ++trial_count) {
		// Generate the radius and the center of the sphere
		SphereT<VectorT> sphere(VectorT::Random(), ud(rng));

		// Generate n + 1 points on the sphere
		std::vector<VectorT> points(VectorT::RowsAtCompileTime + 1);
		auto sampler = get_sphere_surface_sampler(sphere);
		std::generate(points.begin(), points.end(), [&sampler, &rng] () { return sampler.get_sample(rng); });

		// Compute circumsphere
		SphereT<VectorT> circumsphere = get_circumsphere(points.begin(), points.end());

		// Check the circumcircle
		mu_assert(std::fabs(sphere.radius() - circumsphere.radius()) <= 1e-3, "wrong circumsphere radius");
		mu_assert((sphere.center() - circumsphere.center()).squaredNorm() <= 1e-5, "wrong circumsphere center");
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
			// Generate the radius and the center of the sphere
			SphereT<VectorT> seed_sphere(VectorT::Random(), ud(rng));

			// Generate n + 1 points on the support sphere surface
			auto sampler = get_sphere_surface_sampler(seed_sphere);

			std::vector<VectorT> support_points(VectorT::RowsAtCompileTime + 1);
			std::generate(support_points.begin(), support_points.end(), [&sampler, &rng] () { return sampler.get_sample(rng); });			
			SphereT<VectorT> support_sphere = point_collection::get_bounding_sphere(support_points.begin(), support_points.end());

			// Generate points inside the support sphere
			std::vector<VectorT> points(count);
			std::copy(support_points.begin(), support_points.end(), points.begin());
			std::generate(std::next(points.begin(), support_points.size()), points.end(), [&support_sphere, &rng, &gd, &ud] () -> VectorT {
				VectorT P;
				for(Eigen::Index j = 0; j < VectorT::RowsAtCompileTime; ++j)
					P.coeffRef(j) = gd(rng);

				value_type radius = ud(rng);
				P *= support_sphere.radius() * radius / P.norm();
				P += support_sphere.center();
				return P;
			});

			// Get the bounding sphere
			SphereT<VectorT> bounding_sphere = point_collection::get_bounding_sphere(points.begin(), points.end());

			// Check that the bounding sphere and the support sphere are equivalent up to machine precision
			mu_assert(std::fabs(bounding_sphere.radius() - support_sphere.radius()) <= 1e-3, "wrong circumsphere radius");
			mu_assert((bounding_sphere.center() - support_sphere.center()).squaredNorm() <= 1e-5, "wrong circumsphere center");
		}
	}

	// Job done
	return 0;
}



// --- Polygon distance field test --------------------------------------------

int
square_polygon_rasterization_check() {
	std::default_random_engine rng;

	// Setup the distance fields
	Eigen::Matrix<double, Eigen::Dynamic, 2> vertex_array(4, 2);
	vertex_array <<
		 1.,  1.,
		-1.,  1.,
		-1., -1.,
		 1., -1.;

	auto sdf_a = get_polygon2_sdf(vertex_array);
	auto sdf_b = get_box_sdf(Box2d(Eigen::Vector2d(2., 2.)));

	// Check that the bounding spheres match
	auto sphere_a = sdf_a.get_bounding_sphere();
	auto sphere_b = sdf_b.get_bounding_sphere();

	mu_assert(std::fabs(sphere_a.radius() - sphere_b.radius()) <= 1e-3, "wrong bounding sphere radius");
	mu_assert((sphere_a.center() - sphere_b.center()).squaredNorm() <= 1e-5, "wrong bounding sphere center");	

	// Check that the distance functions match for random sample
	auto sampler = get_sphere_volume_sampler(sdf_a.get_bounding_sphere());
	for(unsigned int i = 0; i < 1024; ++i) {
		auto P = sampler.get_sample(rng);

		auto da = sdf_a.dist(P);
		auto db = sdf_b.dist(P);

		mu_assert(std::fabs(da - db) < 1e-3, "wrong signed distance value");
	}

	// Job done
	return 0;
}



// --- Entry point ------------------------------------------------------------

int
all_tests() {
	// Sphere sampling tests
	mu_run_test(sphere_surface_sampling_consistency_check<Eigen::Vector2f>);
	mu_run_test(sphere_surface_sampling_consistency_check<Eigen::Vector2d>);
	mu_run_test(sphere_surface_sampling_consistency_check<Eigen::Vector3f>);
	mu_run_test(sphere_surface_sampling_consistency_check<Eigen::Vector3d>);
	mu_run_test(sphere_surface_sampling_consistency_check<Eigen::Vector4f>);
	mu_run_test(sphere_surface_sampling_consistency_check<Eigen::Vector4d>);
	
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

	// Polygon distance field test
	mu_run_test(square_polygon_rasterization_check);

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



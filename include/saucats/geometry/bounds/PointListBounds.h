#ifndef SAUCATS_GEOMETRY_BOUNDS_POINT_LIST_BOUNDS_H
#define SAUCATS_GEOMETRY_BOUNDS_POINT_LIST_BOUNDS_H

#include <stack>
#include <vector>
#include <memory>
#include <utility>

#include <saucats/geometry/Box.h>
#include <saucats/geometry/Sphere.h>
#include <saucats/utils/ArrayT.h>
#include <saucats/utils/Algorithm.h>



namespace saucats {
	namespace point_collection {
		/*
		 * Compute the minimum volume axis-aligned box bounding of a collection of spheres
		 */

		template <class iterator_type>
		BoxT<typename iterator_type::value_type::vector_type>
		get_bounding_box(iterator_type start_it,
		                 iterator_type end_it) {
			typedef typename iterator_type::vector_type vector_type;

			vector_type min_corner, max_corner;

			bool first_item = true;
			for(; start_it != end_it; ++start_it) {
				if (first_item) {
					min_corner = *start_it;
					max_corner = *end_it;
					first_item = false;
				}
				else {
					min_corner = min_corner.cwiseMin(*start_it);
					max_corner = max_corner.cwiseMax(*end_it);
				}
			}

			return
				BoxT<vector_type>(max_corner - min_corner,
	  		                  (min_corner + max_corner) / 2);
		}



		/*
		 * Compute a minimum volume bounding sphere for a collection of points
		 *
		 * Iterative implementation of Welzl's algorithm, see
		 * "Smallest enclosing disks (balls and ellipsoids)" Emo Welzl 1991
		 */
		template <class iterator_type>
		class PointCollectionBoundingBallT {
		public:
			typedef typename iterator_type::value_type point_type;
			typedef typename point_type::Scalar scalar_type;
			typedef SphereT<point_type> sphere_type;
			typedef std::size_t index_type;
			typedef std::vector<index_type> index_list_type;


			struct Node {
				index_list_type P;
				index_list_type R;
				sphere_type D;
				index_type pivot;
				std::unique_ptr<Node> left;
				std::unique_ptr<Node> right;

				Node(const index_list_type& in_P,
				     const index_list_type& in_R) :
					P(in_P),
					R(in_R),
					pivot(0) { }
			}; // struct Node



			static sphere_type
			get_bounding_sphere(iterator_type start_it,
		                      iterator_type end_it,
			                    scalar_type epsilon) {
				index_list_type P(std::distance(start_it, end_it));
				std::iota(P.begin(), P.end(), 0);

				index_list_type R;
				Node root(P, R);
				traverse(start_it, &root, epsilon);
				return root.D;
			}

		private:
			static sphere_type
			get_boundary(iterator_type start_it,
			             const index_list_type& index_list,
			             scalar_type epsilon) {
				if (index_list.empty())
					return sphere_type();

				std::vector<point_type> point_subset(std::min(index_list.size(), index_list_type::size_type(point_type::RowsAtCompileTime) + 1));
				std::transform(index_list.begin(), std::next(index_list.begin(), point_subset.size()), point_subset.begin(), [start_it](index_type i) { return *std::next(start_it, i); });
			
				sphere_type D = get_circumsphere(point_subset.begin(), point_subset.end());
		
				if (index_list.size() <= index_list_type::size_type(point_type::RowsAtCompileTime) + 1)
					return D;
				
				for(index_type i : index_list)
					if (std::fabs((*std::next(start_it, i) - D.center()).squaredNorm() - D.squared_radius()) > epsilon)
						return sphere_type();
	
				return D;
			}

			static void
			traverse(iterator_type start_it,
			         Node* root,
			         scalar_type epsilon) {
				std::default_random_engine rng;
	
				std::stack<Node*> stack;
				stack.push(root);

				while(!stack.empty()) {
					Node* node = stack.top();
					stack.pop();

					if (node->P.empty() or (node->R.size() >= index_type(point_type::RowsAtCompileTime) + 1))
						node->D = get_boundary(start_it, node->R, epsilon);
					else if (!node->left) {
						std::uniform_int_distribution<> ud(0., 1.);
						node->pivot = *select_randomly(node->P.begin(), node->P.end(), rng);
						index_list_type P;
						std::copy_if(node->P.begin(), node->P.end(), std::back_inserter(P), [node](index_type i) { return i != node->pivot; } );
						node->left = std::make_unique<Node>(P, node->R);
						stack.push(node);
						stack.push(node->left.get());
					}
					else if (!node->right) {
						if (node->left->D.contains(*std::next(start_it, node->pivot)))
							node->D = node->left->D;
						else {
							index_list_type R(node->R);
							R.push_back(node->pivot);
							node->right = std::make_unique<Node>(node->left->P, R);
							stack.push(node);
							stack.push(node->right.get());
						}
					}
					else {
						node->D = node->right->D;
						node->left.reset();
						node->right.reset();
					}
				}
			}
		}; // class PointCollectionBoundingBallT



		template <class iterator_type>
		static SphereT<typename iterator_type::value_type>
		get_bounding_sphere(iterator_type start_it,
		                    iterator_type end_it,
  		                  typename iterator_type::value_type::Scalar epsilon = std::numeric_limits<typename iterator_type::value_type::Scalar>::epsilon()) {
			return PointCollectionBoundingBallT<iterator_type>().get_bounding_sphere(start_it, end_it, epsilon);
		}
	} // namespace point_collection
} // namespace saucats



#endif // SAUCATS_GEOMETRY_BOUNDS_POINT_LIST_BOUNDS_H

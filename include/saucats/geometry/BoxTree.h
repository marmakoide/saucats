#ifndef SAUCATS_GEOMETRY_BOXTREE_H
#define SAUCATS_GEOMETRY_BOXTREE_H

#include <stack>
#include <vector>
#include <list>
#include <map>
#include <iterator>
#include <algorithm>

#include <saucats/geometry/bounds/BoxListBounds.h>
#include <saucats/utils/ArrayT.h>



namespace saucats {
	/*
	 * Implements a generic bound box hierarchy as a binary tree, aka BVH
	 */

	template <class BoxT, class ValueT>
	class BoxTreeT {
	public:
		typedef ValueT value_type;
		typedef BoxT box_type;
		typedef typename BoxT::vector_type vector_type;



		template <typename node_ptr_type>
		class DepthFirstNodeIteratorT {
		public:
			using iterator_category = std::forward_iterator_tag;
			using value_type = node_ptr_type;
			using difference_type = std::ptrdiff_t;
			using pointer = node_ptr_type*;
  	  using reference = node_ptr_type*&;

			DepthFirstNodeIteratorT() { }

			DepthFirstNodeIteratorT(value_type root) {
				m_stack.push(root);
			}

			DepthFirstNodeIteratorT(const DepthFirstNodeIteratorT& other) :
				m_stack(other.m_stack) { }

			DepthFirstNodeIteratorT&
			operator = (const DepthFirstNodeIteratorT& other) {
				m_stack = other.m_stack;
				return *this;
			}
	
			bool
			operator != (const DepthFirstNodeIteratorT& other) const {
				if (m_stack.empty())
					return !other.m_stack.empty();

				if (other.m_stack.empty())
					return true;

				return m_stack.top() != other.m_stack.top();
			}

			DepthFirstNodeIteratorT
			operator ++ () {
				value_type node = m_stack.top();
				m_stack.pop();
				if (!node->is_leaf()) {
					m_stack.push(node->children().first);
					m_stack.push(node->children().second);
				}
				return *this;
			}

			value_type
			operator * () const {
				return m_stack.top();
			}

		private:
			std::stack<value_type> m_stack;
		}; // class DepthFirstNodeIteratorT



		class Node {
		public:
			typedef std::pair<Node*, Node*> node_pair_type;



			Node(const value_type& value = value_type(),
			     const box_type& box = box_type()) :
				m_value(value),
				m_box(box),
				m_children(node_pair_type(nullptr, nullptr)) { }

			inline DepthFirstNodeIteratorT<Node*> begin() {
				return DepthFirstNodeIteratorT<Node*>(this);
			}

			inline DepthFirstNodeIteratorT<Node*> end() {
				return DepthFirstNodeIteratorT<Node*>();
			}

			inline DepthFirstNodeIteratorT<Node const*> begin() const {
				return DepthFirstNodeIteratorT<Node const*>(this);
			}

			inline DepthFirstNodeIteratorT<Node const*> end() const {
				return DepthFirstNodeIteratorT<Node const*>();
			}

			inline value_type& value() {
				return m_value;
			}

			inline const value_type& value() const {
				return m_value;
			}

			inline box_type& box() {
				return m_box;
			}

			inline const box_type& box() const {
				return m_box;
			}

			inline node_pair_type& children() {
				return m_children;
			}

			inline const node_pair_type& children() const {
				return m_children;
			}

			inline bool is_leaf() const {
				return (m_children.first == nullptr) and (m_children.second == nullptr);
			}

		private:
			value_type m_value;
			box_type m_box;
			node_pair_type m_children;
		}; // class Node



		BoxTreeT() { }

		BoxTreeT(const BoxTreeT& other) {
			if (other.m_node_array.empty())
				m_node_array = ArrayT<Node>(0);
			else
				deep_copy(other.root());
		}
 
		BoxTreeT&
		operator = (const BoxTreeT& other) {
			if (other.m_node_array.empty())
				m_node_array = ArrayT<Node>(0);
			else
				deep_copy(other.root());

			return *this;
		}
	
		inline Node* root() {
			return &(m_node_array[0]);
		}

		inline const Node* root() const {
			return &(m_node_array[0]);
		}



		template <class iterator_type, class bounding_func_type>
		static BoxTreeT
		from_kd_construction(iterator_type first, iterator_type last,
		                     const bounding_func_type& bounding_func) {
			typedef std::vector<Node> item_list_type;

			// Build a list of the input items
			item_list_type item_list;
			item_list.reserve(std::distance(first, last));
			std::transform(first, last, std::back_inserter(item_list), [bounding_func](const ValueT& value) -> Node { return Node(value, bounding_func(value)); });

			// Top-down construction of the balltree
			Node* root = new Node();
			std::stack< std::pair<Node*, item_list_type> > stack;
			stack.push(std::make_pair(root, item_list));

			while(!stack.empty()) {
				// Pop a node
				auto node_item_list = stack.top();
				stack.pop();

				Node* node = node_item_list.first;
				item_list_type& item_list = node_item_list.second;

				// Item list with a single member : copy as a leaf
				if (item_list.size() == 1)
					*node = *(item_list.begin());
				// Item list with multiple members
				else {
					// Compute the bounding sphere of the item list
					std::vector<box_type> box_list;
					box_list.reserve(item_list.size());
					std::transform(item_list.begin(), item_list.end(), std::back_inserter(box_list), [](const Node& node) -> box_type { return node.box(); });
					node->box() = box_collection::get_bounding_box(box_list.begin(), box_list.end());

					// Split the item list in two along longest axis
					auto bbox = box_collection::get_bounding_box(box_list.begin(), box_list.end());
					Eigen::Index axis;
					bbox.half_extent().maxCoeff(&axis);
	
					std::nth_element(item_list.begin(), item_list.begin() + item_list.size() / 2, item_list.end(), [axis](const Node& a, const Node& b) -> bool { return a.box().center()[axis] < b.box().center()[axis]; });

					auto median_it = item_list.begin();
					std::advance(median_it, item_list.size() / 2);

					item_list_type left_item_list;
					std::copy(item_list.begin(), median_it, std::back_inserter(left_item_list));				

					item_list_type right_item_list;
					std::copy(median_it, item_list.end(), std::back_inserter(right_item_list));				

					// Creates the two children
					node->children().first = new Node();
					node->children().second = new Node();

					// Attach the two children and push to stack
					stack.push(std::make_pair(node->children().first, left_item_list));
					stack.push(std::make_pair(node->children().second, right_item_list));
				}
			}

			// Job done
			BoxTreeT ret;
			ret.deep_copy(root);
			return ret;
		}

	private:
		void deep_copy(const Node* root) {
			// Setup the copied in a contiguous block of memory
			m_node_array = ArrayT<Node>(std::distance(root->begin(), root->end()));

			// Copy the nodes depth first
			std::size_t node_count = 0;
			std::stack< std::pair<const Node*, std::size_t> > stack;
			stack.push(std::make_pair(root, node_count++));
			while(!stack.empty()) {
				auto node_id = stack.top();
				stack.pop();

				m_node_array[node_id.second].value() = node_id.first->value();
				m_node_array[node_id.second].box() = node_id.first->box();

				if (!node_id.first->is_leaf()) {
					std::size_t left_id = node_count++;
					std::size_t right_id = node_count++;

					m_node_array[node_id.second].children().first = &(m_node_array[left_id]);
					m_node_array[node_id.second].children().second = &(m_node_array[right_id]);

					stack.push(std::make_pair(node_id.first->children().first, left_id));
					stack.push(std::make_pair(node_id.first->children().second, right_id));
				}
			}
		}

		ArrayT<Node> m_node_array;
	}; // class BoxTreeT
} // namespace saucats



#endif // SAUCATS_GEOMETRY_BOXTREE_H

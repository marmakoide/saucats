#ifndef SAUCATS_UTILS_ARRAY_T_H
#define SAUCATS_UTILS_ARRAY_T_H

#include <cstddef>
#include <algorithm>



namespace saucats {
	/*
	 * STL-friendly implementation of a fixed-sized scalar array, with the size
	 * of the array decided at runtime.
	 *
	 * Its advantages over a plain array allocated with 'new' are
	 *   - Works out of the box with most of STL
	 *   - Allocation/deallocation is far more automated
	 *   - As efficient as plain array 
	 */

	template <class T>
	class ArrayT {
	public:
		typedef T         value_type;
		typedef T*        iterator;
		typedef const T*  const_iterator;
		typedef T&        reference;
		typedef const T&  const_reference;
		typedef size_t    size_type;
		typedef ptrdiff_t difference_type;



		inline ArrayT() : m_size(0), m_data(0) { }

		inline ArrayT(size_t size) : m_size(size) {
			allocate();
		}

		inline ArrayT(const ArrayT<T>& other) : m_size(other.size()) {
			allocate();
			std::copy(other.begin(), other.end(), begin());
		}

		inline ~ArrayT() {
			dispose();
		}



		inline ArrayT<T>& operator = (const ArrayT<T>& other) {
			if (size() != other.size()) {
				dispose();
				m_size = other.size();
				allocate();
			}

			std::copy(other.begin(), other.end(), begin());
			return *this;
		}



		inline iterator begin() {
			return m_data;
		}

		inline const_iterator begin() const {
			return m_data;
		}

		inline iterator end() {
			return m_data + m_size;
		}

		inline const_iterator end() const {
			return m_data + m_size;
		}



		inline T* data() {
			return m_data;
		}

		inline const T* data() const {
			return m_data;
		}



		inline bool empty() const {
			return m_size == 0;
		}

		inline size_t size() const {
			return m_size;
		}

		inline T& front() {
			return m_data[0];
		}

		inline const T& front() const {
			return m_data[0];
		}

		inline T& back() {
			return m_data[m_size - 1];
		}

		inline const T& back() const {
			return m_data[m_size - 1];
		}

		inline T& operator () (size_t i) {
			return m_data[i];
		}

		inline const T& operator () (size_t i) const {
			return m_data[i];
		}

		inline T& operator [] (size_t i) {
			return m_data[i];
		}

		inline const T& operator [] (size_t i) const {
			return m_data[i];
		}

		inline void fill(const T& inValue) {
			std::fill(begin(), end(), inValue);
		}

	private:
		void allocate() {
			m_data = new T[m_size];
		}

		void dispose() {
			if (m_data)
				delete[] m_data;
		}



		size_t m_size;
		T* m_data;
	}; // class ArrayT
} // namespace saucats



#endif // SAUCATS_UTILS_ARRAY_T_H

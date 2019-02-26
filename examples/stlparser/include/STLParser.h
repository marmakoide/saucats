#ifndef STLPARSER_H
#define STLPARSER_H

#include <fstream>
#include <iostream>
#include <Eigen/Dense>



template <class scalar_type>
class STLParser {
public:
	typedef Eigen::Matrix<scalar_type, 3, 3> triangle_data;



	/*
	class Exception {
	public:
		Exception(const std::string& message) :
			m_message(message) { }

		inline const std::string message() const {
			return m_message;
		}

	private:
		std::string m_message;
	}; // class Exception
	*/



	template <class iterator_type>	
	static void load(const std::string& path, iterator_type out_it) {
		std::ifstream input(path.c_str(), std::ios::in | std::ios::binary);
		load(input, out_it);
	}



	template <class iterator_type>	
	static void load(std::istream& input, iterator_type out_it) {
		// Read header
		char header_data[80];
		input.read(header_data, 80);

    // Read triangle count
		char triangle_count_data[4];
		std::uint32_t triangle_count;

		input.read(triangle_count_data, 4);
		triangle_count = 0;
		for(int i = 0; i < 4; ++i)
			triangle_count = 256 * triangle_count + triangle_count_data[3-i];
		
    // For each triangle
    for(; triangle_count > 0; triangle_count -= 1) {
			// Read triangle data
			char facet_data[50];
			input.read(facet_data, 50);
			triangle_data T = Eigen::Map<Eigen::Matrix<float, 3, 3, Eigen::RowMajor> >((float*)(facet_data + 12)).cast<scalar_type>();

			// Append the triangle to the output collection
			++out_it;
			*out_it = T;
 		}
	}
}; // class STLParser



#endif // STLPARSER_H

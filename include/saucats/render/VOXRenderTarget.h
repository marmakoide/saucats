#ifndef SAUCATS_RENDER_VOX_RENDER_TARGET_H
#define SAUCATS_RENDER_VOX_RENDER_TARGET_H

#include <map>
#include <sstream>
#include <fstream>
#include <iostream>
#include <optional>



namespace saucats {
	/*
	 * Render to a VOX file, a file format support by the Magica and Goxel voxel editors
	 */

	class VOXRenderTarget {
	public:
		typedef std::uint32_t size_type;



		class VoxelPos {
		public:
			inline VoxelPos(size_type x, size_type y, size_type z) :
				m_x(x),
				m_y(y),
				m_z(z) { }

			inline VoxelPos(const VoxelPos& other) :
				m_x(other.m_x),
				m_y(other.m_y),
				m_z(other.m_z) { }

			bool
			operator == (const VoxelPos& other) const {
				return
					(m_x == other.m_x) &&
					(m_y == other.m_y) &&
					(m_z == other.m_z);
			}

			bool
			operator != (const VoxelPos& other) const {
				return
					(m_x != other.m_x) ||
					(m_y != other.m_y) ||
					(m_z != other.m_z);
			}

			bool
			operator < (const VoxelPos& other) const {
				if (m_x < other.m_x)
					return true;

				if (m_x > other.m_x)
					return false;

				if (m_y < other.m_y)
					return true;

				if (m_y > other.m_y)
					return false;

				if (m_z < other.m_z)
					return true;

				if (m_z > other.m_z)
					return false;

				return false;
			}

			inline size_type x() const {
				return m_x;
			}

			inline size_type y() const {
				return m_y;
			}

			inline size_type z() const {
				return m_z;
			}

		private:
			size_type m_x;
			size_type m_y;
			size_type m_z;
		}; // class VoxelPos



		VOXRenderTarget(const std::string& path,
		                size_type w,
		                size_type h,
		                size_type d) :
			m_w(w),
			m_h(h),
			m_d(d),
			m_path(path) { }

		inline size_type width() const {
			return m_w;
		}

		inline size_type height() const {
			return m_h;		
		}

		inline size_type depth() const {
			return m_d;		
		}

		inline void set_voxel(size_type i, size_type j, size_type k, std::uint8_t id) {
			m_voxel_data[VoxelPos(i, j, k)] = id;
		}

		inline void start_render() {
			m_voxel_data.clear();
			m_file = std::ofstream(m_path, std::ios::out | std::ios::binary);
		}

		inline void end_render() {
			// Write to the file
			write_header(m_file);

			// Write the main chunk
			std::ostringstream children_chunks_data;
			write_size_chunk(children_chunks_data, m_w, m_h, m_d);
			write_xyzi_chunk(children_chunks_data, m_voxel_data);
			write_chunk(m_file, "MAIN", "", children_chunks_data.str());

			// Close the file
			m_file.close();
		}

	private:
		static void
		write_header(std::ostream& out) {
			// Write the magic string
			char magic_string[4] = { 'V', 'O', 'X', ' ' };
			out.write(magic_string, 4);

			// Write the version number
			std::int32_t version_number = 150;
			out.write(reinterpret_cast<char*>(&version_number), 4);
		}

		static void
		write_chunk(std::ostream& out,
                const char* chunk_identifier,
		            const std::string& data,
		            const std::string& children_data) {
			// Write the chunk identifier
			out.write(chunk_identifier, 4);

			// Write the chunk header
			std::uint32_t chunk_size = data.size();
			out.write(reinterpret_cast<char*>(&chunk_size), 4);

			std::uint32_t children_chunks_size = children_data.size();
			out.write(reinterpret_cast<char*>(&children_chunks_size), 4);

			// Write the chunk data
			out.write(data.c_str(), data.size());
			out.write(children_data.c_str(), children_data.size());
		}

		static void
		write_size_chunk(std::ostream& out,
		                 size_type x,
		                 size_type y,
		                 size_type z) {
			// Generate the chunk's data
			std::ostringstream chunk_data;
			chunk_data.write(reinterpret_cast<char*>(&x), 4);
			chunk_data.write(reinterpret_cast<char*>(&y), 4);
			chunk_data.write(reinterpret_cast<char*>(&z), 4);

			// Write the chunk
			write_chunk(out, "SIZE", chunk_data.str(), "");
		}

		static void
		write_xyzi_chunk(std::ostream& out,
		                 const std::map<VoxelPos, std::uint8_t>& voxel_data) {
			// Generate the chunk's data
			std::ostringstream chunk_data;
			
			std::uint32_t voxel_count = voxel_data.size();
			chunk_data.write(reinterpret_cast<char*>(&voxel_count), 4);

			int i = 0;
			for(std::map<VoxelPos, std::uint8_t>::const_iterator it = voxel_data.cbegin(); it != voxel_data.cend(); ++it,++i) {
				std::uint8_t x = (*it).first.x();
				std::uint8_t y = (*it).first.y();
				std::uint8_t z = (*it).first.z();
				std::uint8_t xyzi_data[4] = { x, y, z, (*it).second };
				chunk_data.write(reinterpret_cast<char*>(&xyzi_data), 4);
			}

			// Write the chunk
			write_chunk(out, "XYZI", chunk_data.str(), "");
		}

	private:
		size_type m_w, m_h, m_d;
		std::map<VoxelPos, std::uint8_t> m_voxel_data;
		std::string m_path;
		std::ofstream m_file;
	}; // class VOXRenderTarget
} // namespace saucats



#endif // SAUCATS_RENDER_VOX_RENDER_TARGET_H

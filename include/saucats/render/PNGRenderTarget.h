#ifndef SAUCATS_RENDER_PNG_RENDER_TARGET_H
#define SAUCATS_RENDER_PNG_RENDER_TARGET_H

#define PNG_SKIP_SETJMP_CHECK
#include <cstdio>
#include <cstdint>

#include <png.h>

#include <Eigen/Dense>



namespace saucats {
	class PngIOException {
	public:
		PngIOException(const std::string& message) :
			m_message(message) { }

		inline const std::string message() const {
			return m_message;
		}

	private:
		std::string m_message;
	}; // class PngIOException



	/*
	 * Render to a PNG file
	 */

	class PNGRenderTarget {
	public:
		PNGRenderTarget(const std::string& path,
		                int w, int h) :
			m_w(w),
			m_h(h),
			m_path(path),
			m_pixel_data(0),
			m_pixel_row_data(0),
			m_fp(0) {
			m_pixel_data = new png_byte[3 * m_w * m_h];
			m_pixel_row_data = new png_bytep[m_h];
			for(int i = 0; i < m_h; ++i)
				m_pixel_row_data[i] = m_pixel_data + 3 * m_w * i;

		}

		~PNGRenderTarget() {
			delete[] m_pixel_row_data;
			delete[] m_pixel_data;
		}

		inline int width() const {
			return m_w;
		}

		inline int height() const {
			return m_h;		
		}

		template <class color_type>
		inline void set_pixel(int i, int j, const color_type& color) {
			// Compute the pixel address
			std::uint8_t* pixel = (std::uint8_t*)m_pixel_data;
			pixel += j * 3 * m_w;
			pixel += i * 3;

			// Write the pixel
			pixel[0] = (std::uint8_t)std::floor(255. * color.coeff(0));
			pixel[1] = (std::uint8_t)std::floor(255. * color.coeff(1));
			pixel[2] = (std::uint8_t)std::floor(255. * color.coeff(2));		
		}

		inline void start_render() {
			// Open the file (better have problem before starting a long, expensive rendering)
			m_fp = fopen(m_path.c_str(), "w");
			if (!m_fp)
				throw PngIOException("fopen failed");
		}

		inline void end_render() {
			// Creates libpng structures
			png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
			if (!png_ptr)
				throw PngIOException("png_create_write_struct failed");
	
			png_infop info_ptr = png_create_info_struct(png_ptr);
			if (!info_ptr)
				throw PngIOException("png_create_info_struct failed");

			if (setjmp(png_jmpbuf(png_ptr)))
				throw PngIOException("png_init_io failed");

			png_init_io(png_ptr, m_fp);

			// Write file header
			if (setjmp(png_jmpbuf(png_ptr)))
				throw PngIOException("Error while writing header");

			png_set_IHDR(png_ptr, info_ptr,
			             m_w,
			             m_h,
	  	             8,
			             PNG_COLOR_TYPE_RGB,
			             PNG_INTERLACE_NONE,
			             PNG_COMPRESSION_TYPE_BASE,
			             PNG_FILTER_TYPE_BASE);

			png_write_info(png_ptr, info_ptr);

			// Write image data
			if (setjmp(png_jmpbuf(png_ptr)))
				throw PngIOException("Error while writing images data");

			png_write_image(png_ptr, m_pixel_row_data);

			// Write end of file
			if (setjmp(png_jmpbuf(png_ptr)))
				throw PngIOException("Error while writing end of file");

			png_write_end(png_ptr, NULL);

			// Destroy libpng structures
			if (info_ptr)
				png_destroy_info_struct(png_ptr, &info_ptr);

			if (png_ptr)
				png_destroy_write_struct(&png_ptr, &info_ptr);

			// Close the file
			fclose(m_fp);
		}

	private:
		int m_w, m_h;
		std::string m_path;

		png_bytep m_pixel_data;
		png_bytep* m_pixel_row_data;

		FILE* m_fp;
	}; // class PNGRenderTarget
} // namespace saucats



#endif // SAUCATS_RENDER_PNG_RENDER_TARGET_H

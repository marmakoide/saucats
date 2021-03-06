#ifndef SAUCATS_RENDER_COLOR_RAMP_H
#define SAUCATS_RENDER_COLOR_RAMP_H

#include <saucats/render/ColorMap.h>



namespace saucats {
	/*
	 * Interpolation operators
	 */

	struct NearestNeighbourInterpolation {
		template <class interpolant_type>
		inline interpolant_type operator () (double k,
		    	                               const interpolant_type& a,
		    	                               const interpolant_type& b) const {
			if (k < .5)
				return a;
			return b;
		}
	}; // struct NearestNeighbourInterpolation



	struct LinearInterpolation {
		template <class interpolant_type>		
		inline interpolant_type operator () (double k,
		                                     const interpolant_type& a,
		                                     const interpolant_type& b) const {
			return (1. - k) * a + k * b;
		}
	}; // struct LinearInterpolation



	template <class color_array_type, class interpolation_type>
	class ColorRamp {
	public:
		typedef typename color_array_type::Scalar scalar_type;
		typedef Eigen::Matrix<typename color_array_type::Scalar, color_array_type::ColsAtCompileTime, 1> color_type;

		ColorRamp(const color_array_type& color_array,
		          const interpolation_type& interpolation) :
			m_color_array(color_array),
			m_interpolation(interpolation) { }

		color_type
		get_color(scalar_type u) const {
			scalar_type u_int = 0.;
			scalar_type u_frac = std::modf(std::fabs(u), &u_int);

			scalar_type v = (m_color_array.rows() - 1) * u_frac;
			scalar_type v_int = 0;
			scalar_type v_frac = std::modf(v, &v_int);

			Eigen::Index i = (Eigen::Index)v_int;
			return
				m_interpolation(v_frac,
				                m_color_array.row(i).eval(),
				                m_color_array.row(i+1).eval());
		}

	private:
		color_array_type m_color_array;
		interpolation_type m_interpolation;
	}; // class ColorRamp



	template <class color_array_type, class interpolation_type>
	ColorRamp<color_array_type, interpolation_type>
	get_color_ramp(const color_array_type& color_array,
  	             const interpolation_type& interpolation) {
		return ColorRamp<color_array_type, interpolation_type>(color_array, interpolation);
	}	

	template <class color_array_type>
	ColorRamp<color_array_type, LinearInterpolation>
	get_color_ramp(const color_array_type& color_array) {
		LinearInterpolation interpolation;
		return ColorRamp<color_array_type, LinearInterpolation>(color_array, interpolation);
	}	



	// --- AbsColorRamp ---------------------------------------------------------

	template <class color_ramp_type>
	class AbsColorRamp {
	public:
		typedef typename color_ramp_type::scalar_type scalar_type;
		typedef typename color_ramp_type::color_type color_type;

		AbsColorRamp(const color_ramp_type& color_ramp) :
			m_color_ramp(color_ramp) { }

		inline color_type get_color(scalar_type u) const {
			return m_color_ramp.get_color(std::fabs(u));
		}

	private:
		color_ramp_type m_color_ramp;
	}; // class AbsColorRampMap



	template <class color_ramp_type>
	AbsColorRamp<color_ramp_type>
	get_abs_color_ramp(const color_ramp_type& color_ramp) {
		return AbsColorRamp<color_ramp_type>(color_ramp);
	}



	// --- SignedColorRamp --------------------------------------------------------

	template <class color_ramp_type>
	class SignedColorRamp {
	public:
		typedef typename color_ramp_type::scalar_type scalar_type;
		typedef typename color_ramp_type::color_type color_type;

		SignedColorRamp(const color_ramp_type& A,
		                const color_ramp_type& B) : 
			m_A(A),
			m_B(B) { }

		inline color_type get_color(scalar_type u) const {
			if (u < 0)
				return m_A.get_color(-u);

			return m_B.get_color(u);
		}

	private:
		color_ramp_type m_A;
		color_ramp_type m_B;
	}; // class SignedRampMap



	template <class color_ramp_type>
	SignedColorRamp<color_ramp_type>
	get_signed_color_ramp(const color_ramp_type& A,
	                      const color_ramp_type& B) {
		return SignedColorRamp<color_ramp_type>(A, B);
	}

} // namespace saucats



#endif // SAUCATS_RENDER_COLOR_RAMP_H

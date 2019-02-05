#ifndef SAUCATS_SDF_CIRCLE_SDF_H
#define SAUCATS_SDF_CIRCLE_SDF_H

#include <saucats/geometry/Circle.h>



namespace saucats {
	/*
	 Euclidean signed distance to a sphere
	 */

	template <class VectorT>
	class CircleSDF {
	public:
		typedef typename VectorT::Scalar scalar_type;
		typedef VectorT vector_type;
		typedef CircleT<VectorT> circle_type;
		typedef SphereT<VectorT> sphere_type;



		inline CircleSDF(const circle_type& circle)
			: m_circle(circle) { }

		template <typename InVectorT>
		inline typename InVectorT::Scalar
		dist(const InVectorT& X) const {
			InVectorT Y = X - m_circle.center();
			typename InVectorT::Scalar d = m_circle.normal().dot(Y);
			InVectorT H = Y - d * m_circle.normal();
			Eigen::Matrix<typename InVectorT::Scalar, 2, 1> I(H.norm() - m_circle.radius(), d);
			return I.norm();
		}

		inline sphere_type get_bounding_sphere() const {
			return get_bounding_sphere(m_circle.center(), m_circle.radius());
		}

	private:
		circle_type m_circle;
	}; // class CircleSDF



	template <class VectorT>
	CircleSDF<VectorT>
	get_circle_sdf(const CircleT<VectorT>& circle) {
		return CircleSDF<VectorT>(circle);
	}
} // namespace saucats



#endif // SAUCATS_SDF_CIRCLE_SDF_H

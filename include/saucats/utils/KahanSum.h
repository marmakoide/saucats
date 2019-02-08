#ifndef SAUCATS_UTILS_KAHAN_SUM_H
#define SAUCATS_UTILS_KAHAN_SUM_H



namespace saucats {
	// Convenience object to use Kahan summation
	template <class ValueT>
	class KahanSum {
	public:
		inline KahanSum() :
			m_sum(0),
			m_c(0) { }

		inline KahanSum(ValueT sum) :
			m_sum(sum),
			m_c(0) { }

		inline KahanSum(const KahanSum& other) :
			m_sum(other.m_sum),
			m_c(other.m_c) { }

		KahanSum& operator += (const ValueT& x) {
			ValueT y = x - m_c;
			ValueT t = m_sum + y;
			m_c = (t - m_sum) - y;
			m_sum = t;
			return *this;
		}

		inline KahanSum& operator -= (const ValueT& x) {
			*this += -x;
			return *this;
		}

		inline const ValueT& operator () () const {
			return m_sum;
		}

	private:
		ValueT m_sum, m_c;
	}; // class KahanSum
}



#endif // SAUCATS_UTILS_KAHAN_SUM_H

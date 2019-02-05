#ifndef SAUCATS_UTILS_STOP_WATCH_H
#define SAUCATS_UTILS_STOP_WATCH_H

#include <chrono>
#include <atomic>



/*
 * Multi-thread friendly stopwatch implementation taken from
 * https://codereview.stackexchange.com/questions/196245/extremely-simple-timer-class-in-c
 */

namespace saucats {
	template <typename Clock = std::chrono::high_resolution_clock>
	class StopWatch {
	public:
		StopWatch() :
			m_start_point(Clock::now()) { }

		template <typename Rep = typename Clock::duration::rep, typename Units = typename Clock::duration>
		Rep elapsed_time() const {
			std::atomic_thread_fence(std::memory_order_relaxed);
			auto counted_time = std::chrono::duration_cast<Units>(Clock::now() - m_start_point).count();
			std::atomic_thread_fence(std::memory_order_relaxed);
			return static_cast<Rep>(counted_time);
		}

	private:
		const typename Clock::time_point m_start_point;
	}; // class StopWatch



	using precise_stopwatch = StopWatch<>;
	using system_stopwatch = StopWatch<std::chrono::system_clock>;
	using monotonic_stopwatch = StopWatch<std::chrono::steady_clock>;
} // namespace saucats



#endif // SAUCATS_UTILS_STOP_WATCH_H

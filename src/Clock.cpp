#include "crpropa/Clock.h"

#if defined(WIN32) || defined(_WIN32)

#define USE_WINDOWS_TIMERS
#define WIN32_LEAN_AND_MEAN
#define NOWINRES
#define NOMCX
#define NOIME
#ifdef _XBOX
#include <Xtl.h>
#else
#include <windows.h>
#endif
#include <ctime>

#else
#include <sys/time.h>
#endif

#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#include <stdexcept>
#endif

namespace crpropa {

#ifdef WIN32
class Clock::Impl {
	LARGE_INTEGER clockFrequency;
	DWORD startTick;
	LONGLONG prevElapsedTime;
	LARGE_INTEGER startTime;
public:
	Impl() {
		QueryPerformanceFrequency(&clockFrequency);
		reset();
	}

	void reset() {
		QueryPerformanceCounter(&startTime);
		startTick = GetTickCount();
		prevElapsedTime = 0;
	}

	/// Returns the time in ms since the last call to reset or since
	/// the btClock was created.
	double getSecond() {
		LARGE_INTEGER currentTime;
		QueryPerformanceCounter(&currentTime);
		LONGLONG elapsedTime = currentTime.QuadPart - startTime.QuadPart;

		// Compute the number of millisecond ticks elapsed.
		unsigned long msecTicks = (unsigned long) (1000 * elapsedTime
				/ clockFrequency.QuadPart);

		// Check for unexpected leaps in the Win32 performance counter.
		// (This is caused by unexpected data across the PCI to ISA
		// bridge, aka south bridge.  See Microsoft KB274323.)
		unsigned long elapsedTicks = GetTickCount() - startTick;
		signed long msecOff = (signed long) (msecTicks - elapsedTicks);
		if (msecOff < -100 || msecOff > 100) {
			// Adjust the starting time forwards.
			LONGLONG msecAdjustment = std::min(
					msecOff * clockFrequency.QuadPart / 1000,
					elapsedTime - prevElapsedTime);
			startTime.QuadPart += msecAdjustment;
			elapsedTime -= msecAdjustment;
		}

		// Store the current elapsed time for adjustments next time.
		prevElapsedTime = elapsedTime;

		// Convert to microseconds.
		unsigned long usecTicks = (unsigned long) (1000000 * elapsedTime
				/ clockFrequency.QuadPart);

		return double(usecTicks) / 1000000;
	}
};
#else

class Clock::Impl {
	struct timeval startTime;
public:
	Impl() {
		reset();
	}

	void reset() {
		gettimeofday(&startTime, 0);
	}

	/// Returns the time in since the last call to reset or since
	/// the btClock was created.
	double getTime() {
		struct timeval currentTime;
		gettimeofday(&currentTime, 0);
		double t = double(currentTime.tv_sec - startTime.tv_sec);
		t += double(currentTime.tv_usec - startTime.tv_usec) / 1000000.;
		return t;
	}
};
#endif

Clock::Clock() {
	impl = new Impl;
}

Clock::~Clock() {
	delete impl;
}

void Clock::reset() {
	impl->reset();
}

double Clock::getSecond() {
	return impl->getTime();
}

double Clock::getMillisecond() {
	return impl->getTime() * 1000;
}

#ifdef _OPENMP

// see http://stackoverflow.com/questions/8051108/using-the-openmp-threadprivate-directive-on-static-instances-of-c-stl-types
const static int MAX_THREAD = 256;

struct CLOCK_TLS_ITEM {
	Clock r;
	char padding[(sizeof(Clock) / 64 + 1) * 64 - sizeof(Clock)];
};

Clock &Clock::getInstance() {
#ifdef _MSC_VER
	__declspec(align(64)) static CLOCK_TLS_ITEM tls[MAX_THREAD];
#else
	__attribute__ ((aligned(64))) static CLOCK_TLS_ITEM tls[MAX_THREAD];
#endif
	int i = omp_get_thread_num();
	if (i >= MAX_THREAD)
	throw std::runtime_error("crpropa::Clock: more than MAX_THREAD threads!");
	return tls[i].r;
}
#else
Clock &Clock::getInstance() {
	static Clock r;
	return r;
}
#endif

} /* namespace scs */

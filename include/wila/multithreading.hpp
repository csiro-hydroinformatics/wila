#pragma once

#include <string>
#include <vector>
#include <exception>
#include <mutex>
#include <atomic>
#include <thread>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <boost/lexical_cast.hpp>


// On Pearcey, TIME_UTC is not defined with these include, and the threadpool fails to compile with:
// error: TIME_UTC was not declared in this scope
//           xtime_get(&xt, TIME_UTC);
#ifndef _WIN32
#ifndef TIME_UTC
#define TIME_UTC 1
#endif
#endif

#include "boost/threadpool.hpp"

#ifdef _WIN32
#ifdef _MSC_VER
#include <concurrent_vector.h>
using namespace Concurrency;
#endif
#else
#include <tbb/concurrent_vector.h>
using namespace tbb;
#endif

#ifdef _WIN32
	//_set_output_format(_TWO_DIGIT_EXPONENT);
#endif

using namespace std;

namespace mhcpp
{
	namespace utils
	{
		template <typename T = std::function<void()>>
		class CrossThreadExceptions
		{
		public:
			CrossThreadExceptions(const std::vector<T>& tasks, const std::vector<std::function<void()>>& cleanup)
			{
				this->tasks = tasks;
				this->cleanup = cleanup;
			}
			CrossThreadExceptions(const std::vector<T>& tasks)
			{
				this->tasks = tasks;
			}
			~CrossThreadExceptions()
			{
			}

			void TryExecute(T& f)
			{
				try {
					f();
				}
				catch (const std::exception&) {
					threadExceptions.push_back(current_exception());
				}
			}

			void ExecuteTasks(size_t nThreads)
			{
				threadExceptions.clear();
				boost::threadpool::pool tp(nThreads);
				// Using resize leads to a memory leak on Linux, and eventually cannot create threads.
				//if (!tp.size_controller().resize(nThreads)) {
				//	// this can be false on Linux. But, this is unclear how to prevent this. Since the size_controller 
				//  // itself also catches a thread exception and 'just' return false, let's see if we can ignore it.
				//	throw std::runtime_error(string("Unable to set size of threadpool. Caller requested nThreads=") + boost::lexical_cast<string>(nThreads));
				//}

				for (int i = 0; i < tasks.size(); i++)
				{
					T tryExecuteFunc =
						[=]()
					{
						this->TryExecute(tasks[i]);
					};
					boost::threadpool::schedule(tp, tryExecuteFunc);
				}
				tp.wait();

				Cleanup();

				if (threadExceptions.size() > 0)
				{
					rethrow_exception(threadExceptions[0]);
				}
			}

		private:
			void Cleanup()
			{
				try {
					for (int i = 0; i < cleanup.size(); i++)
					{
						cleanup[i]();
					}
				}
				catch (const std::exception&) {
					if (threadExceptions.empty())
						threadExceptions.push_back(current_exception());
					else
						threadExceptions[0] = current_exception();

				}
			}
			concurrent_vector<exception_ptr> threadExceptions;
			std::vector<T> tasks;
			std::vector<std::function<void()>> cleanup;

		};

	}
}

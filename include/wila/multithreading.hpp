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
	namespace threading
	{
		/**
		* \class	ThreadingOptions
		*
		* \brief	Global, process-wide default options relating to threading. 
		*			Required to work within limits set e.g. by a slurm (cluster) job and not automagically detectable.
		*/
		template <typename T = double>
		class ThreadingOptions
		{
		public:
			static size_t DefaultMaxDegreeOfParallelism;
		};

		/**
		* \brief	Global default for multi-threaded operations. 
		*			If 0, up to algorithms to ; typically trying to detect hardware concurrency, 
		*			but this present field is highly preferable for e.g. shared cluster-based jobs.
		*/
		size_t ThreadingOptions<double>::DefaultMaxDegreeOfParallelism = 0;
	}
}

namespace mhcpp
{
	namespace utils
	{
		/**
		* \class	CrossThreadExceptions
		*
		* \brief	A class to facilitate multithreaded tasks execution, taking care 
		*			of cross-thread exceptions that otherwise would cause program termination.
		*/
		template <typename T = std::function<void()>>
		class CrossThreadExceptions
		{
		public:
			/**
			* \brief	Constructor.
			*
			* \param	tasks the tasks that have to be executed by this object
			* \param	cleanup optional cleanup tasks that may need to be performed after the tasks executions
			*/
			CrossThreadExceptions(const std::vector<T>& tasks, const std::vector<std::function<void()>>& cleanup)
			{
				this->tasks = tasks;
				this->cleanup = cleanup;
			}

			/**
			* \brief	Constructor.
			*
			* \param	tasks the tasks that have to be executed by this object
			*/
			CrossThreadExceptions(const std::vector<T>& tasks)
			{
				this->tasks = tasks;
			}

			CrossThreadExceptions()
			{
			}

			~CrossThreadExceptions()
			{
			}

			/**
			* \brief	Sets the size of the thread pool.
			*/
			void PoolSize(size_t nThreads)
			{
				if (tp.size() != nThreads) {
					if (!tp.size_controller().resize(nThreads)) {
						// this can be false on Linux. But, this is unclear how to prevent this. Since the size_controller 
						// itself also catches a thread exception and 'just' return false, let's see if we can ignore it.
						throw std::runtime_error(string("Unable to set size of threadpool. Caller requested nThreads=") + boost::lexical_cast<string>(nThreads));
					}
				}
			}

			/**
			* \brief	Execute a set of tasks
			*
			* \param	tasks the tasks that have to be executed by this object vis the underlying thread pool
			*/
			void ExecuteTasks(const std::vector<T>& tasks)
			{
				this->tasks = tasks;
				ExecuteTasks();
			}

			void ExecuteTasks(size_t nThreads)
			{
				PoolSize(nThreads);
				ExecuteTasks();
			}

		private:

			void TryExecute(T& f)
			{
				try {
					f();
				}
				catch (const std::exception&) {
					threadExceptions.push_back(current_exception());
				}
			}

			void ExecuteTasks()
			{
				threadExceptions.clear();
				if (tasks.size() == 0)
					return;

				if (tp.size() == 0) {
					throw std::logic_error("No thread available in the thread pool. You must use CrossThreadExceptions::PoolSize to set a pool size.");
				}

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


			boost::threadpool::pool tp;

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

#pragma once

#include <string>
#include <vector>
#include <exception>
#include <mutex>
#include <thread>
#include <iostream>
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
#include <concurrent_vector.h>
#else
#include <tbb/concurrent_vector.h>
#endif

#ifdef _WIN32
using namespace Concurrency;
#else
using namespace tbb;
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
				catch (const std::exception& e) {
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
				catch (const std::exception& e) {
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


		template <typename T>
		class InstanceCounter
		{
		private:
			static std::atomic<int> instances;
		protected:
			InstanceCounter() { instances++; }
		public:
			static int NumInstances() { return InstanceCounter<T>::instances; };
			virtual ~InstanceCounter() { instances--; }
		};

		template <typename T>
		std::atomic<int> InstanceCounter<T>::instances(0);

		template<typename T>
		bool AreSetEqual(const std::vector<T>& a, const std::vector<T>& b)
		{
			std::set<T> sa(a.begin(), a.end());
			std::set<T> sb(b.begin(), b.end());
			return(sa == sb);
		}

		template<typename K, typename V>
		const K& GetKey(const std::pair<K, V>& keyValue)
		{
			return keyValue.first;
		}

		template<typename K, typename V>
		const std::vector<K> GetKeys(const std::map<K, V>& aMap)
		{
			std::vector<K> keys(aMap.size());
			transform(aMap.begin(), aMap.end(), keys.begin(), GetKey<K, V>);
			return keys;
		}

		template<typename ElemType>
		void PrintVec(const std::vector<ElemType>& hist, std::ostream& stream)
		{
			int n = hist.size();
			for (size_t i = 0; i < n; ++i)
				stream << i << ": " << std::to_string(hist[i]) << std::endl;
		}

		template<typename ElemType>
		void PrintHistogram(const std::vector<ElemType>& hist, std::ostream& stream, int nstars = 100, char c = '*')
		{
			ElemType total = std::accumulate(hist.begin(), hist.end(), 0);
			size_t n = hist.size();
			std::vector<std::string> s(n);
			for (size_t i = 0; i < n; ++i)
				s[i] = std::string(hist[i] * nstars / total, c);
			PrintVec<ElemType>(s);
		}

		template<typename ElemType>
		std::vector<double> Normalize(const std::vector<ElemType>& hist)
		{
			size_t n = hist.size();
			std::vector<double> p(n);
			ElemType total = std::accumulate(hist.begin(), hist.end(), 0);
			for (size_t i = 0; i < n; ++i)
				p[i] = (double)hist[i] / total;
			return p;
		}

		template<typename ElemType>
		std::vector<ElemType> RelativeDiff(const std::vector<ElemType>& expected, const std::vector<ElemType>& b)
		{
			std::vector<ElemType> result(expected.size());
			for (size_t i = 0; i < expected.size(); i++)
			{
				result[i] = (std::abs(expected[i] - b[i]) / expected[i]);
			}
			return result;
		}

		template<typename ElemType>
		void PrintValues(const std::vector<ElemType>& hist, std::ostream& stream, bool proportions = false)
		{
			int n = hist.size();
			if (!proportions)
			{
				PrintVec(hist, stream);
			}
			else
			{
				auto p = Normalize(hist);
				PrintVec(p, stream);
			}
		}
	}
}

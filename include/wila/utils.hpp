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
#include <concurrent_vector.h>
#else
#include <tbb/concurrent_vector.h>
#endif

#ifdef _WIN32
using namespace Concurrency;
#else
using namespace tbb;
#endif

#ifdef _WIN32
	_set_output_format(_TWO_DIGIT_EXPONENT);
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


		template <typename T>
		class InstanceCounter
		{
		private:
			static std::atomic<int> instances;
		protected:
			InstanceCounter() { instances++; }
			InstanceCounter(const InstanceCounter& src) { instances++; }
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

		template<typename T>
		bool IsSubset(const std::vector<T>& tested, const std::vector<T>& container)
		{
			std::set<T> sa(tested.begin(), tested.end());
			std::set<T> sb(container.begin(), container.end());
			return(sa <= sb);
		}

		template<typename K, typename V>
		const K& GetKey(const std::pair<K, V>& keyValue)
		{
			return keyValue.first;
		}

		template<typename K = string, typename V = string>
		std::map<K, V> MergeDictionaries(const std::map<K, V>& first, const std::map<K, V>& second)
		{
			std::map<K, V> result(first);
			for (auto& x : second)
			{
				result[x.first] = x.second;
			}
			return result;
		}

		template<typename K = string, typename V = string>
		bool HasKey(const std::map<K, V>& dict, const string& key)
		{
			return (dict.find(key) != dict.end());
		}

		template<typename K = string, typename V = string>
		std::map<K, V> Subset(const std::vector<K>& keys, const std::map<K, V>& dict, bool allowMissing = false)
		{
			std::map<K, V> result;
			for (const auto& k : keys)
			{
				if (!HasKey(dict, k) && !allowMissing)
					throw std::logic_error("std::map subsetting key not found: " + k);
				else
					result.emplace(k, V(dict.at(k)));
			}
			return result;
		}

		template<typename K, typename V>
		const std::vector<K> GetKeys(const std::map<K, V>& aMap)
		{
			std::vector<K> keys(aMap.size());
			transform(aMap.begin(), aMap.end(), keys.begin(), GetKey<K, V>);
			return keys;
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

		template<typename T>
		string ToScientificString(T value)
		{
			std::stringstream ss(stringstream::out);
			ss.precision(6);
			ss << std::scientific << value;
			return ss.str();
		}

		template<typename T>
		string ToString(const T& value)
		{
			return std::to_string(value);
		}

		template<>
		inline string ToString<double>(const double& value)
			// Note: one must use inline, otherwise users of this template only library 
			// would probably have linker issues with this specialization.
		{
			return ToScientificString<double>(value);
		}

		template<>
		inline string ToString<std::string>(const std::string& value)
		{
			return value;
		}

		template<typename ElemType>
		void PrintVecLine(const std::vector<ElemType>& v, std::ostream& stream, const std::string& sep=",")
		{
			int n = v.size();
			for (size_t i = 0; i < (n-1); ++i)
				stream << ToString<ElemType>(v[i]) << sep;
			stream << ToString<ElemType>(v[(n-1)]);
		}

		template<typename ElemType>
		void PrintVec(const std::vector<ElemType>& hist, std::ostream& stream)
		{
			int n = hist.size();
			for (size_t i = 0; i < n; ++i)
				stream << i << ": " << ToString<ElemType>(hist[i]) << std::endl;
		}

		template<typename ElemType>
		void PrintHistogram(const std::vector<ElemType>& hist, std::ostream& stream, int nstars = 100, char c = '*')
		{
			ElemType total = std::accumulate(hist.begin(), hist.end(), 0);
			size_t n = hist.size();
			std::vector<std::string> s(n);
			for (size_t i = 0; i < n; ++i)
				s[i] = ToString<ElemType>(hist[i] * nstars / total, c);
			PrintVec<ElemType>(s);
		}

		template < typename K, typename ElemType >
		void PrintRow(const std::vector<K>& keys, std::map<K, std::vector<ElemType>> & sk, size_t i, std::ostream& stream, const std::string& sep=",")
		{
			int n = keys.size();
			string s;
			for (size_t c = 0; c < (n-1); ++c)
			{
				s = keys[c];
				stream << ToString<ElemType>(sk[s][i]) << sep;
			}
			s = keys[(n-1)];
			stream << ToString<ElemType>(sk[s][i]);
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

		template<typename T>
		T LogVarValue(T value, string name, std::ostream& stream = std::cout)
		{
#ifdef LOG_VALUE
			stream << name << ": " << mhcpp::utils::ToString(value) << std::endl;
#endif
			return value;
		}

	}
}

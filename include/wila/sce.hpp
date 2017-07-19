#pragma once

#include <algorithm>
#include <atomic>
#include <thread>
#include <mutex>
#include "sce.h"
#include "core.hpp"
#include "utils.hpp"
#include "evaluations.hpp"
#include "multithreading.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"

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

namespace mhcpp
{
	namespace logging
	{
		/// <summary>
		/// A facade for logging information from optimisation processes, to avoid coupling to specific frameworks.
		/// </summary>
		template<class TSys>
		class ILoggerMh
		{
		public:
			//virtual void Write(std::vector<IBaseObjectiveScores> scores, const std::map<string, string>& ctags) = 0;
			//virtual void Write(FitnessAssignedScores<double> worstPoint, const std::map<string, string>& ctags) = 0;
			virtual void Write(TSys* newPoint, const std::map<string, string>& ctags) = 0;
			virtual void Write(const string& message, const std::map<string, string>& ctags) = 0;
			virtual void Write(const FitnessAssignedScores<double, TSys>& scores, const std::map<string, string>& tags) = 0;
			virtual void Write(const std::vector<FitnessAssignedScores<double, TSys>>& scores, const std::map<string, string>& ctags) = 0;
			virtual void Write(const std::vector<IObjectiveScores<TSys>>& scores, const std::map<string, string>& tags) = 0;
			virtual void Reset() = 0;

			virtual std::map<string, vector<string>>GetStringData() = 0;
			virtual std::map<string, vector<double>>GetNumericData() = 0;
			virtual int GetLength() = 0;
			virtual ILoggerMh<TSys>* CreateNew() = 0;

			virtual ~ILoggerMh () {}
		};

		class LoggerMhHelper
		{
		public:
			template<class T, class TSys>
			static void Write(const FitnessAssignedScores<T,TSys>& scores, const std::map<string, string>& tags, ILoggerMh<TSys>* logger)
			{
				if (logger != nullptr)
					logger->Write(scores, tags);
			}

			template<class T, class TSys>
			static void Write(const std::vector<FitnessAssignedScores<T, TSys>>& scores, const std::map<string, string>& ctags, ILoggerMh<TSys>* logger)
			{
				if (logger != nullptr)
					logger->Write(scores, ctags);
			}


			template<class TSys>
			static void Write(const std::vector<IObjectiveScores<TSys>>& scores, const std::map<string, string>& tags, ILoggerMh<TSys>* logger)
			{
				if (logger != nullptr)
					logger->Write(scores, tags);
			}

			template<class TSys>
			static void Write(string infoMsg, const std::map<string, string>& tags, ILoggerMh<TSys>* logger)
			{
				if (logger != nullptr)
					logger->Write(infoMsg, tags);
			}

			template<typename K = string, typename V = string>
			static std::map<K, V> MergeDictionaries(const std::map<K, V>& first, const std::map<K, V>& second)
			{
				return mhcpp::utils::MergeDictionaries(first, second);
			}

			template<typename K = string, typename V = string>
			static bool HasKey(const std::map<K, V>& m, const string& key)
			{
				return mhcpp::utils::HasKey(m, key);
			}

			static std::map<string, string> CreateTag(const std::initializer_list<std::tuple<string, string>>& tuples)
			{
				std::map<string, string> result;
				for (auto& x : tuples)
				{
					result[get<0>(x)] = get<1>(x);
				}
				return result;
			}

			static std::map<string, string> CreateTag(const std::tuple<string, string>& x)
			{
				std::map<string, string> result;
				result[get<0>(x)] = get<1>(x);
				return result;
			}

			static std::tuple<string, string> MkTuple(const string& key, const string& value)
			{
				return std::make_tuple(key, value);
			}

		};


		template<class TSys>
		class SimpleLogger : public ILoggerMh < TSys >
		{
		public:

			SimpleLogger()
			{
				// Nothing
			}

			~SimpleLogger()
			{
				// Nothing
			}

			void Write(TSys* newPoint, const std::map<string, string>& ctags)
			{
				Add(LogEntry(newPoint->GetValues(), ctags));
			}
			void Write(const string& message, const std::map<string, string>& ctags)
			{
				Add(LogEntry(message, ctags));
			}
			void Write(const FitnessAssignedScores<double, TSys>& scores, const std::map<string, string>& tags)
			{
				// TODO: add fitness if/when use case required.
				Write(scores.Scores(), tags);
			}

			void Write(const std::vector<FitnessAssignedScores<double, TSys>>& scores, const std::map<string, string>& tags)
			{
				// TODO: add fitness if/when use case required.
				Write(FitnessAssignedScores<double, TSys>::GetScores(scores), tags);
			}

			static vector<map<string, double>> flattenInfo(const std::vector<IObjectiveScores<TSys>>& scores)
			{
				vector<map<string, double>> r;
				for (auto& s : scores)
					r.push_back(flattenInfo(s));
				return r;
			}

			static map<string, double> flattenInfo(const IObjectiveScores<TSys>& scores)
			{
				return LoggerMhHelper::MergeDictionaries<string, double>(
					scores.GetObjectiveValues(),
					scores.GetParameterValues()
					);
			}

			void Write(const IObjectiveScores<TSys>& scores, const std::map<string, string>& tags)
			{
				Add(LogEntry(flattenInfo(scores), tags));
			}

			void Write(const std::vector<IObjectiveScores<TSys>>& scores, const std::map<string, string>& tags)
			{
				vector<map<string, double>> logInfo = flattenInfo(scores);
				Add(LogEntry(logInfo, tags));
			}
			void Reset()
			{
				entries.clear();
			}

			vector<string> NamesStringData()
			{
				std::set<string> keys;
				for (auto& e : entries)
					e.AddStringDataKeys(keys);
				vector<string> res;
				for (auto& k : keys)
					res.push_back(k);
				return res;
			}

			vector<string> NamesNumericData()
			{
				std::set<string> keys;
				for (auto& e : entries)
					e.AddNumericDataKeys(keys);
				vector<string> res;
				for (auto& k : keys)
					res.push_back(k);
				return res;
			}

			std::map<string, vector<string>> GetStringData()
			{
				vector<string> keys = NamesStringData();
				int n = GetLength();
				map<string, vector<string>> result;
				for (auto& k : keys)
					result[k] = vector<string>(n);
				int offset = 0;
				for (auto& e : entries)
				{
					int ne = e.GetLength();
					e.FillTags(keys, result, offset);
					offset += ne;
				}
				return result;
			}

			std::map<string, vector<double>> GetNumericData()
			{
				vector<string> keys = NamesNumericData();
				int n = GetLength();
				map<string, vector<double>> result;
				for (auto& k : keys)
					result[k] = vector<double>(n);
				int offset = 0;
				for (auto& e : entries)
				{
					int ne = e.GetLength();
					e.FillNumerics(keys, result, offset);
					offset += ne;
				}
				return result;
			}

			int GetLength()
			{
				int result = 0;
				for (auto& e : entries)
					result += e.GetLength();
				return result;
			}

			ILoggerMh<TSys>* CreateNew()
			{
				return new SimpleLogger<TSys>();
			}

		private:
			class LogEntry
			{
			public:
				LogEntry(const string& message, const map<string, string>& tags)
				{
					map<string, string> m;
					m["Message"] = message;
					SetTags(LoggerMhHelper::MergeDictionaries<>(m, tags));
				}


				LogEntry(const map<string, double>& data, const map<string, string>& tags)
				{
					this->data.push_back(data);
					SetTags(tags);
				}

				LogEntry(const vector<map<string, double>>& data, const map<string, string>& tags)
				{
					for (auto& d : data)
						this->data.push_back(d);
					SetTags(tags);
				}

				LogEntry(const LogEntry& src)
				{
					this->data = src.data;
					this->tags = src.tags;
				}

				LogEntry(LogEntry&& src)
				{
					std::swap(this->data, src.data);
					std::swap(this->tags, src.tags);
				}

				LogEntry()
				{
					// Nothing
				}

				~LogEntry()
				{
					// Nothing
				}


				LogEntry& operator=(const LogEntry& src)
				{
					if (&src == this){
						return *this;
					}
					this->data = src.data;
					this->tags = src.tags;
					return *this;
				}

				LogEntry& operator=(LogEntry&& src)
				{
					if (&src == this){
						return *this;
					}
					std::swap(this->data, src.data);
					std::swap(this->tags, src.tags);
					return *this;
				}

				void AddStringDataKeys(std::set<string>& keys)
				{
					for (auto& tag : this->tags)
						keys.emplace(tag.first);
				}

				void AddNumericDataKeys(std::set<string>& keys)
				{
					if (data.size() == 0)
						return;
					// HACK, but OK if we are requiring uniform num data keys
					for (auto& kvp : data[0])
						keys.emplace(kvp.first);
				}

				void FillTags(const vector<string>& keys, map<string, vector<string>> & toFill, int offset)
				{
					int n = GetLength();
					for (auto& k : keys)
					{
						if (LoggerMhHelper::HasKey(tags, k))
						{
							vector<string>& v = toFill[k];
							for (size_t i = 0; i < n; i++)
								v[i + offset] = tags[k];
						}
					}
				}

				void FillNumerics(const vector<string>& keys, map<string, vector<double>> & toFill, int offset)
				{
					for (auto& k : keys)
					{
						vector<double>& v = toFill[k];
						for (size_t i = 0; i < data.size(); i++)
						{
							if (!LoggerMhHelper::HasKey(data[i], k))
								throw std::logic_error("Simple log: keys for numeric data must always exist: " + k);
							v[i + offset] = data[i][k];
						}
					}
				}

				size_t GetLength() { return data.size(); }

			private:
				vector<map<string, double>> data;
				map<string, string> tags;
				void SetTags(const map<string, string>& tags)
				{
#ifdef _DEBUG
					if (LoggerMhHelper::HasKey(tags, ""))
						throw std::logic_error("tags keys cannot be empty strings");
#endif
					this->tags = tags;
				}
			};

			concurrent_vector<LogEntry> entries;

			void Add(const LogEntry& entry)
			{
				this->entries.push_back(entry);
			}
		};

	}

	namespace utils {
		template<typename TSys>
		void PrintLoggerTo(mhcpp::logging::ILoggerMh<TSys>* logger, std::ostream& stream)
		{
			int n = logger->GetLength();
			std::map<string, vector<string>> s = logger->GetStringData();
			std::map<string, vector<double>> d = logger->GetNumericData();
			const std::vector<string> sk = GetKeys<string,vector<string>>(s);
			const std::vector<string> dk = GetKeys<string,vector<double>>(d);
			std::string sep(",");
			PrintVecLine<string>(sk, stream, sep);
			stream << sep;
			PrintVecLine<string>(dk, stream, sep);
			stream << std::endl;
			for (size_t i = 0; i < n; ++i)
			{
				PrintRow<string,string>(sk, s, i, stream, sep);
				stream << sep;
				PrintRow<string,double>(dk, d, i, stream, sep);
				stream << std::endl;
			}
		}
	}

}

namespace mhcpp
{
	namespace optimization
	{
		using namespace mhcpp::logging;
		using namespace mhcpp::objectives;

		template<typename S>
		S CreateSceParamsForProblemOfDimension(int n, int nshuffle)
		{
			throw std::logic_error("Not Implemented");
		}

		template<typename S>
		S CreateDefaultSceParams()
		{
			throw std::logic_error("Not Implemented");
		}

		template<>
		inline SceParameters CreateDefaultSceParams<SceParameters>()
		{
			const int n = 4;
			const int nshuffle = 40;

			SceParameters result;
			result.P = n + 2;
			result.Pmin = n + 2;
			result.M = 2 * n + 1;
			result.Q = std::max(result.M - 2, 2);
			result.Alpha = 1;
			result.Beta = result.M;
			result.NumShuffle = nshuffle;
			result.TrapezoidalDensityParameter = 1.0;
			result.ReflectionRatio = -1.0;
			result.ContractionRatio = 0.5;

			return result;
		}

		template<typename S>
		S AdjustSceParamsForProblemOfDimension(const S& sceParams, int n)
		{
			throw std::logic_error("Not Implemented");
		}

		template<>
		inline SceParameters AdjustSceParamsForProblemOfDimension<SceParameters>(const SceParameters& sceParams, int n)
		{
			if (n <= 0)
				throw new std::logic_error("There must be at least one free parameter to calibrate");
			SceParameters result = sceParams;
			result.M = 2 * n + 1;
			result.Q = std::max(result.M - 2, 2);
			result.Beta = result.M;

			return result;
		}

		template<>
		inline SceParameters CreateSceParamsForProblemOfDimension<SceParameters>(int n, int nshuffle)
		{
			if (n <= 0)
				throw new std::logic_error("There must be at least one free parameter to calibrate");
			SceParameters result = CreateDefaultSceParams<SceParameters>();
			result.NumShuffle = nshuffle;
			return AdjustSceParamsForProblemOfDimension(result, n);
		}

		template<typename T>
		class Complex;

		template<typename T>
		class ShuffledComplexEvolution;

		template<typename T>
		class SubComplex
		{
		private:
			size_t nbObjectiveEvaluations = 0;

			void IncrementEvaluations(size_t n)
			{
				nbObjectiveEvaluations += n;
			}

		public:

			SubComplex(const std::vector<IObjectiveScores<T>>& complexPopulation, IObjectiveEvaluator<T>* evaluator, int q, int alpha,
				IRandomNumberGeneratorFactory<> rng, ICandidateFactory<T> * candidateFactory, RngInt<>& discreteRng, 
				IFitnessAssignment<double, T> fitnessAssignment, ILoggerMh<T>* logger = nullptr, const std::map<string, string>& tags = std::map<string, string>(), SceOptions options = SceOptions::RndInSubComplex, double contractionRatio = 0.5, double reflectionRatio = -1.0)
			{
				Init(complexPopulation, evaluator, q, alpha, rng, candidateFactory, fitnessAssignment, logger, tags, options, contractionRatio, reflectionRatio, discreteRng);
			}

			SubComplex(Complex<T>& complex)
			{
				this->complex = &complex;
				Init(complex.scores, complex.evaluator, complex.q, complex.alpha,
					complex.rng, complex.candidateFactory, complex.fitnessAssignment, complex.logger, complex.tags,
					complex.options, complex.ContractionRatio, complex.ReflectionRatio, complex.discreteGenerator);
			}

			~SubComplex()
			{
				// Nothing;
			}

			bool IsCancelledOrFinished()
			{
				if (complex != nullptr)
					return this->complex->IsCancelledOrFinished();
				return false;
			}

			std::map<string, string> createTagConcat(const std::initializer_list<std::tuple<string, string>>& tuples)
			{
				return LoggerMhHelper::MergeDictionaries<>(LoggerMhHelper::CreateTag(tuples), this->tags);
			}

			size_t EvaluationCount()
			{
				return nbObjectiveEvaluations;
			}

			void Evolve()
			{
				int a = 0;
				while (a < alpha && !IsCancelledOrFinished())
				{
					std::vector<FitnessAssignedScores<double, T>> withoutWorstPoint;
					FitnessAssignedScores<double, T> worstPoint = FindWorstPoint(evolved, withoutWorstPoint);
					loggerWrite(worstPoint, createTagConcat({
						LoggerMhHelper::MkTuple("Message", "Worst point in subcomplex"),
						createTagCatComplexNo() }));
					loggerWrite(withoutWorstPoint, createTagConcat({
						LoggerMhHelper::MkTuple("Message", "Subcomplex without worst point"),
						createTagCatComplexNo() }
						));
					T centroid = getCentroid(withoutWorstPoint);

					bool success;
					T reflectedPoint = reflect(worstPoint, centroid, success);
					if (success)
					{
						std::vector<FitnessAssignedScores<double, T>> candidateSubcomplex;
						// TODO: if we are using multi-objective Pareto optimisation, the following two lines should be re-scrutinized. 
						// While the behavior in the C# implementation was sensical, there may be some smarter ways to test the 
						// improvement in the fitness. Also a valid remark in single-objective case.
						FitnessAssignedScores<double, T> fitReflectedPoint = evaluateNewSet(reflectedPoint, withoutWorstPoint, candidateSubcomplex);
						if (fitReflectedPoint.CompareTo(worstPoint) <= 0)
						{
							replaceEvolved(candidateSubcomplex);
							//deleteElements(evolved);
							//evolved = candidateSubcomplex;

							loggerWrite(fitReflectedPoint, createTagConcat({
								LoggerMhHelper::MkTuple("Message", "Reflected point in subcomplex"),
								createTagCatComplexNo() }));
						}
						else
						{
							clear(candidateSubcomplex);
							loggerWrite(fitReflectedPoint,
								createTagConcat({ LoggerMhHelper::MkTuple("Message", "Reflected point in subcomplex - Failed"), createTagCatComplexNo() }));
							candidateSubcomplex = contractionOrRandom(withoutWorstPoint, worstPoint, centroid);
							//if (candidateSubcomplex == nullptr) // this can happen if the feasible region of the parameter space is not convex.
							//	candidateSubcomplex = fitnessAssignment.AssignFitness(bufferComplex);
							replaceEvolved(candidateSubcomplex);
						}
					}
					else
					{
						// 2012-02-02 A change to fit the specs of the Duan 1993 paper, to validate the use for AWRA-L.
						// This change is in line after discussions with Neil Viney
						// TODO: After discussion with Neil Viney (2012-02-03): Allow for a strategy where the generation 
						// of the random point can be based on another hypercube than the complex. Duan documents that, but this may
						// prevent a faster convergence.
						//evolved = contractionOrRandom(withoutWorstPoint, worstPoint, centroid, bufferComplex);

						replaceEvolved(replaceWithRandom(withoutWorstPoint, worstPoint));
#ifdef _DEBUG
						//CheckParameterFeasible(evolved);
#endif

					}
					a++;
				}
			}

			std::vector<IObjectiveScores<T>> WholePopulation()
			{
				auto scores = aggregatePoints(evolved, leftOutFromSubcomplex);
#ifdef _DEBUG
				//CheckParameterFeasible(scores);
#endif
				return scores;
			}

			static FitnessAssignedScores<double, T> FindWorstPoint(std::vector<FitnessAssignedScores<double, T>>& subComplex, std::vector<FitnessAssignedScores<double, T>>& pointRemoved)
			{
				auto tmp = AsPointers<FitnessAssignedScores<double, T>>(subComplex);
				FitnessAssignedScores<double, T>::Sort(tmp);
				auto worst = tmp[tmp.size() - 1];
				pointRemoved.clear();
				for (auto& point : subComplex)
				{
					if (&point != worst)
						pointRemoved.push_back(point);
				}
				return *worst;
			}

		private:
			std::vector<IObjectiveScores<T>> leftOutFromSubcomplex;
			std::vector<FitnessAssignedScores<double, T>> evolved;
			SceOptions options;
			IObjectiveEvaluator<T>* evaluator;
			IFitnessAssignment<double, T> fitnessAssignment;
			double ContractionRatio;
			double ReflectionRatio;
			Complex<T>* complex = nullptr;
			int alpha, q;
			IRandomNumberGeneratorFactory<> rng;
			ICandidateFactory<T> * cf = nullptr;
			RngInt<> * discreteRng = nullptr;
			std::map<string, string> tags;
			ILoggerMh<T>* logger = nullptr;

			void Init(const std::vector<IObjectiveScores<T>>& complexPopulation, IObjectiveEvaluator<T>* evaluator, int q, int alpha,
				IRandomNumberGeneratorFactory<> rng, ICandidateFactory<T> * candidateFactory,
				IFitnessAssignment<double, T> fitnessAssignment, ILoggerMh<T>* logger, const std::map<string, string>& tags, SceOptions options,
				double contractionRatio, double reflectionRatio, RngInt<>& discreteRng)
			{
				this->options = options;
				this->evaluator = evaluator;
				this->fitnessAssignment = fitnessAssignment;
				this->ContractionRatio = contractionRatio;
				this->ReflectionRatio = reflectionRatio;
				this->alpha = alpha;
				this->q = q;
				this->rng = rng;
				this->cf = candidateFactory;
				if(cf==nullptr)
					cf = new UniformRandomSamplingFactory<T>(rng, complexPopulation[0].SystemConfiguration());
				this->discreteRng = &discreteRng;
				this->tags = tags;
				this->logger = logger;
#ifdef _DEBUG
				//CheckParameterFeasible(complexPopulation);
#endif
				evolved = getSubComplex(complexPopulation, leftOutFromSubcomplex);
#ifdef _DEBUG
				//CheckParameterFeasible(leftOutFromSubcomplex);
#endif
			}

			void clear(std::vector<FitnessAssignedScores<double, T>>& vec)
			{
				//for (auto ptr : vec)
				//	delete ptr;
				vec.clear();
			}

			void replaceEvolved(const std::vector<FitnessAssignedScores<double, T>>& candidateSubcomplex)
			{
				clear(evolved);
				evolved = candidateSubcomplex;
			}

			std::vector<FitnessAssignedScores<double, T>> replaceWithRandom(const std::vector<FitnessAssignedScores<double, T>>& withoutWorstPoint, const FitnessAssignedScores<double, T>& worstPoint)
			{
				std::vector<FitnessAssignedScores<double, T>> result;
				if ((options & SceOptions::RndInSubComplex) == SceOptions::RndInSubComplex)
				{
					// TODO: what was SceOptions::ReflectionRandomization about? Needs re-clarification to assess whether worth porting.
					//						if ((options & SceOptions::ReflectionRandomization) == SceOptions::ReflectionRandomization)
					result = generateRandomWithinSubcomplex(withoutWorstPoint, worstPoint);
					//						else
					//							result = generateRandomWithinShuffleBounds(worstPoint, withoutWorstPoint);
				}
				else
				{
					result = addRandomInHypercube(withoutWorstPoint, this->evolved);
				}
				return result;
			}

			T getCentroid(const std::vector<FitnessAssignedScores<double, T>>& population)
			{
				auto tmp = convertAllToHyperCube(population);
				return T::GetCentroid(tmp);
			}

			void loggerWrite(const string& infoMsg, const std::map<string, string>& ctags)
			{
				LoggerMhHelper::Write(infoMsg, tags, this->logger);
			}

			void loggerWrite(const std::vector<IObjectiveScores<T>>& scores, const std::map<string, string>& ctags)
			{
				auto tags = LoggerMhHelper::MergeDictionaries<>(this->tags, ctags);
				LoggerMhHelper::Write(scores, tags, this->logger);
			}

			void loggerWrite(const FitnessAssignedScores<double, T>& scores, const std::map<string, string>& ctags)
			{
				auto tags = LoggerMhHelper::MergeDictionaries<>(this->tags, ctags);
				LoggerMhHelper::Write(scores, tags, this->logger);
			}

			void loggerWrite(const std::vector<FitnessAssignedScores<double, T>>& scores, const std::map<string, string>& ctags)
			{
				auto tags = LoggerMhHelper::MergeDictionaries<>(this->tags, ctags);
				LoggerMhHelper::Write<double,T>(scores, tags, this->logger);
			}

			std::tuple<string, string> createTagCatComplexNo()
			{
				string id = (this->complex == nullptr ? string("NULL") : this->complex->ComplexId);
				return LoggerMhHelper::MkTuple("Category", "Complex No " + id);
			}

			void loggerWrite(const T& point, const std::map<string, string>& ctags)
			{
			    if (logger != nullptr)
			        logger.Write(point, ctags);
			}

			void loggerWrite(const IObjectiveScores<T>& point, const std::map<string, string>& ctags)
			{
				std::vector < IObjectiveScores<T> > pts = { point };
				this->loggerWrite(pts, ctags);
			}

			std::vector<IObjectiveScores<T>> aggregatePoints(const std::vector<FitnessAssignedScores<double, T>>& subComplex, const std::vector<IObjectiveScores<T>>& leftOutFromSubcomplex)
			{
				std::vector<IObjectiveScores<T>> result = FitnessAssignedScores<double, T>::GetScores(subComplex);
				for (size_t i = 0; i < leftOutFromSubcomplex.size(); i++)
				{
					result.push_back(leftOutFromSubcomplex[i]);
				}
				return result;
			}

			FitnessAssignedScores<double, T> evaluateNewSet(T newPoint, const std::vector<FitnessAssignedScores<double, T>>& withoutWorstPoint, std::vector<FitnessAssignedScores<double, T>>& candidateSubcomplex)
			{
				IObjectiveScores<T> scoreNewPoint = evaluator->EvaluateScore(newPoint);
				IncrementEvaluations(1);

				std::vector<IObjectiveScores<T>> scores = FitnessAssignedScores<double, T>::GetScores(withoutWorstPoint);
				scores.push_back(scoreNewPoint);
				auto fitness = fitnessAssignment.AssignFitness(scores);
				candidateSubcomplex.clear();
				for (size_t i = 0; i < fitness.size(); i++)
				{
					candidateSubcomplex.push_back(FitnessAssignedScores<double, T>(fitness[i]));
				}
				//return Array.Find<FitnessAssignedScores<double, T>>(candidateSubcomplex, (x = > (x.Scores() == scoreNewPoint)));
				return (candidateSubcomplex[candidateSubcomplex.size() - 1]);
			}

			T reflect(const FitnessAssignedScores<double, T>& worstPoint, T centroid, bool& success)
			{
				//double ratio = -1.0;
				double ratio = this->ReflectionRatio;
				return performHomothecy(worstPoint, centroid, ratio, success);
			}

			T contract(FitnessAssignedScores<double, T> worstPoint, T centroid, bool& success)
			{
				//double ratio = 0.5;
				double ratio = this->ContractionRatio;
				return performHomothecy(worstPoint, centroid, ratio, success);
			}

			static T performHomothecy(const FitnessAssignedScores<double, T>& worstPoint, T centroid, double ratio, bool& success)
			{
				T candidate = centroid.HomotheticTransform(worstPoint.Scores().SystemConfiguration(), ratio);
				success = candidate.IsFeasible();
				return candidate;
			}

			static std::vector<T> convertAllToHyperCube(const std::vector<FitnessAssignedScores<double, T>>& points)
			{
				return FitnessAssignedScores<double, T>::GetSystemConfigurations(points);
			}

			static std::vector<T> convertAllToHyperCube(const std::vector<IObjectiveScores<T>>& points)
			{
				return IObjectiveScores<T>::GetSystemConfigurations(points);
			}

			static std::vector<IObjectiveScores<T>> convertToScores(const std::vector<FitnessAssignedScores<double, T>>& points)
			{
				return FitnessAssignedScores<double, T>::GetScores(points);
			}

			std::vector<FitnessAssignedScores<double, T>> removePoint(const std::vector<FitnessAssignedScores<double, T>>& subComplex, FitnessAssignedScores<double, T> worstPoint)
			{
				std::vector<FitnessAssignedScores<double, T>> result;
				for (size_t i = 0; i < subComplex.size(); i++)
				{
					auto p = subComplex.at(i);
					if (p != worstPoint)
						result.push_back(p);
				}
				return result;
			}

			std::vector<FitnessAssignedScores<double, T>> contractionOrRandom(const std::vector<FitnessAssignedScores<double, T>>& withoutWorstPoint,
				const FitnessAssignedScores<double, T>& worstPoint, const T& centroid)
			{
				std::vector<FitnessAssignedScores<double, T>> result;
				std::vector<FitnessAssignedScores<double, T>> candidateSubcomplex;
				bool success;
				T contractionPoint = contract(worstPoint, centroid, success);

				if (success)
				{
					FitnessAssignedScores<double, T> trialPoint = evaluateNewSet(contractionPoint, withoutWorstPoint, candidateSubcomplex);
					if (trialPoint.CompareTo(worstPoint) <= 0)
					{
						result = candidateSubcomplex;
						loggerWrite(trialPoint, createTagConcat({
							LoggerMhHelper::MkTuple("Message", "Contracted point in subcomplex"),
							createTagCatComplexNo() }));
						return result;
					}
					else
					{
						clear(candidateSubcomplex);
						loggerWrite(trialPoint, createTagConcat({
							LoggerMhHelper::MkTuple("Message", "Contracted point in subcomplex-Failed"),
							createTagCatComplexNo() }));
					}
				}
				else
				{
					auto msg = "Contracted point unfeasible";
					loggerWrite(msg, createTagConcat({ LoggerMhHelper::MkTuple("Message", msg), createTagCatComplexNo() }));
				}
				// 2012-02-14: The Duan et al 1993 paper specifies to use the complex to generate random points. However, comparison to a Matlab
				// implementation showed a slower rate of convergence. 
				// result = addRandomInHypercube(withoutWorstPoint, bufferComplex);
				result = replaceWithRandom(withoutWorstPoint, worstPoint);
				return result;
			}

/*
			std::vector<FitnessAssignedScores<double, T>> generateRandomWithinShuffleBounds(const FitnessAssignedScores<double, T>& worstPoint, const std::vector<FitnessAssignedScores<double, T>>  withoutWorstPoint)
			{
				auto sbcplx = convertAllToHyperCube(merge(withoutWorstPoint, worstPoint));
				auto wp = worstPoint.Scores().GetSystemConfiguration() as IHyperCube < double >;
				auto newPoint = wp.Clone() as IHyperCube < double >;
				auto varnames = newPoint.GetVariableNames();
				auto rand = hyperCubeOps.GenerateRandomWithinHypercube(sbcplx);
				for (int i = 0; i < varnames.Length; i++)
				{
					auto v = varnames[i];
					auto value = 2 * centroid.GetValue(v) - wp.GetValue(v);
					if (value < wp.GetMinValue(v) || value > wp.GetMaxValue(v))
						newPoint.SetValue(v, rand.GetValue(v));
					else
						newPoint.SetValue(v, value);
				}
				auto newScore = evaluator->EvaluateScore(newPoint);
				loggerWrite(newScore, createTagConcat(
					LoggerMhHelper::MkTuple("Message", "Adding a partially random point"),
					LoggerMhHelper::MkTuple("Category", "Complex No " + complexId)
					));
				return fitnessAssignment.AssignFitness(aggregate(newScore, withoutWorstPoint));
			}
*/
			std::vector<FitnessAssignedScores<double, T>> generateRandomWithinSubcomplex(const std::vector<FitnessAssignedScores<double, T>>& withoutWorstPoint, const FitnessAssignedScores<double, T>& worstPoint)
			{
				// 2012-02-14: The Duan et al 1993 paper specifies to use the complex to generate random points. However, comparison to a Matlab
				// implementation showed a slower rate of convergence. 
				std::vector<FitnessAssignedScores<double, T>> result;
				auto subCplx = merge(withoutWorstPoint, worstPoint);
				result = addRandomInHypercube(withoutWorstPoint, subCplx);
#ifdef _DEBUG
				//CheckParameterFeasible(result);
#endif
				return result;
			}

			std::vector<FitnessAssignedScores<double, T>> addRandomInHypercube(const std::vector<FitnessAssignedScores<double, T>>& withoutWorstPoint, const std::vector<FitnessAssignedScores<double, T>>& popForHypercubeDefn)
			{
				return 	addRandomInHypercube(withoutWorstPoint, convertToScores(popForHypercubeDefn));
			}

			std::vector<FitnessAssignedScores<double, T>> addRandomInHypercube(const std::vector<FitnessAssignedScores<double, T>>& withoutWorstPoint, const std::vector<IObjectiveScores<T>>& popForHypercubeDefn)
			{
				auto tmp = convertAllToHyperCube(popForHypercubeDefn);
				T newPoint = this->cf->CreateRandomCandidate(tmp);
				//if (newPoint == nullptr)
				//{
				string msg = "Random point within hypercube bounds is unfeasible";
				loggerWrite(msg, createTagConcat({ LoggerMhHelper::MkTuple("Message", msg), createTagCatComplexNo() }));
				//return null;
				auto newScore = evaluator->EvaluateScore(newPoint);
				IncrementEvaluations(1);
				loggerWrite(newScore, createTagConcat({
					LoggerMhHelper::MkTuple("Message", "Adding a random point in hypercube"),
					createTagCatComplexNo() }
					));
				std::vector<IObjectiveScores<T>> newSubComplex = aggregate(newScore, withoutWorstPoint);
#ifdef _DEBUG
				//CheckParameterFeasible(newSubComplex);
#endif
				return fitnessAssignment.AssignFitness(newSubComplex);
			}

			std::vector<FitnessAssignedScores<double, T>> getSubComplex(const std::vector<IObjectiveScores<T>>& bufferComplex, std::vector<IObjectiveScores<T>>& leftOutFromSubcomplex)
			{
				auto fitnessPoints = fitnessAssignment.AssignFitness(bufferComplex);
#ifdef _DEBUG
				//CheckParameterFeasible(fitnessPoints);
#endif
				FitnessAssignedScores<double, T>::Sort(fitnessPoints);
#ifdef _DEBUG
				//CheckParameterFeasible(fitnessPoints);
#endif
				std::vector<FitnessAssignedScores<double, T>> leftOut;
				std::vector<FitnessAssignedScores<double, T>> subset = SampleFrom(*discreteRng, fitnessPoints, q, leftOut, false);
#ifdef _DEBUG
				//CheckParameterFeasible(leftOut);
#endif
				leftOutFromSubcomplex.clear();
				for (size_t i = 0; i < leftOut.size(); i++)
				{
					leftOutFromSubcomplex.push_back(leftOut[i].Scores());
				}
#ifdef _DEBUG
				//CheckParameterFeasible(leftOutFromSubcomplex);
#endif
				return subset;
			}

			static std::vector<IObjectiveScores<T>> merge(const std::vector<FitnessAssignedScores<double, T>>& withoutWorstPoint, const FitnessAssignedScores<double, T>& worstPoint)
			{
				std::vector<IObjectiveScores<T>> tmp = convertToScores(withoutWorstPoint);
				tmp.push_back(worstPoint.Scores());
				return tmp;
			}

			std::vector<IObjectiveScores<T>> aggregatePoints(T newPoint, std::vector<IObjectiveScores<T>> withoutWorstPoint)
			{
				IncrementEvaluations(1);
				return aggregate(evaluator->EvaluateScore(newPoint), withoutWorstPoint);
			}

			std::vector<IObjectiveScores<T>> aggregate(const IObjectiveScores<T>& newPoint, std::vector<FitnessAssignedScores<double, T>> withoutWorstPoint)
			{
				std::vector<IObjectiveScores<T>> result = convertToScores(withoutWorstPoint);
				result.push_back(newPoint);
				return result;
			}
		};

		template<typename T>
		class Complex // : IComplex
		{
			std::map<string, string> createTagConcat(const std::initializer_list<std::tuple<string, string>>& tuples)
			{
				return LoggerMhHelper::MergeDictionaries<>(LoggerMhHelper::CreateTag(tuples), this->tags);
			}

			friend SubComplex<T>::SubComplex(Complex&);

			std::vector<IObjectiveScores<T>> scores;
			int q = 10;
			int alpha = 2;
			int beta = 3;
			RngInt<> discreteGenerator;
			IRandomNumberGeneratorFactory<> rng;
			ICandidateFactory<T> * candidateFactory = nullptr;
			IFitnessAssignment<double, T> fitnessAssignment;
			// IHyperCubeOperations* hyperCubeOps;
			ILoggerMh<T>* logger = nullptr; // new Log4netAdapter();
			double ContractionRatio = 0.5;
			double ReflectionRatio = -1;

			std::map<string, string> tags;
			double factorTrapezoidalPDF = -1;
			SceOptions options = SceOptions::None;


			bool ownEvaluator = true;
			bool ownCandidateFactory = true;
			bool IsCancelled = false;

			ITerminationCondition<T, ShuffledComplexEvolution<T>>* terminationCondition;
			bool allowPrematureTermination = false;

			size_t nbObjectiveEvaluations = 0;

			void IncrementEvaluations(size_t n)
			{
				nbObjectiveEvaluations += n;
			}

		protected:

			IObjectiveEvaluator<T>* evaluator = nullptr;

			void Init(const std::vector<IObjectiveScores<T>>& scores, IObjectiveEvaluator<T>* evaluator, bool ownEvaluator, int q, int alpha = 2, int beta = 3, double factorTrapezoidalPDF = -1,
				SceOptions options = SceOptions::None, double reflectionRatio = -1.0, double contractionRatio = 0.5)
			{
				// TODO checks on consistencies.
				this->scores = scores;
				this->q = q;
				this->alpha = alpha;
				this->beta = beta;
				this->factorTrapezoidalPDF = factorTrapezoidalPDF;
				this->options = options;
				this->ReflectionRatio = reflectionRatio;
				this->ContractionRatio = contractionRatio;
				this->evaluator = evaluator;
				this->ownEvaluator = ownEvaluator;
			}

			Complex(const std::vector<IObjectiveScores<T>>& scores, IObjectiveEvaluator<T>* evaluator, bool ownEvaluator, int q, int seed = 0, int alpha = 2, int beta = 3, double factorTrapezoidalPDF = -1,
				SceOptions options = SceOptions::None, double reflectionRatio = -1.0, double contractionRatio = 0.5) :
				rng(seed),
				discreteGenerator(CreateTrapezoidalRng(scores.size(), rng.CreateNewStd(), factorTrapezoidalPDF))
			{
				Init(scores, evaluator, ownEvaluator, q, alpha, beta, factorTrapezoidalPDF, options, reflectionRatio, contractionRatio);
			}

		public:
			string ComplexId;

			Complex(const std::vector<IObjectiveScores<T>>& scores, 
				IObjectiveEvaluator<T>* evaluator, bool ownEvaluator, IRandomNumberGeneratorFactory<> rng, ICandidateFactory<T>* candidateFactory, bool ownCandidateFactory,
				IFitnessAssignment<double, T> fitnessAssignment, ITerminationCondition<T, ShuffledComplexEvolution<T>>& terminationCondition, bool allowPrematureTermination=true,
				ILoggerMh<T>* logger = nullptr, const std::map<string, string>& tags = std::map<string, string>(), int q=10, int alpha=2, int beta=3, double factorTrapezoidalPDF = -1,
				SceOptions options = SceOptions::None, double reflectionRatio = -1.0, double contractionRatio = 0.5) 
				:
				discreteGenerator(CreateTrapezoidalRng(scores.size(), rng.CreateNewStd(), factorTrapezoidalPDF))
			{
				Init(scores, evaluator, ownEvaluator, q, alpha, beta, factorTrapezoidalPDF, options, reflectionRatio, contractionRatio);
				this->fitnessAssignment = fitnessAssignment;
				this->candidateFactory = candidateFactory;
				this->ownCandidateFactory = ownCandidateFactory;
				this->logger = logger;
				this->tags = tags;
				this->rng = rng;
				this->terminationCondition = &terminationCondition;
				this->allowPrematureTermination = allowPrematureTermination;
			}

			virtual ~Complex()
			{
				if (ownEvaluator)
					if (this->evaluator != nullptr)
					{
						delete evaluator;
						evaluator = nullptr;
					}
				if (ownCandidateFactory)
					if (this->candidateFactory != nullptr)
					{
						delete candidateFactory;
						candidateFactory = nullptr;
					}
			}

			void AllowPrematureTermination(bool allow) { allowPrematureTermination = allow; }
			bool AllowsPrematureTermination() { return allowPrematureTermination; }

			bool IsFinished()
			{
				if (!allowPrematureTermination)
					return false;
				if (terminationCondition == nullptr)
					return false;
				if (!terminationCondition->IsThreadSafe())
					return false;
				return terminationCondition->IsFinished();
			}

			bool IsCancelledOrFinished()
			{
				return (IsCancelled || IsFinished());
			}

			size_t EvaluationCount()
			{
				return nbObjectiveEvaluations;
			}

			void Evolve()
			{
				//if (Thread.CurrentThread.Name == nullptr)
				//{
				//	Thread.CurrentThread.Name = ComplexId;
				//}
				int b; // counters for alpha and beta parameters
				b = 0;
				while (b < beta && !IsCancelledOrFinished())
				{
					SubComplex<T> subComplex(*this);
					subComplex.Evolve();
					IncrementEvaluations(subComplex.EvaluationCount());
					this->scores = subComplex.WholePopulation();
#ifdef _DEBUG
					//CheckParameterFeasible(scores);
#endif
					b++;
				}

			}

			const std::vector<IObjectiveScores<T>> GetObjectiveScores()
			{
				return this->scores;
			}

		};

		template<typename T>
		class ShuffledComplexEvolution
			: public IEvolutionEngine<T>,
			public IPopulation<double, T>
			//where T : ICloneableSystemConfiguration
		{
			// needed? friend Complexes<T>(SCE&);

		public:

			typedef ITerminationCondition<T, ShuffledComplexEvolution<T>> TerminationCondition;

			ShuffledComplexEvolution(const IObjectiveEvaluator<T>& evaluator,
				const ICandidateFactorySeed<T>& candidateFactory,
				const TerminationCondition& terminationCondition,
				const SceParameters& sceParameters,
				bool allowPrematureComplexTermination = true,
				IRandomNumberGeneratorFactory<> rng = IRandomNumberGeneratorFactory<>(),
				IFitnessAssignment<double, T> fitnessAssignment = IFitnessAssignment<double, T>(),
				const std::map<string, string>& logTags = std::map<string, string>())
			{
				IObjectiveEvaluator<T>* pEval = evaluator.Clone();
				Init(pEval, candidateFactory, terminationCondition,
					allowPrematureComplexTermination,
					rng,
					fitnessAssignment,
					sceParameters.P,
					sceParameters.Pmin,
					sceParameters.M,
					sceParameters.Q,
					sceParameters.Alpha,
					sceParameters.Beta,
					sceParameters.NumShuffle,
					sceParameters.TrapezoidalDensityParameter,
					logTags,
					SceOptions::None,
					sceParameters.ReflectionRatio,
					sceParameters.ContractionRatio);
			}

			ShuffledComplexEvolution(IObjectiveEvaluator<T>* evaluator,
				const ICandidateFactorySeed<T>& candidateFactory,
				const TerminationCondition& terminationCondition,
				const SceParameters& sceParameters,
				bool allowPrematureComplexTermination = true,
				IRandomNumberGeneratorFactory<> rng = IRandomNumberGeneratorFactory<>(),
				IFitnessAssignment<double, T> fitnessAssignment = IFitnessAssignment<double, T>(),
				const std::map<string, string>& logTags = std::map<string, string>())
			{
				Init(evaluator, candidateFactory, terminationCondition,
					allowPrematureComplexTermination,
					rng,
					fitnessAssignment,
					sceParameters.P,
					sceParameters.Pmin,
					sceParameters.M,
					sceParameters.Q,
					sceParameters.Alpha,
					sceParameters.Beta,
					sceParameters.NumShuffle,
					sceParameters.TrapezoidalDensityParameter,
					logTags,
					SceOptions::None,
					sceParameters.ReflectionRatio,
					sceParameters.ContractionRatio);
			}

			ShuffledComplexEvolution(const ShuffledComplexEvolution& src)
			{
				throw std::logic_error("copy constructor for ShuffledComplexEvolution not yet supported");
				// CopyBasicsFrom(src);
				// this->evaluator = src.evaluator;
				// this->candidatefactory etc.
				// this->populationInitializer = src.populationInitializer;
				// TODO 
				// ILoggerMh<TSys>* Logger;
			}

			ShuffledComplexEvolution(ShuffledComplexEvolution&& src)
			{
				MoveFrom(src);
			}

			ShuffledComplexEvolution& operator=(ShuffledComplexEvolution&& src)
			{
				if (&src == this) {
					return *this;
				}
				MoveFrom(src);
				return *this;
			}

			ShuffledComplexEvolution& operator=(const ShuffledComplexEvolution& src)
			{
				if (&src == this) {
					return *this;
				}
				throw std::logic_error("copy assignment for ShuffledComplexEvolution not yet supported");
				// CopyBasicsFrom(src);
				// this->evaluator = src.evaluator;
				// this->populationInitializer = src.populationInitializer;
				// TODO 
				// ILoggerMh<TSys>* Logger;
				return *this;
			}

			virtual ~ShuffledComplexEvolution()
			{
				if (evaluator != nullptr)
				{
					delete evaluator;
					evaluator = nullptr;
				}
				if (populationInitializer != nullptr)
				{
					delete populationInitializer;
					populationInitializer = nullptr;
				}
				if (candidateFactory != nullptr)
				{
					delete candidateFactory;
					candidateFactory = nullptr;
				}
				if (logger != nullptr)
				{
					delete logger;
					logger = nullptr;
				}
			}

			void SetLogger()
			{
				if (logger != nullptr)
					delete logger;
				logger = new SimpleLogger<T>();
			}

			ILoggerMh<T>* GetLogger()
			{
				return logger;
			}

			void ResetLog()
			{
				if (logger != nullptr)
					logger->Reset();
			}

			void IncrementEvaluations(size_t n)
			{
				nbObjectiveEvaluations += n;
			}

			IOptimizationResults<T> Evolve()
			{
				Reset();

				this->populationInitializer = candidateFactory->Create();

				std::vector<IObjectiveScores<T>> scores = EvaluateScores(initialisePopulation());
				loggerWrite(scores, createSimpleMsg("Initial Population", "Initial Population"));
				IncrementEvaluations(scores.size());
				auto isFinished = terminationCondition.IsFinished();
				if (isFinished)
				{
					logTerminationConditionMet();
					return packageResults(scores);
				}
				CreateComplexes(scores);

				currentShuffle = 1;
				isFinished = terminationCondition.IsFinished();
				if (isFinished) logTerminationConditionMet();
				while (!isFinished && !isCancelled)
				{
					EvolveComplexes();
					for (size_t i = 0; i < complexes.size(); i++)
					{
						IncrementEvaluations(complexes.at(i)->EvaluationCount());
					}
					string shuffleMsg = "Shuffling No " + std::to_string(currentShuffle);
					std::vector<IObjectiveScores<T>> shufflePoints = AggregateComplexes();
					loggerWrite(shufflePoints, createSimpleMsg(shuffleMsg, shuffleMsg));
					SetPopulationToShuffle(shufflePoints);
					loggerWrite(PopulationAtShuffling[0], createSimpleMsg("Best point in shuffle", shuffleMsg));
					ShuffleComplexes();
					currentShuffle++;
					isFinished = terminationCondition.IsFinished();
					if (isFinished) logTerminationConditionMet();
				}
				return packageResults(complexes);
			}

			size_t PopulationSize() const
			{
				return (this->p * this->m);
			}

			size_t NumComplexes()
			{
				return (this->p);
			}

			std::vector<FitnessAssignedScores<double, T>> Population()
			{
				if (complexes.size() == 0) return std::vector<FitnessAssignedScores<double, T>>();
				return sortByFitness(complexes.Aggregate());
			}

			void UseMultiThreading(bool use)
			{
				this->useMultiThreading = use;
			}

			void SetMaxDegreeOfParallelism(size_t maximum)
			{
				maxDegreeOfParallelism = std::min(GetMaxHardwareConcurrency(), std::max((size_t)1, maximum));
			}

			void SetMaxDegreeOfParallelismHardwareMinus(size_t freeCoresRemaining = 1)
			{
				size_t hardwareMax = GetMaxHardwareConcurrency();
				maxDegreeOfParallelism = std::max<size_t>(1, hardwareMax - freeCoresRemaining);
			}

			size_t GetMaxHardwareConcurrency()
			{
				return (size_t)std::thread::hardware_concurrency();
			}

			int GetMaxDegreeOfParallelism()
			{
				return maxDegreeOfParallelism;
			}

			size_t EvaluationCount()
			{
				return nbObjectiveEvaluations;
			}

			void AllowComplexPrematureTermination(bool allow) { allowPrematureComplexTermination = allow; }
			bool AllowsComplexPrematureTermination() { return allowPrematureComplexTermination; }

			// the next section has protected methods solely to facilitate unit testing.
		protected:

			void Reset()
			{
				isCancelled = false;
				nbObjectiveEvaluations = 0;
				terminationCondition.Reset();
				ResetLog();
				if (this->populationInitializer != nullptr) delete this->populationInitializer;
			}

			std::vector<IObjectiveScores<T>> EvaluateScores(const std::vector<T>& population)
			{
				return Evaluations::EvaluateScores(this->evaluator, population, this->useMultiThreading, GetMaxDegreeOfParallelism());
			}

			void CreateComplexes(const std::vector<IObjectiveScores<T>>& scores)
			{
				this->complexes = partition(scores);
			}

			void EvolveComplexes()
			{
				for (int i = 0; i < complexes.size(); i++)
					complexes.at(i)->ComplexId = std::to_string(i);
				if (useMultiThreading && evaluator->IsCloneable())
				{
					int nThreads = this->GetMaxDegreeOfParallelism();

					vector<std::function<void()>> tasks;
					for (int i = 0; i < complexes.size(); i++)
					{
						auto cplx = complexes.at(i);
						std::function<void()> evolveFunc =
						[=]()
						{
							cplx->Evolve();
						};
						tasks.push_back(evolveFunc);
					}

					CrossThreadExceptions<std::function<void()>> cte(tasks);
					cte.ExecuteTasks(nThreads);
				}
				else
				{
					for (int i = 0; i < complexes.size(); i++)
						complexes.at(i)->Evolve();
				}

				// Optionally add some log information.
			}

			std::vector<IObjectiveScores<T>> AggregateComplexes()
			{
				return complexes.Aggregate();
			}

			void SetPopulationToShuffle(const std::vector<IObjectiveScores<T>>& shufflePoints)
			{
				this->PopulationAtShuffling = sortByFitness(shufflePoints);
			}

			void ShuffleComplexes()
			{
				complexes = shuffle(complexes);
			}

			std::vector<T> CreateCandidates(size_t n)
			{
				if (this->populationInitializer == nullptr)
					this->populationInitializer = candidateFactory->Create();
				return populationInitializer->CreateRandomCandidates(n);
			}

		private:

			void Init(IObjectiveEvaluator<T>* evaluator,
				const ICandidateFactorySeed<T>& candidateFactory,
				const TerminationCondition& terminationCondition,
				bool allowPrematureComplexTermination,
				IRandomNumberGeneratorFactory<> rng,
				IFitnessAssignment<double, T> fitnessAssignment,
				int p = 5,
				int pmin = 5,
				int m = 13,
				int q = 7,
				int alpha = 3,
				int beta = 13,
				int numShuffle = 15,
				double trapezoidalPdfParam = 1.8,
				const std::map<string, string>& logTags = std::map<string, string>(),
				SceOptions options = SceOptions::None, double reflectionRatio = -1.0, double contractionRatio = 0.5)
			{
				if (m < 2)
					throw std::logic_error("M is too small");

				if (q > m)
					throw std::logic_error("Q must be less than or equal to M");

				SetMaxDegreeOfParallelismHardwareMinus(1);

				this->evaluator = evaluator;
				this->candidateFactory = candidateFactory.Clone();
				this->populationInitializer = nullptr;
				this->terminationCondition = terminationCondition;
				//if (this->terminationCondition == nullptr)
				//	this->terminationCondition = new MaxShuffleTerminationCondition();
				this->terminationCondition.SetEvolutionEngine(this);
				this->allowPrematureComplexTermination = allowPrematureComplexTermination;
				this->p = p;
				this->pmin = pmin;
				this->m = m;
				this->q = q;
				this->alpha = alpha;
				this->beta = beta;
				this->numShuffle = numShuffle;
				this->rng = rng;
				//if (this->rng == nullptr)
				//	this->rng = new BasicRngFactory(0);
				this->fitnessAssignment = fitnessAssignment;
				//if (this->fitnessAssignment == nullptr)
				//	this->fitnessAssignment = new DefaultFitnessAssignment();
				this->logTags = logTags;
				this->trapezoidalPdfParam = trapezoidalPdfParam;
				this->options = options;
				this->ReflectionRatio = reflectionRatio;
				this->ContractionRatio = contractionRatio;
			}


			class Complexes
			{
			public:
				ShuffledComplexEvolution<T>* sce = nullptr;
				Complexes(const std::vector<FitnessAssignedScores<double, T>>& sortedScores, ShuffledComplexEvolution<T>* sce)// int p, int m, int shuffleCount, const std::map<string, string>& logTags)
				{
					this->sce = sce;
					int p = sce->p;
					int m = sce->m;
					for (int a = 0; a < p; a++)
					{
						std::vector<FitnessAssignedScores<double, T>> sample;
						for (int k = 1; k <= m; k++)
							sample.push_back(sortedScores[a + p * (k - 1)]);
						std::vector<IObjectiveScores<T>> scores = getScores(sample);
						Complex<T>* complex = createComplex(scores);
						complex->ComplexId = std::to_string(sce->currentShuffle) + "_" + std::to_string(a + 1);
						complexes.push_back(complex);
					}
				}

				Complexes() { }

				Complexes(const Complexes& src)
				{
					throw std::logic_error("deep copy construction of Complexes is not supported");
					this->sce = src.sce;
					this->complexes = src.complexes;
				}

				Complexes(Complexes&& src)
				{
					this->sce = std::move(src.sce);
					src.sce = nullptr;
					std::swap(this->complexes, src.complexes);
				}

				Complexes& operator=(const Complexes& src)
				{
					throw std::logic_error("deep copy construction of Complexes is not supported");
					if (&src == this){
						return *this;
					}
					this->sce = src.sce;
					this->complexes = src.complexes;
					return *this;
				}

				Complexes& operator=(Complexes&& src)
				{
					if (&src == this){
						return *this;
					}
					this->sce = std::move(src.sce);
					src.sce = nullptr;
					std::swap(this->complexes, src.complexes);
					return *this;
				}


				~Complexes()
				{
					DisposeComplexes();
				}

				Complex<T>* at(size_t i) { return complexes.at(i); }

				size_t size() const { return complexes.size(); }

				std::vector<IObjectiveScores<T>> Aggregate()
				{
					std::vector<IObjectiveScores<T>> result;
					for (auto& c : complexes)
					{
						auto scores = c->GetObjectiveScores();
#ifdef _DEBUG
						//CheckParameterFeasible(scores);
#endif
						for (auto& s : scores)
						{
							result.push_back(s);
						}
					}
					return result;
				}

			private:

				void DisposeComplexes()
				{
					for (size_t i = 0; i < complexes.size(); i++)
						if (complexes[i] != nullptr)
						{
							delete complexes[i];
							complexes[i] = nullptr;
						}
					this->sce = nullptr;
				}

				Complex<T>* createComplex(std::vector<IObjectiveScores<T>> scores)
				{
					auto loggerTags = LoggerMhHelper::MergeDictionaries<>(sce->logTags,
						LoggerMhHelper::CreateTag(LoggerMhHelper::MkTuple("CurrentShuffle", std::to_string(sce->currentShuffle))));

					// TODO: reconsider how to handle parallel executions.
					bool ownedPtr = (sce->evaluator->IsCloneable());
					IObjectiveEvaluator<T>* evaluator = (ownedPtr ? sce->evaluator->Clone() : sce->evaluator);

					// We should always create a new candidate factory, 
					// to have reproducible outputs whether multi-threaded or not
					ICandidateFactory<T>* candidateFactory = sce->candidateFactory->Create();
					bool ownedCandidateFactory = true;

					return new Complex<T>(scores, evaluator, ownedPtr, sce->rng,
						candidateFactory, ownedCandidateFactory,
						sce->fitnessAssignment, sce->terminationCondition, sce->allowPrematureComplexTermination,
						sce->logger, loggerTags, sce->q, sce->alpha, sce->beta, sce->trapezoidalPdfParam,
						sce->options, sce->ReflectionRatio, sce->ContractionRatio);
				}

				std::vector<Complex<T>*> complexes;
			};

		protected:
			void MoveFrom(ShuffledComplexEvolution& src)
			{
				this->terminationCondition = std::move(src.terminationCondition);
				this->allowPrematureComplexTermination = src.allowPrematureComplexTermination;
				this->terminationCondition.SetEvolutionEngine(this);
				this->p = std::move(src.p);
				this->pmin = std::move(src.pmin);
				this->m = std::move(src.m);
				this->q = std::move(src.q);
				this->alpha = std::move(src.alpha);
				this->beta = std::move(src.beta);
				this->numShuffle = std::move(src.numShuffle);
				this->options = std::move(src.options);
				this->ReflectionRatio = std::move(src.ReflectionRatio);
				this->ContractionRatio = std::move(src.ContractionRatio);
				this->logTags = std::move(src.logTags);
				this->trapezoidalPdfParam = std::move(src.trapezoidalPdfParam);
				this->fitnessAssignment = std::move(src.fitnessAssignment);
				this->rng = std::move(src.rng);
				this->complexes = std::move(src.complexes);
				this->maxDegreeOfParallelism = std::move(src.maxDegreeOfParallelism);

				this->evaluator = std::move(src.evaluator);
				this->populationInitializer = std::move(src.populationInitializer);
				this->candidateFactory = std::move(src.candidateFactory);
				this->logger = std::move(src.logger);

				src.evaluator = nullptr;
				src.populationInitializer = nullptr;
				src.candidateFactory = nullptr;
				src.logger = nullptr;
			}

			void CopyBasicsFrom(const ShuffledComplexEvolution& src)
			{
				this->terminationCondition = src.terminationCondition;
				//if (this->terminationCondition == src.nullptr)
				//	this->terminationCondition = src.new MaxShuffleTerminationCondition();
				this->terminationCondition.SetEvolutionEngine(this);
				this->maxDegreeOfParallelism = src.maxDegreeOfParallelism;
				this->p = src.p;
				this->pmin = src.pmin;
				this->m = src.m;
				this->q = src.q;
				this->alpha = src.alpha;
				this->beta = src.beta;
				this->numShuffle = src.numShuffle;
				this->options = src.options;
				this->ReflectionRatio = src.ReflectionRatio;
				this->ContractionRatio = src.ContractionRatio;
				this->logTags = src.logTags;
				this->trapezoidalPdfParam = src.trapezoidalPdfParam;
				this->fitnessAssignment = src.fitnessAssignment;
				this->rng = src.rng;
				this->complexes = src.complexes;
				this->logger = src.logger->CreateNew();
			}

		private:

			bool useMultiThreading = true;
			Complexes complexes;
			size_t nbObjectiveEvaluations = 0;
			std::map<string, string> logTags;
			ICandidateFactorySeed<T>* candidateFactory = nullptr;
			IObjectiveEvaluator<T>* evaluator = nullptr;
			ICandidateFactory<T>* populationInitializer = nullptr;
			TerminationCondition terminationCondition;
			bool allowPrematureComplexTermination = true;
			IRandomNumberGeneratorFactory<> rng;
			IFitnessAssignment<double, T> fitnessAssignment;
			double trapezoidalPdfParam;

			ILoggerMh<T>* logger = nullptr;

			int pmin = 5;
			int p = 5, m = 27, q = 14, alpha = 3, beta = 27;
			int numShuffle = -1;

			int seed = 0;

			//class IComplex
			//{
			//public:
			//	//virtual std::vector<IObjectiveScores<T>> GetObjectiveScores() = 0;
			//	//virtual void Evolve() = 0;
			//	virtual std::vector<IObjectiveScores<T>> GetObjectiveScores() const { return std::vector<IObjectiveScores<T>>(); }
			//	virtual void Evolve() { ; }
			//	string ComplexId;
			//	bool IsCancelled;
			//};

			//IObjectiveEvaluator<ISystemConfiguration> evaluator;
			//string fullLogFileName = @"c:\tmp\logMoscem.csv";
			//string paretoLog = @"c:\tmp\logParetoMoscem.csv";

			/*
			class MaxShuffleTerminationCondition : ITerminationCondition<T>
			{
			MaxShuffleTerminationCondition()
			{
			}
			ShuffledComplexEvolution<T> algorithm;
			bool IsFinished()
			{
			return algorithm.currentShuffle >= algorithm.numShuffle;
			}

			#region ITerminationCondition Members

			void SetEvolutionEngine(IEvolutionEngine<T> engine)
			{
			this->algorithm = (ShuffledComplexEvolution<T>)engine;
			}

			#endregion
			}

			class MarginalImprovementTerminationCheck : MaxWalltimeCheck, ITerminationCondition<T>
			{
			MarginalImprovementTerminationCheck(double maxHours, double tolerance, int cutoffNoImprovement)
			: base(maxHours)
			{
			this->tolerance = tolerance;
			this->maxConverge = cutoffNoImprovement;
			}

			IPopulation<double> algorithm;
			void SetEvolutionEngine(IEvolutionEngine<T> engine)
			{
			this->algorithm = (IPopulation<double>)engine;
			}


			double oldBest = double.NaN;
			double tolerance = 1e-6;
			int converge = 0;
			int maxConverge = 10;
			bool IsFinished()
			{
			// https://jira.csiro.au/browse/WIRADA-129
			//current SWIFT SCE implementation uses this algorithm to define convergence and it normally guarantees
			// reproducible optimum is found. It needs two parameters, Tolerance (normally of the order of 10e-6)
			// and maxConverge (normally of the order of 10)

			if (this->HasReachedMaxTime())
			return true;
			std::vector<FitnessAssignedScores<double>> currentPopulation = algorithm.Population;
			if (currentPopulation == nullptr)
			return false;
			auto currentBest = currentPopulation.First().FitnessValue;
			if (double.IsNaN(oldBest))
			{
			oldBest = currentBest;
			return false;
			}
			if (Math.Abs(currentBest - oldBest) <= Math.Abs(oldBest * tolerance))
			{
			converge++;
			}
			else
			{
			converge = 0;
			}
			oldBest = currentBest;
			if (converge > maxConverge) return true;
			return false;
			}
			}

			class CoefficientOfVariationTerminationCondition : ITerminationCondition<T>
			{
			ShuffledComplexEvolution<T> algorithm;
			double threshold;
			double maxHours;
			Stopwatch stopWatch;
			// FIXME: consider something where the termination criteria is customizable to an extent.
			// Func<double[], double> statistic;

			CoefficientOfVariationTerminationCondition(double threshold = 2.5e-2, double maxHours = 1.0)
			{
			this->threshold = threshold;
			this->maxHours = maxHours;
			this->stopWatch = new Stopwatch();
			stopWatch.Start();
			}

			void SetEvolutionEngine(IEvolutionEngine<T> engine)
			{
			this->algorithm = (ShuffledComplexEvolution<T>) engine;
			}

			bool IsFinished()
			{
			if (this->HasReachedMaxTime())
			return true;
			if (algorithm.numShuffle >= 0 && algorithm.currentShuffle >= algorithm.numShuffle)
			return true;
			if (algorithm.PopulationAtShuffling == nullptr)
			return false; // start of the algorithm.
			int n = (int)Math.Ceiling(algorithm.PopulationAtShuffling.Length / 2.0);
			auto tmp = algorithm.PopulationAtShuffling.Where((x, i) = > i < n).ToArray();
			auto popToTest = Array.ConvertAll<FitnessAssignedScores<double>, IObjectiveScores<T>>(tmp, (x = > x.Scores()));
			return IsBelowCvThreshold(popToTest);
			}

			double GetMaxParameterCoeffVar(std::vector<IObjectiveScores<T>> population)
			{
			auto pSets = ConvertAllToHyperCube(population);
			auto varNames = pSets[0].GetVariableNames();
			double[] coeffVar = new double[varNames.Length];
			for (int i = 0; i < varNames.Length; i++)
			{
			coeffVar[i] = calcCoeffVar(MetaheuristicsHelper.GetValues(pSets, varNames[i]));
			}
			return MetaheuristicsHelper.GetMaximum(coeffVar);
			}

			double calcCoeffVar(double[] p)
			{
			double sum = p.Sum();
			double mean = sum / p.Length;
			double[] diffsMean = Array.ConvertAll(p, x = > x - mean);
			double sumSqrDiffs = Array.ConvertAll(diffsMean, x = > x * x).Sum();
			double sdev = Math.Sqrt(sumSqrDiffs / (p.Length - 1));
			if (mean == 0)
			if (sdev == 0) return 0;
			else return double.PositiveInfinity;
			else
			return Math.Abs(sdev / mean);

			}

			bool IsBelowCvThreshold(std::vector<IObjectiveScores<T>> population)
			{
			return GetMaxParameterCoeffVar(population) < threshold;
			}

			bool HasReachedMaxTime()
			{
			double hoursElapsed = this->stopWatch.Elapsed.TotalHours;
			if (this->maxHours <= 0)
			return true;
			else if (this->maxHours < hoursElapsed)
			return true;
			else
			return false;
			}

			double RemainingHours
			{
			get
			{
			return this->maxHours - this->stopWatch.Elapsed.TotalHours;
			}
			}
			}

			class FalseTerminationCondition : ITerminationCondition<T>
			{
			void SetEvolutionEngine(IEvolutionEngine<T> engine)
			{
			// Nothing
			}

			bool IsFinished()
			{
			return false;
			}
			}
			
			class MaxWalltimeTerminationCondition : MaxWalltimeCheck, ITerminationCondition<T>
			{
			MaxWalltimeTerminationCondition(double maxHours) : base(maxHours)
			{
			}

			virtual void SetEvolutionEngine(IEvolutionEngine<T> engine)
			{
			// Nothing
			}
			bool IsFinished()
			{
			return this->HasReachedMaxTime();
			}
			}


			static ITerminationCondition<T> CreateMaxShuffleTerminationCondition()
			{
			return new MaxShuffleTerminationCondition();
			}


			static T[] ConvertAllToHyperCube(std::vector<IObjectiveScores<T>> population)
			{
			auto tmp = Array.ConvertAll<IObjectiveScores<T>, T>(population, (x = > x.GetSystemConfiguration()));
			return tmp;
			}

			*/

			size_t currentShuffle = 0;

			size_t maxDegreeOfParallelism = 1;

			bool isCancelled = false;
			//IComplex currentComplex;

			//CancellationTokenSource tokenSource = new CancellationTokenSource();
			SceOptions options = SceOptions::None;

			double ContractionRatio;
			double ReflectionRatio;

			void loggerWrite(const std::vector<IObjectiveScores<T>>& scores, const std::map<string, string>& ctags)
			{
				auto tags = LoggerMhHelper::MergeDictionaries<>(this->logTags, ctags);
				LoggerMhHelper::Write(scores, tags, this->logger);
			}

			void loggerWrite(const FitnessAssignedScores<double, T>& scores, const std::map<string, string>& ctags)
			{
				auto tags = LoggerMhHelper::MergeDictionaries<>(this->logTags, ctags);
				LoggerMhHelper::Write(scores, tags, this->logger);
			}

			void loggerWrite(const string& msg, const std::map<string, string>& ctags)
			{
				auto tags = LoggerMhHelper::MergeDictionaries<>(this->logTags, ctags);
				LoggerMhHelper::Write(msg, tags, this->logger);
			}

			void Cancel()
			{
				isCancelled = true;
				//if (currentComplex != nullptr)
				//	currentComplex.IsCancelled = isCancelled;
				//tokenSource.Cancel();
			}

			static IOptimizationResults<T> packageResults(Complexes& complexes)
			{
				//saveLog( logPopulation, fullLogFileName );
				std::vector<IObjectiveScores<T>> population = complexes.Aggregate();
				//saveLogParetoFront( population );
				return packageResults(population);
			}

			static IOptimizationResults<T> packageResults(const std::vector<IObjectiveScores<T>>& population)
			{
				// cater for cases where we have null references (e.g. if the termination condition was in the middle of the population creation)
				// return new BasicOptimizationResults<T>(population.Where(p = > (p != nullptr)).ToArray());
				return IOptimizationResults<T>(population);
			}

			void logTerminationConditionMet()
			{
				auto tags = createSimpleMsg("Termination condition", "Termination condition");
				loggerWrite(string("Termination condition using ") + typeid(terminationCondition).name() + " is met", tags);
			}

			string GetDescription()
			{
				throw std::logic_error("Not implemented");
			}

			Complexes shuffle(Complexes& complexes)
			{
				std::vector<IObjectiveScores<T>> population = complexes.Aggregate();
				auto newComplexes = partition(population);
				return newComplexes;
			}

			//std::vector<IObjectiveScores<T>> evaluateScores(IObjectiveEvaluator<T>* evaluator, const std::vector<T>& population)
			//{
			//	return Evaluations::EvaluateScores(evaluator, population);
			//}

			std::vector<T> initialisePopulation()
			{
				return populationInitializer->CreateRandomCandidates(this->PopulationSize());
			}

			Complexes partition(const std::vector<FitnessAssignedScores<double, T>>& sortedScores)
			{
				if (currentShuffle > 0)
					if (this->pmin < this->p)
						this->p = this->p - 1;
				Complexes c = Complexes(sortedScores, this);
				return c;
			}

			// https://github.com/jmp75/metaheuristics/issues/3
			// ITerminationCondition<T> createMaxWalltimeCondition(ITerminationCondition<T> terminationCondition)
			// {
			// 	auto t = terminationCondition as CoefficientOfVariationTerminationCondition;
			// 	if (t == nullptr)
			// 		return new FalseTerminationCondition();
			// 	else
			// 		return new MaxWalltimeTerminationCondition(t.RemainingHours);
			// }

			Complexes partition(const std::vector<IObjectiveScores<T>>& scores)
			{
				auto sortedScores = sortByFitness(scores);
				//logPoints( currentShuffle, sortedScores );
				auto complexes = partition(sortedScores);
				return complexes;
			}

			std::vector<FitnessAssignedScores<double, T>> sortByFitness(const std::vector<IObjectiveScores<T>>& scores)
			{
				auto fittedScores = fitnessAssignment.AssignFitness(scores);
				FitnessAssignedScores<double, T>::Sort(fittedScores);
				//std::sort(fittedScores.begin(), fittedScores.end());
				return fittedScores;
			}

			static std::vector<IObjectiveScores<T>> getScores(const std::vector<FitnessAssignedScores<double, T>>& fitnessedScores)
			{
				std::vector<IObjectiveScores<T>> result;
				for (int i = 0; i < fitnessedScores.size(); i++)
					result.push_back(fitnessedScores[i].Scores());
				return result;
			}

			/*
			class ComplexEvolutionEvent : EventArgs, IMonitoringEvent
			{
			IObjectiveScores[] scoresSet;

			ComplexEvolutionEvent( std::vector<IComplex> complexes )
			{
			List<IObjectiveScores> list = new List<IObjectiveScores>( );

			foreach( IComplex complex in complexes )
			{
			foreach( IObjectiveScores scores in complex.GetObjectiveScores( ) )
			{
			list.Add( scores );
			}
			}

			this->scoresSet = list.ToArray( );
			}

			IObjectiveScores[] ScoresSet
			{
			get { return scoresSet; }
			}
			}

			*/

			std::vector<FitnessAssignedScores<double, T>> PopulationAtShuffling;

			std::map<string, string> createSimpleMsg(string message, string category)
			{
				return LoggerMhHelper::CreateTag({
					LoggerMhHelper::MkTuple("Message", message),
					LoggerMhHelper::MkTuple("Category", category) });
			}

			public: 
				int CurrentShuffle() { return currentShuffle; }
				// Method to facilitate testing: given the current state of the optimizer, create a population
				std::vector<T> CreatePopulation() const
				{
					auto cf = candidateFactory->Clone();
					auto gen = cf->Create();
					std::vector<T> result = gen->CreateRandomCandidates(this->PopulationSize());
					delete cf;
					delete gen;
					return result;
				}

		};

		template<typename TSys, typename TEngine = IEvolutionEngine<TSys>>
		class MaxWalltimeCheck : public TerminationCheck<TSys, TEngine>
		{
		public:
			MaxWalltimeCheck(double maxHours)
			{
				this->maxHours = maxHours;
				Reset();
			}

			virtual ~MaxWalltimeCheck()
			{
			}

			bool HasReachedMaxTime() const
			{
				if (maxHours <= 0)
					return false;
				boost::posix_time::ptime currentTime(boost::posix_time::second_clock::local_time());
				double secondsElapsed = (currentTime - startTime).total_seconds();
				double hoursElapsed = secondsElapsed / 3600.0;
				return (hoursElapsed >= maxHours);
			}

			virtual void Reset()
			{
				startTime = boost::posix_time::second_clock::local_time();
			}

			virtual bool IsFinished(TEngine* engine)
			{
				return HasReachedMaxTime();
			}

			virtual bool IsThreadSafe() { return true; };

			virtual TerminationCheck<TSys, TEngine>* Clone() const
			{
				// TOCHECK: is this the behavior we want (think parallel operations)
				auto result = new MaxWalltimeCheck(maxHours);
				result->startTime = this->startTime;
				return result;
			}


		protected:
			double maxHours;
			boost::posix_time::ptime startTime;
		};

		template<typename TSys, typename TEngine = IEvolutionEngine<TSys>>
		class MaxNumberSceShuffles : public TerminationCheck<TSys, TEngine>
		{
		public:
			MaxNumberSceShuffles(int maxShuffles)
			{
				this->maxShuffles = maxShuffles;
				Reset();
			}

			virtual ~MaxNumberSceShuffles()
			{
			}

			virtual void Reset()
			{
				// nothing
			}

			virtual bool IsFinished(TEngine* engine)
			{
				return (engine->CurrentShuffle() >= maxShuffles);
			}

			virtual bool IsThreadSafe() { return true; };

			virtual TerminationCheck<TSys, TEngine>* Clone() const
			{
				auto result = new MaxNumberSceShuffles(this->maxShuffles);
				return result;
			}
		protected:
			int maxShuffles;
		};

		template<typename TSys, typename TEngine = IEvolutionEngine<TSys>>
		class MarginalImprovementTerminationCheck :
			public MaxWalltimeCheck<TSys, TEngine>
		{
		public:
			MarginalImprovementTerminationCheck(double maxHours, double tolerance, int cutoffNoImprovement) :
				MaxWalltimeCheck<TSys, TEngine>(maxHours)
			{
				this->tolerance = tolerance;
				this->maxConverge = cutoffNoImprovement;
			}

			virtual void Reset()
			{
				MaxWalltimeCheck<TSys, TEngine>::Reset();
				converge = 0;
			}

			virtual bool IsFinished(TEngine* engine)
			{
				return isFinished(engine);
			}

			virtual bool IsThreadSafe()
			{
				// When complexes check the SCE population, there are some race conditions, perhaps due to 
				// use of move semantics in some classes. Too tricky to handle, and this check 
				// should only be used in thread safe contexts.
				return false;
			}

			virtual TerminationCheck<TSys, TEngine>* Clone() const
			{
				// TOCHECK: is this the behavior we want (think parallel operations)
				auto result = new MarginalImprovementTerminationCheck(this->maxHours, tolerance, maxConverge);
				result->Reset();
				return result;
			}

		protected:
			bool isFinished(TEngine* engine)
			{
				if (engine == nullptr) throw std::invalid_argument(string("Argument must not be nullptr"));
				IPopulation<double, TSys>* sce = dynamic_cast<IPopulation<double, TSys>*>(engine);
				if (sce == nullptr) throw std::invalid_argument(string("Argument 'engine' is not a ") + typeid(IPopulation<double, TSys>).name());

				if (MaxWalltimeCheck<TSys, TEngine>::HasReachedMaxTime())
					return true;
				auto currentPopulation = sce->Population();
				if (currentPopulation.size() == 0)
					return false;
				double currentBest = currentPopulation[0].FitnessValue();
				if (std::isnan(oldBest))
				{
					oldBest = currentBest;
					return false;
				}
				if (std::abs(currentBest - oldBest) <= std::abs(oldBest * tolerance))
				{
					converge++;
				}
				else
				{
					converge = 0;
				}
				oldBest = currentBest;
				if (converge > maxConverge)
					return true;
				return false;
			}
			//std::function<bool(TEngine*)> CreateNew(MarginalImprovementTerminationCheck& mitc)
			//{
			//	return [&mitc](TEngine* e)
			//	{
			//		return mitc.IsFinished(e);
			//	};
			//}
		private:
			double tolerance;
			int maxConverge;

			double oldBest = std::numeric_limits<double>::quiet_NaN();
			int converge = 0;
		};


		template<typename TSys, typename TEngine = IEvolutionEngine<TSys>>
		class MaxIterationTerminationCheck :
			public TerminationCheck<TSys, TEngine>
		{
		private:
			size_t maxIterations = 0;
		public:
			MaxIterationTerminationCheck(size_t maxIterations) :
				TerminationCheck<TSys, TEngine>()
			{
				this->maxIterations = maxIterations;
			}

			virtual void Reset()
			{
				// Nothing done
			}

			virtual bool IsFinished(TEngine* engine)
			{
				return isFinished(engine);
			}

			virtual bool IsThreadSafe()
			{
				return true;
			}

			virtual TerminationCheck<TSys, TEngine>* Clone() const
			{
				return new MaxIterationTerminationCheck<TSys, TEngine>(this->maxIterations);
			}

		protected:
			bool isFinished(TEngine* engine)
			{
				return (this->maxIterations <= engine->EvaluationCount());
			}
		};
	}
}

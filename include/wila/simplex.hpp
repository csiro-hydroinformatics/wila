#pragma once

#include <algorithm>
#include <atomic>
#include <thread>
#include <mutex>
#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/threadpool.hpp"
//#include "sce.h"
#include "core.hpp"
#include "utils.hpp"
#include "evaluations.hpp"
#include "multithreading.hpp"
#include "logging.hpp"


#ifdef _MSC_VER
#include <concurrent_vector.h>
#else
#include <tbb/concurrent_vector.h>
#endif

#ifdef _MSC_VER
using namespace Concurrency;
#else
using namespace tbb;
#endif

namespace mhcpp
{
	namespace optimization
	{
		using namespace mhcpp;
		using namespace mhcpp::logging;
		using namespace mhcpp::objectives;

		template<typename T>
		class SimplexPopulation
		{
		private:
			size_t nbObjectiveEvaluations = 0;

			void IncrementEvaluations(size_t n)
			{
				nbObjectiveEvaluations += n;
			}

		private:
			std::vector<FitnessAssignedScores<double, T>> evolved;
			SceOptions options;
			IObjectiveEvaluator<T>* evaluator;
			IFitnessAssignment<double, T> fitnessAssignment;
			double ContractionRatio;
			double ReflectionRatio;
			int alpha, q;
			//IRandomNumberGeneratorFactory<> rng;
			ICandidateFactory<T>* cf = nullptr;
			std::map<string, string> tags;
			ILoggerMh<T>* logger = nullptr;

			void Init(const std::vector<IObjectiveScores<T>>& complexPopulation, IObjectiveEvaluator<T>* evaluator, 
				ICandidateFactory<T> * candidateFactory,
				IFitnessAssignment<double, T> fitnessAssignment, ILoggerMh<T>* logger, const std::map<string, string>& tags, SceOptions options,
				double contractionRatio, double reflectionRatio)
			{
				this->options = options;
				this->evaluator = evaluator;
				this->fitnessAssignment = fitnessAssignment;
				this->ContractionRatio = contractionRatio;
				this->ReflectionRatio = reflectionRatio;
				//this->rng = rng;
				this->cf = candidateFactory;
				//if (cf == nullptr)
				//	cf = new UniformRandomSamplingFactory<T>(rng, complexPopulation[0].SystemConfiguration());
				//this->discreteRng = &discreteRng;
				this->tags = tags;
				this->logger = logger;
#ifdef _DEBUG
				//CheckParameterFeasible(complexPopulation);
#endif
				evolved = getSimplexPopulation(complexPopulation);
#ifdef _DEBUG
				//CheckParameterFeasible(leftOutFromSubcomplex);
#endif
			}

		public:

			SimplexPopulation(const std::vector<IObjectiveScores<T>>& complexPopulation, IObjectiveEvaluator<T>* evaluator, 
				ICandidateFactory<T> * candidateFactory, 
				IFitnessAssignment<double, T> fitnessAssignment, ILoggerMh<T>* logger = nullptr, const std::map<string, string>& tags = std::map<string, string>(), SceOptions options = SceOptions::RndInSimplexPopulation)
			{
				Init(complexPopulation, evaluator, q, alpha, candidateFactory, fitnessAssignment, logger, tags, options, contractionRatio, reflectionRatio);
			}

			~SimplexPopulation()
			{
				// Nothing;
			}

			//bool IsCancelledOrFinished()
			//{
			//	if (complex != nullptr)
			//		return this->complex->IsCancelledOrFinished();
			//	return false;
			//}

			std::map<string, string> createTagConcat(const std::initializer_list<std::tuple<string, string>>& tuples)
			{
				return LoggerMhHelper::MergeDictionaries<>(LoggerMhHelper::CreateTag(tuples), this->tags);
			}

			size_t EvaluationCount()
			{
				return nbObjectiveEvaluations;
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

			void clear(std::vector<FitnessAssignedScores<double, T>>& vec)
			{
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
				if ((options & SceOptions::RndInSimplexPopulation) == SceOptions::RndInSimplexPopulation)
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
				LoggerMhHelper::Write<double, T>(scores, tags, this->logger);
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
				std::vector<IObjectiveScores<T>> newSimplexPopulation = aggregate(newScore, withoutWorstPoint);
#ifdef _DEBUG
				//CheckParameterFeasible(newSimplexPopulation);
#endif
				return fitnessAssignment.AssignFitness(newSimplexPopulation);
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
		class SimplexNelderMead
		{
		private:
			size_t nbObjectiveEvaluations = 0;

			void IncrementEvaluations(size_t n)
			{
				nbObjectiveEvaluations += n;
			}
			//std::vector<IObjectiveScores<T>> leftOutFromSimplexNelderMead;
			std::vector<FitnessAssignedScores<double, T>> evolved;
			SceOptions options;
			IObjectiveEvaluator<T>* evaluator;
			IFitnessAssignment<double, T> fitnessAssignment;
			double ContractionRatio;
			double ReflectionRatio;
			IRandomNumberGeneratorFactory<> rng;
			ICandidateFactory<T>* cf = nullptr;
			std::map<string, string> tags;
			ILoggerMh<T>* logger = nullptr;

		public:

			typedef ITerminationCondition<T, SimplexNelderMead<T>> TerminationCondition;

			SimplexNelderMead(const IObjectiveEvaluator<T>& evaluator,
				const ICandidateFactorySeed<T>& candidateFactory,
				const TerminationCondition& terminationCondition,
				const IRandomNumberGeneratorFactory<>& rng,
				IFitnessAssignment<double, T> fitnessAssignment,
				double reflectionRatio = -1.0,
				double expansionRatio = +1.5,
				double contractionRatio = 0.5,
				//ILoggerMh<T>* logger = nullptr,
				const std::map<string, string>& tags = std::map<string, string>(), SceOptions options = SceOptions::RndInSimplexPopulation)
			{
				this->options = options;
				if (!evaluator.IsCloneable()) throw std::logic_error("objective evaluator must be cloneable if calling this simplex constructor");
				this->evaluator = evaluator.Clone();
				this->fitnessAssignment = fitnessAssignment;
				this->ContractionRatio = contractionRatio;
				this->ReflectionRatio = reflectionRatio;
				//this->alpha = alpha;
				//this->q = q;
				this->rng = rng;
				this->cf = candidateFactory;
				if (cf == nullptr)
					cf = new UniformRandomSamplingFactory<T>(rng, simplexPopulation[0].SystemConfiguration());
				this->discreteRng = &discreteRng;
				this->tags = tags;
				this->logger = logger;
				evolved = getSimplexNelderMead(simplexPopulation, leftOutFromSimplexNelderMead);
			}

			~SimplexNelderMead()
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
				while (!IsCancelledOrFinished())
				{
					std::vector<FitnessAssignedScores<double, T>> withoutWorstPoint;
					FitnessAssignedScores<double, T> worstPoint = FindWorstPoint(evolved, withoutWorstPoint);
					loggerWrite(worstPoint, createTagConcat({
						LoggerMhHelper::MkTuple("Message", "Worst point in simplex"),
						createTagCatComplexNo() }));
					loggerWrite(withoutWorstPoint, createTagConcat({
						LoggerMhHelper::MkTuple("Message", "simplex without worst point"),
						createTagCatComplexNo() }
					));
					T centroid = getCentroid(withoutWorstPoint);

					bool success;
					//Compute reflected point x r = x o + alpha(x o - x n + 1) 
					T reflectedPoint = reflect(worstPoint, centroid, success);
					if (success)
					{
						std::vector<FitnessAssignedScores<double, T>> candidateWithReflection;
						// TODO: if we are using multi-objective Pareto optimisation, the following two lines should be re-scrutinized. 
						// While the behavior in the C# implementation was sensical, there may be some smarter ways to test the 
						// improvement in the fitness. Also a valid remark in single-objective case.
						FitnessAssignedScores<double, T> fitReflectedPoint = evaluateNewSet(reflectedPoint, withoutWorstPoint, candidateWithReflection);
						//	If the reflected point is better than the second worst, 
						if (fitReflectedPoint.CompareTo(worstPoint) <= 0)
						{
							FitnessAssignedScores<double, T> bestPoint = withoutWorstPoint[0];
							if (fitReflectedPoint.CompareTo(bestPoint) <= 0)
							{
								//This looks promising - try an expansion in the direction of the reflected point.
								bool expandSuccess;
								T expandedPoint = expand(reflectedPoint, centroid, expandSuccess);
								if (expandSuccess)
								{
									std::vector<FitnessAssignedScores<double, T>> candidateWithExpansion;
									FitnessAssignedScores<double, T> fitExpandedPoint = evaluateNewSet(expandedPoint, withoutWorstPoint, candidateWithExpansion);


									// CAUTION: the following check is pretty obvious for single objective, but 
									//for the case of a fitness derived from multi-objective dominance, since the dominance is relative to a population, 
									//  I am not sure the direct comparison is neceessarily valid. we may need to instead calculate a fitness against the following instead:
									// std::vector<FitnessAssignedScores<double, T>> dummy;
									// FitnessAssignedScores<double, T> dummyFitEval = evaluateNewSet(expandedPoint, candidateWithReflection, tmp);
									// In practice we only need single objective so sill park the above 
									// There should be a check and exception thrown for MO for now, but not readily implemented.
									if (fitExpandedPoint.CompareTo(fitReflectedPoint) <= 0)
									{
										replaceEvolved(candidateWithExpansion);
										loggerWrite(fitReflectedPoint, createTagConcat({
											LoggerMhHelper::MkTuple("Message", "Expansion point in simplex"),
											createTagCatComplexNo() }));
									}
								}
								else
								{
									// Expansion did not outcompete reflection - stick with reflection
									replaceEvolved(candidateWithReflection);
									loggerWrite(fitReflectedPoint, createTagConcat({
										LoggerMhHelper::MkTuple("Message", "Reflected point in simplex"),
										createTagCatComplexNo() }));
								}
							}
							else
							{
								// the reflected point is not better than the best, i.e.f(x 1) = f(x r) < f(x n)
									//	then obtain a new simplex by replacing the worst point x n + 1 
								replaceEvolved(candidateWithReflection);
								loggerWrite(fitReflectedPoint, createTagConcat({
									LoggerMhHelper::MkTuple("Message", "Reflected point in simplex"),
									createTagCatComplexNo() }));
							}
						}
						else
						{
							loggerWrite(fitReflectedPoint,
								createTagConcat({ LoggerMhHelper::MkTuple("Message", "Reflected point in simplex - Failed"), createTagCatComplexNo() }));
							auto candidateContraction = contractionOrRandom(withoutWorstPoint, worstPoint, centroid);
							//if (candidateSimplexNelderMead == nullptr) // this can happen if the feasible region of the parameter space is not convex.
							//	candidateSimplexNelderMead = fitnessAssignment.AssignFitness(bufferComplex);
							replaceEvolved(candidateContraction);
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
				auto scores = aggregatePoints(evolved, leftOutFromSimplexNelderMead);
#ifdef _DEBUG
				//CheckParameterFeasible(scores);
#endif
				return scores;
			}

			static FitnessAssignedScores<double, T> FindWorstPoint(std::vector<FitnessAssignedScores<double, T>>& SimplexNelderMead, std::vector<FitnessAssignedScores<double, T>>& pointRemoved)
			{
				auto tmp = AsPointers<FitnessAssignedScores<double, T>>(SimplexNelderMead);
				FitnessAssignedScores<double, T>::Sort(tmp);
				auto worst = tmp[tmp.size() - 1];
				pointRemoved.clear();
				for (auto& point : SimplexNelderMead)
				{
					if (&point != worst)
						pointRemoved.push_back(point);
				}
				return *worst;
			}

		private:

			void clear(std::vector<FitnessAssignedScores<double, T>>& vec)
			{
				//for (auto ptr : vec)
				//	delete ptr;
				vec.clear();
			}

			void replaceEvolved(const std::vector<FitnessAssignedScores<double, T>>& candidateSimplexNelderMead)
			{
				clear(evolved);
				evolved = candidateSimplexNelderMead;
			}

			std::vector<FitnessAssignedScores<double, T>> replaceWithRandom(const std::vector<FitnessAssignedScores<double, T>>& withoutWorstPoint, const FitnessAssignedScores<double, T>& worstPoint)
			{
				std::vector<FitnessAssignedScores<double, T>> result;
				if ((options & SceOptions::RndInSimplexNelderMead) == SceOptions::RndInSimplexNelderMead)
				{
					// TODO: what was SceOptions::ReflectionRandomization about? Needs re-clarification to assess whether worth porting.
					//						if ((options & SceOptions::ReflectionRandomization) == SceOptions::ReflectionRandomization)
					result = generateRandomWithinSimplexNelderMead(withoutWorstPoint, worstPoint);
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
				LoggerMhHelper::Write<double, T>(scores, tags, this->logger);
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

			std::vector<IObjectiveScores<T>> aggregatePoints(const std::vector<FitnessAssignedScores<double, T>>& SimplexNelderMead, const std::vector<IObjectiveScores<T>>& leftOutFromSimplexNelderMead)
			{
				std::vector<IObjectiveScores<T>> result = FitnessAssignedScores<double, T>::GetScores(SimplexNelderMead);
				for (size_t i = 0; i < leftOutFromSimplexNelderMead.size(); i++)
				{
					result.push_back(leftOutFromSimplexNelderMead[i]);
				}
				return result;
			}

			FitnessAssignedScores<double, T> evaluateNewSet(T newPoint, const std::vector<FitnessAssignedScores<double, T>>& withoutWorstPoint, std::vector<FitnessAssignedScores<double, T>>& candidateSimplexNelderMead)
			{
				IObjectiveScores<T> scoreNewPoint = evaluator.EvaluateScore(newPoint);
				IncrementEvaluations(1);

				std::vector<IObjectiveScores<T>> scores = FitnessAssignedScores<double, T>::GetScores(withoutWorstPoint);
				scores.push_back(scoreNewPoint);
				auto fitness = fitnessAssignment.AssignFitness(scores);
				candidateSimplexNelderMead.clear();
				for (size_t i = 0; i < fitness.size(); i++)
				{
					candidateSimplexNelderMead.push_back(FitnessAssignedScores<double, T>(fitness[i]));
				}
				//return Array.Find<FitnessAssignedScores<double, T>>(candidateSimplexNelderMead, (x = > (x.Scores() == scoreNewPoint)));
				return (candidateSimplexNelderMead[candidateSimplexNelderMead.size() - 1]);
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

			std::vector<FitnessAssignedScores<double, T>> removePoint(const std::vector<FitnessAssignedScores<double, T>>& SimplexNelderMead, FitnessAssignedScores<double, T> worstPoint)
			{
				std::vector<FitnessAssignedScores<double, T>> result;
				for (size_t i = 0; i < SimplexNelderMead.size(); i++)
				{
					auto p = SimplexNelderMead.at(i);
					if (p != worstPoint)
						result.push_back(p);
				}
				return result;
			}

			std::vector<FitnessAssignedScores<double, T>> contractionOrRandom(const std::vector<FitnessAssignedScores<double, T>>& withoutWorstPoint,
				const FitnessAssignedScores<double, T>& worstPoint, const T& centroid)
			{
				std::vector<FitnessAssignedScores<double, T>> result;
				std::vector<FitnessAssignedScores<double, T>> candidateSimplexNelderMead;
				bool success;
				T contractionPoint = contract(worstPoint, centroid, success);

				if (success)
				{
					FitnessAssignedScores<double, T> trialPoint = evaluateNewSet(contractionPoint, withoutWorstPoint, candidateSimplexNelderMead);
					if (trialPoint.CompareTo(worstPoint) <= 0)
					{
						result = candidateSimplexNelderMead;
						loggerWrite(trialPoint, createTagConcat({
							LoggerMhHelper::MkTuple("Message", "Contracted point in SimplexNelderMead"),
							createTagCatComplexNo() }));
						return result;
					}
					else
					{
						clear(candidateSimplexNelderMead);
						loggerWrite(trialPoint, createTagConcat({
							LoggerMhHelper::MkTuple("Message", "Contracted point in SimplexNelderMead-Failed"),
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
							auto newScore = evaluator.EvaluateScore(newPoint);
							loggerWrite(newScore, createTagConcat(
								LoggerMhHelper::MkTuple("Message", "Adding a partially random point"),
								LoggerMhHelper::MkTuple("Category", "Complex No " + complexId)
								));
							return fitnessAssignment.AssignFitness(aggregate(newScore, withoutWorstPoint));
						}
			*/
			std::vector<FitnessAssignedScores<double, T>> generateRandomWithinSimplexNelderMead(const std::vector<FitnessAssignedScores<double, T>>& withoutWorstPoint, const FitnessAssignedScores<double, T>& worstPoint)
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
				auto newScore = evaluator.EvaluateScore(newPoint);
				IncrementEvaluations(1);
				loggerWrite(newScore, createTagConcat({
					LoggerMhHelper::MkTuple("Message", "Adding a random point in hypercube"),
					createTagCatComplexNo() }
				));
				std::vector<IObjectiveScores<T>> newSimplexNelderMead = aggregate(newScore, withoutWorstPoint);
#ifdef _DEBUG
				//CheckParameterFeasible(newSimplexNelderMead);
#endif
				return fitnessAssignment.AssignFitness(newSimplexNelderMead);
			}

			std::vector<FitnessAssignedScores<double, T>> getSimplexNelderMead(const std::vector<IObjectiveScores<T>>& bufferComplex, std::vector<IObjectiveScores<T>>& leftOutFromSimplexNelderMead)
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
				leftOutFromSimplexNelderMead.clear();
				for (size_t i = 0; i < leftOut.size(); i++)
				{
					leftOutFromSimplexNelderMead.push_back(leftOut[i].Scores());
				}
#ifdef _DEBUG
				//CheckParameterFeasible(leftOutFromSimplexNelderMead);
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
				return aggregate(evaluator.EvaluateScore(newPoint), withoutWorstPoint);
			}

			std::vector<IObjectiveScores<T>> aggregate(const IObjectiveScores<T>& newPoint, std::vector<FitnessAssignedScores<double, T>> withoutWorstPoint)
			{
				std::vector<IObjectiveScores<T>> result = convertToScores(withoutWorstPoint);
				result.push_back(newPoint);
				return result;
			}
		};

	}
}

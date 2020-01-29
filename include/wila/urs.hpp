#pragma once

#include <algorithm>
#include <atomic>
#include <thread>
#include <mutex>
#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/threadpool.hpp"
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
		using namespace mhcpp::logging;
		using namespace mhcpp::objectives;


		template<typename T>
		class UniformRandomSamplingOptimizer
			: public IEvolutionEngine<T>,
			public IPopulation<double, T>
			//where T : ICloneableSystemConfiguration
		{
			// needed? friend Complexes<T>(SCE&);

		public:

			typedef ITerminationCondition<T, UniformRandomSamplingOptimizer<T>> TerminationCondition;

			UniformRandomSamplingOptimizer(IObjectiveEvaluator<T>* evaluator,
				const ICandidateFactorySeed<T>& candidateFactory,
				const TerminationCondition& terminationCondition,
				size_t populationSize = 100, 
				IRandomNumberGeneratorFactory<> rng = IRandomNumberGeneratorFactory<>(),
				IFitnessAssignment<double, T> fitnessAssignment = IFitnessAssignment<double, T>(),
				const std::map<string, string>& logTags = std::map<string, string>())
			{
				Init(evaluator, candidateFactory, terminationCondition,
					rng,
					fitnessAssignment,
					populationSize,
					logTags);
			}

			UniformRandomSamplingOptimizer(const UniformRandomSamplingOptimizer& src)
			{
				throw std::logic_error("copy constructor for UniformRandomSamplingOptimizer not yet supported");
			}

			UniformRandomSamplingOptimizer(UniformRandomSamplingOptimizer&& src)
			{
				MoveFrom(src);
			}

			UniformRandomSamplingOptimizer& operator=(UniformRandomSamplingOptimizer&& src)
			{
				if (&src == this) {
					return *this;
				}
				MoveFrom(src);
				return *this;
			}

			UniformRandomSamplingOptimizer& operator=(const UniformRandomSamplingOptimizer& src)
			{
				if (&src == this) {
					return *this;
				}
				throw std::logic_error("copy assignment for UniformRandomSamplingOptimizer not yet supported");
				// CopyBasicsFrom(src);
				// this->evaluator = src.evaluator;
				// this->populationInitializer = src.populationInitializer;
				// TODO 
				// ILoggerMh<TSys>* Logger;
				return *this;
			}

			virtual ~UniformRandomSamplingOptimizer()
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

			std::vector<IObjectiveScores<T>> best;

			std::vector<FitnessAssignedScores<double, T>> sortByFitness(const std::vector<IObjectiveScores<T>>& scores)
			{
				auto fittedScores = fitnessAssignment.AssignFitness(scores);
				FitnessAssignedScores<double, T>::Sort(fittedScores);
				//std::sort(fittedScores.begin(), fittedScores.end());
				return fittedScores;
			}

			std::vector<FitnessAssignedScores<double, T>> Population()
			{
				return sortByFitness(best);
			}

			IOptimizationResults<T> Evolve()
			{
				Reset();

				this->populationInitializer = candidateFactory->Create();

				bool isFinished = false;
				while (!isFinished)
				{
					currentRound = 1;
					auto candidates = populationInitializer->CreateRandomCandidates(this->PopulationSize());
					std::vector<IObjectiveScores<T>> scores = EvaluateScores(candidates);
					loggerWrite(scores, createSimpleMsg("Round " + to_string(currentRound), "Round " + to_string(currentRound)));
					IncrementEvaluations(scores.size());
					for (auto& s : scores)
					{
						best.push_back(s);
					}
					std::vector<FitnessAssignedScores<double, T>> fittedScores = sortByFitness(best);
					best.clear();
					for (size_t i = 0; i < std::min(this->PopulationSize(), fittedScores.size()); i++)
					{
						best.push_back(fittedScores.at(i).Scores());
					}
					currentRound++;
					isFinished = terminationCondition.IsFinished();
					if (isFinished)
					{
						logTerminationConditionMet();
						return packageResults(best);
					}
				}
			}

			size_t PopulationSize() const
			{
				return this->popSize;
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

			CrossThreadExceptions<std::function<void()>> cte;

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
				IRandomNumberGeneratorFactory<> rng,
				IFitnessAssignment<double, T> fitnessAssignment,
				size_t populationSize = 100,
				const std::map<string, string>& logTags = std::map<string, string>())
			{
				this->popSize = populationSize;
				size_t mt = mhcpp::threading::ThreadingOptions<double>::DefaultMaxDegreeOfParallelism;
				if (mt > 0)
					SetMaxDegreeOfParallelism(mt);
				else
					SetMaxDegreeOfParallelismHardwareMinus(1);

				this->evaluator = evaluator;
				this->candidateFactory = candidateFactory.Clone();
				this->populationInitializer = nullptr;
				this->terminationCondition = terminationCondition;
				this->terminationCondition.SetEvolutionEngine(this);
				this->rng = rng;
				this->fitnessAssignment = fitnessAssignment;
				this->logTags = logTags;
			}

		protected:
			void MoveFrom(UniformRandomSamplingOptimizer& src)
			{
				this->terminationCondition = std::move(src.terminationCondition);
				this->terminationCondition.SetEvolutionEngine(this);
				this->logTags = std::move(src.logTags);
				this->fitnessAssignment = std::move(src.fitnessAssignment);
				this->rng = std::move(src.rng);
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

			void CopyBasicsFrom(const UniformRandomSamplingOptimizer& src)
			{
				this->terminationCondition = src.terminationCondition;
				this->terminationCondition.SetEvolutionEngine(this);
				this->maxDegreeOfParallelism = src.maxDegreeOfParallelism;
				this->logTags = src.logTags;
				this->fitnessAssignment = src.fitnessAssignment;
				this->rng = src.rng;
				this->logger = src.logger->CreateNew();
			}

		private:

			bool useMultiThreading = true;
			size_t nbObjectiveEvaluations = 0;
			std::map<string, string> logTags;
			ICandidateFactorySeed<T>* candidateFactory = nullptr;
			IObjectiveEvaluator<T>* evaluator = nullptr;
			ICandidateFactory<T>* populationInitializer = nullptr;
			TerminationCondition terminationCondition;
			IRandomNumberGeneratorFactory<> rng;
			IFitnessAssignment<double, T> fitnessAssignment;

			ILoggerMh<T>* logger = nullptr;

			size_t popSize = 100;
			int seed = 0;
			size_t currentRound = 0;

			size_t maxDegreeOfParallelism = 1;

			bool isCancelled = false;

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

			std::vector<T> initialisePopulation()
			{
				return populationInitializer->CreateRandomCandidates(this->PopulationSize());
			}

			std::map<string, string> createSimpleMsg(string message, string category)
			{
				return LoggerMhHelper::CreateTag({
					LoggerMhHelper::MkTuple("Message", message),
					LoggerMhHelper::MkTuple("Category", category) });
			}

		public:
			int CurrentRound() { return currentRound; }
			// Method to facilitate testing: given the current state of the optimizer, create a population
			//std::vector<T> CreatePopulation() const
			//{
			//	auto cf = candidateFactory->Clone();
			//	auto gen = cf->Create();
			//	std::vector<T> result = gen->CreateRandomCandidates(this->PopulationSize());
			//	delete cf;
			//	delete gen;
			//	return result;
			//}

		};
	}
}

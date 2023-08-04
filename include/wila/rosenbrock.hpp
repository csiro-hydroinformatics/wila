#pragma once

#include <algorithm>
#include <atomic>
#include <thread>
#include <mutex>
#include "boost/date_time/posix_time/posix_time.hpp"
// #include "boost/threadpool.hpp"
#include "core.hpp"
#include "utils.hpp"
#include "evaluations.hpp"
// #include "multithreading.hpp"
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
		class Complex // : IComplex
		{
			std::map<string, string> createTagConcat(const std::initializer_list<std::tuple<string, string>>& tuples)
			{
				return LoggerMhHelper::MergeDictionaries<>(LoggerMhHelper::CreateTag(tuples), this->tags);
			}

			// friend SubComplex<T>::SubComplex(Complex&);
			friend SubComplex<T>;

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
		class RosenbrockOptimizer
			: public IEvolutionEngine<T>
		{

		public:

			typedef ITerminationCondition<T, RosenbrockOptimizer<T>> TerminationCondition;

			RosenbrockOptimizer(const IObjectiveEvaluator<T>& evaluator,
				const T& startingPoint,
				const TerminationCondition& terminationCondition,
				double alpha = 1.4,
				double beta = 0.7,
				IAlgebraProvider algebraprovider = nullptr,
				IFitnessAssignment<double, T> fitnessAssignment = IFitnessAssignment<double, T>(),
				const std::map<string, string>& logTags = std::map<string, string>())
			{
				IObjectiveEvaluator<T>* pEval = evaluator.Clone();
				this->countingEvaluator = new CountingEvaluator(evaluator);
				this->startingPoint = evaluator.EvaluateScore( startingPoint );
				this->terminationCondition = terminationCondition;
				terminationCondition.SetEvolutionEngine( this );
				this->alpha = alpha;
				this->beta = beta;
				this->AlgebraProvider = algebraprovider;
				this->logTags = logTags;
			}
			RosenbrockOptimizer(const RosenbrockOptimizer& src)
			{
				throw std::logic_error("copy constructor for RosenbrockOptimizer not yet supported");
				// CopyBasicsFrom(src);
				// this->evaluator = src.evaluator;
				// this->candidatefactory etc.
				// this->populationInitializer = src.populationInitializer;
				// TODO 
				// ILoggerMh<TSys>* Logger;
			}

			RosenbrockOptimizer(RosenbrockOptimizer&& src)
			{
				MoveFrom(src);
			}

			RosenbrockOptimizer& operator=(RosenbrockOptimizer&& src)
			{
				if (&src == this) {
					return *this;
				}
				MoveFrom(src);
				return *this;
			}

			RosenbrockOptimizer& operator=(const RosenbrockOptimizer& src)
			{
				if (&src == this) {
					return *this;
				}
				throw std::logic_error("copy assignment for RosenbrockOptimizer not yet supported");
				return *this;
			}

			virtual ~RosenbrockOptimizer()
			{
				if (countingEvaluator != nullptr)
				{
					delete countingEvaluator;
					countingEvaluator = nullptr;
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

				if (useMultiThreading && evaluator->IsCloneable())
				{
					int nThreads = this->GetMaxDegreeOfParallelism();
					cte.PoolSize(nThreads);
				}

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
			}

			std::vector<IObjectiveScores<T>> EvaluateScores(const std::vector<T>& population)
			{
				return Evaluations::EvaluateScores(this->evaluator, population, this->useMultiThreading, GetMaxDegreeOfParallelism());
			}

			void CreateComplexes(const std::vector<IObjectiveScores<T>>& scores)
			{
				this->complexes = partition(scores);
			}

			CrossThreadExceptions<std::function<void()>> cte;

			void EvolveComplexes()
			{
				for (int i = 0; i < complexes.size(); i++)
					complexes.at(i)->ComplexId = std::to_string(i);
				if (useMultiThreading && evaluator->IsCloneable())
				{
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
					// call to PoolSize as candidate fix for https://jira.csiro.au/browse/WIRADA-476
					int nThreads = this->GetMaxDegreeOfParallelism();
					cte.PoolSize(nThreads);
					cte.ExecuteTasks(tasks);
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
				double trapezoidalPdfParam = 1.0,
				const std::map<string, string>& logTags = std::map<string, string>(),
				SceOptions options = SceOptions::None, double reflectionRatio = -1.0, double contractionRatio = 0.5)
			{
				if (m < 2)
					throw std::logic_error("M is too small");

				if (q > m)
					throw std::logic_error("Q must be less than or equal to M");

				size_t mt = mhcpp::threading::ThreadingOptions<double>::DefaultMaxDegreeOfParallelism;
				if (mt > 0)
					SetMaxDegreeOfParallelism(mt);
				else
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

		protected:
			void MoveFrom(RosenbrockOptimizer& src)
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

			void CopyBasicsFrom(const RosenbrockOptimizer& src)
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

			size_t nbObjectiveEvaluations = 0;
			std::map<string, string> logTags;
			CountingEvaluator<T>* countingEvaluator = nullptr;
			TerminationCondition terminationCondition;
			bool allowPrematureComplexTermination = true;
			IRandomNumberGeneratorFactory<> rng;
			IFitnessAssignment<double, T> fitnessAssignment;
			double trapezoidalPdfParam;

			ILoggerMh<T>* logger = nullptr;

			size_t maxDegreeOfParallelism = 1;

			bool isCancelled = false;

			//CancellationTokenSource tokenSource = new CancellationTokenSource();

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
		class PopulationStdDevTerminationCheck :
			public MaxWalltimeCheck<TSys, TEngine>
		{
		public:
			/**
			* \brief	Termination criteria based on the spread of parameter values in the a population
			*
			* \param 	maxRelativeStdDev	threshold of maximum relative standard deviation ( max_over_parameters(std dev / (max-min))) below which convergence is considered achieved
			* \param 	values	maxHours  A fallback maximum wall time allowed for runtime
			*/
			PopulationStdDevTerminationCheck(double maxRelativeStdDev, double maxHours = 0.1) :
				MaxWalltimeCheck<TSys, TEngine>(maxHours)
			{
				this->maxRelativeStdDev = maxRelativeStdDev;
			}

			virtual void Reset()
			{
				MaxWalltimeCheck<TSys, TEngine>::Reset();
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
				auto result = new PopulationStdDevTerminationCheck(this->maxRelativeStdDev, this->maxHours);
				result->Reset();
				return result;
			}

		protected:
			bool isFinished(TEngine* engine)
			{
				if (engine == nullptr) throw std::invalid_argument(string("Argument must not be nullptr"));
				IPopulation<double, TSys>* optimAlgorithm = dynamic_cast<IPopulation<double, TSys>*>(engine);
				if (optimAlgorithm == nullptr) throw std::invalid_argument(string("Argument 'engine' is not a ") + typeid(IPopulation<double, TSys>).name());

				if (MaxWalltimeCheck<TSys, TEngine>::HasReachedMaxTime())
					return true;
				auto currentPopulation = 
					FitnessAssignedScores<double, TSys>::GetSystemConfigurations(optimAlgorithm->Population());
				if (currentPopulation.size() < 2)
					return false;
				std::map<string, double> rsdevs = mhcpp::GetRelativeSdev(currentPopulation);
				rsdevs = mhcpp::utils::RemoveNotFinite(rsdevs);

				auto v = mhcpp::utils::GetValues(rsdevs);
				double maxsd = *std::max_element(v.begin(), v.end());
				return (maxsd <= this->maxRelativeStdDev);
			}
		private:
			double maxRelativeStdDev;
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

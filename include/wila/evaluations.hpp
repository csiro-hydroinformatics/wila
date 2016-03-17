#pragma once

#include <vector>
#include <boost/threadpool.hpp>
#include <thread>
#include "core.hpp"
#include "boost/threadpool.hpp"

namespace mhcpp
{
	namespace objectives
	{
		class Evaluations
		{
		private:
			Evaluations() {}

			template<typename T>
			class EvaluateScoresTask
			{
			public:
				EvaluateScoresTask(const IObjectiveEvaluator<T>& evaluator, std::vector<T>& population) : population(population)
				{
					this->evaluator = evaluator.Clone();
				}

				std::vector<T>& population;
				IObjectiveEvaluator<T>* evaluator = nullptr;
				std::vector<IObjectiveScores<T>> result;
				~EvaluateScoresTask()
				{
					if (evaluator != nullptr)
					{
						delete evaluator;
						evaluator = nullptr;
					}
				}
				void EvaluateScores()
				{
					for (size_t i = 0; i < population.size(); i++)
						result.push_back(evaluator->EvaluateScore(population[i]));
				}
			};

		public:
			~Evaluations() {}

			template<typename T>
			static vector<vector<T>> MakeBins(vector<T> population, int numBins)
			{
				int div = population.size() / numBins;
				int remainder = population.size() % numBins;
				vector<vector<T>> result;
				int offset = 0;
				for (int i = 0; i < numBins; i++)
				{
					int len = (i < remainder ? div + 1 : div);
					result.push_back(vector<T>(len));
					std::copy(population.begin() + offset, population.begin() + offset + len, result[i].begin());
					offset += len;
				}
				return result;
			}

			template<typename T>
			static std::vector<IObjectiveScores<T>> EvaluateScores(IObjectiveEvaluator<T>* evaluator, const std::vector<T>& population, bool useMultiThreading, int nThreads)
			{
				std::vector<IObjectiveScores<T>> result;

				if (useMultiThreading && evaluator->IsCloneable())
				{
					int nParallel = nThreads;
					vector<vector<T>> subPop = MakeBins(population, nParallel);

					vector<EvaluateScoresTask<T>*> taskPkgs;
					vector<std::function<void()>> tasks;
					for (size_t i = 0; i < subPop.size(); i++)
					{
						auto eval = new EvaluateScoresTask<T> (*evaluator, subPop.at(i));
						taskPkgs.push_back(eval);
						std::function<void()> evalFunc =
							[=]()
						{
							eval->EvaluateScores();
						};
						tasks.push_back(evalFunc);
					}

					CrossThreadExceptions<std::function<void()>> cte(tasks);
					vector<exception_ptr> exceptions;
					try {
						cte.ExecuteTasks(nParallel);
						for (size_t i = 0; i < taskPkgs.size(); i++)
						{
							auto r = taskPkgs[i]->result;
							for (size_t j = 0; j < r.size(); j++)
								result.push_back(r[j]);
						}
					}
					catch (const std::exception& e) {
						exceptions.push_back(current_exception());
					}
					for (size_t i = 0; i < taskPkgs.size(); i++)
						delete taskPkgs[i];

					if(exceptions.size() > 0)
						rethrow_exception(exceptions[0]);
				}
				else
				{
					for (size_t i = 0; i < population.size(); i++)
						result.push_back(evaluator->EvaluateScore(population[i]));
				}

				return result;
			}
		};
	}
}

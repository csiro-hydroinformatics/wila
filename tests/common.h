#pragma once

#include <random>

#include "wila/core.hpp"
#include "wila/sce.hpp"
#include "wila/urs.hpp"

using namespace mhcpp;
using namespace mhcpp::optimization;

using Hc = HyperCube<double>;
using Sce = ShuffledComplexEvolution<Hc>;
using SceTc = Sce::TerminationCondition;
using Urs = UniformRandomSamplingOptimizer<Hc>;
using UrsTc = Urs::TerminationCondition;


void BuildTestHc(Hc& goal);
Hc CreateTestHc();
CandidateFactorySeed<Hc> CreateCandidateFactorySeed(unsigned int seed, const Hc& hc);

Hc createTestHc(double a, double b, double aMin = 1, double bMin = 3, double aMax = 2, double bMax = 4);
bool assertHyperCube(const Hc& hc, double a, double b, double tolerance = 1.0e-9);
bool isIn(const std::string& s, const std::vector<std::string>& b);
bool allIn(const std::vector<std::string>& a, const std::vector<std::string>& b);
bool sameSets(const std::vector<std::string>& a, const std::vector<std::string>& b);
bool sameValues(const std::vector<std::string>& keys, const Hc& a, const Hc& b);
bool sameMinValues(const std::vector<std::string>& keys, const Hc& a, const Hc& b);
bool sameMaxValues(const std::vector<std::string>& keys, const Hc& a, const Hc& b);
bool assertEqual(const Hc& a, const Hc& b);
bool assertValuesNotEqual(const Hc& a, const Hc& b);

template <typename T=Sce>
typename T::TerminationCondition CreateMaxNumShuffle(int maxCount)
{
	MaxNumberSceShuffles<Hc, T> c(maxCount);
	typename T::TerminationCondition terminationCondition(c);
	return terminationCondition;
}

template <typename T = Sce>
typename T::TerminationCondition CreateCounterTermination(int maxCount)
{
	CounterTestFinished<Hc, T> c(maxCount);
	typename T::TerminationCondition terminationCondition(c);
	return terminationCondition;
}

SceTc CreateMaxIterationTermination(int maxIterations);
SceTc CreateStdDevTermination(double maxRelativeStdDev, double maxHours = 0.1);
Sce CreateSceQuadraticGoal(Hc& goal, const ITerminationCondition<HyperCube < double >, ShuffledComplexEvolution<Hc > >&terminationCondition);
Urs CreateUrsQuadraticGoal(Hc& goal, const ITerminationCondition<Hc, Urs>& terminationCondition);

template < typename THC = Hc >
typename ShuffledComplexEvolution<THC>::TerminationCondition CreateWallClockTermination(double seconds)
// The 'typename' keyword above is necessary at least for visual studio
// for the return type to be recognised as a type.
{
	MaxWalltimeCheck<THC, ShuffledComplexEvolution<THC>> c(seconds / 3600);
	return typename ShuffledComplexEvolution<THC>::TerminationCondition(c);
}

class ThrowsException : public IObjectiveEvaluator < Hc >
{
public:
	ThrowsException(const ThrowsException& src);
	ThrowsException(const Hc& goal);
	~ThrowsException();

	IObjectiveScores<Hc> EvaluateScore(const Hc& systemConfiguration);
	bool IsCloneable() const;
	IObjectiveEvaluator<Hc> * Clone() const;
private:
	Hc goal;
};

Sce CreateSceQuadraticGoalThrowsException(Hc& goal, const ITerminationCondition<HyperCube < double >, Sce>&terminationCondition);
Sce* CreateSceQuadraticGoalPtr(Hc& goal, const ITerminationCondition<HyperCube < double >, Sce>&terminationCondition);

template<typename T>
bool requireEqual(const vector<T>& a, const vector<T>& b)
{
	if (a.size() != b.size()) return false;
	for (size_t i = 0; i < a.size(); i++)
	{
		if (a[i] != b[i]) return false;
	}
	return true;
}

template<typename T>
bool requireRelativeEqual(const vector<T>& expected, const vector<T>& b, T tolerance)
{
	if (expected.size() != b.size()) return false;
	for (size_t i = 0; i < expected.size(); i++)
	{
		if ((std::abs(expected[i] - b[i]) / expected[i]) > tolerance) return false;
	}
	return true;
}


typedef default_random_engine test_random_engine;

class TestRandomEngine : 
	public test_random_engine,
	public InstanceCounter<TestRandomEngine>
{
public:
	TestRandomEngine() : test_random_engine()
	{
	}

	TestRandomEngine(unsigned int seed) : test_random_engine(seed)
	{
	}

	// The cast  to (const test_random_engine&) is important - will not compile otherwise.
	TestRandomEngine(const TestRandomEngine& eng) : test_random_engine((const test_random_engine&)eng), InstanceCounter<TestRandomEngine>(eng)
	{
	}

	TestRandomEngine(TestRandomEngine& eng) : test_random_engine((test_random_engine&)eng)
	{	
	}

};

template<typename T>
UniformRandomSamplingFactory<T> createTestUnifrand(int seed = 0)
{
	IRandomNumberGeneratorFactory<> rng(seed);
	Hc goal = createTestHc(1.5, 3.4);
	UniformRandomSamplingFactory<T> prand(rng, goal);
	return prand;
}

template<typename T>
std::vector < IObjectiveScores<T> >createTestScores(int m, int seed = 0)
{
	std::vector < IObjectiveScores<T> > scores;
	auto prand = createTestUnifrand<T>();
	Hc goal = createTestHc(1.5, 3.4);
	TopologicalDistance<T> evaluator(goal);

	for (size_t i = 0; i < m; i++)
		scores.push_back(evaluator.EvaluateScore(prand.CreateRandomCandidate()));

	return scores;
}


template<typename T>
class ThreeDimArray // surely some utilities out there? quicker for now, but look out for .
{
public:
	int a, b, c;
	ThreeDimArray(int a, int b, int c)
	{
		this->a = a;
		this->b = b;
		this->c = c;

		Values = new T**[a];
		for (int i = 0; i < a; i++)
		{
			Values[i] = new T*[b];
			for (int j = 0; j < b; j++)
			{
				Values[i][j] = new T[c];
			}
		}
	}
	T*** Values;

	bool IsPairwiseDistinct()
	{
		for (int i_1 = 0; i_1 < a; i_1++)
		{
			for (int i_2 = 0; i_2 < a; i_2++)
			{
				if (i_2 == i_1) continue;
				for (int j_1 = 0; j_1 < b; j_1++)
					for (int j_2 = 0; j_2 < b; j_2++)
					{
						if (j_2 == j_1) continue;
						if (!AreDistinct(Values[i_1][j_1], Values[i_2][j_2]))
							return false;
					}
			}
		}
		return true;
	}

	string FailedPairwiseDistinct()
	{
		for (int i_1 = 0; i_1 < a; i_1++)
		{
			for (int i_2 = 0; i_2 < a; i_2++)
			{
				if (i_2 == i_1) continue;
				for (int j_1 = 0; j_1 < b; j_1++)
					for (int j_2 = 0; j_2 < b; j_2++)
					{
						if (j_2 == j_1) continue;
						if (!AreDistinct(Values[i_1][j_1], Values[i_2][j_2]))
						{
							return
								string("Sequences are not distinct for [") +
								boost::lexical_cast<string>(i_1) +
								string("][") +
								boost::lexical_cast<string>(j_1) +
								string("]") +

								string(" and [") +
								boost::lexical_cast<string>(i_2) +
								string("][") +
								boost::lexical_cast<string>(j_2) +
								string("]");
						}
					}
			}
		}
		return "";
	}

	bool AreDistinct(T* v1, T* v2)
	{
		bool result = false;
		for (int i = 0; i < c; i++)
		{
			if (v2[i] != v1[i])
				result = true;
		}
		return result;
	}

	bool Equals(ThreeDimArray<T>& other)
	{
		bool diffDim = (other.a != a) || (other.b != b) || (other.c != c);
		if (diffDim) return false;
		for (int i = 0; i < a; i++)
			for (int j = 0; j < b; j++)
				for (int k = 0; k < c; k++)
				{
					if (Values[i][j][k] != other.Values[i][j][k])
						return false;
				}
		return true;

	}
	~ThreeDimArray()
	{
		for (int i = 0; i < a; i++)
		{
			for (int j = 0; j < b; j++)
			{
				delete[] Values[i][j];
			}
			delete[] Values[i];
		}
		delete[] Values;
	}
};


ThreeDimArray<int> * CreateRandValues(int seed = 0, int numFactories = 2, int numEngines = 2, int numDraw = 10);

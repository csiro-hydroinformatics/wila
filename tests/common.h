#pragma once

#include <random>

#include "wila/core.hpp"
#include "wila/sce.hpp"

using namespace mhcpp;
using namespace mhcpp::optimization;

HyperCube<double> createTestHc(double a, double b, double aMin = 1, double bMin = 3, double aMax = 2, double bMax = 4);
bool assertHyperCube(const HyperCube<double>& hc, double a, double b, double tolerance = 1.0e-9);
bool isIn(const std::string& s, const std::vector<std::string>& b);
bool allIn(const std::vector<std::string>& a, const std::vector<std::string>& b);
bool sameSets(const std::vector<std::string>& a, const std::vector<std::string>& b);
bool sameValues(const std::vector<std::string>& keys, const HyperCube<double>& a, const HyperCube<double>& b);
bool sameMinValues(const std::vector<std::string>& keys, const HyperCube<double>& a, const HyperCube<double>& b);
bool sameMaxValues(const std::vector<std::string>& keys, const HyperCube<double>& a, const HyperCube<double>& b);
bool assertEqual(const HyperCube<double>& a, const HyperCube<double>& b);
bool assertValuesNotEqual(const HyperCube<double>& a, const HyperCube<double>& b);

ShuffledComplexEvolution<HyperCube<double>>::TerminationCondition CreateCounterTermination(int maxCount);
ShuffledComplexEvolution<HyperCube<double>>::TerminationCondition CreateMaxIterationTermination(int maxIterations);
ShuffledComplexEvolution<HyperCube<double>> CreateQuadraticGoal(HyperCube<double>& goal, const ITerminationCondition<HyperCube < double >, ShuffledComplexEvolution<HyperCube<double> > >&terminationCondition);

template < typename THC = HyperCube<double> >
typename ShuffledComplexEvolution<THC>::TerminationCondition CreateWallClockTermination(double seconds)
// The 'typename' keyword above is necessary at least for visual studio
// for the return type to be recognised as a type.
{
	MaxWalltimeCheck<THC, ShuffledComplexEvolution<THC>> c(seconds / 3600);
	return typename ShuffledComplexEvolution<THC>::TerminationCondition(c);
}

class ThrowsException : public IObjectiveEvaluator < HyperCube<double> >
{
public:
	ThrowsException(const ThrowsException& src);
	ThrowsException(const HyperCube<double>& goal);
	~ThrowsException();

	IObjectiveScores<HyperCube<double>> EvaluateScore(const HyperCube<double>& systemConfiguration);
	bool IsCloneable() const;
	IObjectiveEvaluator<HyperCube<double>> * Clone() const;
private:
	HyperCube<double> goal;
};

ShuffledComplexEvolution<HyperCube<double>> CreateQuadraticGoalThrowsException(HyperCube<double>& goal, const ITerminationCondition<HyperCube < double >, ShuffledComplexEvolution<HyperCube<double>>>&terminationCondition);
ShuffledComplexEvolution<HyperCube<double>>* CreateQuadraticGoalPtr(HyperCube<double>& goal, const ITerminationCondition<HyperCube < double >, ShuffledComplexEvolution<HyperCube<double>>>&terminationCondition);

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
	HyperCube<double> goal = createTestHc(1.5, 3.4);
	UniformRandomSamplingFactory<T> prand(rng, goal);
	return prand;
}

template<typename T>
std::vector < IObjectiveScores<T> >createTestScores(int m, int seed = 0)
{
	std::vector < IObjectiveScores<T> > scores;
	auto prand = createTestUnifrand<T>();
	HyperCube<double> goal = createTestHc(1.5, 3.4);
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

#pragma once

#include "wila/core.hpp"
#include "wila/sce.hpp"

using namespace mhcpp;

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

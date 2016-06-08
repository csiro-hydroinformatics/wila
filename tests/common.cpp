#include "common.h"

using namespace mhcpp;
using namespace mhcpp::optimization;

HyperCube<double> createTestHc(double a, double b, double aMin, double bMin, double aMax, double bMax) {
	HyperCube<double> hc;
	hc.Define("a", aMin, aMax, a);
	hc.Define("b", bMin, bMax, b);
	return hc;
}

bool assertHyperCube(const HyperCube<double>& hc, double a, double b, double tolerance)
{
	return
		(hc.Dimensions() == 2) &&
		(std::abs(hc.GetValue("a") - a) < tolerance) &&
		(std::abs(hc.GetValue("b") - b) < tolerance);
}

bool isIn(const std::string& s, const std::vector<std::string>& b)
{
	return (std::find(b.begin(), b.end(), s) != b.end());
}

bool allIn(const std::vector<std::string>& a, const std::vector<std::string>& b)
{
	for (auto s : a)
	{
		if (!isIn(s, b))
			return false;
	}
	return true;
}

bool sameSets(const std::vector<std::string>& a, const std::vector<std::string>& b)
{
	if (a.size() != b.size()) return false;
	if (!allIn(a, b)) return false;
	if (!allIn(b, a)) return false;
	return true;
}

bool sameValues(const std::vector<std::string>& keys, const HyperCube<double>& a, const HyperCube<double>& b)
{
	for (auto& s : keys)
		if (a.GetValue(s) != b.GetValue(s))
			return false;
	return true;
}

bool sameMinValues(const std::vector<std::string>& keys, const HyperCube<double>& a, const HyperCube<double>& b)
{
	for (auto& s : keys)
		if (a.GetMinValue(s) != b.GetMinValue(s))
			return false;
	return true;
}

bool sameMaxValues(const std::vector<std::string>& keys, const HyperCube<double>& a, const HyperCube<double>& b)
{
	for (auto& s : keys)
		if (a.GetMaxValue(s) != b.GetMaxValue(s))
			return false;
	return true;
}

bool assertEqual(const HyperCube<double>& a, const HyperCube<double>& b)
{
	return
		(a.Dimensions() == b.Dimensions()) &&
		(sameSets(a.GetVariableNames(), b.GetVariableNames())) &&
		(sameValues(a.GetVariableNames(), a, b)) &&
		(sameMinValues(a.GetVariableNames(), a, b)) &&
		(sameMaxValues(a.GetVariableNames(), a, b));
}

bool assertValuesNotEqual(const HyperCube<double>& a, const HyperCube<double>& b)
{
	return
		(a.Dimensions() == b.Dimensions()) &&
		(sameSets(a.GetVariableNames(), b.GetVariableNames())) &&
		(!sameValues(a.GetVariableNames(), a, b));
}

ShuffledComplexEvolution<HyperCube<double>>::TerminationCondition CreateCounterTermination(int maxCount)
{
	CounterTestFinished<HyperCube<double>, ShuffledComplexEvolution<HyperCube<double>>> c(maxCount);
	ShuffledComplexEvolution<HyperCube<double>>::TerminationCondition terminationCondition(c);
	return terminationCondition;
}

ShuffledComplexEvolution<HyperCube<double>>::TerminationCondition CreateMaxIterationTermination(int maxIterations)
{
	MaxIterationTerminationCheck<HyperCube<double>, ShuffledComplexEvolution<HyperCube<double>>> c(maxIterations);
	return ShuffledComplexEvolution<HyperCube<double>>::TerminationCondition(c);
}

ShuffledComplexEvolution<HyperCube<double>> CreateQuadraticGoal(HyperCube<double>& goal, const ITerminationCondition<HyperCube < double >, ShuffledComplexEvolution<HyperCube<double>>>&terminationCondition)
{
	SceParameters sceParams = CreateSceParamsForProblemOfDimension(5, 20);
	// TODO: check above
	sceParams.P = 5;
	sceParams.Pmin = 3;

	goal.Define("a", 1, 2, 1);
	goal.Define("b", 3, 4, 3);
	HyperCube<double> hc = goal;
	TopologicalDistance<HyperCube < double > >  * evaluator = new TopologicalDistance<HyperCube < double > >(goal);
	CandidateFactorySeed<HyperCube < double >> seeding(0, hc);

	ShuffledComplexEvolution<HyperCube<double>> opt(evaluator, seeding, terminationCondition, sceParams);
	return opt;
}

ThrowsException::ThrowsException(const ThrowsException& src)
{
	this->goal = src.goal;
}
ThrowsException::ThrowsException(const HyperCube<double>& goal) { this->goal = goal; }
ThrowsException::~ThrowsException() {}

IObjectiveScores<HyperCube<double>> ThrowsException::EvaluateScore(const HyperCube<double>& systemConfiguration)
{
	throw std::runtime_error("An error occured in ThrowsException::EvaluateScore");
}

bool ThrowsException::IsCloneable() const
{
	return true;
}

IObjectiveEvaluator<HyperCube<double>> * ThrowsException::Clone() const
{
	return new ThrowsException(*this);
}


ShuffledComplexEvolution<HyperCube<double>> CreateQuadraticGoalThrowsException(HyperCube<double>& goal, const ITerminationCondition<HyperCube < double >, ShuffledComplexEvolution<HyperCube<double>>>&terminationCondition)
{
	SceParameters sceParams = CreateSceParamsForProblemOfDimension(5, 20);
	// TODO: check above
	sceParams.P = 5;
	sceParams.Pmin = 3;
	goal.Define("a", 1, 2, 1);
	goal.Define("b", 3, 4, 3);
	HyperCube<double> hc = goal;
	ThrowsException* evaluator = new ThrowsException(goal);
	CandidateFactorySeed<HyperCube < double >> seeding(0, hc);

	ShuffledComplexEvolution<HyperCube<double>> opt(evaluator, seeding, terminationCondition, sceParams);
	return opt;
}

ShuffledComplexEvolution<HyperCube<double>>* CreateQuadraticGoalPtr(HyperCube<double>& goal, const ITerminationCondition<HyperCube < double >, ShuffledComplexEvolution<HyperCube<double>>>&terminationCondition)
{
	SceParameters sceParams = CreateSceParamsForProblemOfDimension(5, 20);
	// TODO: check above
	sceParams.P = 5;
	sceParams.Pmin = 3;

	goal.Define("a", 1, 2, 1);
	goal.Define("b", 3, 4, 3);
	HyperCube<double> hc = goal;
	TopologicalDistance<HyperCube < double > >  * evaluator = new TopologicalDistance<HyperCube < double > >(goal);
	CandidateFactorySeed<HyperCube < double >> seeding(0, hc);

	return new ShuffledComplexEvolution<HyperCube<double>>(evaluator, seeding, terminationCondition, sceParams);
}

ThreeDimArray<int> * CreateRandValues(int seed, int numFactories, int numEngines, int numDraw)
{
	ThreeDimArray<int> * result = new ThreeDimArray<int>(numFactories, numEngines, numDraw);
	IRandomNumberGeneratorFactory<> rngf(seed);

	for (int i = 0; i < numFactories; i++)
	{
		auto rngf_2 = rngf.CreateNew();
		for (int j = 0; j < numEngines; j++)
		{
			auto e = rngf_2.CreateNewStd();
			for (int k = 0; k < numDraw; k++)
			{
				result->Values[i][j][k] = e();
			}
		}
	}
	return result;
}

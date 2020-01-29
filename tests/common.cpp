#include "common.h"

using namespace mhcpp;
using namespace mhcpp::optimization;

Hc createTestHc(double a, double b, double aMin, double bMin, double aMax, double bMax) {
	Hc hc;
	hc.Define("a", aMin, aMax, a);
	hc.Define("b", bMin, bMax, b);
	return hc;
}

bool assertHyperCube(const Hc& hc, double a, double b, double tolerance)
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

bool sameValues(const std::vector<std::string>& keys, const Hc& a, const Hc& b)
{
	for (auto& s : keys)
		if (a.GetValue(s) != b.GetValue(s))
			return false;
	return true;
}

bool sameMinValues(const std::vector<std::string>& keys, const Hc& a, const Hc& b)
{
	for (auto& s : keys)
		if (a.GetMinValue(s) != b.GetMinValue(s))
			return false;
	return true;
}

bool sameMaxValues(const std::vector<std::string>& keys, const Hc& a, const Hc& b)
{
	for (auto& s : keys)
		if (a.GetMaxValue(s) != b.GetMaxValue(s))
			return false;
	return true;
}

bool assertEqual(const Hc& a, const Hc& b)
{
	return
		(a.Dimensions() == b.Dimensions()) &&
		(sameSets(a.GetVariableNames(), b.GetVariableNames())) &&
		(sameValues(a.GetVariableNames(), a, b)) &&
		(sameMinValues(a.GetVariableNames(), a, b)) &&
		(sameMaxValues(a.GetVariableNames(), a, b));
}

bool assertValuesNotEqual(const Hc& a, const Hc& b)
{
	return
		(a.Dimensions() == b.Dimensions()) &&
		(sameSets(a.GetVariableNames(), b.GetVariableNames())) &&
		(!sameValues(a.GetVariableNames(), a, b));
}

//SceTc CreateCounterTermination(int maxCount)
//{
//	CounterTestFinished<Hc, Sce> c(maxCount);
//	SceTc terminationCondition(c);
//	return terminationCondition;
//}

//SceTc CreateMaxNumShuffle(int maxCount)
//{
//	MaxNumberSceShuffles<Hc, Sce> c(maxCount);
//	SceTc terminationCondition(c);
//	return terminationCondition;
//}

SceTc CreateMaxIterationTermination(int maxIterations)
{
	MaxIterationTerminationCheck<Hc, Sce> c(maxIterations);
	return SceTc(c);
}

SceTc CreateStdDevTermination(double maxRelativeStdDev, double maxHours)
{
	PopulationStdDevTerminationCheck<Hc, Sce> c(maxRelativeStdDev, maxHours);
	return SceTc(c);
}

void BuildTestHc(Hc& goal)
{
	goal.Define("a", 1, 2, 1);
	goal.Define("b", 3, 4, 3);
}

Hc CreateTestHc()
{
	Hc goal;
	BuildTestHc(goal);
	return goal;
}

CandidateFactorySeed<Hc> CreateCandidateFactorySeed(unsigned int seed, const Hc& hc)
{
	CandidateFactorySeed<Hc> seeding(seed, hc);
	return seeding;
}

Sce CreateSceQuadraticGoal(Hc& goal, const ITerminationCondition<Hc, Sce>& terminationCondition)
{
	SceParameters sceParams = CreateSceParamsForProblemOfDimension<SceParameters>(5, 20);
	// TODO: check above
	sceParams.P = 5;
	sceParams.Pmin = 3;
	BuildTestHc(goal);
	Hc hc = goal;
	TopologicalDistance<Hc>* evaluator = new TopologicalDistance<Hc>(goal);
	CandidateFactorySeed<Hc> seeding(0, hc);

	Sce opt(evaluator, seeding, terminationCondition, sceParams);
	return opt;
}

Urs CreateUrsQuadraticGoal(Hc& goal, const ITerminationCondition<Hc, Urs>& terminationCondition)
{
	BuildTestHc(goal);
	Hc hc = goal;
	TopologicalDistance<Hc>* evaluator = new TopologicalDistance<Hc>(goal);
	CandidateFactorySeed<Hc> seeding(0, hc);

	Urs opt(evaluator, seeding, terminationCondition);
	return opt;
}

ThrowsException::ThrowsException(const ThrowsException& src)
{
	this->goal = src.goal;
}
ThrowsException::ThrowsException(const Hc& goal) { this->goal = goal; }
ThrowsException::~ThrowsException() {}

IObjectiveScores<Hc> ThrowsException::EvaluateScore(const Hc& systemConfiguration)
{
	throw std::runtime_error("An error occured in ThrowsException::EvaluateScore");
}

bool ThrowsException::IsCloneable() const
{
	return true;
}

IObjectiveEvaluator<Hc> * ThrowsException::Clone() const
{
	return new ThrowsException(*this);
}

Sce CreateSceQuadraticGoalThrowsException(Hc& goal, const ITerminationCondition<Hc, Sce>&terminationCondition)
{
	SceParameters sceParams = CreateSceParamsForProblemOfDimension<SceParameters>(5, 20);
	// TODO: check above
	sceParams.P = 5;
	sceParams.Pmin = 3;
	goal.Define("a", 1, 2, 1);
	goal.Define("b", 3, 4, 3);
	Hc hc = goal;
	ThrowsException* evaluator = new ThrowsException(goal);
	CandidateFactorySeed<Hc> seeding(0, hc);

	Sce opt(evaluator, seeding, terminationCondition, sceParams);
	return opt;
}

Sce* CreateSceQuadraticGoalPtr(Hc& goal, const ITerminationCondition<Hc, Sce>&terminationCondition)
{
	SceParameters sceParams = CreateSceParamsForProblemOfDimension<SceParameters>(5, 20);
	// TODO: check above
	sceParams.P = 5;
	sceParams.Pmin = 3;

	goal.Define("a", 1, 2, 1);
	goal.Define("b", 3, 4, 3);
	Hc hc = goal;
	TopologicalDistance<Hc >  * evaluator = new TopologicalDistance<Hc >(goal);
	CandidateFactorySeed<Hc> seeding(0, hc);

	return new Sce(evaluator, seeding, terminationCondition, sceParams);
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

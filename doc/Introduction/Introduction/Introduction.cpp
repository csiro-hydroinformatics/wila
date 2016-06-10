
#include "wila/core.hpp"
#include "wila/sce.hpp"


using namespace mhcpp;
using namespace mhcpp::optimization;

int main()
{
	// If you have to remember *one* thing from this tutorial, it is the three essential players 
	// 1- the tunable parameters of the optimisation problem, here a hypercube
	// 2- the object evaluating the score(s) for a given "system configuration" (IObjectiveEvaluator<T>)
	// 3- the optimisation algorithm, implementing IEvolutionEngine<T>
	// It is recommended that you declare variables typed as interfaces, not concrete classes, 
	// whenever possible...

	// wila uses templates to make the most of C++ strengths in compile time checking
	// the following two lines are just convenience aliasing of variable types.
	using HC = HyperCube < double >;
	using SceType = ShuffledComplexEvolution<HC>;


	// Let's define a simple objective , the traditionall "L2" geometric distance from a point
	HC goal;
	goal.Define("a", 1, 2, 1.234);
	goal.Define("b", 3, 4, 3.1415);
	goal.Define("c", -10, 10, -4);
	goal.Define("d", -100, 123, -3.1415);
	TopologicalDistance<HC> measure(goal);

	// The optimizer we will use requires a population; let's define one. 
	// By default this will seed using a uniform random distribution on the feasible parameter space.
	CandidateFactorySeed<HC> seeding(0, goal);



	SceParameters sceParams = CreateSceParamsForProblemOfDimension(goal.Dimensions(), 20);
	sceParams.P = 5;
	sceParams.Pmin = 3;


	// A termination criteria for the optimization. 
	// for the sake of the example:
	MaxNumberSceShuffles<HC, SceType> c(10);
	auto terminationCondition = /*typename*/ SceType::TerminationCondition(c);
	
	
	ShuffledComplexEvolution<HC> opt(measure, seeding, terminationCondition, sceParams);

	// wila is designed to make use of multiple core machines if available and desired. 
	size_t nCores = opt.GetMaxHardwareConcurrency();
	size_t nThreads = std::min((size_t)3, nCores - 1);
	opt.SetMaxDegreeOfParallelism(nThreads);


	IObjectiveEvaluator<HC>& evaluator = measure;
	IEvolutionEngine<HC>& finder = opt;


	IOptimizationResults<HC> results = finder.Evolve();

	results.PrintTo(std::cout);

	return 0;
}


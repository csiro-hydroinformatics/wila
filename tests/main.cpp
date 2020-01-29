#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#include <set>
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include <iterator>
#include "catch/catch.hpp"


//#define WILA_USE_VLD
#ifdef WILA_USE_VLD
#include<vld.h>
#endif

//#define LOG_VALUE

#include "common.h"

#include <wila/interop_c_cpp.hpp>
#include "wila/urs.hpp"

#define REQUIRE_WITHIN_ABSOLUTE_TOLERANCE( expected, actual, delta) REQUIRE( (abs(expected - actual) < delta) )

void checkLOneTolerance(Hc point, const std::map<string, double>& expected, double delta)
{
	for (const std::pair<string, double>& p : expected)
	{
		double actual = point.GetValue(p.first);
		double expected = p.second;
		REQUIRE_WITHIN_ABSOLUTE_TOLERANCE(expected, actual, delta);
	}
}

using namespace std;
using namespace mhcpp;
using namespace mhcpp::random;
using namespace mhcpp::optimization;
using namespace mhcpp::utils;

void setDefaults() {
	mhcpp::threading::ThreadingOptions<>::DefaultMaxDegreeOfParallelism = 4;
}

TEST_CASE("Number formating is scientific by default", "[utils]") {
	REQUIRE(ToString(1.23456789) == "1.234568e+00");
	REQUIRE(ToString(1.23456711) == "1.234567e+00");
	REQUIRE(ToString(1.23) == "1.230000e+00");
	REQUIRE(ToString(12.3456789) == "1.234568e+01");
	REQUIRE(ToString(1.23456789e33) == "1.234568e+33");
	REQUIRE(ToString(1.23456789e-33) == "1.234568e-33");
	REQUIRE(ToString(-1.23456789e-33) == "-1.234568e-33");

	REQUIRE(ToString(123456789) == "123456789");
	REQUIRE(ToString(123456789123456789ul) == "123456789123456789");
}

TEST_CASE("Basic hypercubes", "[sysconfig]") {

	GIVEN("A 2 dimensional hypercube")
	{
		Hc hc = createTestHc(1.5, 3.3);
		vector<string> keys;
		keys.push_back("a");
		keys.push_back("b");

		WHEN("Read values after initial definition") {
			REQUIRE(hc.Dimensions() == 2);
			REQUIRE(mhcpp::utils::AreSetEqual(keys, hc.GetVariableNames()));
			REQUIRE(hc.GetValue("a") == 1.5);
			REQUIRE(hc.GetMinValue("a") == 1);
			REQUIRE(hc.GetMaxValue("a") == 2);
			REQUIRE(hc.GetValue("b") == 3.3);
			REQUIRE(hc.GetMinValue("b") == 3);
			REQUIRE(hc.GetMaxValue("b") == 4);
		}

		WHEN("Change definition values after initial definition") {

			hc.SetValue("a", 1.1);
			hc.SetMinValue("a", 0.5);
			hc.SetMaxValue("a", 1.99);
			hc.SetValue("b", 5);
			hc.SetMinValue("b", 0);
			hc.SetMaxValue("b", 10);

			REQUIRE(hc.GetValue("a") == 1.1);
			REQUIRE(hc.GetValue("b") == 5);
			REQUIRE(hc.GetMinValue("a") == 0.5);
			REQUIRE(hc.GetMaxValue("a") == 1.99);
			REQUIRE(hc.GetMinValue("b") == 0);
			REQUIRE(hc.GetMaxValue("b") == 10);
		}
	}
}

TEST_CASE("Calculation of the spread of a population of hypercubes", "[sysconfig]")
{
	std::vector<Hc> points;

	double aMin = 1.1;
	double bMin = 3.3;
	double aMax = 2.2;
	double bMax = 4.5;
	std::function<Hc(double, double)> f = [&](double a, double b)
	{
		return createTestHc(a, b, aMin, bMin, aMax, bMax);
	};
	points.push_back(f(1.2, 2.54));
	points.push_back(f(1.27, 2.56));
	points.push_back(f(1.32, 2.87));

	auto sdevs = mhcpp::GetStdDev(points);
	//> sd(c(1.2, 1.27, 1.32))
	//[1] 0.06027714
	//> sd(c(2.54, 2.56, 2.87))
	//[1] 0.1850225
	REQUIRE_WITHIN_ABSOLUTE_TOLERANCE(0.06027714, sdevs["a"], 1e-7);
	REQUIRE_WITHIN_ABSOLUTE_TOLERANCE(0.1850225, sdevs["b"], 1e-7);

	auto rsdevs = mhcpp::GetRelativeSdev(points);
	//> sd(c(1.2, 1.27, 1.32)) / (2.2 - 1.1)
	//[1] 0.0547974
	REQUIRE_WITHIN_ABSOLUTE_TOLERANCE(0.0547974, rsdevs["a"], 1e-7);
	//> sd(c(2.54, 2.56, 2.87)) / (4.5 - 3.3)
	//[1] 0.1541854
	REQUIRE_WITHIN_ABSOLUTE_TOLERANCE(0.1541854, rsdevs["b"], 1e-7);


	std::map<string, double> dict;
	dict["a"] = 1.0;
	dict["b"] = std::numeric_limits<double>::quiet_NaN();
	dict["bb"] = -std::numeric_limits<double>::quiet_NaN();
	dict["c"] = -std::numeric_limits<double>::infinity();
	dict["d"] = std::numeric_limits<double>::infinity();
	dict["e"] = 11.0;
	dict["f"] = 123.0;
	rsdevs = mhcpp::utils::RemoveNotFinite(dict);

	auto keys = mhcpp::utils::GetKeys(rsdevs);

	REQUIRE( requireEqual(keys, std::vector<string>({ "a","e", "f" })) );

}

TEST_CASE("Calculation of a centroid", "[sysconfig]") 
{
	GIVEN("A population of hypercubes")
	{
		std::vector<Hc> points;
		points.push_back(createTestHc(1.1, 2.1));
		points.push_back(createTestHc(1.2, 2.2));
		auto c = Hc::GetCentroid(points);
		WHEN("Population of two points") {
			THEN("Expected barycentre") {
				REQUIRE(assertHyperCube(c, 1.15, 2.15));
			}
		}
		points.push_back(createTestHc(1.9, 2.6));
		c = Hc::GetCentroid(points);
		WHEN("Population of three points") {
			THEN("Expected barycentre") {
				REQUIRE(assertHyperCube(c, 1.4, 2.3));
			}
		}
	}
}

TEST_CASE("Basic objective evaluator", "[objectives]") {

	GIVEN("Single-objective calculator, L2 distance")
	{
		Hc hc = createTestHc(1.5, 3.3);
		Hc goal = createTestHc(1, 3);

		TopologicalDistance<HyperCube < double > > evaluator(goal);

		WHEN("Evaluating distance") {
			IObjectiveScores<Hc> scores = evaluator.EvaluateScore(hc);
			THEN("Gets one objective value with expected value"){
				REQUIRE(scores.ObjectiveCount() == 1);
				REQUIRE(scores.ObjectiveName(0) != "");
				REQUIRE(scores.Value(0) == std::sqrt(0.25 + 0.09));
			}
		}
	}
}

TEST_CASE("Sort objective scores", "[objectives]") {

	GIVEN("Single-objective calculator, L2 (Euclidian) distance")
	{
		Hc goal = createTestHc(0,0);

		//IObjectiveEvaluator<HyperCube < double > >* evaluator = new TopologicalDistance<HyperCube < double > >(goal);
		TopologicalDistance<HyperCube < double > > evaluator(goal);

		vector<Hc> points({
			createTestHc(1, 2),
			createTestHc(10, 10),
			createTestHc(3, 3),
			createTestHc(1,1)
		});
		string scoreName = "L2 distance";
		vector<IObjectiveScores<HyperCube < double > >> scores = IObjectiveEvaluator<HyperCube < double > >::EvaluateScores(evaluator, points);
		IObjectiveScores<HyperCube < double > >::Sort(scores, scoreName);

		//int n = scores.size();
		//for (size_t i = 0; i < n; ++i)
		//	std::cout << i << ": " << scores[i].ToString() << std::endl;

		REQUIRE(scores.size() == 4);
		REQUIRE(scores[0].Value(0) == std::sqrt(1 + 1));
		REQUIRE(scores[1].Value(0) == std::sqrt(1 + 4));
		REQUIRE(scores[2].Value(0) == std::sqrt(9 + 9));
		REQUIRE(scores[3].Value(0) == std::sqrt(100 + 100));
	}
}

TEST_CASE("RNG basics", "[rng]") {
	// These test cases were largely written for chasing up the bug
	// https://jira.csiro.au/browse/WIRADA-341


	// https://jira.csiro.au/browse/WIRADA-392
	// Windows output, run with #define LOG_VALUE   as of 2016-08-09
//first number sampled : 3499211612
//firstInt : 581869302
//secondInt : 2412496532
//x1[0] : 2991312382
//x2[0] : 1068398491
//nextSeed = factory.PeekNext() : 2991312382
//nextInitSeed = factory.PeekNext() : 3062119789

	GIVEN("An mhcpp::random::default_wila_random_engine")
	{
		unsigned int seed = 123;
		mhcpp::random::default_wila_random_engine stdRng;
		WHEN("The first number is sampled from it") {
			THEN("It is not the seed that was used")
			{
				unsigned int myValue = LogVarValue<unsigned int>(stdRng(), "first number sampled");
				REQUIRE(myValue == 3499211612);
				REQUIRE_FALSE(myValue == seed);
			}
		}
		WHEN("another engine is created, seeded with the sampled value out of the first") {
			THEN("The first value sampled from each of the two engine differs")
			{
				unsigned int seedNewRng = stdRng();				
				mhcpp::random::default_wila_random_engine stdRng2(seedNewRng);
				unsigned int firstInt =  LogVarValue<unsigned int>(stdRng() , "firstInt");
				unsigned int secondInt = LogVarValue<unsigned int>(stdRng2(), "secondInt");
				REQUIRE_FALSE(firstInt == seedNewRng);
				REQUIRE_FALSE(firstInt == secondInt);
				REQUIRE(firstInt == 581869302);
				REQUIRE(secondInt == 2412496532);
			}
		}
	}
	GIVEN("A random integer generator IRandomNumberGeneratorFactory<>")
	{
		IRandomNumberGeneratorFactory<> factory(123);
		WHEN("Another is created with a different seed") {
			IRandomNumberGeneratorFactory<> f2(456);
			THEN("The sequence of random numbers is different")
			{
				REQUIRE_FALSE(factory.Equals(f2));
				auto x1 = factory.Next(10);
				LogVarValue<unsigned int>(x1[0], "x1[0]");
				auto x2 = f2.Next(10);
				LogVarValue<unsigned int>(x2[0], "x2[0]");
				REQUIRE(x1[0] != x1[1]);
				for (size_t i = 0; i < 10; i++)
					REQUIRE(x1[i] != x2[i]);
				REQUIRE(x1[0] == 2991312382);
				REQUIRE(x2[0] == 1068398491);
			}
		}
		WHEN("A copy is created by assignment") {
			IRandomNumberGeneratorFactory<> f2 = factory;
			THEN("The sequence of random numbers is the same")
			{
				REQUIRE(factory.Equals(f2));
				auto x1 = factory.Next(10);
				auto x2 = f2.Next(10);
				REQUIRE(x1[0] != x1[1]);
				for (size_t i = 0; i < 10; i++)
					REQUIRE(x1[i] == x2[i]);
			}
		}
		WHEN("A new factory is generated from the first one via CreateNew") {
			IRandomNumberGeneratorFactory<> f2;
			THEN("The initial factory changes seed after CreateNew, if we peek")
			{
				auto nextSeed = LogVarValue<unsigned int>(factory.PeekNext(), "nextSeed=factory.PeekNext()");
				REQUIRE(nextSeed == factory.PeekNext());
				f2 = factory.CreateNew();
				auto nextInitSeed = LogVarValue<unsigned int>(factory.PeekNext(), "nextInitSeed=factory.PeekNext()");
				REQUIRE_FALSE(nextSeed == nextInitSeed);
				REQUIRE(nextSeed == 2991312382);
				REQUIRE(nextInitSeed == 3062119789);
				AND_THEN("The first seed produced by each of these rng factories differs")
				{
					REQUIRE(factory.PeekNext() != f2.PeekNext());
					REQUIRE(factory.Next() != f2.Next());
					REQUIRE_FALSE(factory.Equals(f2));
					auto x1 = factory.Next(10);
					auto x2 = f2.Next(10);
					REQUIRE(x1[0] != x1[1]);
					for (size_t i = 0; i < 10; i++)
						REQUIRE_FALSE(x1[i] == x2[i]);
				}
			}
		}
	}
}

TEST_CASE("RNG memory management", "[memory]") {
	auto count = [&]() {return TestRandomEngine::NumInstances(); };
	int initialCount = count();
	using rngtype = IRandomNumberGeneratorFactory<TestRandomEngine>;
	rngtype* factory = new rngtype(123);
	REQUIRE((initialCount + 1) == count());
	rngtype* f2 = new rngtype(456);
	REQUIRE((initialCount + 2) == count());
	rngtype f3 = factory->CreateNew();
	REQUIRE((initialCount + 3) == count());
	
	mhcpp::random::uniform_real_distribution_double dist(0, 1);
	auto rng = f3.CreateVariateGenerator<mhcpp::random::uniform_real_distribution_double>(dist, 333);
	REQUIRE((initialCount + 4) == count());
	delete factory;
	REQUIRE((initialCount + 3) == count());
	delete f2;
	REQUIRE((initialCount + 2) == count());
}

TEST_CASE("RNG factories and resulting RNG engines", "[rng]") {
	GIVEN("RNG factories spawned hierarchically")
	{
		int seed = 456;
		int numFactories = 2;
		int numEngines = 3;
		int numDraw = 10;

		vector<vector<IRandomNumberGeneratorFactory<>>> rngs(numFactories);
		vector<vector<mhcpp::random::default_wila_random_engine>> rngEngines(numFactories);
			
		IRandomNumberGeneratorFactory<> rngf(seed);
		for (int i = 0; i < numFactories; i++)
		{
			auto rngf_2 = rngf.CreateNew();
			REQUIRE_FALSE(rngf_2.Equals(rngf));
			for (int j = 0; j < numEngines; j++)
			{
				rngs[i].push_back(rngf_2);
				auto e = rngf_2.CreateNewStd();
				rngEngines[i].push_back(e);
			}
		}

		for (int i = 0; i < numFactories; i++)
			for (int j = 0; j < numEngines; j++)
				for (int i_1 = 0; i_1 < numFactories; i_1++)
					for (int j_1 = 0; j_1 < numEngines; j_1++)
					{
						if ((i != i_1) && (j != j_1))
						{
							string indicesMsg = 
								string("[") +
								boost::lexical_cast<string>(i) +
								string("][") +
								boost::lexical_cast<string>(j) +
								string("]") +
								string(" and [") +
								boost::lexical_cast<string>(i_1) +
								string("][") +
								boost::lexical_cast<string>(j_1) +
								string("]");
							if (rngs[i][j].Equals(rngs[i_1][j_1]))
							{
								string msg =
									string("RNG factories are equal for ") + indicesMsg;
								FAIL(msg);
							}
							if (rngEngines[i][j] == rngEngines[i_1][j_1])
							{
								string msg =
									string("RNG engines are equal for ") + indicesMsg;
								FAIL(msg);
							}
						}
					}
	}
}

TEST_CASE("RNG sequences", "[rng]") {
	GIVEN("A random integer generator IRandomNumberGeneratorFactory<>")
	{
		auto v1 = CreateRandValues();
		THEN("The sequence of random numbers is different")
		{
			bool b = v1->IsPairwiseDistinct();
			string msg = "";
			if (!b)
				WARN(v1->FailedPairwiseDistinct());
			REQUIRE(b);
		}

		auto v2 = CreateRandValues();
		WHEN("A different sequence is create with the same default seed and all inputs equal") 
		{
			THEN("these two distinct randomized sequences are identical in their content")
			{
				REQUIRE(v1->Equals(*v2));
			}
		}
		delete v2;
		auto v3 = CreateRandValues(3);
		WHEN("A sequence is create with a different seed and all other inputs equal")
		{
			THEN("these two distinct randomized sequences differ in their content")
			{
				REQUIRE_FALSE(v1->Equals(*v3));
			}
		}
		delete v1;
		auto v4 = CreateRandValues(3);
		WHEN("Two sequences are created but with the same but non default seed")
		{
			THEN("these two distinct randomized sequences are identical in their content")
			{
				REQUIRE(v3->Equals(*v4));
			}
		}
		delete v3;
		delete v4;
	}
}

TEST_CASE("trapezoidal, discrete RNG to sample from a population of points, as used to create the sub-complexes", "[rng]")
{
	mhcpp::random::default_wila_random_engine generator(234);
	const int ncandidates = 10;
	RngInt<> rng = CreateTrapezoidalRng(ncandidates, generator);

	const int nrolls = 100000; // number of experiments
	auto p = SampleFrom(rng, nrolls);
	//std::cout << "a discrete_distribution:" << std::endl;

	//PrintHistogram(p, std::cout);
	//PrintValues(p, std::cout, true);
	//PrintValues(p, std::cout, false);

	int n = ncandidates;
	std::vector<double> expectedProportions(n);
	double total = n*(n + 1) / 2; // sum ints from 1 to n;
	for (size_t i = 1; i <= n; i++)
		expectedProportions[i - 1] = (n + 1 - i) / total;

	// PrintVec(RelativeDiff(expectedProportions, Normalize(p)), std::cout);
	// Build on windows with MSVCC, I get a max relative deviation of <3%, and ~4.2% on linux/gcc5
	REQUIRE(requireRelativeEqual(expectedProportions, Normalize(p), 4.5 * 1e-2));

}

TEST_CASE("URS RNG basics", "[rng]") {
	Hc hc;
	hc.Define("a", 1, 2, 1.5);
	hc.Define("b", 3, 4, 3.3);
	GIVEN("An uniform random sampler")
	{
		auto rng = UniformRandomSamplingFactory<Hc>(IRandomNumberGeneratorFactory<>(), hc);

		WHEN("Creating a random point with default template") {
			Hc p = rng.CreateRandomCandidate();
			THEN("Feasible parameter space is the same as the template, but values are different from template")
			{
				REQUIRE(p.GetMinValue("a") == 1.0);
				REQUIRE(p.GetMaxValue("a") == 2.0);
				REQUIRE(p.GetMinValue("b") == 3.0);
				REQUIRE(p.GetMaxValue("b") == 4.0);

				REQUIRE(assertValuesNotEqual(p, hc));

				REQUIRE(p.GetValue("a") >= 1.0);
				REQUIRE(p.GetValue("b") >= 3.0);

				REQUIRE(p.GetValue("a") <= 2.0);
				REQUIRE(p.GetValue("b") <= 4.0);
			}
		}

		WHEN("A copy of the URS factory is made, and samples taken from each") {
			auto rngCopy = rng;
			THEN("The two URS objects are equal") {
				REQUIRE(rng.Equals(rngCopy));
			}
			auto p1 = rng.CreateRandomCandidate();
			auto p2 = rngCopy.CreateRandomCandidate();
			AND_THEN("The sampled values are the same as returned by the original one") {
				REQUIRE(rng.Equals(rngCopy));
				REQUIRE(assertEqual(p1, p2));
				REQUIRE(assertValuesNotEqual(p2, hc));
			}
			AND_WHEN("A new factory is created from these existing ones") {
				auto urs1 = rng.CreateNew();
				auto urs2 = rngCopy.CreateNew();
				auto p1_2 = urs1->CreateRandomCandidate();
				auto p2_2 = urs2->CreateRandomCandidate();
				THEN("The sampled values from these two new factories are identical") {
					REQUIRE(assertEqual(p1_2, p2_2));
				}
				AND_THEN("But differ from ones from the original factories") {
					REQUIRE(assertValuesNotEqual(p1_2, p1));
					REQUIRE(assertValuesNotEqual(p2_2, p2));
				}
				delete urs1;
				delete urs2;
			}
		}
	}
}

TEST_CASE("Complex for SCE, single objective", "[optimizer]") {

	setDefaults();
	using T = HyperCube < double > ;
	int hcStart = T::NumInstances();
	int m = 20;
	int q = 10, alpha = 2, beta = 3;
	IRandomNumberGeneratorFactory<> rng(2);
	auto unif = createTestUnifrand<T>(421);

	std::vector < IObjectiveScores<T> > scores = createTestScores<T>(m, 123);
	REQUIRE((hcStart + m + 1) == T::NumInstances());
	IFitnessAssignment<double, T> fitnessAssignment;
	//IHyperCubeOperations* hyperCubeOperations
	//ILoggerMh* logger = nullptr;
	T goal = createTestHc(1.5, 3.4);
	TopologicalDistance<T> evaluator(goal);

	WHEN("looking for the point with the worst fitness value")
	{
		std::vector<FitnessAssignedScores<double, T>> fvec;
		std::vector<FitnessAssignedScores<double, T>> subpopulation;
		for (size_t i = 0; i < scores.size(); i++)
		{
			fvec.push_back(FitnessAssignedScores<double, T>(scores[i], i));
		}
		auto worst = SubComplex<T>::FindWorstPoint(fvec, subpopulation);
		THEN("The score of the worst point found is indeed the worst fitness we assigned in the population")
		{
			REQUIRE(worst.FitnessValue() == (double)(scores.size() - 1));
			REQUIRE(subpopulation.size() == (scores.size() - 1));
			for (size_t i = 0; i < subpopulation.size(); i++)
				REQUIRE(subpopulation[i].FitnessValue() != worst.FitnessValue());
		}
	}
	WHEN("Building and running a subcomplex") {
		auto discreteRng = CreateTrapezoidalRng(scores.size(), rng.CreateNewStd());
		int beforeSubCplx = T::NumInstances();
		SubComplex<T>* scplx = new SubComplex<T>
			(scores, &evaluator, q, alpha, rng, &unif, discreteRng,
			fitnessAssignment);
		THEN("The subcomplex evolution completes without exception")
		{
			// Calling Evolve works without exceptions.
			int beforeEvolve = T::NumInstances();
			REQUIRE(beforeEvolve == (beforeSubCplx + scores.size()));
			scplx->Evolve();
			int afterEvolve = T::NumInstances();
			REQUIRE(beforeEvolve == afterEvolve);
			AND_THEN("Retrieving the final population does return a vector of expected size = m")
			{
				std::vector<IObjectiveScores<T>> finalPop;
				REQUIRE_NOTHROW(finalPop = scplx->WholePopulation());
				REQUIRE(finalPop.size() == m);
				// TODO: further tests.
			}
		}
		delete scplx;
		REQUIRE(beforeSubCplx == T::NumInstances());
	}
	WHEN("Building and running a complex")
	{
		auto unif = createTestUnifrand<T>(421);
		ShuffledComplexEvolution<Hc>::TerminationCondition terminationCondition;
		terminationCondition.RequireEngine = false;

		//Complex<T> cplx_noargs;

		Complex<T>* cplx = new Complex<T>
			(scores, &evaluator, false, rng, &unif, false,
			fitnessAssignment, terminationCondition, true, nullptr, std::map<string, string>(), q, alpha, beta);
		THEN("The complex evolution completes without exception")
		{
			REQUIRE_NOTHROW(cplx->Evolve());
			// TODO: further tests.
			delete cplx;
		}
	}
	scores.clear();
	// at this point only remaining parameterizers should be goal, unif and evaluator
	REQUIRE((hcStart + 3) >= T::NumInstances());

}

TEST_CASE("SCE basic port", "[optimizer]") {
	setDefaults();
	GIVEN("A 2D Hypercube")
	{
		auto terminationCondition = CreateCounterTermination(100);
		Hc goal;
		ShuffledComplexEvolution<Hc> opt = CreateSceQuadraticGoal(goal, terminationCondition);
		WHEN("optimizing single-thread") {
			opt.UseMultiThreading(false);
			auto results = opt.Evolve();
			REQUIRE(results.size() > 0);
			//results.PrintTo(std::cout);
			auto first = results[0];
			REQUIRE(first.ObjectiveCount() == 1);
		}
		WHEN("optimizing multi-threaded") {
			opt.UseMultiThreading(true);
			auto results = opt.Evolve();
			REQUIRE(results.size() > 0);
			//results.PrintTo(std::cout);
			auto first = results[0];
			REQUIRE(first.ObjectiveCount() == 1);
		}
	}
}

TEST_CASE("exception handling is thread-safe", "[optimizer]") {
	setDefaults();
	GIVEN("An objective calculation that triggers an std::exception")
	{
		auto terminationCondition = CreateCounterTermination(100);
		Hc goal;
		ShuffledComplexEvolution<Hc> opt = CreateSceQuadraticGoalThrowsException(goal, terminationCondition);
		WHEN("optimizing multi-threaded") {
			opt.UseMultiThreading(true);
			string msg;
			try {
				auto results = opt.Evolve();
			}
			catch (std::exception& e)
			{
				msg = e.what();
			}
		}
	}
}

TEST_CASE("Termination conditions", "[optimizer]") {
	setDefaults();

	// TODO (not thread safe now)
	//GIVEN("A marginal improvement termination condition")
	GIVEN("A runtime length termination condition")
	{
		auto terminationCondition = CreateWallClockTermination<>(5.0);
		Hc goal;
		ShuffledComplexEvolution<Hc> opt = CreateSceQuadraticGoal(goal, terminationCondition);

		WHEN("") {
			auto results = opt.Evolve();
			REQUIRE(results.size() > 0);
			//results.PrintTo(std::cout);
			auto first = results[0];
			REQUIRE(first.ObjectiveCount() == 1);
		}
	}

	GIVEN("A max nb iteration termination condition")
	{
		auto terminationCondition = CreateMaxIterationTermination(1000);
		Hc goal;

		ShuffledComplexEvolution<Hc> opt = CreateSceQuadraticGoal(goal, terminationCondition);
		WHEN("Complexes are allowed to check on termination criterion, and the criterion is thread safe") {
			REQUIRE(terminationCondition.IsThreadSafe());
			opt.AllowComplexPrematureTermination(true);
			auto results = opt.Evolve();
			REQUIRE(results.size() > 0);
			THEN("The number of evaluations is very close to the maximum nb of iterations") {
				REQUIRE(opt.EvaluationCount() > 1000);
				REQUIRE(opt.EvaluationCount() < 1020);
			}
		}
		WHEN("Complexes are not allowed to check on termination criterion, and the criterion is thread safe") {
			REQUIRE(terminationCondition.IsThreadSafe());
			opt.AllowComplexPrematureTermination(false);
			auto results = opt.Evolve();
			REQUIRE(results.size() > 0);
			THEN("The number of evaluations is more than the maximum nb of iterations, but not too much") {
				REQUIRE(opt.EvaluationCount() > 1001);
				REQUIRE(opt.EvaluationCount() < 1100);
			}
		}
	}

	GIVEN("A maximum parameter population std dev as a termination condition")
	{
		auto terminationCondition = CreateStdDevTermination(1.0/100);
		Hc goal;

		ShuffledComplexEvolution<Hc> opt = CreateSceQuadraticGoal(goal, terminationCondition);
		WHEN("Optimizing") {
			opt.AllowComplexPrematureTermination(false);
			auto results = opt.Evolve();
			REQUIRE(results.size() > 2);
			THEN("The number of evaluations is more than TODO, but not too much") {
				REQUIRE(opt.EvaluationCount() > 33); // initial pop size.
				//REQUIRE(opt.EvaluationCount() < 1100);
			}
		}
	}
}

TEST_CASE("Memory management", "[memory]") {
	setDefaults();
	auto terminationCondition = CreateWallClockTermination<>(5.0);
	using T = HyperCube < double > ;
	T goal;
	int hcBeforeCreation = T::NumInstances();
	ShuffledComplexEvolution<T>* opt = CreateSceQuadraticGoalPtr(goal, terminationCondition);
	int hcAfterCreation = T::NumInstances();
	//REQUIRE(hcStart == hcEnd);
	int hcBeforeEvolve = T::NumInstances();
	opt->SetLogger();
	IOptimizationResults<T> results = opt->Evolve();
	int hcAfterEvolve = T::NumInstances();
	results = IOptimizationResults<T>();
	int hcAfterDeleteResults = T::NumInstances();
	delete opt;
	int hcAfterDeleteOpt = T::NumInstances();
	REQUIRE(hcAfterDeleteOpt == hcBeforeCreation);
}

TEST_CASE("Multiple evolve runs in sequence", "[optimizer]") {
	// This test checks that an segfault issue reported running on Linux is not present.

	setDefaults();
	auto terminationCondition = CreateWallClockTermination<>(14.4);
	using T = HyperCube < double >;
	T goal;

	int hcInitialCount = T::NumInstances();
	int iterNum;
	for (int i = 0; i <= 7; i++) {
		iterNum=i+1;
		REQUIRE(iterNum == iterNum); // hack to check via cmd line in verbose mode.
		int hcNewLoopCount = T::NumInstances();
		REQUIRE(hcInitialCount == hcNewLoopCount);
		ShuffledComplexEvolution<Hc> opt = CreateSceQuadraticGoal(goal, terminationCondition);
		size_t nCores = opt.GetMaxHardwareConcurrency();
		size_t nThreads = std::min((size_t)3, nCores - 1);
		opt.SetMaxDegreeOfParallelism(nThreads);
		auto results = opt.Evolve();
		REQUIRE(results.size() > 0);
		auto first = results[0];
		REQUIRE(first.ObjectiveCount() == 1);
	}
}

std::map<string, double> mkPt (double a, double b) {
	return std::map<string, double>({ { "a", a },{ "b", b } });
};

TEST_CASE("SCE behavior is identical across OSes", "[rng]") {
	setDefaults();
	GIVEN("A 2D Hypercube")
	{
		auto terminationCondition = CreateMaxNumShuffle(1);
		Hc goal;
		ShuffledComplexEvolution<Hc> opt = CreateSceQuadraticGoal(goal, terminationCondition);
		double delta = 1e-5;

		WHEN("Initializing the population of points") {
			auto points = opt.CreatePopulation();
			// mhcpp::utils::PrintTo<Hc>(scores, std::cout);
			// as of 2016-08-09 on Windows:
			// > testwila.cmd "SCE behavior is identical across OSes"
			// 0: a:1.947638, b: 3.883420,
			// 10: a: 1.448429, b: 3.749307,
			// 54: a: 1.895343, b: 3.158786,
			checkLOneTolerance(points[0], mkPt(1.947638, 3.883420), delta);
			checkLOneTolerance(points[10], mkPt(1.448429, 3.749307), delta);
			checkLOneTolerance(points[54], mkPt(1.895343, 3.158786), delta);
		}
		WHEN("optimizing single-thread") {
			opt.UseMultiThreading(false);
			opt.SetLogger();
			auto results = opt.Evolve();
			//results.PrintTo(std::cout);
			// as of 2016-08-09 on Windows:
			//0: L2 distance:0.098551, a:1.074382, b:3.064649,
			//10: L2 distance:1.207733, a:1.968343, b:3.721755,
			//54: L2 distance:1.311178, a:1.874004, b:3.977397,
			std::vector<Hc> points = IObjectiveScores<Hc>::GetSystemConfigurations(results);
#ifdef LOG_VALUE
			mhcpp::utils::PrintTo< Hc >(points, std::cout);
			mhcpp::logging::ILoggerMh< Hc >* logger = opt.GetLogger();
			mhcpp::utils::PrintLoggerTo< Hc >(logger, std::cout);
#endif
			checkLOneTolerance(points[0], mkPt(1.074382, 3.064649), delta);
			checkLOneTolerance(points[10], mkPt(1.968343, 3.721755), delta);
			checkLOneTolerance(points[54], mkPt(1.874004, 3.977397), delta);
		}
		WHEN("optimizing multi-threaded") {
			opt.UseMultiThreading(true);
			opt.SetLogger();
			auto results = opt.Evolve();
			auto points = IObjectiveScores<Hc>::GetSystemConfigurations(results);
#ifdef LOG_VALUE
			mhcpp::utils::PrintTo< Hc >(points, std::cout);
			mhcpp::logging::ILoggerMh< Hc >* logger = opt.GetLogger();
			mhcpp::utils::PrintLoggerTo< Hc >(logger, std::cout);
#endif
			checkLOneTolerance(points[0], mkPt(1.074382, 3.064649), delta);
			checkLOneTolerance(points[10], mkPt(1.968343, 3.721755), delta);
			checkLOneTolerance(points[54], mkPt(1.874004, 3.977397), delta);
		}
	}
}

TEST_CASE("Distribution sampling is deterministic and identical across platforms", "[rng]") {
	mhcpp::random::uniform_real_distribution_double dist(0, 1);
	mhcpp::random::default_wila_random_engine eng(342);
	auto e1 = LogVarValue<unsigned int>(eng(), "e1=eng()");
	auto v = LogVarValue<double>(dist(eng), "v=dist(eng)");
//e1=eng(): 279320595
//v=dist(eng): 0.985974
	REQUIRE(e1 == 279320595);
	REQUIRE_WITHIN_ABSOLUTE_TOLERANCE(.985974, v, 1e-5);
}

TEST_CASE("Candidate factory seed deterministic across platforms", "[rng]") {
	GIVEN("a IRandomNumberGeneratorFactory")
	{
//> testwila.cmd "Candidate factory seed deterministic across platforms"
//peek1=f.PeekNext(): 2357136044
//peek2=f.PeekNext(): 2357136044
//peek3=f.PeekNext(): 2546248239
//e1=engine(): 2714253906
//e2=engine(): 1225636010
//e3=engineSeeded(): 3008354540
//e4=engineSeeded(): 440739714
//seedVg=f.PeekNext(): 3071714933
//v=vg(): 0.759871
		IRandomNumberGeneratorFactory<> f(0);
		auto peek1 = LogVarValue<unsigned int>(f.PeekNext(), "peek1=f.PeekNext()");
		auto peek2 = LogVarValue<unsigned int>(f.PeekNext(), "peek2=f.PeekNext()");
		IRandomNumberGeneratorFactory<> f2 = f.CreateNew();
		auto peek3 = LogVarValue<unsigned int>(f.PeekNext(), "peek3=f.PeekNext()");
		auto engine = f.CreateNewEngine();
		auto e1 = LogVarValue<unsigned int>(engine(), "e1=engine()");
		auto e2 = LogVarValue<unsigned int>(engine(), "e2=engine()");
		auto engineSeeded = IRandomNumberGeneratorFactory<>::CreateNewEngine(666);
		auto e3 = LogVarValue<unsigned int>(engineSeeded(), "e3=engineSeeded()");
		auto e4 = LogVarValue<unsigned int>(engineSeeded(), "e4=engineSeeded()");
		mhcpp::random::uniform_real_distribution_double dist(0, 1);
		auto seedVg = LogVarValue<unsigned int>(f.PeekNext(), "seedVg=f.PeekNext()");
		auto vg = f.CreateVariateGenerator<>(dist);
		auto v = LogVarValue<double>(vg(), "v=vg()");
		REQUIRE(peek1 == 2357136044);
		REQUIRE(peek1 == peek2);
		REQUIRE(peek3 == 2546248239);
		REQUIRE(e1 == 2714253906);
		REQUIRE(e2 == 1225636010);
		REQUIRE(e3 == 3008354540);
		REQUIRE(e4 == 440739714);
		REQUIRE(seedVg == 3071714933);		
		REQUIRE_WITHIN_ABSOLUTE_TOLERANCE(0.759871, v, 1e-5);
	}
	GIVEN("A candidate factory seed")
	{
		Hc hc = CreateTestHc();
		auto c = CreateCandidateFactorySeed(0, hc);
		double delta = 1e-5;

		WHEN("Initializing the population of points") {
			auto cf = c.Create();
			auto points = cf->CreateRandomCandidates(3);
			//mhcpp::utils::PrintTo<Hc>(points, std::cout);
			delete cf;
			// as of 2016-08-09 on Windows:
//> testwila.cmd "Candidate factory seed deterministic across platforms"
//0: a:1.947638, b:3.883420,
//1: a:1.061278, b:3.673165,
//2: a:1.173725, b:3.378660,			
			checkLOneTolerance(points[0], mkPt(1.947638, 3.883420), delta);
			checkLOneTolerance(points[1], mkPt(1.061278, 3.673165), delta);
			checkLOneTolerance(points[2], mkPt(1.173725, 3.378660), delta);
		}
	}
}

TEST_CASE("Test templated generic utilities for maps and vectors", "[utils]")
{
	std::map < string, string > dict = {
		{string("a"), string("A")},
		{string("b"), string("B")},
		{string("c"), string("C")}
	};
	vector<string> sub = { "b" };
	auto res = mhcpp::utils::Subset(sub, dict);
	REQUIRE(res.size() == 1);
	REQUIRE(mhcpp::utils::HasKey(res, "b"));
}

TEST_CASE("Optimizer log data interop creation/disposal")
{

	// A minimal unit test after refactoring log data creation/disposal out of RPP and swift.
	SimpleLogger<Hc> logger;

	Hc goal = createTestHc(1, 3);
	TopologicalDistance<Hc> evaluator(goal);

	IObjectiveScores<Hc> scores = evaluator.EvaluateScore(createTestHc(1.5, 3.3));
	logger.Write(scores, std::map<string, string>({ { "key1", "value1" } }));

	scores = evaluator.EvaluateScore(createTestHc(2.5, 3.3));
	logger.Write(scores, std::map<string, string>({ { "key1", "value2" } }));
	scores = evaluator.EvaluateScore(createTestHc(2.5, 3.4));
	logger.Write(scores, std::map<string, string>({ { "key2", "value1" } }));

	logger.Write("Some message", std::map<string, string>({ { "key3", "attribute of the message" } }));

	auto log = mhcpp::interop::get_optimizer_log_data<Hc>(logger);

	REQUIRE(log->LogLength == 3);
	REQUIRE(log->StringDataCount == 4);
	REQUIRE(log->NumericDataCount == 3);

	cinterop::disposal::dispose_of(*log);
	delete log;

}

TEST_CASE("URS port for perf analysis", "[optimizer]") {
	setDefaults();
	using mhcpp::optimization::UniformRandomSamplingOptimizer;
	GIVEN("A 2D Hypercube")
	{
		auto terminationCondition = CreateCounterTermination<Urs>(10);
		Hc goal;
		UniformRandomSamplingOptimizer<Hc> opt = CreateUrsQuadraticGoal(goal, terminationCondition);
		WHEN("optimizing single-thread") {
			opt.UseMultiThreading(false);
			auto results = opt.Evolve();
			REQUIRE(results.size() > 0);
			//results.PrintTo(std::cout);
			auto first = results[0];
			REQUIRE(first.ObjectiveCount() == 1);
		}
		WHEN("optimizing multi-threaded") {
			opt.UseMultiThreading(true);
			auto results = opt.Evolve();
			REQUIRE(results.size() > 0);
			//results.PrintTo(std::cout);
			auto first = results[0];
			REQUIRE(first.ObjectiveCount() == 1);
		}
	}
}

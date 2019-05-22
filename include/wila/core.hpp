#pragma once

#include <string>
#include <vector>
#include <atomic>
#include <map>
#include <random>
#include <functional>
#include "random.hpp"
#include "utils.hpp"
#include "core.h"
#include "constructs.hpp"

using namespace std;

namespace mhcpp
{
	using namespace mhcpp::random;
	using namespace mhcpp::utils;

	template<typename T = double>
	class IHyperCubeSetBounds : public IHyperCube < T > //where T : IComparable
	{
	public:
		virtual ~IHyperCubeSetBounds() {}
		virtual void SetMinValue(const string& variableName, T value) = 0;
		virtual void SetMaxValue(const string& variableName, T value) = 0;
		virtual void SetMinMaxValue(const string& variableName, T min, T max, T value) = 0;
	};

	template<typename TSysConfig>
	class ICandidateFactory
	{
	public:
		virtual ~ICandidateFactory() {}
		virtual TSysConfig CreateRandomCandidate() = 0;
		virtual TSysConfig CreateRandomCandidate(const TSysConfig& point) = 0;
		virtual TSysConfig CreateRandomCandidate(const std::vector<TSysConfig>& points) = 0;
		virtual ICandidateFactory<TSysConfig>* CreateNew() = 0;
		vector<TSysConfig> CreateRandomCandidates(size_t n)
		{
			vector<TSysConfig> v(n);
			for (size_t i = 0; i < n; i++)
				v[i] = CreateRandomCandidate();
			return v;
		}

	};

	template<typename TSysConfig>
	class ICandidateFactorySeed
	{
	public:
		virtual ICandidateFactory<TSysConfig>* Create() const = 0;
		virtual ICandidateFactorySeed<TSysConfig>* Clone() const = 0;
		virtual ~ICandidateFactorySeed() {/*Nothing*/};
	};

	/// <summary>
	/// A generic interface for one or more objective scores derived from the evaluation of a candidate system configuration.
	/// </summary>
	/// <typeparam name="TSysConf">The type of the system configuration</typeparam>
	template<typename TSysConf>
	class IObjectiveScores : //: public IBaseObjectiveScores //where TSysConf : ISystemConfiguration
		public SetOfScores
	{
	public:
		IObjectiveScores(const TSysConf& sysConfig, const string& name, double value, bool maximizable = false)
		{
			this->sys = sysConfig;
			this->objectives.push_back(ObjectiveValue(name, value, maximizable));
		}

		IObjectiveScores() {}

		virtual ~IObjectiveScores() {}

		IObjectiveScores(const IObjectiveScores<TSysConf>& src)
		{
			this->sys = src.sys;
			this->objectives = src.objectives;
		}

		IObjectiveScores(IObjectiveScores<TSysConf>&& src)
		{
			std::swap(this->sys, src.sys);
			std::swap(this->objectives, src.objectives);
		}


		IObjectiveScores<TSysConf>& operator=(const IObjectiveScores<TSysConf> &src)
		{
			if (&src == this){
				return *this;
			}
			this->sys = src.sys;
			this->objectives = src.objectives;
			return *this;
		}

		IObjectiveScores<TSysConf>& operator=(IObjectiveScores<TSysConf>&& src)
		{
			if (&src == this){
				return *this;
			}
			std::swap(this->sys, src.sys);
			std::swap(this->objectives, src.objectives);
			return *this;
		}

		/// <summary>
		/// Gets the system configuration that led to these scores.
		/// </summary>
		virtual TSysConf SystemConfiguration() const { return sys; }

		/// <summary>
		/// Gets the number of objective scores in this instance.
		/// </summary>
		virtual size_t ObjectiveCount() const { return this->objectives.size(); }

		virtual std::string ObjectiveName(size_t i) const { return this->objectives[i].Name; }

		virtual vector<string> ObjectiveNames() const {
			vector<string> result;
			size_t n = ObjectiveCount();
			for (size_t i = 0; i < n; i++)
				result.push_back(ObjectiveName(i));
			return result;
		}

		virtual double Value(size_t i) const { return objectives[i].Value; } //= 0;

		virtual bool Maximizable(size_t i) const { return objectives[i].Maximizable; } //= 0;

		string ToString() const
		{
			string s;
			for (size_t i = 0; i < this->ObjectiveCount(); i++)
			{
				s += this->objectives[i].ToString();
				s += ", ";
			}
			s += SystemConfiguration().ToString();
			return s;
		}

		map<string, double> GetObjectiveValues() const
		{
			map<string, double> s;
			for (size_t i = 0; i < this->ObjectiveCount(); i++)
				s[this->ObjectiveName(i)] = this->Value(i);
			return s;
		}

		map<string, double> GetParameterValues() const
		{
			return SystemConfiguration().GetValues();
		}

		static std::vector<TSysConf> GetSystemConfigurations(const std::vector<IObjectiveScores<TSysConf>>& points)
		{
			std::vector<TSysConf> result;
			for (size_t i = 0; i < points.size(); i++)
			{
				result.push_back(points[i].SystemConfiguration());
			}
			return result;
		}

		static void Sort(std::vector<IObjectiveScores<TSysConf>>& points, const string& scoreName)
		{
			if (points.size() == 0)
				return; // to be confirmed, but can make sense.
			auto names = points.at(0).ObjectiveNames();
			int index = -1;
			for (size_t i = 0; i < names.size(); i++)
				if (names[i] == scoreName) {
					index = i; break;
				}
			if (index < 0) throw std::logic_error(string("Score not found: ") + scoreName);

			std::function<bool(const IObjectiveScores<TSysConf>& elem1, const IObjectiveScores<TSysConf>& elem2)> comparator =
			[&index](const IObjectiveScores<TSysConf>& elem1, const IObjectiveScores<TSysConf>& elem2)
				{
					return IObjectiveScores<TSysConf>::BetterThan(elem1, elem2, index);
				};

			std::stable_sort(points.begin(), points.end(), comparator);
		}

	private:
		static bool BetterThan(const IObjectiveScores<TSysConf>& elem1, const IObjectiveScores<TSysConf>& elem2, int index)
		{
			bool lessThan = (elem1.Value(index) < elem2.Value(index));
			bool moreThan = (elem1.Value(index) > elem2.Value(index));
			if(elem1.Maximizable(index))
				return moreThan;
			return lessThan;
		}

		class ObjectiveValue
		{
		public:
			ObjectiveValue(){}

			ObjectiveValue(string name, double value, bool maximizable) :
				Name(name), Value(value), Maximizable(maximizable)
			{
			}

			ObjectiveValue(const ObjectiveValue& src) :
				Name(src.Name), Value(src.Value), Maximizable(src.Maximizable)
			{
			}
			string Name;
			double Value;
			bool Maximizable;
			string ToString() const
			{
				return Name + ":" + mhcpp::utils::ToString(Value);
			}
		};

		TSysConf sys;
		std::vector<ObjectiveValue> objectives;

	};

	namespace utils {
		template<typename T>
		void PrintTo(const std::vector<T>& scores, std::ostream& stream)
		{
			int n = scores.size();
			for (size_t i = 0; i < n; ++i)
				stream << i << ": " << scores[i].ToString() << std::endl;
		}
	}

	template<typename T>
	class IOptimizationResults // : public std::vector < IObjectiveScores<T> >
		//where T : ISystemConfiguration
	{
		std::vector<IObjectiveScores<T>> scores;

	public:
		typedef typename std::vector<IObjectiveScores<T>>::const_iterator const_iterator;
		IOptimizationResults() {}

		IOptimizationResults(const std::vector<IObjectiveScores<T>>& scores)
		{
			this->scores = scores;
		}

		IOptimizationResults(const IOptimizationResults& src)
		{
			this->scores = src.scores;
		}

		IOptimizationResults(const IOptimizationResults&& src)
		{
			this->scores = std::move(src.scores);
		}

		IOptimizationResults& operator=(const IOptimizationResults& src)
		{
			if (&src == this) {
				return *this;
			}
			this->scores = src.scores;
			return *this;
		}

		IOptimizationResults& operator=(const IOptimizationResults&& src)
		{
			if (&src == this) {
				return *this;
			}
			this->scores = std::move(src.scores);
			return *this;
		}

		operator std::vector<IObjectiveScores<T>>()
		{
			return this->scores;
		}

		size_t size() const
		{
			return scores.size();
		}

		const_iterator begin() const
		{
			return scores.begin();
		}

		const_iterator end() const
		{
			return scores.end();
		}

		const IObjectiveScores<T>& operator[](size_t index) const
		{
			return scores[index];
		}

		void PrintTo(std::ostream& stream)
		{
			mhcpp::utils::PrintTo<IObjectiveScores<T>>(scores, stream);
		}
	};

	template<typename TSysConfig>
	class UniformRandomSamplingFactory : 
		public ICandidateFactory < TSysConfig >
	{
	public:
		UniformRandomSamplingFactory(const IRandomNumberGeneratorFactory<>& rng, const TSysConfig& t)
		{
			this->rng = rng;
			unsigned int samplerSeed = UniformRandomSamplingFactory::CreateSamplerSeed(rng);
			//if (!t.SupportsThreadSafloning)
			//	throw new ArgumentException("This URS factory requires cloneable and thread-safe system configurations");
			this->t = t;
			//this->hcOps = CreateIHyperCubeOperations();
			SetSampler(samplerSeed);
		}

		UniformRandomSamplingFactory(const UniformRandomSamplingFactory& src)
		{
			this->rng = src.rng;
			this->t = src.t;
			unsigned int samplerSeed = CreateSamplerSeed(rng);
			SetSampler(samplerSeed);
		}

		UniformRandomSamplingFactory(UniformRandomSamplingFactory&& src)
		{
			this->rng = std::move(src.rng);
			this->t = std::move(src.t);
			unsigned int samplerSeed = CreateSamplerSeed(rng);
			SetSampler(samplerSeed);
		}

		UniformRandomSamplingFactory& operator=(const UniformRandomSamplingFactory &src)
		{
			if (&src == this) {
				return *this;
			}
			this->rng = src.rng;
			this->t = src.t;
			return *this;
		}

		UniformRandomSamplingFactory& operator=(UniformRandomSamplingFactory&& src)
		{
			if (&src == this) {
				return *this;
			}
			this->rng = std::move(src.rng);
			this->t = std::move(src.t);
			return *this;
		}

		~UniformRandomSamplingFactory() { }

		bool Equals(const UniformRandomSamplingFactory& other) const
		{
			if (&other == this) {
				return true;
			}
			bool samplerRngEquals = sampler.RngEngineEquals(other.sampler);
			bool rngEquals = this->rng.Equals(other.rng);
			return rngEquals && samplerRngEquals;
		}

		//IHyperCubeOperations CreateNew(const IRandomNumberGeneratorFactory& rng)
		//{
		//	return new HyperCubeOperations(rng.CreateFactory());
		//}

		TSysConfig CreateRandomCandidate()
		{
			return CreateRandomCandidate(t);
		}

		TSysConfig CreateRandomCandidate(const TSysConfig& bounds)
		{
			TSysConfig rt(t);
			for (auto& vname : rt.GetVariableNames())
			{
				double min = bounds.GetMinValue(vname);
				double max = bounds.GetMaxValue(vname);
				double d = Urand();
				//if (d > 1 || d < 0)
				//	throw std::logic_error("[0-1] uniform distribution, but got a sample out of this interval");
				if (abs(min - max) <= (std::numeric_limits<double>::epsilon()*100.) )
				{
					rt.SetValue(vname, min);
				} else {
					rt.SetValue(vname, min + d * (max - min));
				}
			}
			return rt;
		}

		TSysConfig CreateRandomCandidate(const vector<TSysConfig>& population)
		{
			if (population.size() == 0)
				throw std::logic_error("There must be at least one point in the population to define a feasible space for the random point");
			TSysConfig bounds(t);
			for (auto& vname : bounds.GetVariableNames())
			{
				double min = std::numeric_limits<double>::lowest();		// should this be std::numeric_limits<T>::lowest() ?
				double max = std::numeric_limits<double>::max();		// should this be std::numeric_limits<T>::max() ?
				for (size_t i = 0; i < population.size(); i++)
				{
					double x = population[i].GetValue(vname);
					min = std::max(min, x);
					max = std::min(max, x);
				}
				bounds.SetMaxValue(vname, max);
				bounds.SetMinValue(vname, min);
			}
			auto result = CreateRandomCandidate(bounds);
#ifdef _DEBUG
			//CheckParameterFeasible(result);
#endif
			return result;
		}

		ICandidateFactory<TSysConfig>* CreateNew()
		{
			return new UniformRandomSamplingFactory<TSysConfig>(this->rng.CreateNew(), this->t);
		}

	private:

		static unsigned int CreateSamplerSeed(const IRandomNumberGeneratorFactory<>& rng)
		{
			auto tmp_rng = rng;
			unsigned int s = tmp_rng.Next();
			unsigned int samplerSeed = s ^ 0xFFFFFFFF;
			return samplerSeed;
		}

		void SetSampler(unsigned int seed)
		{
			mhcpp::random::uniform_real_distribution_double dist(0, 1);
			sampler = rng.CreateVariateGenerator<mhcpp::random::uniform_real_distribution_double>(dist, seed);
		}

		double Urand()
		{
			return sampler();
		}

		RngReal<> sampler;
		IRandomNumberGeneratorFactory<> rng;
		TSysConfig t;
	};

	template<typename TSysConfig, typename TCandidateFactory = UniformRandomSamplingFactory<TSysConfig>>
	class CandidateFactorySeed : public ICandidateFactorySeed<TSysConfig>
	{
	private:
		unsigned int seed;
		TSysConfig templateConfig;
	public:
		CandidateFactorySeed()
		{
			this->seed = 0;
		}

		~CandidateFactorySeed()
		{
			// Nothing
		}

		CandidateFactorySeed(const CandidateFactorySeed& src)
		{
			this->seed = src.seed;
			this->templateConfig = src.templateConfig;
		}

		CandidateFactorySeed(unsigned int seed, const TSysConfig& templateConfig)
		{
			this->seed = seed;
			this->templateConfig = templateConfig;
		}

		ICandidateFactory<TSysConfig>* Create() const
		{
			auto rng = IRandomNumberGeneratorFactory<>(seed);
			return new TCandidateFactory(rng, templateConfig);
		}
		
		ICandidateFactorySeed<TSysConfig>* Clone() const
		{
			return new CandidateFactorySeed<TSysConfig, TCandidateFactory>(*this);
		}
	};


	template<typename T>
	class IEvolutionEngine
	{
	public:
		virtual IOptimizationResults<T> Evolve() = 0;
	};

	template<typename T, typename TEngine = IEvolutionEngine<T>>
	class TerminationCheck
	{
	protected:
		TerminationCheck() {}

	public:
		virtual void Reset() { }
		virtual bool IsFinished(TEngine* engine) = 0;
		virtual bool IsThreadSafe() = 0;
		virtual TerminationCheck* Clone() const = 0;
		virtual ~TerminationCheck()
		{
		}
	};

	template<typename TSys, typename TEngine = IEvolutionEngine<TSys>>
	class CounterTestFinished : public TerminationCheck<TSys, TEngine>
	{
		int counter = 0;
		int maxChecks = 0;
	public:
		CounterTestFinished(int maxChecks)
		{
			this->maxChecks = maxChecks;
		}
		void Reset() { counter = 0; }

		bool IsFinished(TEngine* engine)
		{
			counter++;
			return (counter >= maxChecks);
		}

		bool IsThreadSafe() { return true; }; // TODO: probably not, unless counter is made atomic...

		TerminationCheck<TSys, TEngine>* Clone() const
		{
			// TOCHECK: is this the behavior we want (think parallel operations)
			return new CounterTestFinished(maxChecks);
		}
	};

	template<typename T, typename TEngine = IEvolutionEngine<T>>
	class ITerminationCondition
	{
	private:
		class AlwaysFinished : public TerminationCheck < T, TEngine >
		{
		public:
			bool IsFinished(TEngine* engine) { return true; }
			bool IsThreadSafe() { return true; };
			TerminationCheck < T, TEngine >* Clone() const { return new AlwaysFinished(); };

			~AlwaysFinished() {}
		};

		void DeleteCheck()
		{
			if (Check != nullptr) {
				delete Check;
				Check = nullptr;
			}
		}

	public:
		ITerminationCondition()
		{
			this->Check = new AlwaysFinished();
		}

		ITerminationCondition(const TerminationCheck<T,TEngine>& isFinishedFunc)
		{
			DeleteCheck();
			this->Check = isFinishedFunc.Clone();
		}

		ITerminationCondition(const ITerminationCondition& src)
		{
			engine = src.engine;
			Check = src.Check->Clone();
		}

		ITerminationCondition(ITerminationCondition&& src)
		{
			engine = std::move(src.engine);
			Check = std::move(src.Check);
			src.engine = nullptr;
			src.Check = nullptr;
		}

		ITerminationCondition& operator=(ITerminationCondition&& src)
		{
			if (&src == this){
				return *this;
			}
			engine = std::move(src.engine);
			Check = std::move(src.Check);
			src.engine = nullptr;
			src.Check = nullptr;
			return *this;
		}

		ITerminationCondition& operator=(const ITerminationCondition& src)
		{
			if (&src == this){
				return *this;
			}
			engine = src.engine;
			Check = src.Check->Clone();
			return *this;
		}

		virtual ~ITerminationCondition()
		{
			DeleteCheck();
		}

		void SetEvolutionEngine(TEngine* engine) { this->engine = engine; };

		bool RequireEngine = true;

		bool IsFinished()
		{
			if (engine == nullptr && RequireEngine)
				throw std::logic_error("The optimization engine is not yet set - cannot check for termination criterion");
			return Check->IsFinished(engine);
		}

		bool IsThreadSafe()
		{
			return Check->IsThreadSafe();
 		}

		void Reset()
		{
			Check->Reset();
		}
	private:
		TerminationCheck<T, TEngine> * Check = nullptr;
		//std::function<bool(TEngine*)> Check;
		TEngine* engine = nullptr;
	};

	/// <summary>
	/// Capture a fitness score derived from a candidate system configuration and its objective scores.
	/// </summary>
	/// <typeparam name="T">The type of fitness used to compare system configuration.</typeparam>
	template<typename T, typename TSys>
 	class FitnessAssignedScores //: IComparable<FitnessAssignedScores<T>> where T : IComparable
	{
	public:
		/// <summary>
		/// Creates a FitnessAssignedScores, a union of a candidate system configuration and its objective scores, and an overall fitness score.
		/// </summary>
		/// <param name="scores">Objective scores</param>
		/// <param name="fitnessValue">Fitness value, derived from the scores and context information such as a candidate population.</param>
		FitnessAssignedScores(const IObjectiveScores<TSys>& scores, T fitnessValue)
		{
			this->scores = scores;
			this->fitnessValue = fitnessValue;
		}

		FitnessAssignedScores()
		{
			this->fitnessValue = T();
		}

		FitnessAssignedScores(const FitnessAssignedScores<T, TSys>& src)
		{
			this->scores = src.scores;
			this->fitnessValue = src.fitnessValue;
		}

		FitnessAssignedScores(FitnessAssignedScores<T, TSys>&& src)
		{
			std::swap(this->scores, src.scores);
			std::swap(this->fitnessValue, src.fitnessValue);
		}


		FitnessAssignedScores<T, TSys>& operator=(const FitnessAssignedScores<T, TSys> &src)
		{
			if (&src == this){
				return *this;
			}
			this->scores = src.scores;
			this->fitnessValue = src.fitnessValue;
			return *this;
		}

		FitnessAssignedScores<T, TSys>& operator=(FitnessAssignedScores<T, TSys>&& src)
		{
			if (&src == this){
				return *this;
			}
			std::swap(this->scores, src.scores);
			std::swap(this->fitnessValue, src.fitnessValue);
			return *this;
		}

		/// <summary>
		/// Gets the objective scores
		/// </summary>
		IObjectiveScores<TSys> Scores() const { return scores; }

		static std::vector<IObjectiveScores<TSys>> GetScores(const std::vector<FitnessAssignedScores<double, TSys>>& points)
		{
			std::vector<IObjectiveScores<TSys>> result;
			for (size_t i = 0; i < points.size(); i++)
			{
				result.push_back(points[i].Scores());
			}
			return result;
		}

		static std::vector<TSys> GetSystemConfigurations(const std::vector<FitnessAssignedScores<double, TSys>>& points)
		{
			std::vector<TSys> result;
			for (size_t i = 0; i < points.size(); i++)
			{
				result.push_back(points[i].Scores().SystemConfiguration());
			}
			return result;
		}

		/// <summary>
		/// Gets the fitness value that has been assigned to the candidate system configuration and its objective scores
		/// </summary>
		T FitnessValue() const { return fitnessValue; }

		/// <summary>
		/// Compares two FitnessAssignedScores<T>.
		/// </summary>
		/// <param name="other">Object to compare with this object</param>
		/// <returns>an integer as necessary to implement IComparable</returns>
		int CompareTo(const FitnessAssignedScores<T, TSys>& other) const
		{
			int result;
			if (FitnessValue() < other.FitnessValue())
				result = -1;
			else if (FitnessValue() == other.FitnessValue())
				result = 0;
			else
				result = 1;
			return result;
		}

		int CompareTo(const FitnessAssignedScores<T, TSys>* other)
		{
			return CompareTo(*other);
		}

		static bool BetterThan(const FitnessAssignedScores<T, TSys>& elem1, const FitnessAssignedScores<T, TSys>& elem2)
		{
			if (elem1.CompareTo(elem2) < 0)
				return true;
			return false;
		}

		static bool BetterThanPtr(const FitnessAssignedScores<T, TSys>* elem1, const FitnessAssignedScores<T, TSys>* elem2)
		{
			return BetterThan(*elem1, *elem2);
		}

		static void Sort(std::vector<FitnessAssignedScores<T, TSys>*>& points)
		{
			std::stable_sort(points.begin(), points.end(), FitnessAssignedScores<T, TSys>::BetterThanPtr);
		}

		static void Sort(std::vector<const FitnessAssignedScores<T, TSys>*>& points)
		{
			std::stable_sort(points.begin(), points.end(), FitnessAssignedScores<T, TSys>::BetterThanPtr);
		}

		static void Sort(std::vector<FitnessAssignedScores<T, TSys>>& points)
		{
			std::stable_sort(points.begin(), points.end(), FitnessAssignedScores<T, TSys>::BetterThan);
		}

		string ToString() const
		{
			return FitnessValue().ToString() + ", " + Scores().ToString();
		}

	private:

		IObjectiveScores<TSys> scores;
		T fitnessValue;

	};

	template<typename TVal, typename TSys>
	class IFitnessAssignment
	{
	public:
		std::vector<FitnessAssignedScores<TVal, TSys>> AssignFitness(const std::vector<IObjectiveScores<TSys>>& scores)
		{
			std::vector<FitnessAssignedScores<TVal, TSys>> result;
			// for (auto& s : scores) // <== causes a debug assertion failure with VS - "vector iterator not incrementable"
			for (size_t i = 0; i < scores.size(); i++)
			{
				if (scores[i].ObjectiveCount() != 1)
					throw std::logic_error("Fitness score currently must be derived exactly from one objective");
				TVal scoreVal = scores[i].Value(0);
				result.push_back(FitnessAssignedScores<TVal, TSys>(scores[i], scores[i].Maximizable(0) ? -scoreVal : scoreVal));
			}
			return result;
		}
	};

	template<typename TSysConf>
	class IObjectiveEvaluator
	{
	public:
		virtual ~IObjectiveEvaluator() {}
		/// <summary>
		/// Evaluate the objective values for a candidate system configuration
		/// </summary>
		/// <param name="systemConfiguration">candidate system configuration</param>
		/// <returns>An object with one or more objective scores</returns>
		virtual IObjectiveScores<TSysConf> EvaluateScore(const TSysConf& systemConfiguration) = 0;
		virtual bool IsCloneable() const { return false; }
		virtual IObjectiveEvaluator * Clone() const 
		{ 
			throw std::logic_error(string("Clone operation is not supported by default for ") + typeid(IObjectiveEvaluator<TSysConf>).name());
		}

		static std::vector<IObjectiveScores<TSysConf>> EvaluateScores(IObjectiveEvaluator<TSysConf>& evaluator, const vector<TSysConf>& systemConfigurations)
		{
			std::vector<IObjectiveScores<TSysConf>> result;
			for (size_t i = 0; i < systemConfigurations.size(); i++)
				result.push_back(evaluator.EvaluateScore(systemConfigurations.at(i)));
			return result;
		}

	};

	template<typename T, typename TSys>
	class IPopulation
	{
	public:
		virtual std::vector<FitnessAssignedScores<T, TSys>> Population() = 0;
	};

	template<typename THyperCube>
	THyperCube GetCentroid(const std::vector<THyperCube>& points)
	{
		// TODO: surely some vector libraries to reuse (Boost?)
		if (points.size() == 0) throw std::logic_error("Cannot take centroid of empty set of points");
		vector<string> names = points[0].GetVariableNames();
		vector<double> coords(names.size());
		coords.assign(coords.size(), 0);
		for (auto& p : points)
			for (size_t i = 0; i < coords.size(); i++)
				coords[i] += p.GetValue(names[i]);
		for (size_t i = 0; i < coords.size(); i++)
			coords[i] /= points.size();
		THyperCube centroid = points[0];
		for (size_t i = 0; i < coords.size(); i++)
		{
			centroid.SetValue(names[i], coords[i]);
		}
		return centroid;
	}

	template<typename T=double>
	double Sdev(const std::vector<T>& v) {
		T sum = std::accumulate(std::begin(v), std::end(v), 0.0);
		T m = sum / v.size();

		T accum = 0.0;
		std::for_each(std::begin(v), std::end(v), [&](const T d) {
			accum += (d - m) * (d - m);
		});
		return std::sqrt(accum / (v.size() - 1));
	}

	template<typename THyperCube>
	std::map<string,double> GetStdDev(const std::vector<THyperCube>& points)
	{
		size_t n = points.size();
		if (n <= 1) throw std::logic_error("Cannot take standard deviation of empty or unary set of points");
		vector<string> names = points[0].GetVariableNames();
		vector<double> sdevs(names.size());
		vector<double> v(n);
		std::map<string, double> out;
		for (size_t k = 0; k < names.size(); k++)
		{
			for (size_t i = 0; i < n; i++)
				v[i] = points[i].GetValue(names[k]);
			out[names[k]] = Sdev(v);
		}
		return out;
	}

	template<typename THyperCube>
	std::map<string, double> GetRelativeSdev(const std::vector<THyperCube>& points)
	{
		std::map<string, double> sdevs = GetStdDev<THyperCube>(points);
		auto p = points[0];
		vector<string> names = p.GetVariableNames();
		std::map<string, double> out;
		for (size_t i = 0; i < names.size(); i++)
		{
			string k = names[i];
			double delta = std::abs(p.GetMaxValue(k) - p.GetMinValue(k));
			out[k] = sdevs[k] / delta; // Yes, may lead to nan...
		}
		return out;
	}

	template<typename T>
	T Reflect(T point, T reference, T factor)
	{
		return reference + ((point - reference) * factor);
	}

	template<typename THyperCube, typename T>
	THyperCube HomotheticTransform(const THyperCube& centre, const THyperCube& from, double factor)
	{
		THyperCube result(from);
		auto varnames = centre.GetVariableNames();
		for (auto& v : varnames)
		{
			T mx = from.GetMaxValue(v);
			T mn = from.GetMinValue(v);
			T epsilon;
			if (mx != 0.0)
			{
				// https://stackoverflow.com/questions/9999221/double-precision-decimal-places
				// TODO check if macros are available DBL_DIG? Other std::numeric_limits things?
				epsilon = std::abs((T)1e-13 * mx);
			}
			else
			{
				// HACK: use an absolute value tolerance even for zero valued bounds
				// This is potentially problematic but less likely to be an issue than letting the reflection go ahead unchecked.
				// TODO: rethink. 
				epsilon = (T)1e-13;
			}
			if (std::abs(mx - mn) > epsilon) {
				double newVal = mhcpp::Reflect(from.GetValue(v), centre.GetValue(v), factor);
				result.SetValue(v, newVal);
			}
		}
		return result;
	}

	template<typename T>
	class HyperCube : public IHyperCubeSetBounds < T > //where T : IComparable
		, public InstanceCounter<HyperCube<T>>
	{
	public:
		HyperCube() {}

		virtual ~HyperCube() {}

		HyperCube(const HyperCube& src)
		{
			this->def = src.def;
		}
		vector<string> GetVariableNames() const {
			return mhcpp::utils::GetKeys(def);
		}
		void Define(string name, T min, T max, T value) {
			def[name] = mhcpp::parameters::ParameterInterval<double>(name, min, max, value);
		}
		size_t Dimensions() const { return def.size(); }
		T GetValue(const string& variableName) const { return def.at(variableName).Value; }
		T GetMaxValue(const string& variableName) const { return def.at(variableName).Max; }
		T GetMinValue(const string& variableName) const { return def.at(variableName).Min; }
		void SetValue(const string& variableName, T value)    { def[variableName].Value = value; }
		void SetMinValue(const string& variableName, T value) { def[variableName].Min = value; }
		void SetMaxValue(const string& variableName, T value) { def[variableName].Max = value; }
		void SetMinMaxValue(const string& variableName, T min, T max, T value) {
			SetMinValue(variableName, min);
			SetMaxValue(variableName, max);
			SetValue(variableName, value);
		}
		//IHyperCube<T> HomotheticTransform(IHyperCube<T> point, double factor) {}
		string GetConfigurationDescription() const { return ""; }
		void ApplyConfiguration(void* system) {}

		static HyperCube<T> GetCentroid(const std::vector<HyperCube<T>>& points)
		{
			return mhcpp::GetCentroid<HyperCube<T>>(points);
		}

		HyperCube<T> HomotheticTransform(const HyperCube<T>& from, T factor)
		{
			return mhcpp::HomotheticTransform<HyperCube<T>, T>(*this, from, factor);
		}

		virtual bool IsFeasible() const
		{
			auto varnames = GetVariableNames();
			for (auto& v : varnames)
			{
				T mx = this->def.at(v).Max;
				T mn = this->def.at(v).Min;
				if (std::abs(mx - mn) >= (std::numeric_limits<T>::epsilon() * 100))
				{
					if (!this->def.at(v).IsFeasible()) return false;
				}
			}
			return true;
		}


	private:
		std::map<string, mhcpp::parameters::ParameterInterval<double>> def;
	};

	template<typename TSysConf>
	class TopologicalDistance : public IObjectiveEvaluator < TSysConf >
	{
	public:
		TopologicalDistance(const TopologicalDistance& src) 
		{ 
			this->goal = src.goal; 
		}
		TopologicalDistance(const TSysConf& goal) { this->goal = goal; }
		~TopologicalDistance() {}

		/// <summary>
		/// Evaluate the objective values for a candidate system configuration
		/// </summary>
		/// <param name="systemConfiguration">candidate system configuration</param>
		/// <returns>An object with one or more objective scores</returns>
		IObjectiveScores<TSysConf> EvaluateScore(const TSysConf& systemConfiguration)
		{
			double sumsqr = 0;
			vector<string> varNames = goal.GetVariableNames();
			for (auto& v : varNames)
			{
				auto d = goal.GetValue(v) - systemConfiguration.GetValue(v);
				sumsqr += d*d;
			}
			return IObjectiveScores<TSysConf>(systemConfiguration, "L2 distance", std::sqrt(sumsqr));
		}

		bool IsCloneable() const 
		{ 
			return true;
		}

		IObjectiveEvaluator<TSysConf> * Clone() const
		{
			return new TopologicalDistance(*this);
		}

	private:
		TSysConf goal;
	};

	namespace utils{

		template<typename TSysConf = HyperCube<double>>
		bool CheckParameterFeasible(const IObjectiveScores<TSysConf>& s)
		{
			if (!s.SystemConfiguration().IsFeasible())
			{
				string ps = s.ToString();
				std::cout << ps;
				return false;
			}
			return true;
		}

		template<typename TSysConf = HyperCube<double>>
		bool CheckParameterFeasible(const TSysConf& p)
		{
			if (!p.IsFeasible())
			{
				string ps = p.ToString();
				std::cout << ps;
				return false;
			}
			return true;
		}

		template<typename TSysConf = HyperCube<double>>
		bool CheckParameterFeasible(const std::vector<IObjectiveScores<TSysConf>>& vec)
		{
			for (size_t i = 0; i < vec.size(); i++)
			{
				if (!CheckParameterFeasible(vec[i])) return false;
			}
			return true;
		}

		template<typename TSysConf = HyperCube<double>>
		bool CheckParameterFeasible(const FitnessAssignedScores<double, TSysConf>& p)
		{
			return CheckParameterFeasible(p.Scores());
		}

		template<typename TSysConf = HyperCube<double>>
		bool CheckParameterFeasible(const std::vector<FitnessAssignedScores<double, TSysConf>>& vec)
		{
			for (size_t i = 0; i < vec.size(); i++)
			{
				if (!CheckParameterFeasible(vec[i])) return false;
			}
			return true;
		}
	}

}


#pragma once

#include <random>
#include <numeric>
#include <set>
#include <boost/random.hpp>
#include <boost/random/uniform_real_distribution.hpp>

namespace mhcpp
{
	namespace random
	{
		using namespace std;

		// In order to solve WIRADA-392 we cannot use the distributions in std::
		// as they behave differently across compilers.
		typedef boost::random::uniform_real_distribution<double> uniform_real_distribution_double;
		typedef boost::random::discrete_distribution<int> discrete_distribution_int;
		typedef boost::random::normal_distribution<double> normal_distribution_double;

		// using typedef for GCC backward compatibility - replace with using statements later.
		typedef boost::mt19937 default_wila_random_engine;

		template<typename RNG, typename Distribution>
		class VariateGenerator
		{
		public:
			typedef typename Distribution::result_type result_type;

			VariateGenerator() { }

			VariateGenerator(const Distribution& d) : _dist(d) { }

			/**
			* Constructs a @c VariateGenerator object with the associated
			* \uniform_random_number_generator eng and the associated
			* \random_distribution d.
			*
			* Throws: If and what the copy constructor of RNG or
			* Distribution throws.
			*/

			VariateGenerator(const RNG& e, const Distribution& d)
				: _eng(e), _dist(d) { }

			VariateGenerator(const VariateGenerator& vg)
				: _eng(vg._eng), _dist(vg._dist) { }

			VariateGenerator(VariateGenerator&& vg)
			{
				this->_eng = std::move(vg._eng);
				this->_dist = std::move(vg._dist);
			}

		private:
			template<typename RNG1>
			void Clean(typename enable_if<std::is_pointer<RNG1>::value, RNG1>::type e)
			{
				delete e;
			}

			template<typename RNG1>
			void Clean(typename enable_if<!std::is_pointer<RNG1>::value, RNG1>::type e) {}
		public:

			~VariateGenerator()
			{
				Clean<RNG>(this->_eng);
			}

			VariateGenerator& operator=(const VariateGenerator& vg)
			{
				if (&vg == this) {
					return *this;
				}
				this->_eng = vg._eng;
				this->_dist = vg._dist;
				return *this;
			}

			VariateGenerator& operator=(VariateGenerator&& vg)
			{
				if (&vg == this) {
					return *this;
				}
				this->_eng = std::move(vg._eng);
				this->_dist = std::move(vg._dist);
				return *this;
			}

			bool RngEngineEquals(const VariateGenerator& other) const
			{
				if (&other == this) {
					return true;
				}
				return _eng == (other._eng);
			}

			result_type operator()() { return _dist(_eng); }

			const Distribution& distribution()
			{
				return _dist;
			}

			vector<result_type> Generate(size_t n, std::function<double(double)>& f)
			{
				vector<result_type> v(n);
				auto g = [&]()
				{
					result_type x = (*this)();
					return f(x);
				};
				std::generate(v.begin(), v.end(), g);
				return v;
			}

			vector<result_type> Generate(size_t n)
			{
				auto identity = [&](double x) { return x; };
				return Generate(n, identity);
			}

		private:
			RNG _eng;
			Distribution _dist;
		};

#if defined(__GNUC__) && !defined(__clang__)
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#if GCC_VERSION < 40800
#define WILA_USE_TYPEDEF_RNG
#endif
#endif

#ifdef WILA_USE_TYPEDEF_RNG
		typedef VariateGenerator < default_wila_random_engine, mhcpp::random::discrete_distribution_int > RngInt;
		typedef VariateGenerator < default_wila_random_engine, mhcpp::random::uniform_real_distribution_double > RngReal;
		typedef VariateGenerator < default_wila_random_engine, mhcpp::random::normal_distribution_double > RngNorm;
#else
		/*template <typename RNG = default_wila_random_engine>
		typedef VariateGenerator < RNG, mhcpp::random::discrete_distribution_int > RngInt;
*/
		template <typename RNG = default_wila_random_engine>
		using RngInt = VariateGenerator < RNG, mhcpp::random::discrete_distribution_int >;

		template <typename RNG = default_wila_random_engine>
		using RngReal = VariateGenerator < RNG, mhcpp::random::uniform_real_distribution_double >;

		template <typename RNG = default_wila_random_engine>
		using RngNorm = VariateGenerator < RNG, mhcpp::random::normal_distribution_double >;
#endif

		template<typename RNG = default_wila_random_engine>
		class IRandomNumberGeneratorFactory
		{

		private:
			/** \brief	A random number generator factory. */

			template<typename RngEngine = default_wila_random_engine>
			class RandomNumberGeneratorFactory
			{
			private:
				RngEngine seedEngine;
			public:
				typedef RngEngine engine_type; // No typename needed here. See http://stackoverflow.com/questions/6489351/nested-name-specifier
				//	http://stackoverflow.com/questions/495021/why-can-templates-only-be-implemented-in-the-header-file

				RandomNumberGeneratorFactory(unsigned int seed) : seedEngine(seed)
				{
				}

				RandomNumberGeneratorFactory() : seedEngine(0)
				{
				}

				RandomNumberGeneratorFactory(const RandomNumberGeneratorFactory& src)
				{
					this->seedEngine = src.seedEngine;
				}

				RandomNumberGeneratorFactory(RandomNumberGeneratorFactory&& src)
				{
					this->seedEngine = std::move(src.seedEngine);
				}

				RandomNumberGeneratorFactory& operator = (const RandomNumberGeneratorFactory& src)
				{
					if (&src == this) {
						return *this;
					}
					this->seedEngine = src.seedEngine;
					return *this;
				}

				RandomNumberGeneratorFactory& operator = (RandomNumberGeneratorFactory&& src)
				{
					if (&src == this) {
						return *this;
					}
					this->seedEngine = std::move(src.seedEngine);
					return *this;
				}

				virtual ~RandomNumberGeneratorFactory()
				{
				}

				bool Equals(const RandomNumberGeneratorFactory& other) const
				{
					return seedEngine == (other.seedEngine);
				}

				RngEngine * CreateNewEngine()
				{
					return new RngEngine(seedEngine());
				}

				RandomNumberGeneratorFactory<RngEngine> * CreateNewFactory()
				{
					return new RandomNumberGeneratorFactory<RngEngine>(seedEngine());
				}

				unsigned int operator()() {
					return seedEngine();
				}

				template<class DistributionType = mhcpp::random::uniform_real_distribution_double>
				static VariateGenerator<RngEngine*, DistributionType> * CreateVariateGenerator(RandomNumberGeneratorFactory<RngEngine>& rngf, DistributionType& dist)
				{
					return new VariateGenerator<RngEngine*, DistributionType>(rngf.CreateNewEngine(), dist);
				}

			};

			RandomNumberGeneratorFactory<RNG> rng;

		public:
			IRandomNumberGeneratorFactory() : rng(0)
			{
			}

			IRandomNumberGeneratorFactory(unsigned int seed) : rng(seed)
			{
			}

			IRandomNumberGeneratorFactory(const IRandomNumberGeneratorFactory& src)
			{
				this->rng = src.rng;
			}

			IRandomNumberGeneratorFactory(IRandomNumberGeneratorFactory&& src)
			{
				this->rng = std::move(src.rng);
			}

			IRandomNumberGeneratorFactory& operator = (const IRandomNumberGeneratorFactory& src)
			{
				if (&src == this) {
					return *this;
				}
				this->rng = src.rng;
				return *this;
			}

			IRandomNumberGeneratorFactory& operator = (IRandomNumberGeneratorFactory&& src)
			{
				if (&src == this) {
					return *this;
				}
				this->rng = std::move(src.rng);
				return *this;
			}

			bool Equals(const IRandomNumberGeneratorFactory& other) const
			{
				return rng.Equals(other.rng);
			}

			unsigned int Next()
			{
				return rng();
			}

			unsigned int PeekNext() const
			{
				RandomNumberGeneratorFactory<RNG> rngTmp(rng);
				return rngTmp();
			}

			std::vector<unsigned int> Next(size_t size)
			{
				std::vector<unsigned int> result;
				for (size_t i = 0; i < size; i++)
				{
					result.push_back(Next());
				}
				return result;
			}

			IRandomNumberGeneratorFactory CreateNew()
			{
				return IRandomNumberGeneratorFactory(Next());
			}

			default_wila_random_engine CreateNewStd()
			{
				unsigned int seed = Next();
				return CreateNewStd(seed);
			}

			static default_wila_random_engine CreateNewStd(unsigned int seed)
			{
				return default_wila_random_engine(seed);
			}

			RNG CreateNewEngine()
			{
				unsigned int seed = Next();
				return CreateNewEngine(seed);
			}

			static RNG CreateNewEngine(unsigned int seed)
			{
				return RNG(seed);
			}

			// http://stackoverflow.com/questions/7166799/specialize-template-function-for-template-class
			template<class DistributionType = mhcpp::random::uniform_real_distribution_double>
			VariateGenerator<RNG, DistributionType> CreateVariateGenerator(DistributionType& dist, unsigned int seed) const
			{
				return VariateGenerator<RNG, DistributionType>(CreateNewEngine(seed), dist);
			}

			template<class DistributionType = mhcpp::random::uniform_real_distribution_double>
			VariateGenerator<RNG, DistributionType> CreateVariateGenerator(DistributionType& dist)
			{
				return VariateGenerator<RNG, DistributionType>(CreateNewEngine(), dist);
			}
		};

		template<typename RNG = default_wila_random_engine>
		RngInt<RNG> CreateTrapezoidalRng(size_t n, const RNG& generator, double trapezoidalFactor = -1)
		{
			std::vector<double> weights(n);
			double m = (double)n;
			double sumWeights = n*(n + 1) / 2.0;
			double avgWeight = sumWeights / m;

			if ((trapezoidalFactor <= 0) || (trapezoidalFactor >= 2))
			{
				// default as per the original SCE paper. Note that we do not need
				// to normalize: std::discrete_distribution takes care of it
				for (size_t i = 1; i <= n; i++)
					weights[i - 1] = (double)(n + 1 - i);
			}
			else
			{
				// y = ax + b
				double b = trapezoidalFactor; // y(0) = b
				double a = 2 * (1 - b) / (m - 1);  // because y(n-1) = (2-b) = a * (n-1) + b
				for (size_t i = 0; i < n; i++)
					weights[i] = a * i + b;
			}
			// Would expect to be able to do:
			// mhcpp::random::discrete_distribution_int distribution(weights.begin(), weights.end());
			// but missing constructor in MS implementation. Using workaround derived from
			// http://stackoverflow.com/questions/21959404/initialising-stddiscrete-distribution-in-vs2013
			std::size_t i(0);
			mhcpp::random::discrete_distribution_int distribution(weights.size(),
				0.0, // dummy!
				1.0, // dummy!
				[&weights, &i](double)
			{
				auto w = weights[i];
				++i;
				return w;
			});
			return VariateGenerator<default_wila_random_engine, mhcpp::random::discrete_distribution_int>(generator, distribution);
		}

		template<typename X>
		static std::vector<const X*> AsPointers(const std::vector<X>& vec)
		{
			std::vector<const X*> result;
			for (size_t i = 0; i < vec.size(); i++)
			{
				const X& e = vec[i];
				result.push_back(&e);
			}
			return result;
		}

		template<typename RNG = default_wila_random_engine>
		vector<int> SampleFrom(RngInt<RNG>& drng, size_t nsampling)
		{
			size_t size = drng.distribution().max() - drng.distribution().min() + 1;
			std::vector<int>p(size);

			for (size_t i = 0; i < nsampling; ++i) {
				int number = drng();
				++p[number];
			}
			return p;
		}

		template<typename ElemType>
		vector<ElemType> SampleFrom(RngInt<>& drng, const std::vector<ElemType>& population, size_t n,
			vector<ElemType>& leftOut, bool replace = false)
		{
			if (!replace && population.size() <= n)
				throw std::logic_error("If elements are sampled once, the output size must be less than the population sampled from");

			std::set<int> selected;
			std::set<int> notSelected;
			for (size_t i = 0; i < population.size(); i++)
				notSelected.emplace(i);
			std::vector<ElemType> result(n);
			if (replace)
			{
				for (size_t i = 0; i < n; i++)
				{
					result[i] = population[drng()];
					selected.emplace(i);
					notSelected.erase(i);
				}
			}
			else
			{
				auto src = AsPointers(population);
				int counter = 0;
				while (counter < n)
				{
					int i = drng();
					if (!(src[i] == nullptr))
					{
						result[counter] = *(src[i]);
						src[i] = nullptr;
						selected.emplace(i);
						notSelected.erase(i);
						counter++;
					}
#ifdef _DEBUG
					//CheckParameterFeasible(population);
#endif
				}
				leftOut.clear();
				for (auto index : notSelected)
					leftOut.push_back(population[index]);
#ifdef _DEBUG
				//CheckParameterFeasible(population);
				//CheckParameterFeasible(leftOut);
#endif
			}
			return result;
		}

	}
}

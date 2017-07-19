#pragma once

#include <string>

#include <Rcpp.h>
#include <wila/sce.h>
#include <wila/interop_c_structs.h>
#include <wila/interop_c_cpp.hpp>


#define SCE_ALPHA  "alpha"
#define SCE_BETA "beta"
#define SCE_CONTRACTION_RATIO "contraction_ratio"
#define SCE_M "m"
#define SCE_MAX_NUM_SHUFFLE "max_num_shuffle"
#define SCE_P "p"
#define SCE_PMIN "pmin"
#define SCE_Q "q"
#define SCE_REFLECTION_RATIO "reflection_ratio"
#define SCE_TRAPEZOIDAL_DENSITY_FACTOR "trapezoidal_density_factor"

#define PARAMETERIZER_NAME_ITEM_NAME "Name"
#define PARAMETERIZER_MIN_ITEM_NAME "Min"
#define PARAMETERIZER_MAX_ITEM_NAME "Max"
#define PARAMETERIZER_VALUE_ITEM_NAME "Value"
#define	STRINGSASFACTORS_ITEM_NAME "stringsAsFactors"

using namespace Rcpp;

namespace mhcpp {
	namespace interop {
		namespace r {
			template <typename V=NumericVector>
			SceParameters ToSceParameters(const V& sceParams)
			{
				// TODO: checks presence of param names.
				// std::vector<string> names = (std::vector<string>)(sceParams.names());
				SceParameters s;
				s.Alpha = (int)sceParams[SCE_ALPHA];
				s.Beta = (int)sceParams[SCE_BETA];
				s.ContractionRatio = sceParams[SCE_CONTRACTION_RATIO];
				s.M = (int)sceParams[SCE_M];
				s.NumShuffle = (int)sceParams[SCE_MAX_NUM_SHUFFLE];
				s.P = (int)sceParams[SCE_P];
				s.Pmin = (int)sceParams[SCE_PMIN];
				s.Q = (int)sceParams[SCE_Q];
				s.ReflectionRatio = sceParams[SCE_REFLECTION_RATIO];
				s.TrapezoidalDensityParameter = sceParams[SCE_TRAPEZOIDAL_DENSITY_FACTOR];

				return s;
			}

			template <typename D = DataFrame>
			void GetColumns(const D& parameterSpecs, CharacterVector& pnames, NumericVector& minima, NumericVector& maxima, NumericVector& values)
			{
				pnames = Rcpp::as<CharacterVector>(parameterSpecs[PARAMETERIZER_NAME_ITEM_NAME]);
				minima = Rcpp::as<NumericVector>(parameterSpecs[PARAMETERIZER_MIN_ITEM_NAME]);
				maxima = Rcpp::as<NumericVector>(parameterSpecs[PARAMETERIZER_MAX_ITEM_NAME]);
				values = Rcpp::as<NumericVector>(parameterSpecs[PARAMETERIZER_VALUE_ITEM_NAME]);
			}

			template<typename P = mhcpp::interop::ParameterSetDefinition<>>
			DataFrame HypercubeToDataFrame(P& pdef)
			{
				return Rcpp::DataFrame::create(
					Named(PARAMETERIZER_NAME_ITEM_NAME) = pdef.Names,
					Named(PARAMETERIZER_MIN_ITEM_NAME) = pdef.Mins,
					Named(PARAMETERIZER_MAX_ITEM_NAME) = pdef.Maxs,
					Named(PARAMETERIZER_VALUE_ITEM_NAME) = pdef.Values,
					Named(STRINGSASFACTORS_ITEM_NAME) = false);
			}


		}
	}
}


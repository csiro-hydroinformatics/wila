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

			template<typename P = mhcpp::interop::ParameterSetDefinition<>>
			P HypercubeFromDataFrame(const DataFrame& df)
			{
				using namespace mhcpp::interop;
				vector<string> names = Rcpp::as<vector<string>>(Rcpp::as<CharacterVector>(df[PARAMETERIZER_NAME_ITEM_NAME]));
				vector<double> mins = Rcpp::as<vector<double>>(Rcpp::as<NumericVector>(df[PARAMETERIZER_MIN_ITEM_NAME]));
				vector<double> maxs = Rcpp::as<vector<double>>(Rcpp::as<NumericVector>(df[PARAMETERIZER_MAX_ITEM_NAME]));
				vector<double> values = Rcpp::as<vector<double>>(Rcpp::as<NumericVector>(df[PARAMETERIZER_VALUE_ITEM_NAME]));
				return P(names, mins, maxs, values);
			}

			template<typename P = OptimizerLogData>
			DataFrame OptimizerLogToDataFrame(P& logData)
			{

				auto to_vec = [](char** values, int size) 
				{ 
					return cinterop::utils::to_cpp_string_vector<string>(values, size, false);
				};

				DataFrame df = DataFrame::create();
				int numRows = logData.LogLength;
				if (numRows == 0) return df;

				for (int i = 0; i < logData.StringDataCount; i++)
					df.push_back(to_vec(logData.StringData[i], numRows));

				for (int i = 0; i < logData.NumericDataCount; i++)
				{
					double* numArray = logData.NumericData[i];
					vector<double> vecArray(numRows);
					std::copy(numArray, numArray + numRows, vecArray.begin());
					df.push_back(vecArray);
				}

				vector<string> names = to_vec(logData.NamesStringData, logData.StringDataCount);
				vector<string> tmp = to_vec(logData.NamesNumericData, logData.NumericDataCount);
				for (size_t i = 0; i < tmp.size(); i++)
					names.push_back(tmp[i]);
				df.names() = names;

				return df;
			}

			template <typename S=SceParameters>
			S CreateDefaultSceParameters()
			{
				const int n = 4;
				const int nshuffle = 40;

				SceParameters result;
				result.P = n + 2;
				result.Pmin = n + 2;
				result.M = 2 * n + 1;
				result.Q = std::max(result.M - 2, 2);
				result.Alpha = 1;
				result.Beta = result.M;
				result.NumShuffle = nshuffle;
				result.TrapezoidalDensityParameter = 1.0;
				result.ReflectionRatio = -1.0;
				result.ContractionRatio = 0.5;

				return result;
			}

			template <typename V = NumericVector>
			V FromSceParameters(const SceParameters& s)
			{
				V sceParams = V::create(
					Named(SCE_ALPHA) = s.Alpha,
					Named(SCE_BETA) = s.Beta,
					Named(SCE_CONTRACTION_RATIO) = s.ContractionRatio,
					Named(SCE_M) = s.M,
					Named(SCE_MAX_NUM_SHUFFLE) = s.NumShuffle,
					Named(SCE_P) = s.P,
					Named(SCE_PMIN) = s.Pmin,
					Named(SCE_Q) = s.Q,
					Named(SCE_REFLECTION_RATIO) = s.ReflectionRatio,
					Named(SCE_TRAPEZOIDAL_DENSITY_FACTOR) = s.TrapezoidalDensityParameter
				);
				return sceParams;
			}


		}
	}
}


#pragma once

#include <stdexcept>
#include <algorithm>
#include "../include/wila/sce.h"

namespace mhcpp
{
	namespace optimization
	{
		SceParameters CreateSceParamsForProblemOfDimension(int n, int nshuffle)
		{
			if (n <= 0)
				throw new std::logic_error("There must be at least one free parameter to calibrate");
			SceParameters result = CreateDefaultSceParams();
			result.NumShuffle = nshuffle;
			return AdjustSceParamsForProblemOfDimension(result, n);
		}

		SceParameters AdjustSceParamsForProblemOfDimension(const SceParameters& sceParams, int n)
		{
			if (n <= 0)
				throw new std::logic_error("There must be at least one free parameter to calibrate");
			SceParameters result = sceParams;
			result.M = 2 * n + 1;
			result.Q = std::max(result.M - 2, 2);
			result.Beta = result.M;

			return result;
		}

		SceParameters CreateDefaultSceParams()
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
			result.TrapezoidalDensityParameter = 1.9;
			result.ReflectionRatio = -1.0;
			result.ContractionRatio = 0.5;

			return result;
		}
	}
}

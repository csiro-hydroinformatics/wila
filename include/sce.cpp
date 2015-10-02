#pragma once

#include <stdexcept>
#include <algorithm>
#include "sce.h"

namespace mhcpp
{
	namespace optimization
	{
		SceParameters SceParameters::CreateForProblemOfDimension(int n, int nshuffle)
		{
			if (n <= 0)
				throw new std::logic_error("There must be at least one free parameter to calibrate");
			SceParameters result;
			result.M = 2 * n + 1;
			result.Q = std::max(result.M - 2, 2);
			result.Alpha = 1;
			result.Beta = result.M;
			result.NumShuffle = nshuffle;
			return result;
		}
	}
}

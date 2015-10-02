#pragma once

namespace mhcpp
{
	namespace optimization
	{
		class SceParameters
		{
		public:
			SceParameters()
			{
				this->TrapezoidalDensityParameter = 1.9;
				this->ReflectionRatio = -1.0;
				this->ContractionRatio = 0.5;
			}
			/// <summary>
			/// Number of geometrical transformation for each subcomplex
			/// </summary>
			int Alpha;
			/// <summary>
			/// Number of evolution steps taken by each sub-complex before shuffling occurs
			/// </summary>
			int Beta;
			/// <summary>
			/// Number of complexes
			/// </summary>
			int P = 5;
			/// <summary>
			/// Minimum number of complexes (populations of points)
			/// </summary>
			int Pmin = 3;
			/// <summary>
			/// Number of points per complex
			/// </summary>
			int M;
			/// <summary>
			/// Number of points per SUB-complex
			/// </summary>
			int Q;
			int NumShuffle;
			double TrapezoidalDensityParameter;

			/// <summary>
			/// The homothetic ratio used in the reflection phase of the complex evolution: default -1.0
			/// </summary>
			double ReflectionRatio;

			/// <summary>
			/// The homothetic ratio used in the contraction phase of the complex evolution: default 0.5
			/// </summary>
			double ContractionRatio;

			static SceParameters CreateForProblemOfDimension(int n, int nshuffle);
		};

		enum SceOptions
		{
			None = 0x00,
			ReflectionRandomization = 0x01,
			RndInSubComplex = 0x02,
			FutureOption_1 = 0x04,
			FutureOption_2 = 0x08
		};
	}
}

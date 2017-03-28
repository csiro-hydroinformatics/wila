#pragma once

#ifdef _WIN32
#define EXPORT_DLL_LIB __declspec(dllexport)
#else
#define EXPORT_DLL_LIB // nothing
#endif

typedef struct _SceParameters
{
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
	int P;
	/// <summary>
	/// Minimum number of complexes (populations of points)
	/// </summary>
	int Pmin;
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

} SceParameters;

namespace mhcpp
{
	namespace optimization
	{
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

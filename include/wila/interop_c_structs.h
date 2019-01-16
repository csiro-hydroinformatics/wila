#pragma once

#include <cstddef> // size_t
#include <string.h>

#ifndef STRDUP
#ifdef _WIN32
#define STRDUP _strdup
#else
#define STRDUP strdup
#endif
#endif // !STRDUP

typedef struct _parameter_interval
{
	char* name;
	double value;
	double minimum;
	double maximum;
} parameter_interval;

typedef struct _hypercube_parameter_set
{
	size_t size;
	parameter_interval* parameters;
} hypercube_parameter_set;


struct OptimizerLogData
{
	int LogLength;
	int StringDataCount;
	int NumericDataCount;
	char** NamesNumericData;
	char** NamesStringData;
	double** NumericData;
	char*** StringData;
};





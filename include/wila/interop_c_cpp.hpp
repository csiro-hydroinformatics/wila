#pragma once

#include <cinterop/common_c_interop.h>
#include <cinterop/c_cpp_interop.hpp>
#include <cinterop/object_lifetimes.hpp>
#include "interop_c_structs.h"
#include "core.hpp"
#include "logging.hpp"

namespace cinterop
{
	namespace utils {

		template<>
		inline named_values_vector to_named_values_vector<mhcpp::SetOfScores>(const mhcpp::SetOfScores& x)
		{
			size_t n = x.ObjectiveCount();
			named_values_vector vv;
			vv.size = n;
			vv.names = new char*[n];
			vv.values = new double[n];
			for (size_t i = 0; i < n; i++)
			{
				vv.names[i] = STRDUP(x.ObjectiveName(i).c_str());
				vv.values[i] = x.Value(i);
			}
			return vv;
		}

		template<typename T=double>
		hypercube_parameter_set to_hypercube_parameter_set(const mhcpp::IHyperCube<T>& h)
		{
			hypercube_parameter_set hh;
			vector<string> v = h.GetVariableNames();
			size_t n = v.size();
			hh.size = n;
			hh.parameters = new parameter_interval[n];
			for (size_t i = 0; i < n; i++)
			{
				auto p = &hh.parameters[i];
				p->name = STRDUP(v[i].c_str());
				p->minimum = h.GetMinValue(v[i]);
				p->maximum = h.GetMaxValue(v[i]);
				p->value = h.GetValue(v[i]);
			}
			return hh;
		}

	}

	namespace disposal {
		template<>
		inline void dispose_of<hypercube_parameter_set>(hypercube_parameter_set& d)
		{
			for (size_t i = 0; i < d.size; i++)
			{
				auto p = d.parameters[i];
				if (p.name != nullptr)
				{
					delete[] p.name;
					p.name = nullptr;
				}
			}
			delete[] d.parameters;
		}

		template<>
		inline void dispose_of<OptimizerLogData>(OptimizerLogData& d)
		{
			for (size_t i = 0; i < d.NumericDataCount; i++)
				delete[] d.NumericData[i];
			delete[] d.NumericData;

			cinterop::disposal::free_c_ptr_array<char>(d.NamesNumericData, d.NumericDataCount);
			cinterop::disposal::free_c_ptr_array<char>(d.NamesStringData, d.StringDataCount);
			for (size_t i = 0; i < d.StringDataCount; i++)
				cinterop::disposal::free_c_ptr_array<char>(d.StringData[i], d.LogLength);
			delete[] d.StringData;

		}


	}

}

namespace mhcpp
{
	namespace interop {

		struct ParameterDefinition
		{
			string Name;
			double Min, Max, Value;
		};

		template<typename P = ParameterDefinition>
		class ParameterSetDefinition
		{
		public:
			ParameterSetDefinition() {
			}

			ParameterSetDefinition(
				const vector<string>& names,
				const vector<double>& mins,
				const vector<double>& maxs,
				const vector<double>& values
			) : Names(names), Mins(mins), Maxs(maxs), Values(values)
			{
			}

			ParameterSetDefinition(const std::vector<P>& definitions)
			{
				for (size_t i = 0; i < definitions.size(); i++)
				{
					auto& pd = definitions[i];
					Names.push_back(pd.Name);
					Mins.push_back(pd.Min);
					Maxs.push_back(pd.Max);
					Values.push_back(pd.Value);
				}
			}

			ParameterSetDefinition(const hypercube_parameter_set& ps)
			{
				for (size_t i = 0; i < ps.size; i++)
				{
					const auto& pd = ps.parameters[i];
					Names.push_back(string(pd.name));
					Mins.push_back(pd.minimum);
					Maxs.push_back(pd.maximum);
					Values.push_back(pd.value);
				}
			}

			hypercube_parameter_set* AsInteropStructPtr()
			{
				hypercube_parameter_set* result = new hypercube_parameter_set();
				result->size = this->Names.size();
				result->parameters = new parameter_interval[result->size];
				for (size_t i = 0; i < result->size; i++)
				{
					auto p = &result->parameters[i];
					p->name = STRDUP(this->Names[i].c_str());
					p->minimum = this->Mins[i];
					p->maximum = this->Maxs[i];
					p->value = this->Values[i];
				}
				return result;
			}

			std::vector<string> Names;
			std::vector<double> Mins, Maxs, Values;
		};

		template<typename T = double>
		void transfer_values_from_hypercube_parameter_set(mhcpp::IHyperCubeSetBounds<T>& h, const hypercube_parameter_set& hh, bool setBounds=true, bool strict=true)
		{
			ParameterSetDefinition<> pdef(hh);
			size_t n = pdef.Names.size();

			auto hasItem = [&](const string& s, const vector<string>& v)
			{
				return (std::find(v.begin(), v.end(), s) != v.end());
			};

			auto hNames = h.GetVariableNames();

			for (size_t i = 0; i < n; i++)
			{
				string name = pdef.Names[i];
				if (strict && !hasItem(name, hNames))
					throw std::logic_error(string("Parameter name ") + name + string(" not found in the hypercube to set"));
			}

			for (size_t i = 0; i < n; i++)
			{
				string name = pdef.Names[i];
				if (!strict)
					if (!hasItem(name, hNames))
							continue; // otherwise below will throw an exception, but first loop would have caught the condition anyway.

				if (setBounds)
					h.SetMinMaxValue(name, pdef.Mins[i], pdef.Maxs[i], pdef.Values[i]);
				else
					h.SetValue(name, pdef.Values[i]);
			}
		}


		template<typename T = double>
		T** vector_vector_to_c_array(const vector<vector<T>>& values)
		{
			int* size = nullptr;
			return cinterop::utils::vector_to_c_array<vector<T>, T*>(values, cinterop::utils::vector_identity_to_c_array<T>, size);
		}

		template<typename S = string>
		char*** str_vector_vector_to_char_arrays(const vector<vector<S>>& data, int* size)
		{
			std::function<char**(const vector<S>& x)> conv = [](const vector<S>& x) {return cinterop::utils::to_ansi_char_array<S>(x); };
			return cinterop::utils::vector_to_c_array<vector<S>, char**>(data, conv, size);
		}

		template<typename P>
		OptimizerLogData* get_optimizer_log_data(mhcpp::logging::ILoggerMh<P>& logger)
		{
			OptimizerLogData* result = new OptimizerLogData();
			result->LogLength = logger.GetLength();
			std::map<string, vector<string>> strData = logger.GetStringData();
			std::map<string, vector<double>> numData = logger.GetNumericData();

			vector<string> numKeys = mhcpp::utils::GetKeys<string, vector<double>>(numData);
			vector<string> strKeys = mhcpp::utils::GetKeys<string, vector<string>>(strData);
			vector<vector<double>> numValues = mhcpp::utils::GetValues<string, vector<double>>(numData, numKeys);
			vector<vector<string>> strValues = mhcpp::utils::GetValues<string, vector<string>>(strData, strKeys);
			result->NumericDataCount = numKeys.size();
			result->StringDataCount = strKeys.size();
			result->NamesNumericData = cinterop::utils::to_ansi_char_array<string>(numKeys);
			result->NamesStringData = cinterop::utils::to_ansi_char_array<string>(strKeys);

			result->NumericData = mhcpp::interop::vector_vector_to_c_array<double>(numValues);
			result->StringData = mhcpp::interop::str_vector_vector_to_char_arrays(strValues, &result->StringDataCount);

			return result;
		}

	}
}
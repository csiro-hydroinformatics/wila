#pragma once

#include "interop_c_structs.h"
#include <cinterop/common_c_interop.h>
#include <cinterop/c_cpp_interop.hpp>
#include <cinterop/object_lifetimes.hpp>
#include "core.hpp"

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
			ParameterSetDefinition() {}
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
					const auto& pd = parameters[i];
					Names.push_back(string(pd.name));
					Mins.push_back(pd.minimum);
					Maxs.push_back(pd.maximum);
					Values.push_back(pd.value);
				}
			}

			std::vector<string> Names;
			std::vector<double> Mins, Maxs, Values;
		};
	}
}
#pragma once

#include <string>
#include <vector>

using namespace std;

namespace mhcpp
{
	namespace parameters
	{
		template<typename RealType>
		class ParameterInterval
		{
		public:
			ParameterInterval() {}
			ParameterInterval(string name, RealType min, RealType max, RealType value) :
				Name(name), Max(max), Min(min), Value(value)
			{
			}

			ParameterInterval& operator=(const ParameterInterval& src)
			{
				if (&src == this) {
					return *this;
				}
				this->Name =  src.Name;
				this->Min =   src.Min;
				this->Max =   src.Max;
				this->Value = src.Value;
				return *this;
			}
			ParameterInterval& operator=(ParameterInterval&& src)
			{
				if (&src == this) {
					return *this;
				}
				std::swap(this->Name, src.Name);
				this->Min = src.Min;
				this->Max = src.Max;
				this->Value = src.Value;
				return *this;
			}
			ParameterInterval(const ParameterInterval& src) { *this = src; }
			ParameterInterval(ParameterInterval&& src) { *this = src; }

			string Name;
			RealType Min, Max, Value;
			bool IsFeasible() const { return (Value >= Min) && (Value <= Max); }
		};
	}
}


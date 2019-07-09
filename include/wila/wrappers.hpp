#pragma once

#include <string>
#include <vector>
#include <atomic>
#include <map>
#include <random>
#include <functional>
#include "core.hpp"


namespace mhcpp
{
	namespace types
	{
		using std::string;
		template<typename Type, typename BaseType>
		static Type* As(BaseType* ptr)
		{
			return dynamic_cast<Type*>(ptr);
		}
		template<typename Type, typename BaseType>
		static bool Is(BaseType* ptr)
		{
			return (As<Type>(ptr) != nullptr);
		}

		template<typename Type, typename BaseType>
		static Type* StrictlyAs(BaseType* ptr, const string& errormsg = string(""), const string& nullargmsg = string(""))
		{
			if (ptr == nullptr)
			{
				string msgnull;
				if (nullargmsg == "")
					msgnull = "StrictlyAs: input pointer must not be null";
				else
					msgnull = nullargmsg;
				throw std::logic_error(msgnull);
			}
			Type* res = dynamic_cast<Type*>(ptr);
			if (res == nullptr)
			{
				string msg;
				if (errormsg == "")
					msg = string("StrictlyAs: input pointer ") +
					typeid(BaseType*).name() +
					" points to object of type " +
					typeid(*ptr).name() +
					" which cannot be addressed by pointer of type " +
					typeid(Type*).name();
				else
					msg = errormsg;
				throw std::logic_error(msg);
			}
			return res;
		}
	}

	namespace wrappers
	{
		template<typename HyperCubeParameterizer>
		class HypercubeWrapper :
			public HyperCubeParameterizer,
			public mhcpp::IHyperCubeSetBounds<double>
		{
		private:
			HyperCubeParameterizer* def = nullptr;
		public:
			HypercubeWrapper()
			{
				this->def = nullptr;
			}

			HypercubeWrapper(const HyperCubeParameterizer& def)
			{
				this->def = def.Clone();
			}

			HypercubeWrapper(const HypercubeWrapper& src)
			{
				if (src.def != nullptr)
					this->def = src.def->Clone();
			}

			HypercubeWrapper(HypercubeWrapper&& src)
			{
				this->def = src.def;
				src.def = nullptr;
			}

			virtual ~HypercubeWrapper()
			{
				if (this->def != nullptr)
				{
					delete def;
					def = nullptr;
				}
			}

			HypercubeWrapper& operator= (const HypercubeWrapper& src)
			{
				if (&src == this) {
					return *this;
				}
				if (src.def != nullptr)
					this->def = src.def->Clone();
				return *this;
			}

			HypercubeWrapper& operator= (HypercubeWrapper&& src)
			{
				if (&src == this) {
					return *this;
				}
				this->def = src.def;
				src.def = nullptr;
				return *this;
			}

			HyperCubeParameterizer * InnerParameterizer() const { return def; }

			vector<string> GetVariableNames() const {
				//return def->GetParameterNames();
				return def->GetVariableNames();
			}

			static HypercubeWrapper WrapIfRequired(HyperCubeParameterizer * hcPtr)
			{
				if (hcPtr == nullptr)
					throw std::logic_error("WrapIfRequired: pointer to hypercube must not be null");
				if (mhcpp::types::Is<HypercubeWrapper, HyperCubeParameterizer>(hcPtr))
					return *(mhcpp::types::As<HypercubeWrapper, HyperCubeParameterizer>(hcPtr));
				else
					return HypercubeWrapper(*hcPtr);
			}

			static HypercubeWrapper* WrapIfRequiredNewPtr(HyperCubeParameterizer * hcPtr)
			{
				if (hcPtr == nullptr)
					throw std::logic_error("WrapIfRequired: pointer to hypercube must not be null");
				if (mhcpp::types::Is<HypercubeWrapper, HyperCubeParameterizer>(hcPtr))
					return (mhcpp::types::As<HypercubeWrapper, HyperCubeParameterizer>(hcPtr));
				else
					return new HypercubeWrapper(*hcPtr);
			}

			static HyperCubeParameterizer* UnwrapIfRequiredNewPtr(HyperCubeParameterizer * hcPtr)
			{
				HypercubeWrapper* w = mhcpp::types::As<HypercubeWrapper, HyperCubeParameterizer>(hcPtr);
				if (w != nullptr)
					return(w->InnerParameterizer()->Clone());
				else
					return hcPtr->Clone();
			}

			static HyperCubeParameterizer* UnwrapIfRequired(HyperCubeParameterizer * hcPtr)
			{
				HypercubeWrapper* w = mhcpp::types::As<HypercubeWrapper, HyperCubeParameterizer>(hcPtr);
				if (w != nullptr)
					return(w->InnerParameterizer());
				else
					return hcPtr;
			}

			size_t Dimensions() const { return def->GetNumParameters(); }
			double GetValue(const string& variableName) const { return def->GetValue(variableName); }
			double GetMaxValue(const string& variableName) const { return def->GetMaxValue(variableName); }
			double GetMinValue(const string& variableName) const { return def->GetMinValue(variableName); }

			void SetValue(const string& variableName, double value) { def->SetValue(variableName, value); }
			void SetMinValue(const string& variableName, double value) { def->SetMinValue(variableName, value); }
			void SetMaxValue(const string& variableName, double value) { def->SetMaxValue(variableName, value); }
			void SetMinMaxValue(const string& variableName, double min, double max, double value) {
				SetMinValue(variableName, min);
				SetMaxValue(variableName, max);
				SetValue(variableName, value);
			}
			string GetConfigurationDescription() const { return ""; }
			void ApplyConfiguration(void* system) {/*Nothing?*/ }

			bool SupportsThreadSafeCloning() const { return def->SupportsThreadSafeCloning(); }
			//std::vector<std::string> GetParameterNames() const { return def->GetParameterNames(); }
			bool HasKey(const string &paramName) const { return def->HasKey(paramName); }
			size_t GetNumParameters() const { return def->GetNumParameters(); }
			HyperCubeParameterizer* Clone() const
			{
				return new HypercubeWrapper(*def);
			}

			static HypercubeWrapper GetCentroid(const std::vector<HypercubeWrapper>& points)
			{
				return mhcpp::GetCentroid<HypercubeWrapper>(points);
			}

			HypercubeWrapper HomotheticTransform(const HypercubeWrapper& from, double factor)
			{
				return mhcpp::HomotheticTransform<HypercubeWrapper, double>(*this, from, factor);
			}

			virtual bool IsFeasible() const
			{
				return this->def->IsWithinBounds();
			}

			virtual string ToString() const
			{
				string s;
				return this->def->ToString();
				//auto vnames = GetVariableNames();
				//for (auto& v : vnames)
				//	s += v + ":" + mhcpp::utils::ToString(GetValue(v)) + ", ";
				//return s;
			}

		};

		template<typename HyperCubeParameterizer, typename ObjectiveEvaluator>
		class ObjectiveCalculatorWrapper :
			public IObjectiveEvaluator<HypercubeWrapper<HyperCubeParameterizer>>
		{
		private:
			bool owner = true;
		public:

			typedef HypercubeWrapper<HyperCubeParameterizer> HcWrapper;

			ObjectiveCalculatorWrapper()
			{
			}

			ObjectiveCalculatorWrapper(const ObjectiveCalculatorWrapper& src)
			{
				if (src.def != nullptr)
					this->def = src.def->Clone();
			}

			ObjectiveCalculatorWrapper(ObjectiveEvaluator * objc, bool owner = true)
			{
				this->owner = owner;
				this->def = objc;
			}

			ObjectiveCalculatorWrapper(ObjectiveCalculatorWrapper&& src)
			{
				*this = std::move(src);
			}

			virtual ~ObjectiveCalculatorWrapper()
			{
				if (!owner)
					return;
				if (this->def != nullptr)
				{
					delete def;
					def = nullptr;
				}
			}

			ObjectiveCalculatorWrapper& operator= (ObjectiveCalculatorWrapper&& src)
			{
				if (&src == this) {
					return *this;
				}
				this->def = src.def;
				src.def = nullptr;
				return *this;
			}

			ObjectiveCalculatorWrapper& operator= (const ObjectiveCalculatorWrapper& src)
			{
				if (&src == this) {
					return *this;
				}
				if (src.def != nullptr)
					this->def = src.def->Clone();
				return *this;
			}

			bool IsCloneable() const { return true; }

			IObjectiveEvaluator<HcWrapper> * Clone() const
			{
				return new ObjectiveCalculatorWrapper(*this);
			}

			IObjectiveScores<HcWrapper> EvaluateScore(const HcWrapper& parameterizer)
			{
				auto p = parameterizer.InnerParameterizer();
				if (p == nullptr) throw std::logic_error("Parameterizer is a null pointer - cannot evaluate score");
				double value = def->EvaluateScore(p);
				IObjectiveScores<HcWrapper> res(parameterizer, def->GetName(), value, def->IsMaximizable());
				return res;
			}



		private:
			ObjectiveEvaluator * def = nullptr;

		};




	}
}


#pragma once

namespace mhcpp
{

	using std::string;

	class SetOfScores
	{
	public:
		virtual size_t ObjectiveCount() const = 0;
		virtual std::string ObjectiveName(size_t i) const = 0;
		virtual double Value(size_t i) const = 0;
	};

	class ISystemConfiguration
	{
	public:
		virtual ~ISystemConfiguration() {}
		/// <summary>
		/// Gets an alphanumeric description for this system configuration
		/// </summary>
		virtual string GetConfigurationDescription() const = 0;

		/// <summary>
		/// Apply this system configuration to a compatible system, usually a 'model' in the broad sense of the term.
		/// </summary>
		/// <param name="system">A compatible system, usually a 'model' in the broad sense of the term</param>
		/// <exception cref="ArgumentException">thrown if this system configuration cannot be meaningfully applied to the specified system</exception>
		virtual void ApplyConfiguration(void* system) = 0;
	};

	/// <summary>
	/// Interface for system configurations that are a set of numeric parameters, each with min and max feasible values.
	/// </summary>
	/// <typeparam name="T">A comparable type; typically a float or double, but possibly integer or more esoteric type</typeparam>
	template<typename T = double>
	class IHyperCube : public ISystemConfiguration //where T : IComparable
	{
	public:
		virtual ~IHyperCube() {}
		/// <summary>
		/// Gets the names of the variables defined for this hypercube.
		/// </summary>
		/// <returns></returns>
		virtual vector<string> GetVariableNames() const = 0;

		/// <summary>
		/// Gets the number of dimensions in this hypercube
		/// </summary>
		virtual size_t Dimensions() const = 0;

		/// <summary>
		/// Gets the value for a variable
		/// </summary>
		/// <param name="variableName"></param>
		/// <returns></returns>
		virtual T GetValue(const string& variableName) const = 0;

		/// <summary>
		/// Gets the maximum feasible value for a variable
		/// </summary>
		/// <param name="variableName"></param>
		/// <returns></returns>
		virtual T GetMaxValue(const string& variableName) const = 0;

		/// <summary>
		/// Gets the minimum feasible value for a variable
		/// </summary>
		/// <param name="variableName"></param>
		/// <returns></returns>
		virtual T GetMinValue(const string& variableName) const = 0;

		/// <summary>
		/// Sets the value of one of the variables in the hypercube
		/// </summary>
		/// <param name="variableName"></param>
		/// <param name="value"></param>
		virtual void SetValue(const string& variableName, T value) = 0;

		virtual string ToString() const
		{
			string s;
			auto vnames = GetVariableNames();
			for (auto& v : vnames)
				s += v + ":" + std::to_string(GetValue(v)) + ", ";
			return s;
		}

		virtual std::map<string, double> GetValues() const
		{
			std::map<string, double> s;
			auto vnames = GetVariableNames();
			for (auto& v : vnames)
				s[v] = GetValue(v);
			return s;
		}


	};
}
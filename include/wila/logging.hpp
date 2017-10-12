#pragma once

#include "core.hpp"

namespace mhcpp
{
	namespace logging
	{
		/// <summary>
		/// A facade for logging information from optimisation processes, to avoid coupling to specific frameworks.
		/// </summary>
		template<class TSys>
		class ILoggerMh
		{
		public:
			//virtual void Write(std::vector<IBaseObjectiveScores> scores, const std::map<string, string>& ctags) = 0;
			//virtual void Write(FitnessAssignedScores<double> worstPoint, const std::map<string, string>& ctags) = 0;
			virtual void Write(TSys* newPoint, const std::map<string, string>& ctags) = 0;
			virtual void Write(const string& message, const std::map<string, string>& ctags) = 0;
			virtual void Write(const FitnessAssignedScores<double, TSys>& scores, const std::map<string, string>& tags) = 0;
			virtual void Write(const std::vector<FitnessAssignedScores<double, TSys>>& scores, const std::map<string, string>& ctags) = 0;
			virtual void Write(const std::vector<IObjectiveScores<TSys>>& scores, const std::map<string, string>& tags) = 0;
			virtual void Reset() = 0;

			virtual std::map<string, vector<string>>GetStringData() = 0;
			virtual std::map<string, vector<double>>GetNumericData() = 0;
			virtual int GetLength() = 0;
			virtual ILoggerMh<TSys>* CreateNew() = 0;

			virtual ~ILoggerMh() {}
		};

	}
}

#include "common.h"

using namespace mhcpp;

HyperCube<double> createTestHc(double a, double b, double aMin, double bMin, double aMax, double bMax) {
	HyperCube<double> hc;
	hc.Define("a", aMin, aMax, a);
	hc.Define("b", bMin, bMax, b);
	return hc;
}

bool assertHyperCube(const HyperCube<double>& hc, double a, double b, double tolerance)
{
	return
		(hc.Dimensions() == 2) &&
		(std::abs(hc.GetValue("a") - a) < tolerance) &&
		(std::abs(hc.GetValue("b") - b) < tolerance);
}

bool isIn(const std::string& s, const std::vector<std::string>& b)
{
	return (std::find(b.begin(), b.end(), s) != b.end());
}

bool allIn(const std::vector<std::string>& a, const std::vector<std::string>& b)
{
	for (auto s : a)
	{
		if (!isIn(s, b))
			return false;
	}
	return true;
}

bool sameSets(const std::vector<std::string>& a, const std::vector<std::string>& b)
{
	if (a.size() != b.size()) return false;
	if (!allIn(a, b)) return false;
	if (!allIn(b, a)) return false;
	return true;
}

bool sameValues(const std::vector<std::string>& keys, const HyperCube<double>& a, const HyperCube<double>& b)
{
	for (auto& s : keys)
		if (a.GetValue(s) != b.GetValue(s))
			return false;
	return true;
}

bool sameMinValues(const std::vector<std::string>& keys, const HyperCube<double>& a, const HyperCube<double>& b)
{
	for (auto& s : keys)
		if (a.GetMinValue(s) != b.GetMinValue(s))
			return false;
	return true;
}

bool sameMaxValues(const std::vector<std::string>& keys, const HyperCube<double>& a, const HyperCube<double>& b)
{
	for (auto& s : keys)
		if (a.GetMaxValue(s) != b.GetMaxValue(s))
			return false;
	return true;
}

bool assertEqual(const HyperCube<double>& a, const HyperCube<double>& b)
{
	return
		(a.Dimensions() == b.Dimensions()) &&
		(sameSets(a.GetVariableNames(), b.GetVariableNames())) &&
		(sameValues(a.GetVariableNames(), a, b)) &&
		(sameMinValues(a.GetVariableNames(), a, b)) &&
		(sameMaxValues(a.GetVariableNames(), a, b));
}

bool assertValuesNotEqual(const HyperCube<double>& a, const HyperCube<double>& b)
{
	return
		(a.Dimensions() == b.Dimensions()) &&
		(sameSets(a.GetVariableNames(), b.GetVariableNames())) &&
		(!sameValues(a.GetVariableNames(), a, b));
}

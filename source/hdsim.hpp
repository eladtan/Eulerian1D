#ifndef HDSIM_HPP
#define HDSIM_HPP 1

#include "Primitive.hpp"
#include "Interpolation.hpp"
#include "ideal_gas.hpp"
#include "ExactRS.hpp"
#include "Extensive.hpp"
#include "SourceTerm.hpp"
//#include "Boundary.hpp"
#include "Geometry.hpp"
#include <vector>

using namespace std;

template <class T> vector<T> unique(vector<T> const& v)
{
	std::size_t n = v.size();
	vector<T> res;
	res.reserve(n);
	if (n == 0)
		return res;
	res.push_back(v[0]);
	for (typename vector<T>::const_iterator it = v.begin() + 1; it != v.end(); ++it)
		if (*it == *(it - 1))
			continue;
		else
			res.push_back(*it);
	return res;
}


template <class T> void RemoveVector
(vector<T> &v, vector<size_t> &indeces)
{
	if (indeces.empty())
		return;
	sort(indeces.begin(), indeces.end());
	vector<T> result;
	result.reserve(v.size() - indeces.size());
	int counter = 0;
	for (std::size_t i = 0; i < static_cast<std::size_t>(indeces.back()); ++i)
	{
		if (std::size_t(indeces[std::size_t(counter)]) == i)
			++counter;
		else
			result.push_back(v[i]);
	}
	for (std::size_t i = static_cast<std::size_t>(indeces.back()) + 1; i < v.size(); ++i)
		result.push_back(v[i]);
	v = result;
}

class hdsim
{
private:
	const double cfl_;
	vector<Primitive> cells_;
	vector<double> edges_;
	Interpolation const& interpolation_;
	IdealGas const& eos_;
	RiemannSolver const& rs_;
	double time_;
	size_t cycle_;
	double TotalEcool_;
	vector<pair<Primitive, Primitive> > interp_values_;
	vector<Extensive> fluxes_;
	vector<Extensive> extensives_;
	SourceTerm const& source_;
	Geometry const& geo_;
	const double AMR_ratio_;
	//	BoundarySolution const* BoundarySolution_;
	double dt_suggest_;
	void AMR(void);
public:
	hdsim(double cfl, vector<Primitive> const& cells, vector<double> const& edges, Interpolation const& interp,
		IdealGas const& eos, RiemannSolver const& rs, SourceTerm const& source, Geometry const& geo, const double AMR_ratio = 0);
	~hdsim();
	void TimeAdvance2();
	void TimeAdvance();
	double GetTime()const;
	double GetEcool()const;
	double& GetEcool();
	vector<Primitive>const& GetCells()const;
	vector<Primitive>& GetCells();
	vector<Extensive>const& GetExtensives()const;
	vector<Extensive>& GetExtensives();
	vector<double> const& GetEdges()const;
	vector<double>& GetEdges();
	size_t GetCycle()const;
	void SetTime(double t);
	void SetEcool(double E);
	void SetCycle(size_t cyc);
	void ReCalcCells(vector<Extensive> & extensives);
	void ReCalcExtensives(vector<Primitive> const& cells);
	void SuggestTimeStep(double dt);
};
#endif //HDSIM_HPP

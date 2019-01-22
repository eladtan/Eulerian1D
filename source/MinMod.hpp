#ifndef MINMOD_HPP
#define MINMOD_HPP 1

#include "Interpolation.hpp"
#include "Boundary.hpp"
#include <vector>

using namespace std;

class MinMod : public Interpolation
{
private:
	Boundary const& boundary_;
public:
	MinMod(Boundary const& boundary);
	~MinMod();

	void GetInterpolatedValues(vector<Primitive> const& cells, vector<double> const& edges, vector<pair<Primitive,
		Primitive> > & values,double time)const;
};

#endif

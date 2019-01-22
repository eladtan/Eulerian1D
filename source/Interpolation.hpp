#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP 1

#include "Primitive.hpp"
#include <vector>

using namespace std;

class Interpolation
{
public:
	virtual void GetInterpolatedValues(vector<Primitive> const& cells, vector<double> const& edges, vector<pair<Primitive,
		Primitive> > & values, double time) const = 0;
};

#endif

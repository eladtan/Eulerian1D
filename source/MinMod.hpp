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
	bool SR_;
public:
	MinMod(Boundary const& boundary,bool SR = false);
	~MinMod();

	void GetInterpolatedValues(vector<Primitive> const& cells, vector<double> const& edges, vector<pair<Primitive,
		Primitive> > & values, double time
#ifdef RICH_MPI
		, std::array<Primitive, NGHOSTCELLS * 2> const& ghost_cells,
		std::array<double, 2 * NGHOSTCELLS> const& ghost_edges
#endif
	)const;
};

#endif

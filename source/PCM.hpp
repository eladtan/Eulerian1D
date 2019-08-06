#ifndef PCM_HPP
#define PCM_HPP 1

#include "Interpolation.hpp"
#include "Boundary.hpp"
#include <vector>

using namespace std;

class PCM : public Interpolation
{
private:
	Boundary const& boundary_;
public:
	PCM(Boundary const& boundary);
	~PCM();

	void GetInterpolatedValues(vector<Primitive> const& cells, vector<double> const& edges, vector<pair<Primitive,
		Primitive> > & values, double time
#ifdef RICH_MPI
		, std::array<Primitive, NGHOSTCELLS * 2> const& ghost_cells,
		std::array<double, 2 * NGHOSTCELLS> const& ghost_edges
#endif
	)const;
};

#endif

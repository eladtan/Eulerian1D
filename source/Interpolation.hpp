#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP 1

#include "Primitive.hpp"
#include <vector>
#ifdef RICH_MPI
#include "mpi_comm.hpp"
#endif

using namespace std;

class Interpolation
{
public:
	virtual void GetInterpolatedValues(vector<Primitive> const& cells, vector<double> const& edges, vector<pair<Primitive,
		Primitive> > & values, double time
#ifdef RICH_MPI
		, std::array<Primitive, NGHOSTCELLS * 2> const& ghost_cells,
		std::array<double, 2 * NGHOSTCELLS> const& ghost_edges
#endif
	) const = 0;
};

#endif

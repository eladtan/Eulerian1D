#include "PCM.hpp"
#include <cmath>

PCM::PCM(Boundary const& boundary) :boundary_(boundary)
{}


PCM::~PCM()
{}

void PCM::GetInterpolatedValues(vector<Primitive> const & cells, vector<double> const & edges,
	vector<pair<Primitive, Primitive> >& values, double time) const
{
	size_t N = edges.size();
	values.resize(N);
	// Do bulk edges
	for (size_t i = 0; i < N - 1; ++i)
	{
		values[i].second = cells[i];
		values[i + 1].first = cells[i];
	}
	// Do boundaries
	vector<Primitive> left = boundary_.GetBoundaryValues(cells, edges, 0, time);
	vector<Primitive> right = boundary_.GetBoundaryValues(cells, edges, edges.size() - 1, time);
	values[0].first = left[0];
	values[N - 1].second = right[2];
}

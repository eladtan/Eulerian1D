#include "Boundary.hpp"
#include <algorithm>
#include<cmath>
#include <iostream>
#include <array>

/*BoundarySolution::~BoundarySolution() {}

VacuumInFlow::VacuumInFlow(bool calc_left, bool calc_right) :calc_left_(calc_left), calc_right_(calc_right) {}

pair<RSsolution, RSsolution> VacuumInFlow::GetBoundaryValues(vector<Primitive> const& cells)const
{
	RSsolution left, right;
	left.pressure = 0;
	right.pressure = 0;
	left.velocity = cells[0].velocity;
	right.velocity = cells.back().velocity+3*sqrt(5.*cells.back().pressure/cells.back().density/3.0);
	//right.velocity = 0;
	return pair<RSsolution, RSsolution>(left, right);
}

std::pair<bool, bool> VacuumInFlow::ShouldCalc()const
{
	return pair<bool, bool>(calc_left_, calc_right_);
}
*/
Boundary::~Boundary()
{}

vector<Primitive> RigidWall::GetBoundaryValues(vector<Primitive> const & cells,
	vector<double> const & edges, size_t index,double time, bool SR) const
{
	vector<Primitive> res(3);
	if (index == 0)
	{
		std::array<Primitive, 2> new_cells = { cells[0], cells[1] };
		if (SR)
		{
			for (size_t i = 0; i < 2; ++i)
			{
				double gamma = 1.0 / std::sqrt(1.0 -
					new_cells[i].velocity*new_cells[i].velocity);
				new_cells[i].velocity *= gamma;
			}
		}
		Primitive slope = (new_cells[1] - new_cells[0]) / (0.5*(edges[2] - edges[0]));
		slope.density = 0;
		slope.entropy = 0;
		//slope.velocity = 0;
		slope.pressure = 0;
		res[1] = new_cells[0] - slope*(0.5*(edges[1] - edges[0]));
		res[0] = res[1];
		res[0].velocity = -res[1].velocity;
		res[2] = new_cells[0]+slope*(0.5*(edges[1] - edges[0]));
	}
	else
	{
		size_t N = edges.size();
		std::array<Primitive, 2> new_cells = { cells[N-3], cells[N-2] };
		if (SR)
		{
			for (size_t i = 0; i < 2; ++i)
			{
				double gamma = 1.0 / std::sqrt(1.0 -
					new_cells[i].velocity*new_cells[i].velocity);
				new_cells[i].velocity *= gamma;
			}
		}
		Primitive slope = (new_cells[1] - new_cells[0]) / (0.5*(edges[N-1] - edges[N-3]));
		slope.density = 0;
		slope.entropy = 0;
		//slope.velocity = 0;
		slope.pressure = 0;
		res[1] = new_cells[1] + slope*(0.5*(edges[N-1] - edges[N-2]));
		res[2] = res[1];
		res[2].velocity = -res[1].velocity;
		res[0] = new_cells[1] - slope*(0.5*(edges[N-1] - edges[N-2]));
	}
	if (SR)
	{
		for (size_t i = 0; i < 3; ++i)
		{
			double factor = 1.0 / std::sqrt(1 +
				res[i].velocity*res[i].velocity);
			res[i].velocity *= factor;
		}
	}
	return res;
}

vector<Primitive> Ratchet::GetBoundaryValues(vector<Primitive> const& cells, vector<double> const&
		edges, size_t index, double time, bool SR)const
{
	if (index == 0)
		throw("Not implemented");
	if (cells.back().velocity < 0)
	{
		RigidWall br;
		return br.GetBoundaryValues(cells, edges, index, time,SR);
	}
	else
	{
		FreeFlow br;
		return br.GetBoundaryValues(cells, edges, index, time,SR);
	}
}


vector<Primitive> RigidWall1::GetBoundaryValues(vector<Primitive> const & cells,
	vector<double> const & edges, size_t index, double time, bool SR) const
{
	vector<Primitive> res(3);
	if (index == 0)
	{
		std::array<Primitive, 2> new_cells = { cells[0], cells[1] };
		if (SR)
		{
			for (size_t i = 0; i < 2; ++i)
			{
				double gamma = 1.0 / std::sqrt(1.0 -
					new_cells[i].velocity*new_cells[i].velocity);
				new_cells[i].velocity *= gamma;
			}
		}
		res[1] = new_cells[0];
		res[0] = res[1];
		res[0].velocity = -res[1].velocity;
		res[2] = new_cells[0];
	}
	else
	{
		size_t N = edges.size();
		std::array<Primitive, 2> new_cells = { cells[N-3], cells[N-2] };
		if (SR)
		{
			for (size_t i = 0; i < 2; ++i)
			{
				double gamma = 1.0 / std::sqrt(1.0 -
					new_cells[i].velocity*new_cells[i].velocity);
				new_cells[i].velocity *= gamma;
			}
		}
		res[1] = new_cells[1] ;
		res[2] = res[1];
		res[2].velocity = -res[1].velocity;
		res[0] = new_cells[1];
	}
	if (SR)
	{
		for (size_t i = 0; i < 3; ++i)
		{
			double factor = 1.0 / std::sqrt(1 +
				res[i].velocity*res[i].velocity);
			res[i].velocity *= factor;
		}
	}
	return res;
}


vector<Primitive> FreeFlow::GetBoundaryValues(vector<Primitive> const & cells, vector<double> const & edges,
	size_t index, double time,bool SR) const
{
	vector<Primitive> res(3);
	if (index == 0)
	{
		std::array<Primitive, 2> new_cells = { cells[0], cells[1] };
		if (SR)
		{
			for (size_t i = 0; i < 2; ++i)
			{
				double gamma = 1.0 / std::sqrt(1.0 -
					new_cells[i].velocity*new_cells[i].velocity);
				new_cells[i].velocity *= gamma;
			}
		}
		res[1] = new_cells[0];
		res[0] = res[1];
		res[2] = new_cells[0];
	}
	else
	{
		size_t N = edges.size();
		std::array<Primitive, 2> new_cells = { cells[N-3], cells[N-2] };
		if (SR)
		{
			for (size_t i = 0; i < 2; ++i)
			{
				double gamma = 1.0 / std::sqrt(1.0 -
					new_cells[i].velocity*new_cells[i].velocity);
				new_cells[i].velocity *= gamma;
			}
		}
		res[1] = new_cells[1];
		res[2] = res[1];
		res[0] = new_cells[1];
	/*	const Primitive sl = (cells[i] - cells[i - 1]) / (0.5*(edges[i + 1] - edges[i - 1]));
		const Primitive sr = (cells[i + 1] - cells[i]) / (0.5*(edges[i + 2] - edges[i]));
		const Primitive sc = (cells[i + 1] - cells[i - 1]) / (0.5*(edges[i + 2] + edges[i + 1] - edges[i - 1]*/
	}
	if (SR)
	{
		for (size_t i = 0; i < 3; ++i)
		{
			double factor = 1.0 / std::sqrt(1 +
				res[i].velocity*res[i].velocity);
			res[i].velocity *= factor;
		}
	}
	return res;
}

vector<Primitive> FreeFlow2::GetBoundaryValues(vector<Primitive> const & cells, vector<double> const & edges,
	size_t index, double time,bool SR) const
{
	vector<Primitive> res(3);
	Primitive slope;
	if (index == 0)
	{
		std::array<Primitive, 3> new_cells = { cells[0], cells[1], cells[2] };
		if (SR)
		{
			for (size_t i = 0; i < 3; ++i)
			{
				double gamma = 1.0 / std::sqrt(1.0 -
					new_cells[i].velocity*new_cells[i].velocity);
				new_cells[i].velocity *= gamma;
			}
		}
		size_t i = 1;
		const Primitive sl = (new_cells[i] - new_cells[i - 1]) / (0.5*(edges[i + 1] - edges[i-1]));
		const Primitive sr = (new_cells[i+1] - new_cells[i]) / (0.5*(edges[i + 2] - edges[i]));
		const Primitive sc = (new_cells[i+1] - new_cells[i - 1]) / (0.5*(edges[i + 2]+edges[i + 1] - edges[i - 1]
			- edges[i]));
		if (sl.density*sr.density < 0)
			slope.density = 0;
		else
			slope.density = std::min(std::fabs(sl.density), std::min(std::fabs(sr.density), std::fabs(sc.density))) * (sl.density > 0 ?
				1 : -1);
		if (sl.pressure*sr.pressure < 0)
			slope.pressure = 0;
		else
			slope.pressure = std::min(std::fabs(sl.pressure), std::min(std::fabs(sr.pressure), std::fabs(sc.pressure))) * (sl.pressure > 0 ?
				1 : -1);
		if (sl.velocity*sr.velocity < 0)
			slope.velocity = 0;
		else
			slope.velocity = std::min(std::fabs(sl.velocity), std::min(std::fabs(sr.velocity), std::fabs(sc.velocity))) * (sl.velocity > 0 ?
				1 : -1);
				
		double dx1 = edges[2] - edges[1];
		double dx0 = edges[1] - edges[0];


		res[1] = new_cells[0] - slope*0.5*dx0;
		res[0] = res[1];
		res[2] = new_cells[0] + slope*0.5*dx0;
	}
	else
	{
		size_t N = edges.size();
		std::array<Primitive, 3> new_cells = { cells[N - 4], cells[N -3], cells[N -2] };
		if (SR)
		{
			for (size_t i = 0; i < 3; ++i)
			{
				double gamma = 1.0 / std::sqrt(1.0 -
					new_cells[i].velocity*new_cells[i].velocity);
				new_cells[i].velocity *= gamma;
			}
		}
		size_t i = edges.size()-3;
		const Primitive sl = (new_cells[1] - new_cells[0]) / (0.5*(edges[i + 1] - edges[i - 1]));
		const Primitive sr = (new_cells[2] - new_cells[1]) / (0.5*(edges[i + 2] - edges[i]));
		const Primitive sc = (new_cells[2] - new_cells[0]) / (0.5*(edges[i + 2] + edges[i + 1] - edges[i - 1]
			- edges[i]));
		if (sl.density*sr.density < 0)
			slope.density = 0;
		else
			slope.density = std::min(std::fabs(sl.density), std::min(std::fabs(sr.density), std::fabs(sc.density))) * (sl.density > 0 ?
				1 : -1);
		if (sl.pressure*sr.pressure < 0)
			slope.pressure = 0;
		else
			slope.pressure = std::min(std::fabs(sl.pressure), std::min(std::fabs(sr.pressure), std::fabs(sc.pressure))) * (sl.pressure > 0 ?
				1 : -1);
		if (sl.velocity*sr.velocity < 0)
			slope.velocity = 0;
		else
			slope.velocity = std::min(std::fabs(sl.velocity), std::min(std::fabs(sr.velocity), std::fabs(sc.velocity))) * (sl.velocity > 0 ?
				1 : -1);
		double dxl = edges[i+1] - edges[i];
		double dxr = edges[i+2] - edges[i+1];
		//slope.density = -2 * (cells[i].density + slope.density*0.5*dxl - cells[i + 1].density) / dxr;
		//slope.pressure = -2 * (cells[i].pressure + slope.pressure*0.5*dxl - cells[i + 1].pressure) / dxr;
		//slope.velocity = -2 * (cells[i].velocity + slope.velocity*0.5*dxl - cells[i + 1].velocity) / dxr;

		res[0] = new_cells[2] - slope*0.5*dxr;
		res[1] = new_cells[2] + slope*0.5*dxr;
		res[2] = res[1];
	}
	if (SR)
	{
		for (size_t i = 0; i < 3; ++i)
		{
			double factor = 1.0 / std::sqrt(1 +
				res[i].velocity*res[i].velocity);
			res[i].velocity *= factor;
		}
	}
	return res;
}


vector<Primitive> Periodic::GetBoundaryValues(vector<Primitive> const & cells, vector<double> const & edges, size_t index,
	double time, bool SR) const
{
	vector<Primitive> res(3);
	size_t N = edges.size();
	double L = edges[N - 1] - edges[0];
	Primitive slope0,slopeN;
	std::array<Primitive, 4> new_cells = { cells[0], cells[1], cells[N - 3], cells[N - 2] };
	if (SR)
	{
		for (size_t i = 0; i < 4; ++i)
		{
			double gamma = 1.0 / std::sqrt(1.0 -
				new_cells[i].velocity*new_cells[i].velocity);
			new_cells[i].velocity *= gamma;
		}
	}
	Primitive sr = (new_cells[1] - new_cells[0]) / (0.5*(edges[2] - edges[0]));
	Primitive sl = (new_cells[0] - new_cells[3]) / (0.5*(edges[1] - edges[N - 2] + L));
	Primitive sc = (new_cells[1] - new_cells[3]) / (0.5*(edges[2] + edges[1] - edges[0] - edges[N - 2] + L));
	if (sl.density*sr.density < 0)
		slope0.density = 0;
	else
		slope0.density = min(fabs(sl.density), min(fabs(sr.density), fabs(sc.density))) * (sl.density > 0 ?
			1 : -1);
	if (sl.pressure*sr.pressure < 0)
		slope0.pressure = 0;
	else
		slope0.pressure = min(fabs(sl.pressure), min(fabs(sr.pressure), fabs(sc.pressure))) * (sl.pressure > 0 ?
			1 : -1);
	if (sl.velocity*sr.velocity < 0)
		slope0.velocity = 0;
	else
		slope0.velocity = min(fabs(sl.velocity), min(fabs(sr.velocity), fabs(sc.velocity))) * (sl.velocity > 0 ?
			1 : -1);
	sr = sl;
	sl = (new_cells[3] - new_cells[2]) / (0.5*(edges[N - 1] - edges[N - 3]));
	sc = (new_cells[0] - new_cells[2]) / (0.5*(edges[N - 1] + edges[0] - edges[N - 3] - edges[N - 2] + L));
	if (sl.density*sr.density < 0)
		slopeN.density = 0;
	else
		slopeN.density = min(fabs(sl.density), min(fabs(sr.density), fabs(sc.density))) * (sl.density > 0 ?
			1 : -1);
	if (sl.pressure*sr.pressure < 0)
		slopeN.pressure = 0;
	else
		slopeN.pressure = min(fabs(sl.pressure), min(fabs(sr.pressure), fabs(sc.pressure))) * (sl.pressure > 0 ?
			1 : -1);
	if (sl.velocity*sr.velocity < 0)
		slopeN.velocity = 0;
	else
		slopeN.velocity = min(fabs(sl.velocity), min(fabs(sr.velocity), fabs(sc.velocity))) * (sl.velocity > 0 ?
			1 : -1);

	if (index == 0)
	{
		res[1] = new_cells[0] - slope0*(0.5*(edges[1] - edges[0]));
		res[0] = new_cells[3]+slopeN*(0.5*(edges[N-1] - edges[N-2]));
		res[2] = new_cells[0] + slope0*(0.5*(edges[1] - edges[0]));
	}
	else
	{
		res[1] = new_cells[3] + slopeN*(0.5*(edges[N - 1] - edges[N - 2]));
		res[2] = new_cells[0] - slope0*(0.5*(edges[1] - edges[0]));
		res[0] = new_cells[3] - slopeN*(0.5*(edges[N - 1] - edges[N - 2]));
	}
	if (SR)
	{
		for (size_t i = 0; i < 3; ++i)
		{
			double factor = 1.0 / std::sqrt(1 +
				res[i].velocity*res[i].velocity);
			res[i].velocity *= factor;
		}
	}
	return res;
}

SeveralBoundary::SeveralBoundary(Boundary const & left, Boundary const & right):left_(left),right_(right)
{}

vector<Primitive> SeveralBoundary::GetBoundaryValues(vector<Primitive> const & cells, vector<double> const & edges, size_t index, 
	double time,bool SR) const
{
	if(index==0)
		return left_.GetBoundaryValues(cells, edges, index,time,SR);
	else
		return right_.GetBoundaryValues(cells, edges, index,time,SR);
}

ConstantPrimitive::ConstantPrimitive(Primitive outer):outer_(outer)
{}

vector<Primitive> ConstantPrimitive::GetBoundaryValues(vector<Primitive> const & cells, vector<double> const & edges, size_t index, 
	double time,bool SR) const
{
	vector<Primitive> res(3);
	size_t N = edges.size();
	Primitive left, center, right;
	std::array<Primitive, 5> new_cells = { cells[0], cells[1], cells[N - 3], 
		cells[N-2], outer_ };
	if (SR)
	{
		for (size_t i = 0; i < 5; ++i)
		{
			double gamma = 1.0 / std::sqrt(1.0 -
				new_cells[i].velocity*new_cells[i].velocity);
			new_cells[i].velocity *= gamma;
		}
	}
	if (index == 0)
	{
		left = new_cells[4];
		center = new_cells[0];
		right = new_cells[1];
	}
	else
	{
		left = new_cells[2];
		center = new_cells[3];
		right = new_cells[4];
	}
	double dxl, dxr, dxc;
	if (index == 0)
	{
		dxl = edges[1] - edges[0];
		dxr = 0.5*(edges[2] - edges[0]);
		dxc = 0.5*edges[2] + edges[1] - 1.5*edges[0];
	}
	else
	{
		dxl = 0.5*(edges[N-1] - edges[N-3]);
		dxr = edges[N-1] - edges[N-2];
		dxc = 1.5*edges[N - 1] - edges[N - 2] - 0.5*edges[N - 3];
	}
	Primitive sr = (right - center) / dxr;
	Primitive sl = (center - left) / dxl;
	Primitive sc = (right - left) / dxc;
	Primitive slope0;
	if (sl.density*sr.density < 0)
		slope0.density = 0;
	else
		slope0.density = min(fabs(sl.density), min(fabs(sr.density), fabs(sc.density))) * (sl.density > 0 ?
			1 : -1);
	if (sl.pressure*sr.pressure < 0)
		slope0.pressure = 0;
	else
		slope0.pressure = min(fabs(sl.pressure), min(fabs(sr.pressure), fabs(sc.pressure))) * (sl.pressure > 0 ?
			1 : -1);
	if (sl.velocity*sr.velocity < 0)
		slope0.velocity = 0;
	else
		slope0.velocity = min(fabs(sl.velocity), min(fabs(sr.velocity), fabs(sc.velocity))) * (sl.velocity > 0 ?
			1 : -1);
	if (index == 0)
	{
		res[1] = new_cells[0] - slope0*(0.5*(edges[1] - edges[0]));
		res[0] = new_cells[4];
		res[2] = new_cells[0] + slope0*(0.5*(edges[1] - edges[0]));
	}
	else
	{
		res[1] = new_cells[3] + slope0*(0.5*(edges[N - 1] - edges[N - 2]));
		res[2] = new_cells[4];
		res[0] = new_cells[3] - slope0*(0.5*(edges[N - 1] - edges[N - 2]));
	}
	if (SR)
	{
		for (size_t i = 0; i < 3; ++i)
		{
			double factor = 1.0 / std::sqrt(1 +
				res[i].velocity*res[i].velocity);
			res[i].velocity *= factor;
		}
	}
	return res;
}

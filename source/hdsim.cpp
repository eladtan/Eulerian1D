#include "hdsim.hpp"
#include <algorithm>
#include <iostream>
#include "RiemannSolver.hpp"
#include "universal_error.hpp"
#ifdef RICH_MPI
#include "mpi_comm.hpp"
#endif
#include <cassert>
namespace
{
	Extensive CalcExtensive(Geometry const& geo, std::vector<double> const&
		edges, Primitive const& cell, size_t index)
	{
		Extensive res;
		double vol = geo.GetVolume(edges, index);
		res.mass = vol * cell.density;
		res.momentum = res.mass*cell.velocity;
		res.et = res.mass*cell.energy;
		res.entropy = res.mass*cell.entropy;
		res.energy = res.et + 0.5*res.momentum*res.momentum / res.mass;
		return res;
	}

	Primitive CalcPrimitive(double vol, Extensive const& extensive,
		IdealGas const& eos)
	{
		Primitive cell;
		cell.density = extensive.mass / vol;
		cell.velocity = extensive.momentum / extensive.mass;
		cell.energy = extensive.et / extensive.mass;
		cell.pressure = eos.de2p(cell.density, cell.energy);
		cell.entropy = eos.dp2s(cell.density, cell.pressure);
		return cell;
	}

	Primitive GetSlope(Primitive const& left, Primitive const& center, Primitive const& right,
		double e0, double e1, double e2, double e3)
	{
		Primitive slope;
		const Primitive sl = (center - left) / (0.5*(e2 - e0));
		const Primitive sr = (right - center) / (0.5*(e3 - e1));
		const Primitive sc = (right - left) / (0.5*(e3 + e2 - e0 - e1));
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
		return slope;
	}

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

	std::vector<double> GetVGrid(std::vector<std::pair<Primitive, Primitive> > const& interp_values,double time,std::vector<double> const& edges)
	{
		size_t N = edges.size();
		std::vector<double> res(N, 0);
		double alpha = 5e-5;
#ifdef RICH_MPI
		double vedge= interp_values[0].first.velocity;
		MPI_Bcast(&vedge, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		std::array<double, NGHOSTCELLS * 2> ghost_edges = SendRecvEdges(edges);
		if ((edges[0] - ghost_edges[NGHOSTCELLS - 1]) > alpha*edges[0])
			res[0] = vedge * 1.05;
		else
			res[0] = vedge;
#else
		res[0] = vedge;
#endif
		for (size_t i = 1; i < N; ++i)
		{
			if ((edges[i] - edges[i-1]) > alpha*edges[i])
				res[i] = vedge * 1.05;
			else
				res[i] = vedge;
		}
		return res;
	}
}

hdsim::hdsim(double cfl, vector<Primitive> const& cells, vector<double> const& edges, Interpolation const& interp,
	IdealGas const& eos, RiemannSolver const& rs, SourceTerm const& source, Geometry const& geo,
	const double AMR_ratio) :cfl_(cfl),
	cells_(cells), edges_(edges), interpolation_(interp), eos_(eos), rs_(rs), time_(0), cycle_(0), TotalEcool_(0),
	extensives_(vector<Extensive>()), source_(source), geo_(geo), AMR_ratio_(AMR_ratio), dt_suggest_(-1)
{
#ifdef RICH_MPI
	size_t Ntotal = cells.size();
	int rank = 0, ws = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &ws);
	size_t index_lower = rank * Ntotal / ws;
	size_t index_upper = (rank + 1)*Ntotal / ws;
	std::vector<Primitive> cells_temp(cells.begin() + index_lower,
		cells.begin() + index_upper);
	std::vector<Primitive> cells_temp2(cells.begin(), cells.begin() + 1);
	cells_ = cells_temp;
	std::vector<double> edge_temp(edges.begin() + index_lower, edges.begin() +
		index_upper + 1);
	edges_ = edge_temp;
#endif
	size_t N = cells_.size();
	extensives_.resize(N);
	for (size_t i = 0; i < N; ++i)
	{
		double vol = geo_.GetVolume(edges_, i);;
		extensives_[i].mass = cells_[i].density*vol;
		extensives_[i].momentum = extensives_[i].mass*cells_[i].velocity;
		extensives_[i].energy = 0.5*extensives_[i].momentum*extensives_[i].momentum / extensives_[i].mass +
			eos_.dp2e(cells_[i].density, cells_[i].pressure)*extensives_[i].mass;
		extensives_[i].entropy = eos_.dp2s(cells_[i].density, cells_[i].pressure)*extensives_[i].mass;
		cells_[i].entropy = eos_.dp2s(cells_[i].density, cells_[i].pressure);
		cells_[i].energy = eos_.dp2e(cells_[i].density, cells_[i].pressure);
		cells_[i].LastCool = time_;
		extensives_[i].et = cells_[i].energy*extensives_[i].mass;
	}
}


hdsim::~hdsim()
{}

namespace
{
	void MoveGrid(std::vector<double> const& vgrid, std::vector<double> &edges, double dt)
	{
		size_t N = edges.size();
		for (size_t i = 0; i < N; ++i)
			edges[i] += vgrid[i] * dt;
	}

	double GetTimeStep(vector<Primitive> const& cells, vector<double> const& edges, IdealGas const& eos, double cfl,
		SourceTerm const&source, double &dt_suggest, std::vector<double> const & vgrid)
	{
		double force_inverse_dt = source.GetInverseTimeStep(edges);
		double dt_1 = 0;
		size_t N = cells.size();
		for (size_t i = 0; i < N; ++i)
		{
			dt_1 = std::max(dt_1, 2 * (std::max(std::abs(vgrid[i]-vgrid[i+1]),
				std::max(std::abs(cells[i].velocity - vgrid[i]), std::abs(cells[i].velocity - vgrid[i+1]))) 
				+ eos.dp2c(cells[i].density, cells[i].pressure)) / (edges[i + 1] - edges[i]));
		}
		dt_1 = max(dt_1, force_inverse_dt);
		if (dt_suggest > 0)
		{
			dt_1 = std::max(dt_1, 1.0 / dt_suggest);
			dt_suggest = -1;
		}
#ifdef RICH_MPI
		MPI_Allreduce(MPI_IN_PLACE, &dt_1, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
		return cfl / dt_1;
	}

	void GetFluxes(vector<pair<Primitive, Primitive> > const& interp_values, RiemannSolver const&rs, std::vector<Extensive> &res, IdealGas const& eos, 
		std::vector<double> const & vgrid)
	{
		size_t N = interp_values.size();
		res.resize(N);
		for (int i = 0; i < N; ++i)
			res[i] = rs.SolveRS(interp_values[i].first, interp_values[i].second, eos, vgrid[i]);
	}

	void UpdateExtensives(vector<Extensive> &cells, std::vector<Extensive> &fluxes, double dt, Geometry const& geo,
		std::vector<double> const& edges ,vector<pair<Primitive, Primitive> > const &interp_values, std::vector<double> const& vgrid)
	{
		size_t N = cells.size();
		for (size_t i = 0; i < N; ++i)
		{
			double v = cells[i].momentum / cells[i].mass;

			//double oldEk = cells[i].momentum*cells[i].momentum / cells[i].mass;
			double newP = cells[i].momentum + fluxes[i].momentum * geo.GetArea(edges[i])*dt;
			double newM = cells[i].mass + fluxes[i].mass* geo.GetArea(edges[i])*dt;
			//double newEk0 = newP * newP*0.5 / newM;
			newP = cells[i].momentum - fluxes[i + 1].momentum * geo.GetArea(edges[i + 1])*dt;
			newM = cells[i].mass - fluxes[i + 1].mass* geo.GetArea(edges[i + 1])*dt;
			//double newEk1 = newP * newP*0.5 / newM;

			double dP = -(fluxes[i + 1].momentum * geo.GetArea(edges[i + 1]) - fluxes[i].momentum * geo.GetArea(edges[i]))*dt;
			cells[i].momentum += dP;
			double dE = -(fluxes[i + 1].energy* geo.GetArea(edges[i + 1]) - fluxes[i].energy* geo.GetArea(edges[i]))*dt;
			cells[i].energy += dE;
			double mdot = -(fluxes[i + 1].mass* geo.GetArea(edges[i + 1]) - fluxes[i].mass* geo.GetArea(edges[i]))*dt;
			cells[i].mass += mdot;
			//	double newEk = cells[i].momentum*cells[i].momentum / cells[i].mass;
			cells[i].et += dE - v * dP + 0.5*v*v*mdot;
			//cells[i].et += dE - newEk1 - newEk0 + oldEk;
			cells[i].entropy += -(fluxes[i + 1].entropy* geo.GetArea(edges[i + 1]) - fluxes[i].entropy* geo.GetArea(edges[i]))*dt;
			if (!(cells[i].mass > 0) || !(cells[i].entropy > 0))
			{
				UniversalError eo("Bad extensive update");
				eo.AddEntry("Cell index", i);
#ifdef RICH_MPI
				int rank = 0;
				MPI_Comm_rank(MPI_COMM_WORLD, &rank);
				eo.AddEntry("rank", rank);
#endif
				eo.AddEntry("mass", cells[i].mass);
				eo.AddEntry("Entropy", cells[i].entropy);
				eo.AddEntry("dt", dt);
				eo.AddEntry("Left edge", edges[i]);
				eo.AddEntry("Right edge", edges[i+1]);
				eo.AddEntry("Left density", interp_values[i].first.density);
				eo.AddEntry("Left pressure", interp_values[i].first.pressure);
				eo.AddEntry("Left velocity", interp_values[i].first.velocity);
				eo.AddEntry("Left entropy", interp_values[i].first.entropy);
				eo.AddEntry("Right density", interp_values[i].second.density);
				eo.AddEntry("Right pressure", interp_values[i].second.pressure);
				eo.AddEntry("Right velocity", interp_values[i].second.velocity);
				eo.AddEntry("Right entropy", interp_values[i].second.entropy);
				eo.AddEntry("vgrid", vgrid[i]);
				eo.AddEntry("Left density", interp_values[i + 1].first.density);
				eo.AddEntry("Left pressure", interp_values[i + 1].first.pressure);
				eo.AddEntry("Left velocity", interp_values[i + 1].first.velocity);
				eo.AddEntry("Left entropy", interp_values[i + 1].first.entropy);
				eo.AddEntry("Right density", interp_values[i + 1].second.density);
				eo.AddEntry("Right pressure", interp_values[i + 1].second.pressure);
				eo.AddEntry("Right velocity", interp_values[i + 1].second.velocity);
				eo.AddEntry("Right entropy", interp_values[i + 1].second.entropy);
				eo.AddEntry("vgrid", vgrid[i + 1]);
				throw eo;
			}
		}
	}

	bool ShouldUseEntropy(double dv, double et)
	{
		double ek = dv * dv;
		if (et < 0)
			return true;
		if (et > 0.001*ek)
			return false;
		return true;
	}

	void UpdateCells(vector<Extensive> &extensive, vector<double> const& edges, IdealGas const& eos,
		vector<Primitive> &cells, Geometry const& geo, vector<pair<Primitive, Primitive> > const &interp_values,std::vector<double> const& vgrid)
	{
		size_t N = cells.size();
		for (int i = 0; i < N; ++i)
		{
			double vol = geo.GetVolume(edges, i);
			cells[i].density = extensive[i].mass / vol;
			cells[i].velocity = extensive[i].momentum / extensive[i].mass;
			cells[i].energy = extensive[i].et / extensive[i].mass;
			double et = cells[i].energy;
			bool entropy_use = false;
			if (i == 0)
				entropy_use = ShouldUseEntropy(std::abs(cells[i].velocity - cells[i + 1].velocity), et);
			else
				if (i == (N - 1))
					entropy_use = ShouldUseEntropy(std::abs(cells[i].velocity - cells[i - 1].velocity), et);
				else
					entropy_use = ShouldUseEntropy(std::max(std::abs(cells[i].velocity - cells[i - 1].velocity), std::abs(cells[i].velocity - cells[i + 1].velocity)), et);
			if (entropy_use)
			{
				cells[i].entropy = extensive[i].entropy / extensive[i].mass;
				cells[i].pressure = eos.sd2p(cells[i].entropy , cells[i].density);
				cells[i].energy = eos.dp2e(cells[i].density, cells[i].pressure);
				et = extensive[i].mass*cells[i].energy;
				extensive[i].energy = 0.5*extensive[i].momentum*extensive[i].momentum / extensive[i].mass + et;
				extensive[i].et = et;
			}
			else
			{
				double et2 = (extensive[i].energy - 0.5*extensive[i].momentum*extensive[i].momentum / extensive[i].mass)/extensive[i].mass;
				if (et2 > 0.95*et && et2 < 1.05*et)
				{
					et = et2;
					cells[i].energy = et;
					extensive[i].et = et * extensive[i].mass;
				}			
				cells[i].pressure = eos.de2p(cells[i].density, et);
				cells[i].entropy = eos.dp2s(cells[i].density, cells[i].pressure);
				extensive[i].entropy = cells[i].entropy*extensive[i].mass;
			}
			if (!(cells[i].pressure > 0) || !(cells[i].entropy > 0))
			{
				UniversalError eo("Bad cell update");
				eo.AddEntry("Cell index", i);
#ifdef RICH_MPI
				int rank = 0;
				MPI_Comm_rank(MPI_COMM_WORLD, &rank);
				eo.AddEntry("rank", rank);
#endif
				eo.AddEntry("Left density", interp_values[i].first.density);
				eo.AddEntry("Left pressure", interp_values[i].first.pressure);
				eo.AddEntry("Left velocity", interp_values[i].first.velocity);
				eo.AddEntry("Left entropy", interp_values[i].first.entropy);
				eo.AddEntry("Right density", interp_values[i].second.density);
				eo.AddEntry("Right pressure", interp_values[i].second.pressure);
				eo.AddEntry("Right velocity", interp_values[i].second.velocity);
				eo.AddEntry("Right entropy", interp_values[i].second.entropy);
				eo.AddEntry("vgrid", vgrid[i]);
				eo.AddEntry("Left density", interp_values[i+1].first.density);
				eo.AddEntry("Left pressure", interp_values[i + 1].first.pressure);
				eo.AddEntry("Left velocity", interp_values[i + 1].first.velocity);
				eo.AddEntry("Left entropy", interp_values[i + 1].first.entropy);
				eo.AddEntry("Right density", interp_values[i + 1].second.density);
				eo.AddEntry("Right pressure", interp_values[i + 1].second.pressure);
				eo.AddEntry("Right velocity", interp_values[i + 1].second.velocity);
				eo.AddEntry("Right entropy", interp_values[i + 1].second.entropy);
				eo.AddEntry("vgrid", vgrid[i + 1]);
				throw eo;
			}
			
		}
	}
}

void hdsim::ReCalcCells(vector<Extensive> &extensive)
{
	size_t N = cells_.size();
	for (size_t i = 0; i < N; ++i)
	{
		double vol = geo_.GetVolume(edges_, i);
		cells_[i].density = extensive[i].mass / vol;
		cells_[i].velocity = extensive[i].momentum / extensive[i].mass;
		const double et = extensive[i].et / extensive[i].mass;
		cells_[i].pressure = eos_.de2p(cells_[i].density, et);
		cells_[i].energy = et;
		cells_[i].entropy = eos_.dp2s(cells_[i].density, cells_[i].pressure);
		extensive[i].entropy = cells_[i].entropy*extensive[i].mass;
	}
}

void hdsim::ReCalcExtensives(vector<Primitive> const& cells)
{
	size_t N = cells.size();
	for (size_t i = 0; i < N; ++i)
	{
		double vol = geo_.GetVolume(edges_, i);
		extensives_[i].mass = vol * cells[i].density;
		extensives_[i].momentum = extensives_[i].mass*cells[i].velocity;
		extensives_[i].et = extensives_[i].mass*cells[i].energy;
		extensives_[i].entropy = extensives_[i].mass*cells[i].entropy;
		extensives_[i].energy = extensives_[i].et + 0.5*extensives_[i].momentum*extensives_[i].momentum / extensives_[i].mass;
	}
}

void hdsim::SuggestTimeStep(double dt)
{
	dt_suggest_ = dt;
}

void hdsim::TimeAdvance()
{
#ifdef RICH_MPI
	std::array<Primitive, NGHOSTCELLS * 2> ghost_cells = SendRecvPrimitive(cells_);
	std::array<double, NGHOSTCELLS * 2> ghost_edges = SendRecvEdges(edges_);
#endif
	interpolation_.GetInterpolatedValues(cells_, edges_, interp_values_, time_
#ifdef RICH_MPI
		, ghost_cells, ghost_edges
#endif
	);
	std::vector<double> vgrid(edges_.size(), 0);
	GetFluxes(interp_values_, rs_, fluxes_, eos_,vgrid);
	double dt = GetTimeStep(cells_, edges_, eos_, cfl_, source_, dt_suggest_,vgrid);


	/*if (BoundarySolution_ != 0)
	{
		pair<RSsolution,RSsolution> bvalues = BoundarySolution_->GetBoundaryValues(cells_);
		if (BoundarySolution_->ShouldCalc().first)
			rs_values_[0] = bvalues.first;
		if (BoundarySolution_->ShouldCalc().second)
			rs_values_.back() = bvalues.second;
	}*/

	UpdateExtensives(extensives_, fluxes_, dt, geo_, edges_,interp_values_,vgrid);
	source_.CalcForce(edges_, cells_, time_, extensives_, dt);
	UpdateCells(extensives_, edges_, eos_, cells_, geo_, interp_values_, vgrid);
	time_ += dt;
	++cycle_;
	AMR(
#ifdef RICH_MPI
		ghost_cells , ghost_edges
#endif
	);
}

void hdsim::TimeAdvance2()
{
#ifdef RICH_MPI
	std::array<Primitive, NGHOSTCELLS * 2> ghost_cells = SendRecvPrimitive(cells_);
	std::array<double, NGHOSTCELLS * 2> ghost_edges = SendRecvEdges(edges_);
#endif
	interpolation_.GetInterpolatedValues(cells_, edges_, interp_values_, time_
#ifdef RICH_MPI
		, ghost_cells, ghost_edges
#endif
	);
	std::vector<double> vgrid = GetVGrid(interp_values_,time_,edges_);

	GetFluxes(interp_values_, rs_, fluxes_, eos_, vgrid);

	/*if (BoundarySolution_ != 0)
	{
		pair<RSsolution, RSsolution> bvalues = BoundarySolution_->GetBoundaryValues(cells_);
		if (BoundarySolution_->ShouldCalc().first)
			rs_values_[0] = bvalues.first;
		if (BoundarySolution_->ShouldCalc().second)
			rs_values_.back() = bvalues.second;
	}*/
	double dt = GetTimeStep(cells_, edges_, eos_, cfl_, source_, dt_suggest_, vgrid);
	if (cycle_ == 0)
		dt *= 0.005;

	vector<Extensive> old_extensive(extensives_);
	std::vector<double> oldedges(edges_);
	UpdateExtensives(extensives_, fluxes_, dt, geo_, edges_, interp_values_, vgrid);
	MoveGrid(vgrid, edges_, dt);
	source_.CalcForce(edges_, cells_, time_, extensives_, dt);
	UpdateCells(extensives_, edges_, eos_, cells_, geo_, interp_values_, vgrid);
	time_ += dt;
#ifdef RICH_MPI
	ghost_cells = SendRecvPrimitive(cells_);
	ghost_edges = SendRecvEdges(edges_);
#endif
	interpolation_.GetInterpolatedValues(cells_, edges_, interp_values_, time_
#ifdef RICH_MPI
		, ghost_cells, ghost_edges
#endif
	);
	GetFluxes(interp_values_, rs_, fluxes_, eos_, vgrid);
	/*if (BoundarySolution_ != 0)
	{
		pair<RSsolution, RSsolution> bvalues = BoundarySolution_->GetBoundaryValues(cells_);
		if (BoundarySolution_->ShouldCalc().first)
			rs_values_[0] = bvalues.first;
		if (BoundarySolution_->ShouldCalc().second)
		{
			rs_values_.back() = bvalues.second;
		}
	}
	*/
	UpdateExtensives(extensives_, fluxes_, dt, geo_, edges_, interp_values_, vgrid);
	source_.CalcForce(edges_, cells_, time_, extensives_, dt);
	size_t N = extensives_.size();
	for (size_t i = 0; i < N; ++i)
	{
		extensives_[i] += old_extensive[i];
		extensives_[i] *= 0.5;
	}
#ifdef RICH_MPI
	RedistributeExtensives(extensives_, edges_, cells_);
#endif
	UpdateCells(extensives_, edges_, eos_, cells_, geo_, interp_values_, vgrid);
	++cycle_;
	AMR(
#ifdef RICH_MPI
		ghost_cells, ghost_edges
#endif
	);
}

void hdsim::AMR(
#ifdef RICH_MPI
	std::array<Primitive, NGHOSTCELLS * 2> ghost_cells,
	std::array<double, NGHOSTCELLS * 2> ghost_edges
#endif
)
{
	if (!(AMR_ratio_ > 0))
		return;
	bool planar = (geo_.GetArea(1) < (1 + 1e-7) && geo_.GetArea(1) > (1 - 1e-7));
	size_t Nlevels = 2;
#ifdef RICH_MPI
	Nlevels = NGHOSTCELLS;
#endif
	double min_size = std::pow(2.0, -0.33*static_cast<double>(Nlevels));
	double pratio = 1.01;
	double dratio = 1.01;
	size_t N = cells_.size();
	// remove smooth small cells
	std::vector<size_t> edge_remove;
	for (size_t i = 0; i < N; ++i)
	{
		// was left cell removed?
		if (edge_remove.size() > 0 && (edge_remove.back() == i || edge_remove.back() == i - 1))
			continue;
		double dx = edges_[i + 1] - edges_[i];
		if (i < Nlevels || i >= (N - Nlevels))
		{
			if (dx < (AMR_ratio_* min_size*0.2*(planar ? 1.0 : edges_[i])))
			{
				if (i == 0)
					edge_remove.push_back(1);
				else
					if (i == (N - 1))
						edge_remove.push_back(N - 1);
					else
					{
						if ((edges_[i + 2] - edges_[i + 1]) > (edges_[i] - edges_[i - 1]))
							edge_remove.push_back(i);
						else
							edge_remove.push_back(i + 1);
					}
			}
			continue;
		}
		// do we have too small neighbors
		if ((edges_[i] - edges_[i - 1]) < 0.5*dx || (edges_[i + 2] - edges_[i + 1]) < 0.5*dx)
			continue;
		// Are we too small?
		if (dx < (AMR_ratio_* min_size*(planar ? 1.0 : edges_[i])))
		{
			double pratio2 = 1.1*std::pow((AMR_ratio_ * min_size*(planar ? 1.0 : edges_[i])) / dx, 0.9);;
			double Tratio2 = 1.2*std::pow((AMR_ratio_ * min_size*(planar ? 1.0 : edges_[i])) / dx, 0.1);;
			double dratio2 = 1.2*std::pow((AMR_ratio_ * min_size*(planar ? 1.0 : edges_[i])) / dx, 1.1);
			bool smooth_left = (cells_[i].density < cells_[i - 1].density * dratio2) && (cells_[i].density * dratio2 > cells_[i - 1].density)
				&& (cells_[i].pressure < cells_[i - 1].pressure*pratio2) && (cells_[i].pressure*pratio > cells_[i - 1].pressure);
			bool smooth_right = (cells_[i].density < cells_[i + 1].density * dratio2) && (cells_[i].density * dratio2 > cells_[i + 1].density)
				&& (cells_[i].pressure < cells_[i + 1].pressure*pratio2) && (cells_[i].pressure*pratio2 > cells_[i + 1].pressure);
			bool Tratio_l = (cells_[i].density*cells_[i - 1].pressure) < (Tratio2*cells_[i - 1].density*cells_[i].pressure)
				&& (Tratio2*cells_[i].density*cells_[i - 1].pressure) > (cells_[i - 1].density*cells_[i].pressure);
			bool Tratio_r = (cells_[i].density*cells_[i + 1].pressure) < (Tratio2*cells_[i + 1].density*cells_[i].pressure)
				&& (Tratio2*cells_[i].density*cells_[i + 1].pressure) > (cells_[i + 1].density*cells_[i].pressure);
			if ((smooth_left && smooth_right && Tratio_l && Tratio_r) || (dx < AMR_ratio_ * min_size*0.2*(planar ? 1.0 : edges_[i])))
			{
				if ((edges_[i + 2] - edges_[i + 1]) > (edges_[i] - edges_[i - 1]))
					edge_remove.push_back(i);
				else
					edge_remove.push_back(i + 1);
			}
			else
			{
				if (smooth_right && Tratio_r && ((edges_[i + 2] - edges_[i + 1]) < 1.5*AMR_ratio_*(planar ? 1.0 : edges_[i + 1])))
					edge_remove.push_back(i + 1);
				if (smooth_left && Tratio_l && ((edges_[i] - edges_[i - 1]) < 1.5*AMR_ratio_*(planar ? 1.0 : edges_[i - 1])))
					edge_remove.push_back(i);
			}
		}
		else
		{
			if (dx < (AMR_ratio_*0.6*(planar ? 1.0 : edges_[i])))
			{
				double Tratio2 = 1.1;
				bool smooth = true;
				for (size_t j = 0; j < (2 * Nlevels - 1); ++j)
				{
					if (!((cells_[i - Nlevels + j + 1].density < cells_[i - Nlevels + j].density * dratio) && (cells_[i - Nlevels + j + 1].density * dratio > cells_[i - Nlevels + j].density)
						&& (cells_[i - Nlevels + j + 1].pressure < cells_[i - Nlevels + j].pressure*pratio) && (cells_[i - Nlevels + j + 1].pressure*pratio > cells_[i - Nlevels + j].pressure)))
					{
						smooth = false;
						break;
					}
				}
				if (smooth)
				{
					if ((edges_[i + 2] - edges_[i + 1]) > (edges_[i] - edges_[i - 1]))
						edge_remove.push_back(i);
					else
						edge_remove.push_back(i + 1);
				}
			}
		}
	}
	edge_remove = unique(edge_remove);
	size_t Nremove = edge_remove.size();
	for (size_t i = 0; i < Nremove; ++i)
		extensives_[edge_remove[i] - 1] += extensives_[edge_remove[i]];
	// Remove old cells
	if (Nremove > 0)
	{
		RemoveVector(extensives_, edge_remove);
		RemoveVector(cells_, edge_remove);
		RemoveVector(edges_, edge_remove);
		ReCalcCells(extensives_);
	}
	// split cells ahead of non-smooth
	pratio = 1.05;
	dratio = 1.05;
	N = cells_.size();
	std::vector<size_t> edge_split;
	size_t Nstart = Nlevels;
	size_t Nend = N - Nlevels;
	size_t shift = 0;
	std::vector<Primitive> temp_cells;
	std::vector<double> temp_edges;
#ifdef RICH_MPI
	int rank = 0, ws = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &ws);
	temp_cells.resize(N + NGHOSTCELLS * 2);
	std::copy(cells_.begin(), cells_.end(), temp_cells.begin() + NGHOSTCELLS);
	std::copy(ghost_cells.begin(), ghost_cells.begin() + NGHOSTCELLS, temp_cells.begin());
	std::copy(ghost_cells.begin() + NGHOSTCELLS, ghost_cells.end(), temp_cells.begin() + NGHOSTCELLS + N);
	temp_edges.resize(N + NGHOSTCELLS * 2 + 1);
	std::copy(edges_.begin(), edges_.end(), temp_edges.begin() + NGHOSTCELLS);
	std::copy(ghost_edges.begin(), ghost_edges.begin() + NGHOSTCELLS, temp_edges.begin());
	std::copy(ghost_edges.begin() + NGHOSTCELLS, ghost_edges.end(), temp_edges.begin() + NGHOSTCELLS + N + 1);
	if (rank < (ws - 1))
		Nend = N;
	if (rank > 0)
		Nstart = 0;
	shift = NGHOSTCELLS;
#else
	temp_cells = cells_;
	temp_edges = edges_;
#endif
	for (size_t i = Nstart; i < Nend; ++i)
	{
		double dx = temp_edges[i + 1 + shift] - temp_edges[i + shift];
		if (dx > (AMR_ratio_*min_size * 2 * (planar ? 1.0 : temp_edges[i + shift])))
		{
			double Tratio2 = 1.2;
			// Are we smooth?
			bool smooth = false;
			for (size_t j = 0; j < (Nlevels - 1); ++j)
			{
				bool local_smooth = (temp_cells[i - Nlevels + 1 + j + shift].density < temp_cells[i - Nlevels + j + shift].density
					* dratio) && (temp_cells[i - Nlevels + 1 + j + shift].density * dratio >
						temp_cells[i - Nlevels + j + shift].density)
					&& (temp_cells[i - Nlevels + 1 + j + shift].pressure < temp_cells[i - Nlevels + j + shift].pressure
						*pratio) && (temp_cells[i - Nlevels + 1 + j + shift].pressure*pratio >
							temp_cells[i - Nlevels + j + shift].pressure);
				double l_reduce = std::pow(2.0, -0.33*static_cast<double>(j));
				if (!local_smooth && dx > (AMR_ratio_* 1.4*l_reduce*(planar ? 1.0 : edges_[i])))
				{
					smooth = true;
					break;
				}
			}
			if (!smooth)
			{
				for (size_t j = 1; j < (Nlevels - 1); ++j)
				{
					bool local_smooth = (temp_cells[i + 1 + j + shift].density < temp_cells[i + j + shift].density
						* dratio) && (temp_cells[i + 1 + j + shift].density * dratio >
							temp_cells[i + j + shift].density)
						&& (temp_cells[i + 1 + j + shift].pressure < temp_cells[i + j + shift].pressure
							*pratio) && (temp_cells[i + 1 + j + shift].pressure*pratio >
								temp_cells[i + j + shift].pressure);
					double l_reduce = std::pow(2.0, -0.33*static_cast<double>(Nlevels-j-1));
					if (!local_smooth && dx > (AMR_ratio_*1.4*l_reduce*(planar ? 1.0 : edges_[i])))
					{
						smooth = true;
						break;
					}
				}
			}
			if (!smooth)
				continue;
			bool smooth_left = (temp_cells[i + shift].density < temp_cells[i - 1 + shift].density
				* dratio) && (temp_cells[i + shift].density * dratio >
					temp_cells[i - 1 + shift].density)
				&& (temp_cells[i + shift].pressure < temp_cells[i - 1 + shift].pressure
					*pratio) && (temp_cells[i + shift].pressure*pratio >
						temp_cells[i - 1 + shift].pressure);
			bool smooth_right = (temp_cells[i + shift].density <
				temp_cells[i + 1 + shift].density * dratio) && (temp_cells[i + shift].density
					* dratio > temp_cells[i + 1 + shift].density)
				&& (temp_cells[i + shift].pressure < temp_cells[i + 1 + shift].pressure
					*pratio) && (temp_cells[i + shift].pressure*pratio >
						temp_cells[i + 1 + shift].pressure);
			if (smooth_left && smooth_right)
			{
				edge_split.push_back(i);
			}
		}
	}
	if (edge_split.size() > 0)
	{
		size_t Nadd = edge_split.size();
		size_t counter = 0;
		std::vector<double> new_edges(edges_.size() + Nadd);
		std::vector<Primitive> new_cells(N + Nadd);
		std::vector<Extensive> new_extensive(N + Nadd);
		new_edges[0] = edges_[0];
		for (size_t i = 0; i < N; ++i)
		{
			new_edges[i + counter + 1] = edges_[i + 1];
			new_cells[i + counter] = cells_[i];
			if (edge_split[counter] == i)
			{
				Primitive slope = GetSlope(temp_cells[i + shift - 1], temp_cells[i + shift],
					temp_cells[i + shift + 1], temp_edges[i + shift - 1], temp_edges[i + shift],
					temp_edges[i + shift + 1], temp_edges[i + shift + 2]);
				new_cells[i + counter] = cells_[i] - (0.25*
					(temp_edges[i + shift + 1] - temp_edges[i + shift]))*slope;
				new_cells[i + counter].energy = eos_.dp2e(new_cells[i + counter].density,
					new_cells[i + counter].pressure);
				new_cells[i + counter].entropy = eos_.dp2s(new_cells[i + counter].density,
					new_cells[i + counter].pressure);
				new_edges[i + counter + 1] = 0.5*(edges_[i] + edges_[i + 1]);
				Extensive old = extensives_[i];
				new_extensive[i + counter] = CalcExtensive(geo_, new_edges, new_cells[i + counter], i + counter);
				++counter;
				new_extensive[i + counter] = old;
				new_extensive[i + counter] -= new_extensive[i + counter - 1];
				new_edges[i + counter + 1] = edges_[i + 1];
				new_cells[i + counter] = CalcPrimitive(geo_.GetVolume(new_edges, i + counter),
					new_extensive[i + counter], eos_);
				new_cells[i + counter].entropy = eos_.dp2s(new_cells[i + counter].density,
					new_cells[i + counter].pressure);
				new_extensive[i + counter].entropy = new_extensive[i + counter].mass*
					new_cells[i + counter].entropy;
			}
			else
				new_extensive[i+counter] = CalcExtensive(geo_, new_edges, new_cells[i + counter], i + counter);
			if (counter == Nadd && i < (N - 1))
			{
				std::copy(edges_.begin() + i + 2, edges_.end(), new_edges.begin() + i + counter + 2);
				std::copy(cells_.begin() + i + 1, cells_.end(), new_cells.begin() + i + counter + 1);
				//std::copy(cells_.begin() + i + 1, cells_.end(), new_extensive.begin() + i + counter + 1);
				break;
			}
		}
		cells_ = new_cells;
		edges_ = new_edges;
		extensives_.resize(cells_.size());
		ReCalcExtensives(cells_);
	}
}


void hdsim::SetCycle(size_t cyc)
{
	cycle_ = cyc;
}

double hdsim::GetTime() const
{
	return time_;
}

double hdsim::GetEcool() const
{
	return TotalEcool_;
}

double& hdsim::GetEcool()
{
	return TotalEcool_;
}

vector<Primitive> const & hdsim::GetCells() const
{
	return cells_;
}

vector<Primitive>& hdsim::GetCells()
{
	return cells_;
}

vector<Extensive> const & hdsim::GetExtensives() const
{
	return extensives_;
}

vector<Extensive> & hdsim::GetExtensives()
{
	return extensives_;
}

vector<double> const & hdsim::GetEdges() const
{
	return edges_;
}


vector<double>& hdsim::GetEdges()
{
	return edges_;
}


size_t hdsim::GetCycle() const
{
	return cycle_;
}

void hdsim::SetTime(double t)
{
	time_ = t;
}

void hdsim::SetEcool(double E)
{
	TotalEcool_ = E;
}

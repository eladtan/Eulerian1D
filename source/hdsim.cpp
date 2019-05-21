#include "hdsim.hpp"
#include <algorithm>
#include <iostream>
#include "RiemannSolver.hpp"
#ifdef RICH_MPI
#include "mpi_comm.hpp"
#endif
#include <cassert>
namespace
{
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

	double GetVGrid(std::vector<std::pair<Primitive, Primitive> > const& interp_values)
	{
		return interp_values[0].first.velocity*0.99;
	}
}

hdsim::hdsim(double cfl, vector<Primitive> const& cells, vector<double> const& edges, Interpolation const& interp,
	IdealGas const& eos, RiemannSolver const& rs,SourceTerm const& source, Geometry const& geo, 
	const double AMR_ratio):cfl_(cfl),
	cells_(cells),edges_(edges),interpolation_(interp),eos_(eos),rs_(rs),time_(0),cycle_(0),TotalEcool_(0),
	extensives_(vector<Extensive>()),source_(source),geo_(geo), AMR_ratio_(AMR_ratio), dt_suggest_(-1)
{
#ifdef RICH_MPI
	size_t Ntotal = cells.size();
	int rank = 0,ws=0;
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
		double vol = geo_.GetVolume(edges_,i);;
		extensives_[i].mass = cells_[i].density*vol;
		extensives_[i].momentum = extensives_[i].mass*cells_[i].velocity;
		extensives_[i].energy = 0.5*extensives_[i].momentum*extensives_[i].momentum / extensives_[i].mass +
			eos_.dp2e(cells_[i].density, cells_[i].pressure)*extensives_[i].mass;
		extensives_[i].entropy = eos_.dp2s(cells_[i].density, cells_[i].pressure)*extensives_[i].mass;
		cells_[i].energy = eos_.dp2e(cells_[i].density, cells_[i].pressure);
		cells_[i].LastCool = time_;
		extensives_[i].et = cells_[i].energy*extensives_[i].mass;
	}
}


hdsim::~hdsim()
{}

namespace
{
	void MoveGrid(double vgrid, std::vector<double> &edges,double dt)
	{
		size_t N = edges.size();
		for (size_t i = 0; i < N; ++i)
			edges[i] += vgrid * dt;
	}

	double GetTimeStep(vector<Primitive> const& cells, vector<double> const& edges, IdealGas const& eos, double cfl,
		SourceTerm const&source,double &dt_suggest,double vgrid = 0)
	{
		double force_inverse_dt = source.GetInverseTimeStep(edges);
		double dt_1 = eos.dp2c(cells[0].density, cells[0].pressure) / (edges[1] - edges[0]);
		size_t N = cells.size();
		for (size_t i = 1; i < N; ++i)
		{
			dt_1 = std::max(dt_1,2*(std::abs(cells[i].velocity-vgrid)+ eos.dp2c(cells[i].density, cells[i].pressure)) / (edges[i + 1] - edges[i]));
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
		return cfl/dt_1;
	}

	void GetFluxes(vector<pair<Primitive,Primitive> > const& interp_values, RiemannSolver const&rs,std::vector<Extensive> &res,IdealGas const& eos,double vgrid=0)
	{
		size_t N = interp_values.size();
		res.resize(N);
		for (int i = 0; i < N; ++i)
			res[i] = rs.SolveRS(interp_values[i].first, interp_values[i].second,eos,vgrid);
	}

	void UpdateExtensives(vector<Extensive> &cells, std::vector<Extensive> &fluxes,double dt,Geometry const& geo,
		std::vector<double> const& edges)
	{
		size_t N = cells.size();
		for (size_t i = 0; i < N; ++i)
		{
			double v = cells[i].momentum / cells[i].mass;
			double dP = -(fluxes[i + 1].momentum * geo.GetArea(edges[i+1]) - fluxes[i].momentum * geo.GetArea(edges[i]));
			double oldEk = cells[i].momentum*cells[i].momentum / cells[i].mass;
			cells[i].momentum += dP*dt;
			double dE = -(fluxes[i + 1].energy* geo.GetArea(edges[i + 1]) - fluxes[i].energy* geo.GetArea(edges[i]))*dt;
			cells[i].energy += dE;
			cells[i].mass += -(fluxes[i + 1].mass* geo.GetArea(edges[i + 1]) - fluxes[i].mass* geo.GetArea(edges[i]))*dt;
			double newEk = cells[i].momentum*cells[i].momentum / cells[i].mass;
			cells[i].et += dE - 0.5*(newEk - oldEk);	
			cells[i].entropy += -(fluxes[i + 1].entropy* geo.GetArea(edges[i + 1]) - fluxes[i].entropy* geo.GetArea(edges[i]))*dt;
		}
	}

	bool ShouldUseEntropy(Primitive const& cell, size_t index,double et,double EntropyEt)
	{
		double ek = cell.velocity*cell.velocity;
		if (et < 0)
			return true;
		if (et > 0.001*ek)
			return false;
		return true;
	}

	void UpdateCells(vector<Extensive> &extensive, vector<double> const& edges, IdealGas const& eos,
		vector<Primitive> &cells, Geometry const& geo)
	{
		size_t N = cells.size();
		for (int i = 0; i < N; ++i)
		{
			double vol = geo.GetVolume(edges, i);
			cells[i].density = extensive[i].mass / vol;
			cells[i].velocity = extensive[i].momentum / extensive[i].mass;
			cells[i].energy = extensive[i].et / extensive[i].mass;
			double et = cells[i].energy;
			if (ShouldUseEntropy(cells[i], i,et,eos.dp2e(cells[i].density,eos.sd2p(cells[i].entropy,cells[i].density))))
				cells[i].pressure = eos.sd2p(extensive[i].entropy/extensive[i].mass, cells[i].density);
			else
				cells[i].pressure = eos.de2p(cells[i].density, et);
			cells[i].pressure = std::max(cells[i].pressure, cells[i].density*cells[i].velocity*cells[i].velocity*0.0001);
			assert(cells[i].pressure > 0);
			et = extensive[i].mass*eos.dp2e(cells[i].density, cells[i].pressure);
			extensive[i].energy = 0.5*extensive[i].momentum*extensive[i].momentum / extensive[i].mass +	et;
			extensive[i].et = et;
			cells[i].entropy = eos.dp2s(cells[i].density, cells[i].pressure);
			assert(cells[i].entropy > 0);
			extensive[i].entropy = cells[i].entropy*extensive[i].mass;
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
		extensives_[i].mass = vol*cells[i].density;
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
	interpolation_.GetInterpolatedValues(cells_, edges_, interp_values_,time_);
	GetFluxes(interp_values_, rs_, fluxes_,eos_);
	double dt = GetTimeStep(cells_, edges_, eos_, cfl_, source_,dt_suggest_);


	/*if (BoundarySolution_ != 0)
	{
		pair<RSsolution,RSsolution> bvalues = BoundarySolution_->GetBoundaryValues(cells_);
		if (BoundarySolution_->ShouldCalc().first)
			rs_values_[0] = bvalues.first;
		if (BoundarySolution_->ShouldCalc().second)
			rs_values_.back() = bvalues.second;
	}*/

	UpdateExtensives(extensives_, fluxes_, dt,geo_,edges_);
	source_.CalcForce(edges_, cells_, time_, extensives_,dt);
	UpdateCells(extensives_, edges_, eos_, cells_,geo_);
	time_ += dt;
	++cycle_;
	AMR();
}

void hdsim::TimeAdvance2()
{
	interpolation_.GetInterpolatedValues(cells_, edges_, interp_values_, time_);
	double vgrid = GetVGrid(interp_values_);
#ifdef RICH_MPI
	double torecv = 0;
	MPI_Bcast(&vgrid,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	vgrid = torecv;
#endif
	GetFluxes(interp_values_, rs_, fluxes_, eos_,vgrid);

	/*if (BoundarySolution_ != 0)
	{
		pair<RSsolution, RSsolution> bvalues = BoundarySolution_->GetBoundaryValues(cells_);
		if (BoundarySolution_->ShouldCalc().first)
			rs_values_[0] = bvalues.first;
		if (BoundarySolution_->ShouldCalc().second)
			rs_values_.back() = bvalues.second;
	}*/
	double dt = GetTimeStep(cells_, edges_, eos_, cfl_, source_,dt_suggest_,vgrid);
	if (cycle_ == 0)
		dt *= 0.005;

	vector<Extensive> old_extensive(extensives_);
	std::vector<double> oldedges(edges_);
	UpdateExtensives(extensives_, fluxes_, 0.5*dt,geo_,edges_);
	MoveGrid(vgrid,edges_,0.5*dt);
	source_.CalcForce(edges_, cells_, time_, extensives_, 0.5*dt);	
	UpdateCells(extensives_, edges_, eos_, cells_, geo_);
	time_ += 0.5*dt;

	interpolation_.GetInterpolatedValues(cells_, edges_, interp_values_, time_);
	GetFluxes(interp_values_, rs_, fluxes_,eos_,vgrid);
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
	extensives_ = old_extensive;
	UpdateExtensives(extensives_, fluxes_, dt,geo_,edges_);
	source_.CalcForce(edges_, cells_, time_, extensives_, dt);
	edges_ = oldedges;
	MoveGrid(vgrid,edges_,dt);
#ifdef RICH_MPI
	RedistributeExtensives(extensives_,edges_,cells_);
#endif
	UpdateCells(extensives_, edges_, eos_, cells_,geo_);
	time_ += 0.5*dt;
	++cycle_;
	//AMR();
}

void hdsim::AMR(void)
{
	
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

#ifndef SOURCETERM_HPP
#define SOURCETERM_HPP 1
#define NOMINMAX
#include <algorithm>
#include "Primitive.hpp"
#include "Extensive.hpp"
#include "Geometry.hpp"
#include <vector>

using namespace std;

class SourceTerm
{
public:
	virtual void CalcForce(vector<double> const& edges, vector<Primitive> const& cells, double time,
		vector<Extensive> &extensives,double dt)const=0;

	virtual double GetInverseTimeStep(vector<double> const& edges)const;

	virtual ~SourceTerm();
};


class ZeroForce : public SourceTerm
{

	void CalcForce(vector<double> const& /*edges*/, vector<Primitive> const& cells, double /*time*/,
		vector<Extensive> & extensives,double /*dt*/)const
	{
		return;
	}
};

class SphericalForce : public SourceTerm
{
private:
	Spherical geo_;
public:
	SphericalForce() :geo_(Spherical()) {}

	void CalcForce(vector<double> const& edges, vector<Primitive> const& cells, double /*time*/,
		vector<Extensive> & extensives, double dt)const
	{
		size_t Ncells = cells.size();
		for (size_t i = 0; i < Ncells; ++i)
		{
			double dP = (geo_.GetArea(edges[i + 1]) - geo_.GetArea(edges[i]))*cells[i].pressure*dt;
			//double dP = geo_.GetVolume(edges, i) * 2 * cells[i].pressure * dt / (0.5*(edges[i + 1] + edges[i]));
			extensives[i].momentum += dP;
			extensives[i].et -= dP*cells[i].velocity;
		}
		return;
	}
};

class CylindricalForce : public SourceTerm
{
private:
	Cylindrical geo_;
public:
	CylindricalForce() :geo_(Cylindrical()) {}

	void CalcForce(vector<double> const& edges, vector<Primitive> const& cells, double /*time*/,
		vector<Extensive> & extensives, double dt)const
	{
		size_t Ncells = cells.size();
		for (size_t i = 0; i < Ncells; ++i)
		{
			double dP = (geo_.GetArea(edges[i + 1]) - geo_.GetArea(edges[i]))*cells[i].pressure*dt;
			//double dP = geo_.GetVolume(edges, i) * 2 * cells[i].pressure * dt / (0.5*(edges[i + 1] + edges[i]));
			extensives[i].momentum += dP;
			extensives[i].et -= dP * cells[i].velocity;
		}
		return;
	}
};

class PointGravity : public SourceTerm
{
private:
	const double mass_;
	const double L_;
	Geometry const& geo_;
	mutable double dt_1_;
public:
	PointGravity(double M,Geometry const& geo,double L = 0) :mass_(M),geo_(geo),L_(L),dt_1_(0) {}

	void CalcForce(vector<double> const& edges, vector<Primitive> const& cells, double /*time*/,
		vector<Extensive> & extensives, double dt)const
	{
		size_t Ncells = cells.size();
		dt_1_ = 0;
		for (size_t i = 0; i < Ncells; ++i)
		{
			double r = 0.5*(edges[i + 1] + edges[i]);
			double g = -mass_ / (r*r)+L_*L_/(r*r*r);
			extensives[i].momentum += g * cells[i].density*dt*geo_.GetVolume(edges,i);
			extensives[i].energy += g * cells[i].density*dt*geo_.GetVolume(edges,i)*cells[i].velocity;
			dt_1_ = std::max(dt_1_, std::sqrt(std::abs(g) / (edges[i + 1] - edges[i])));
			dt_1_ = std::max(dt_1_, std::abs(g / std::max(1e-20,cells[i].velocity)));
		}
		dt_1_ *= 8;
		return;
	}

	double GetInverseTimeStep(vector<double> const& edges)const
	{
		return dt_1_;
	}
};

class SeveralSources : public SourceTerm
{
private:
	vector<SourceTerm const*> sources_;
public:
	SeveralSources(vector<SourceTerm const*> sources) :sources_(sources) {}

	void CalcForce(vector<double> const& edges, vector<Primitive> const& cells, double time,
		vector<Extensive> & extensives, double dt)const
	{
		for (size_t i = 0; i < sources_.size(); ++i)
			sources_[i]->CalcForce(edges, cells, time, extensives, dt);
	}

	double GetInverseTimeStep(vector<double> const& edges)const
	{
		double res = 0;
		for (size_t i = 0; i < sources_.size(); ++i)
			res = std::max(res,sources_[i]->GetInverseTimeStep(edges));
		return res;
	}
};

#endif //SOURCETERM_HPP
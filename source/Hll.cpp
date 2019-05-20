#include "Hll.hpp"
#include <algorithm>
#include <cmath>

namespace
{
	class WaveSpeeds
	{
	public:

		WaveSpeeds(double left_i,
			double center_i,
			double right_i, double ps_i) :
			left(left_i),
			center(center_i),
			right(right_i), ps(ps_i) {}

		WaveSpeeds& operator=(WaveSpeeds const& ws)
		{
			left = ws.left;
			center = ws.center;
			right = ws.right;
			ps = ws.ps;
			return *this;
		}

		double left;
		double center;
		double right;
		double ps;
	};

	WaveSpeeds estimate_wave_speeds(Primitive const& left, Primitive const& right,IdealGas const& eos, double &pstar)
	{
		double cl = eos.dp2c(left.density, left.pressure);
		double cr = eos.dp2c(right.density, right.pressure);
		double sl = std::min(left.velocity - cl, right.velocity - cr);
		double sr = std::max(left.velocity + cl, right.velocity + cr);
		if (pstar > 0)
		{
			sl = left.velocity - cl * (pstar > left.pressure ? std::sqrt(0.8*(pstar / left.pressure - 1) + 1) : 1);
			sr = right.velocity + cr * (pstar > right.pressure ? std::sqrt(0.8*(pstar / right.pressure - 1) + 1) : 1);
		}
		Extensive Fl;
		Fl.mass = left.velocity*left.density;
		Fl.momentum = left.velocity*Fl.mass + left.pressure;
		Fl.energy = left.velocity*(left.pressure + left.energy*left.density + 0.5*Fl.mass*left.velocity);
		Extensive Fr;
		Fr.mass = right.velocity*right.density;
		Fr.momentum = right.velocity*Fr.mass + right.pressure;
		Fr.energy = right.velocity*(right.pressure + right.energy*right.density + 0.5*Fr.mass*right.velocity);
		Extensive Ull;
		double denom = 1.0 / (sr - sl);
		Ull.mass = (sr*right.density - sl * left.density + Fl.mass - Fr.mass)*denom;
		Ull.momentum = (sr*right.density*right.velocity - sl * left.density*left.velocity + Fl.momentum - Fr.momentum)*denom;
		Ull.energy =(sr*right.density*(right.energy + 0.5*right.velocity*right.velocity) - sl*left.density*(0.5*left.velocity*left.velocity + left.energy) 
			+ Fl.energy - Fr.energy)*denom;
		pstar = eos.de2p(Ull.mass, std::max(0.0,(Ull.energy - Ull.momentum*Ull.momentum*0.5 )/ Ull.mass));
		WaveSpeeds ws(sl, 0, sr, pstar);
		return ws;
	}

}

Hll::Hll(bool iter):iter_(iter)
{}


Hll::~Hll()
{}


Extensive Hll::SolveRS(Primitive const& left, Primitive const& right, IdealGas const& eos) const
{
	double cl = eos.dp2c(left.density, left.pressure);
	double cr = eos.dp2c(right.density, right.pressure);
	double sl = std::min(left.velocity - cl,right.velocity-cr);
	double sr = std::max(left.velocity + cl, right.velocity + cr);
	if (iter_)
	{
		double pstar = 0;
		WaveSpeeds ws = estimate_wave_speeds(left, right,eos, pstar);
		ws = estimate_wave_speeds(left, right, eos, pstar);
		ws = estimate_wave_speeds(left, right, eos, pstar);
		ws = estimate_wave_speeds(left, right, eos, pstar);
		ws = estimate_wave_speeds(left, right, eos, pstar);
		ws = estimate_wave_speeds(left, right, eos, pstar);
		sl = ws.left;
		sr = ws.right;
	}
	if (sl > 0)
	{
		Extensive res;
		res.mass = left.velocity*left.density;
		res.momentum = left.velocity*res.mass + left.pressure;
		res.energy = left.velocity*(left.pressure + left.energy*left.density + 0.5*res.mass*left.velocity);
		res.entropy = res.mass*left.entropy;
		return res;
	}
	if (sr < 0)
	{
		Extensive res;
		res.mass = right.velocity*right.density;
		res.momentum = right.velocity*res.mass + right.pressure;
		res.energy = right.velocity*(right.pressure + right.energy*right.density + 0.5*res.mass*right.velocity);
		res.entropy = res.mass*right.entropy;
		return res;
	}
	Extensive Fl;
	Fl.mass = left.velocity*left.density;
	Fl.momentum = left.velocity*Fl.mass + left.pressure;
	Fl.energy = left.velocity*(left.pressure + left.energy*left.density + 0.5*Fl.mass*left.velocity);
	Extensive Fr;
	Fr.mass = right.velocity*right.density;
	Fr.momentum = right.velocity*Fr.mass + right.pressure;
	Fr.energy = right.velocity*(right.pressure + right.energy*right.density + 0.5*Fr.mass*right.velocity);
	Extensive res;
	double denom = 1.0 / (sr - sl);
	res.mass = (sr*Fl.mass - sl * Fr.mass + sr * sl*(right.density - left.density))*denom;
	res.momentum = (sr*Fl.momentum - sl * Fr.momentum + sr * sl*(right.density*right.velocity - left.density*left.velocity))*denom;
	res.energy = (sr*Fl.energy - sl * Fr.energy + sr * sl*(right.density*(right.energy+0.5*right.velocity*right.velocity) - left.density*(0.5*left.velocity*left.velocity+left.energy)))*denom;
	res.entropy = res.mass*((res.mass > 0) ? left.entropy : right.entropy);
	return res;
}

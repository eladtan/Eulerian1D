#include "Hllc.hpp"
#include <algorithm>
#include <cmath>

namespace
{
	void PrimitiveToConserved(Primitive const& cell, Extensive &res)
	{
		res.mass = cell.density;
		res.momentum = cell.velocity;
		res.momentum *= res.mass;
		res.energy = res.mass*0.5*cell.velocity*cell.velocity + cell.energy*res.mass;
	}

	void starred_state(Primitive const& state, double sk, double ss, Extensive &res)
	{
		const double dk = state.density;
		const double pk = state.pressure;
		const double uk = state.velocity;
		const double ds = dk * (sk - uk) / (sk - ss);
		const double ek = state.density*(0.5*state.velocity*state.velocity + state.energy);
		res.mass = ds;
		res.momentum = ds * ss;
		res.energy = ek * ds / dk + ds * (ss - uk)*(ss + pk / dk / (sk - uk));
	}

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

	double Hll_pstar(Primitive const& left, Primitive const& right, IdealGas const& eos)
	{
		double cl = eos.dp2c(left.density, left.pressure);
		double cr = eos.dp2c(right.density, right.pressure);
		double sl = std::min(left.velocity - cl, right.velocity - cr);
		double sr = std::max(left.velocity + cl, right.velocity + cr);
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
		Ull.energy = (sr*right.density*(right.energy + 0.5*right.velocity*right.velocity) - sl * left.density*(0.5*left.velocity*left.velocity + left.energy)
			+ Fl.energy - Fr.energy)*denom;
		double pstar = eos.de2p(Ull.mass, std::max(0.0,(Ull.energy - Ull.momentum*Ull.momentum*0.5) / Ull.mass));
		return pstar;
	}

	WaveSpeeds estimate_wave_speeds(Primitive const& left, Primitive const& right, IdealGas const &eos, double pstar)
	{
		double cl = 0, cr = 0;
		const double dl = left.density;
		const double pl = left.pressure;
		const double vl = left.velocity;
		cl = eos.dp2c(dl, pl);
		const double dr = right.density;
		const double pr = right.pressure;
		const double vr = right.velocity;
		cr = eos.dp2c(dr, pr);
		const double sl = vl - cl * (pstar > pl ? std::sqrt(0.8*(pstar / pl - 1) + 1) : 1);
		const double sr = vr + cr * (pstar > pr ? std::sqrt(0.8*(pstar / pr - 1) + 1) : 1);
		const double denom = 1.0 / (dl*(sl - vl) - dr * (sr - vr));
		const double ss = (pr - pl + dl * vl*(sl - vl) - dr * vr*(sr - vr)) *denom;
		const double ps = std::max(0.0, dl * (sl - vl)*(pr - dr * (vr - vl)*(sr - vr)) *denom - pl * dr*(sr - vr) *denom);
		return WaveSpeeds(sl, ss, sr, ps);
	}
}

Hllc::Hllc(bool iter) :iter_(iter)
{}


Extensive Hllc::SolveRS(Primitive const& left, Primitive const& right, IdealGas const& eos) const
{
	double pstar = Hll_pstar(left, right, eos);
	pstar = std::max(pstar, 0.0);
	pstar = 0;

	WaveSpeeds ws2 = estimate_wave_speeds(left, right, eos, pstar);

	/*double old_ps = ws2.ps;
	ws2 = estimate_wave_speeds(left, right, eos, ws2.ps);
	size_t counter = 0;
	while (ws2.ps > 1.01 * old_ps || old_ps > 1.01 * ws2.ps)
	{
		old_ps = ws2.ps;
		ws2 = estimate_wave_speeds(left, right, eos, ws2.ps);
		++counter;
	}*/
	Extensive Fl;
	Fl.mass = left.velocity*left.density;
	Fl.momentum = left.velocity*Fl.mass + left.pressure;
	Fl.energy = left.velocity*(left.pressure + left.energy*left.density + 0.5*Fl.mass*left.velocity);
	Fl.entropy = left.entropy*Fl.mass;
	Extensive Fr;
	Fr.mass = right.velocity*right.density;
	Fr.momentum = right.velocity*Fr.mass + right.pressure;
	Fr.energy = right.velocity*(right.pressure + right.energy*right.density + 0.5*Fr.mass*right.velocity);
	Fr.entropy = right.entropy*Fr.mass;

	Extensive ul, ur, usl, usr;
	PrimitiveToConserved(left, ul);
	PrimitiveToConserved(right, ur);
	starred_state(left, ws2.left, ws2.center, usl);
	starred_state(right, ws2.right, ws2.center, usr);


	if (ws2.left > 0)
		return Fl;
	else
		if (ws2.left <= 0 && ws2.center >= 0)
		{
			Fl.mass += ws2.left*(usl.mass - ul.mass);
			Fl.momentum += ws2.left*(usl.momentum - ul.momentum);
			Fl.energy += ws2.left*(usl.energy - ul.energy);
			Fl.entropy = Fl.mass*left.entropy;
			return Fl;
		}
		else
			if (ws2.center < 0 && ws2.right >= 0)
			{
				Fr.mass += ws2.right*(usr.mass - ur.mass);
				Fr.momentum += ws2.right*(usr.momentum - ur.momentum);
				Fr.energy += ws2.right*(usr.energy - ur.energy);
				Fr.entropy = Fr.mass*right.entropy;
				return Fr;
			}
			else
				if (ws2.right < 0)
					return Fr;
				else
					throw;
}

Hllc::~Hllc()
{
}

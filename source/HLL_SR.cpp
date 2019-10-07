#include "Hll_SR.hpp"
#include <algorithm>
#include <cmath>
#include <cassert>
#include "universal_error.hpp"

namespace
{
	Extensive SR_Primitive2Flux(Primitive const& p, double edge_vel)
	{
		const double g_2 = 1.0 / (1 - p.velocity * p.velocity);
		double gamma = std::sqrt(g_2);
		Extensive res;
		res.mass = p.density*gamma*(p.velocity - edge_vel);
		res.momentum = p.pressure + p.density*(p.energy + 1)*g_2*p.velocity*
			(p.velocity - edge_vel);
		res.energy = p.density*p.energy*g_2*p.velocity + p.density*p.velocity*
			(g_2 - gamma) - edge_vel * (p.density*p.energy*g_2 - p.pressure + p.density*
			(g_2 - gamma));
		res.entropy = 0;
		res.et = 0;
		return res;
	}

	class WaveSpeeds
	{
	public:

		WaveSpeeds(double left_i,
			double right_i) :
			left(left_i),
			right(right_i) {}

		WaveSpeeds& operator=(WaveSpeeds const& ws)
		{
			left = ws.left;
			right = ws.right;
			return *this;
		}

		double left;
		double right;
	};

	WaveSpeeds estimate_wave_speeds(Primitive const& left, Primitive const& right, IdealGas const& eos)
	{
		double cl = 0, cr = 0;
#ifdef RICH_DEBUG
		try
		{
#endif
			cl = eos.dp2c(left.density, left.pressure);
#ifdef RICH_DEBUG
		}
		catch (UniversalError &eo)
		{
			eo.AddEntry("Error in HLLC_SR left sound speed", 0);
			throw eo;
		}
#endif
		const double sig_left = cl * cl*(1 - left.velocity * left.velocity)
			/ (1 - cl * cl);
		double disct = std::sqrt(sig_left*(1 - left.velocity*left.velocity + sig_left));
		const double lamda_minus_left = (left.velocity - disct) / (1 + sig_left);
		const double lamda_plus_left = (left.velocity + disct) / (1 + sig_left);
#ifdef RICH_DEBUG
		try
		{
#endif
			cr = eos.dp2c(right.density, right.pressure);
#ifdef RICH_DEBUG
		}
		catch (UniversalError &eo)
		{
			eo.AddEntry("Error in HLLC_SR right sound speed", 0);
			throw eo;
		}
#endif
		const double sig_right = cr * cr*(1 - right.velocity * right.velocity) / (1 - cr * cr);
		disct = std::sqrt(sig_right*(1 - right.velocity*right.velocity + sig_right));
		const double lamda_minus_right = (right.velocity - disct) / (1 + sig_right);
		const double lamda_plus_right = (right.velocity + disct) / (1 + sig_right);
		const double sl = std::min(lamda_minus_right, lamda_minus_left);
		const double sr = std::max(lamda_plus_left, lamda_plus_right);
		//double sl = -1;
		//double sr = 1;
		assert(sr > sl);
		return WaveSpeeds(sl, sr);
	}

}

Hll_SR::Hll_SR()
{}


Hll_SR::~Hll_SR()
{}


Extensive Hll_SR::SolveRS(Primitive const& left, Primitive const& right, IdealGas const& eos, double vface) const
{
	WaveSpeeds ws = estimate_wave_speeds(left, right, eos);
	double sr = ws.right, sl = ws.left;
	if (ws.left > vface)
	{
		Extensive res = SR_Primitive2Flux(left, vface);
		return res;
	}
	if (ws.right < vface)
	{
		Extensive res = SR_Primitive2Flux(left, vface);
		return res;
	}
	Extensive Fl = SR_Primitive2Flux(left, 0);
	double gl = std::sqrt(1.0 / (1 - left.velocity*left.velocity));
	Extensive Fr = SR_Primitive2Flux(right, 0);
	double gr = std::sqrt(1.0 / (1 - right.velocity*right.velocity));
	Extensive res;
	double denom = 1.0 / (ws.right - ws.left);
	Extensive Ull;
	Ull.mass = (sr*right.density*gr - sl * left.density*gl + Fl.mass - Fr.mass)*denom;
	Ull.momentum = (sr*(right.energy + 1)*right.density*right.velocity*gr*gr
		- sl * (left.energy + 1)*left.density*left.velocity*gl*gl + Fl.momentum - Fr.momentum)*denom;
	Ull.energy = (sr*(right.energy*right.density*gr*gr - right.pressure +
		right.density*(gr*gr - gr)) - sl * (left.energy*left.density*gl*gl - left.pressure +
			left.density*(gl*gl - gl)) + Fl.energy - Fr.energy)*denom;
	res.mass = (sr*Fl.mass - sl * Fr.mass + sr * sl*
		(right.density*gr - left.density*gl))*denom;
	res.momentum = (sr*Fl.momentum - sl * Fr.momentum + sr * sl*
		((right.energy + 1)*right.density*right.velocity*gr*gr -
		(left.energy + 1)*left.density*left.velocity*gl*gl))*denom;
	res.energy = (sr*Fl.energy - sl * Fr.energy + sr * sl*(
		right.energy*right.density*gr*gr - right.pressure + right.density*(gr*gr - gr)
		- (left.energy*left.density*gl*gl - left.pressure +
			left.density*(gl*gl - gl))))*denom;
	res.mass -= vface * Ull.mass;
	res.momentum -= vface * Ull.momentum;
	res.energy -= vface * Ull.energy;
	res.entropy = res.mass*((res.mass > 0) ? left.entropy : right.entropy);
	return res;
}

#include "ExactRS.hpp"
#include "universal_error.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <limits>
#if defined(_MSC_VER)
/* Microsoft C/C++-compatible compiler */
#include <intrin.h>
#else
#include <x86intrin.h>
#endif

namespace
{
	double fastsqrt(double x)
	{
		if (x<std::numeric_limits<float>::min() || x > std::numeric_limits<float>::max())
			return std::sqrt(x);
		double res = static_cast<double>(_mm_cvtss_f32(_mm_rsqrt_ss(_mm_set_ss(static_cast<float>(x)))));
		return x*res*(1.5 - 0.5*res*res*x);
	}

	double CalcFrarefraction(Primitive const& cell, double p,double gamma)
	{
		double cs = fastsqrt(gamma*cell.pressure/ cell.density);
		return 2 * cs*(std::pow(p / cell.pressure, (gamma - 1) / (2 * gamma)) - 1) / (gamma - 1);
	}

	double CalcFshock(Primitive const& cell, double p,double gamma)
	{
		double A = 2 / ((gamma + 1) *cell.density);
		double B = (gamma - 1)*cell.pressure / (gamma + 1);
		return (p - cell.pressure)*fastsqrt(A / (p + B));
	}

	double dCalcFshock(Primitive const& cell, double p, double gamma)
	{
		double A = 2 / ((gamma + 1)*cell.density);
		double B = (gamma - 1)*cell.pressure / (gamma + 1);
		return fastsqrt(A / (B + p))*(1 - (p - cell.pressure) / (2 * (B + p)));
	}

	double dCalcFrarefraction(Primitive const& cell, double p, double gamma)
	{	
		double cs = fastsqrt(gamma*cell.pressure/ cell.density);
		return std::pow(cell.pressure / p, (gamma + 1) / (2 * gamma)) / (cell.density*cs);
	}

	double GetFirstGuess(Primitive const & left, Primitive const & right,double gamma)
	{
		double Pmax = std::max(left.pressure, right.pressure);
		double Pmin = std::min(left.pressure, right.pressure);
		double Q = Pmax / Pmin;
		double csl = fastsqrt(gamma*left.pressure / left.density);
		double csr = fastsqrt(gamma*right.pressure / right.density);
		double res;
		// First try Ppvrs
		double pvrs = 0.5*(Pmax + Pmin) + 0.125*(left.velocity - right.velocity)*(left.density + right.density)*
			(csl + csr);
		if (Q<2 && pvrs >=Pmin && pvrs <= Pmax)
			return pvrs;
		if (pvrs <= Pmin) // Use Two rarefactions
		{
			double z = (gamma - 1) / (2 * gamma);
			res = std::pow((csl + csr - (gamma - 1)*(right.velocity - left.velocity)*0.5) / (csl*
				std::pow(left.pressure, -z)	+ csr*std::pow(right.pressure, -z)), 1 / z);
		}
		else // Two shocks
		{
			double Al = 2 / ((gamma + 1)*left.density);
			double Ar = 2 / ((gamma + 1)*right.density);
			double Bl = (gamma - 1)*left.pressure / (gamma + 1);
			double Br = (gamma - 1)*right.pressure / (gamma + 1);
			res = std::max(0.0, res);
			double gl = fastsqrt(Al / (res + Bl));
			double gr = fastsqrt(Ar / (res + Br));
			res = (gl*left.pressure + gr*right.pressure + left.velocity - right.velocity) / (gl + gr);
			if (res < Pmin)
				res = pvrs;
		}
		return res;
	}

	double GetValue(Primitive const & left, Primitive const & right,double p,double gamma)
	{
		double res = right.velocity-left.velocity;
		if (p > left.pressure)
			res += CalcFshock(left, p, gamma);
		else
			res += CalcFrarefraction(left, p,gamma);
		if (p > right.pressure)
			res += CalcFshock(right, p, gamma);
		else
			res += CalcFrarefraction(right, p, gamma);
		return res;
	}

	double GetdValue(Primitive const & left, Primitive const & right, double p, double gamma)
	{
		double res = 0;
		if (p > left.pressure)
			res += dCalcFshock(left, p, gamma);
		else
			res += dCalcFrarefraction(left, p, gamma);
		if (p > right.pressure)
			res += dCalcFshock(right, p, gamma);
		else
			res += dCalcFrarefraction(right, p, gamma);
		return res;
	}

	double Bisection(Primitive const& left, Primitive const& right,double gamma,
		double max_scale,double minp,double guess)
	{
		/*
		double a = std::min(left.pressure, right.pressure)*1e-10;
		double b = (left.pressure + right.pressure + left.density*left.velocity*left.velocity +
			right.density*right.velocity*right.velocity)*10;*/
		double a = guess*0.1;
		double b = guess * 10;
		double c = guess;
		double dp = 0;
		double valuec = GetValue(left, right, b , gamma);
		double valuea = GetValue(left, right, a, gamma);
		if (valuea*valuec > 0)
		{
			UniversalError eo("Same sign in RS");
			eo.AddEntry("Left density", left.density);
			eo.AddEntry("Left pressure", left.pressure);
			eo.AddEntry("Left velocity", left.velocity);
			eo.AddEntry("Right density", right.density);
			eo.AddEntry("Right pressure", right.pressure);
			eo.AddEntry("Right velocity", right.velocity);
			eo.AddEntry("Init pressure in bisection", guess);
			throw eo;
		}
		int counter = 0;
		while (((b-a) > 1e-10*(minp + c)) || (std::abs(valuec) > 1e-7*max_scale))
		{
			c = (a + b)*0.5;
			valuec = GetValue(left, right, c, gamma);
			if (valuec*valuea > 0)
			{
				a = c;
				valuea = valuec;
			}
			else
				b = c;
			++counter;
			if (counter > 300)
			{
				UniversalError eo("Too many iterations in RS");
				eo.AddEntry("Left density", left.density);
				eo.AddEntry("Left pressure", left.pressure);
				eo.AddEntry("Left velocity", left.velocity);
				eo.AddEntry("Right density", right.density);
				eo.AddEntry("Right pressure", right.pressure);
				eo.AddEntry("Right velocity", right.velocity);
				throw eo;
			}
		}
		return c;
	}

	std::pair<double, double> Get_p_u_star(Primitive const& left, Primitive const& right, double gamma)
	{
		const double eps = 1e-7;
		// Is there a vaccum?
		double dv = right.velocity - left.velocity;
		double soundspeeds = 2 * (fastsqrt(gamma*left.pressure / left.density) + fastsqrt(gamma*right.pressure / right.density)) / (gamma - 1);
		if (dv >= soundspeeds)
		{
			return std::pair<double, double>(0., 0.);
		}
		std::pair<double, double> res;
		res.first = GetFirstGuess(left, right, gamma);
		res.second = 0;
		if (res.first < 0)
		{
			res.first = 0;
			return res;
		}
		double value = GetValue(left, right, res.first, gamma);
		double Cs = std::max(fastsqrt(gamma*left.pressure / left.density), fastsqrt(gamma*right.pressure / right.density));
		double max_scale = std::max(Cs, std::max(std::abs(left.velocity), std::abs(right.velocity)));
		double dp = 0;
		double minp = std::min(left.pressure, right.pressure);
		double p = res.first;
		size_t counter = 0;
		while ((abs(dp) > eps*(minp + res.first)) || (std::abs(value) > eps*max_scale))
		{
			p = res.first;
			dp = value / GetdValue(left, right, res.first, gamma);
			res.first -= std::max(std::min(dp, res.first*0.5), -0.5*res.first);
			value = GetValue(left, right, res.first, gamma);
			++counter;
			if (counter > 30)
			{
				res.first = Bisection(left, right, gamma, max_scale, minp, res.first);
				break;
			}
		}
		double fr = (res.first > right.pressure) ? CalcFshock(right, res.first, gamma) : CalcFrarefraction(
			right, res.first, gamma);
		double fl = (res.first > left.pressure) ? CalcFshock(left, res.first, gamma) : CalcFrarefraction(
			left, res.first, gamma);
		res.second = 0.5*(left.velocity + right.velocity) + 0.5*(fr - fl);
		return res;
	}
}

ExactRS::ExactRS(double gamma):gamma_(gamma)
{}

ExactRS::~ExactRS()
{}

Extensive ExactRS::SolveRS(Primitive const& left, Primitive const& right, IdealGas const& eos,double vface) const
{
	// Solve RS
	std::pair<double, double> p_u_star = Get_p_u_star(left, right, gamma_);
	bool left_shock = p_u_star.first > left.pressure;
	bool right_shock = p_u_star.first > right.pressure;
	
	Extensive res;
	// Did we form a vacuum? 
	if (p_u_star.first == 0.)
		return res;
	// Find branch on
	if (p_u_star.second > vface)
	{
		double al = std::sqrt(left.pressure*gamma_ / left.density);
		if (left_shock)
		{
			double Sl = left.velocity - al * std::sqrt((gamma_ + 1)*p_u_star.first / (2 * gamma_*left.pressure + (gamma_ - 1) / (2 * gamma_)));
			if (Sl > vface)
			{
				res.mass = (left.velocity-vface)*left.density;
				res.momentum = (left.velocity-vface)*left.density*left.velocity + left.pressure;
				res.energy = (left.velocity-vface)*(left.pressure / (gamma_ - 1) + 0.5*left.density*left.velocity*left.velocity)
					+left.velocity*left.pressure;
				res.entropy = res.mass*left.entropy;
			}
			else
			{
				double rho_star = left.density*((p_u_star.first / left.pressure + (gamma_ - 1) / (gamma_ + 1)) / ((gamma_ - 1)*p_u_star.first / ((gamma_ + 1)*left.pressure) + 1));
				res.mass = (p_u_star.second-vface)*rho_star;
				res.momentum = p_u_star.second*res.mass + p_u_star.first;
				res.energy = (p_u_star.second-vface)*(p_u_star.first / (gamma_ - 1) + 0.5*p_u_star.second*rho_star*p_u_star.second)
					+ p_u_star.first*p_u_star.second;
				res.entropy = res.mass*left.entropy;
			}
		}
		else
		{
			double Shl = left.velocity - al;
			if (Shl > vface)
			{
				res.mass = (left.velocity-vface)*left.density;
				res.momentum = left.velocity*res.mass + left.pressure;
				res.energy = (left.velocity-vface)*(left.pressure / (gamma_ - 1)+ 0.5*left.density*left.velocity*left.velocity)
					+left.velocity*left.pressure;
				res.entropy = res.mass*left.entropy;
			}
			else
			{
				double al_star = al * std::pow(p_u_star.first / left.pressure, (gamma_ - 1) / (2 * gamma_));
				double Stl = p_u_star.second - al_star;
				if (Stl < vface)
				{
					double rho_fan = left.density*std::pow(p_u_star.first / left.pressure, 1.0 / gamma_);
					res.mass = (p_u_star.second - vface)*rho_fan;
					res.momentum = p_u_star.second*res.mass + p_u_star.first;
					res.energy = (p_u_star.second - vface)*(p_u_star.first / (gamma_ - 1) + rho_fan * p_u_star.second*p_u_star.second*0.5)
						+ p_u_star.first*p_u_star.second;
					res.entropy = res.mass*left.entropy;
				}
				else
				{
					double rho_fan = left.density*std::pow(2.0 / (gamma_ + 1) + (gamma_ - 1)*(left.velocity - vface) / ((gamma_ + 1)*al), 2.0 / (gamma_ - 1));
					double v = 2.0*(al + (gamma_ - 1)*left.velocity*0.5 + vface) / (gamma_ + 1);
					double p = left.pressure*std::pow(rho_fan / left.density, gamma_);
					res.mass = (v - vface)*rho_fan;
					res.momentum = v * res.mass + p;
					res.energy = (v - vface)*(p / (gamma_ - 1) + 0.5*rho_fan*v*v)+p*v;
					res.entropy = res.mass*left.entropy;
				}
			}
		}
	}
	else
	{
		double ar = std::sqrt(right.pressure*gamma_ / right.density);
		if (right_shock)
		{
			double Sr = right.velocity + ar * std::sqrt((gamma_ + 1)*p_u_star.first / (2 * gamma_*right.pressure + (gamma_ - 1) / (2 * gamma_)));
			if (Sr < vface)
			{
				res.mass = (right.velocity-vface)*right.density;
				res.momentum = right.velocity*res.mass + right.pressure;
				res.energy = (right.velocity-vface)*(right.pressure / (gamma_ - 1) + 0.5*right.density*right.velocity*right.velocity)
					+right.velocity*right.pressure;
				res.entropy = res.mass*right.entropy;
			}
			else
			{
				double rho_star = right.density*((p_u_star.first / right.pressure + (gamma_ - 1) / (gamma_ + 1)) / ((gamma_ - 1)*p_u_star.first / ((gamma_ + 1)*right.pressure) + 1));
				res.mass = (p_u_star.second-vface)*rho_star;
				res.momentum = p_u_star.second*res.mass + p_u_star.first;
				res.energy = (p_u_star.second-vface)*(p_u_star.first / (gamma_ - 1) + 0.5*rho_star*p_u_star.second*p_u_star.second)
					+p_u_star.first*p_u_star.second;;
				res.entropy = res.mass*right.entropy;
			}
		}
		else
		{
			double ar_star = ar * std::pow(p_u_star.first / right.pressure, (gamma_ - 1) / (2 * gamma_));
			double Str = p_u_star.second + ar_star;
			if (Str > vface)
			{
				double rho_fan = right.density*std::pow(p_u_star.first / right.pressure, 1.0 / gamma_);
				res.mass = (p_u_star.second-vface)*rho_fan;
				res.momentum = p_u_star.second*res.mass + p_u_star.first;
				res.energy = (p_u_star.second-vface)*(p_u_star.first / (gamma_ - 1)+0.5*rho_fan*p_u_star.second*p_u_star.second)
					+p_u_star.first*p_u_star.second;
				res.entropy = res.mass*right.entropy;
			}
			else
			{
				double Shr = right.velocity + ar;
				if (Shr > vface)
				{
					double rho_fan = right.density*std::pow(2.0 / (gamma_ + 1) - (gamma_ - 1)*(right.velocity - vface) / ((gamma_ + 1)*ar), 2.0 / (gamma_ - 1));
					double v = 2.0*(-ar + (gamma_ - 1)*right.velocity*0.5 + vface) / (gamma_ + 1);
					double p = right.pressure*std::pow(rho_fan / right.density, gamma_);
					res.mass = (v - vface) * rho_fan;
					res.momentum = v * res.mass + p;
					res.energy = (v - vface) * (p / (gamma_ - 1) + 0.5*rho_fan*v*v)+v*p;
					res.entropy = res.mass*right.entropy;
				}
				else
				{
					res.mass = (right.velocity - vface)*right.density;
					res.momentum = right.velocity*res.mass + right.pressure;
					res.energy = (right.velocity - vface)*(right.pressure / (gamma_ - 1) + 0.5*right.density*right.velocity*right.velocity)
						+right.velocity*right.pressure;
					res.entropy = res.mass*right.entropy;
				}
			}
		}
	}
	if (!(std::isfinite(res.mass)))
	{
		std::cout << "Bad ExactRS output, rhol = " << left.density << " rhor = " << right.density << " pl = " << left.pressure << " pr = " << right.pressure <<
			" ul = " << left.velocity << " ur = " << right.velocity << " rhostar = " << res.mass << " pstar = " << p_u_star.first << " ustar = " << p_u_star.second << std::endl;
		throw;
	}
	return res;
}

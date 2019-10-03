#include <cmath>
#include "ideal_gas.hpp"
#include "universal_error.hpp"

IdealGas::IdealGas(double AdiabaticIndex, bool SR):
  g_(AdiabaticIndex), SR_(SR) {}

double IdealGas::getAdiabaticIndex(void) const
{
  return g_;
}

double IdealGas::dp2e(double d, double p) const
{
	if(SR_)
		return p * g_ / (d*(g_ - 1));
	else
		return p/d/(g_-1);
}

double IdealGas::de2p(double d, double e) const
{
	if (e < 0)
		throw UniversalError("Negative thermal energy");
	if(SR_)
		return e * (g_ - 1)*d / g_;
	else
		return (g_-1)*e*d;
}

double IdealGas::dp2c(double d, double p) const
{
	if (d < 0 || p < 0)
	{
		UniversalError eo("Imaginary Cs");
		eo.AddEntry("Density", d);
		eo.AddEntry("Pressure", p);
		throw eo;
	}
	if(SR_)
		return std::sqrt(g_*p / (d + p * g_ / ((g_ - 1))));
	else
		return sqrt(g_*p/d);
}

double IdealGas::de2c(double d, double e) const
{
  double p = de2p(d, e);
  if(SR_)
	  return std::sqrt(g_*p / (d*(1 + e)));
  else
	  return sqrt(g_*p/d);
}

double IdealGas::dp2s(double d, double p) const
{
  return p*pow(d,-g_);
}

double IdealGas::sd2p(double s, double d) const
{
	if (d < 0 || s < 0)
	{
		UniversalError eo("Imaginary pressure");
		eo.AddEntry("Density", d);
		eo.AddEntry("Entropy", s);
		throw eo;
	}
  return s*pow(d,g_);
}

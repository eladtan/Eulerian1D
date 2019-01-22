#ifndef IDEAL_GAS_HPP
#define IDEAL_GAS_HPP 1

class IdealGas
{
private:

  double g_;

public:

  explicit IdealGas(double AdiabaticIndex);

  double getAdiabaticIndex(void) const;

  double dp2e(double d, double p) const;

  double de2p(double d, double e) const;

  double dp2c(double d, double p) const;

  double de2c(double d, double e) const;

  double dp2s(double d, double p) const;

  double sd2p(double s, double d) const;
};

#endif // IDEAL_GAS_HPP

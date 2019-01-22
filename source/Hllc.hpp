#ifndef HLLC_HPP
#define HLLC_HPP 1
#include <RiemannSolver.hpp>

class Hllc : public RiemannSolver
{
private:
	bool iter_;
public:
	Hllc(bool iter);

	Extensive SolveRS(Primitive const& left, Primitive const& right, IdealGas const& eos) const;

	~Hllc();
};

#endif

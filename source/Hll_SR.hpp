#ifndef HLL_SR_HPP
#define HLL_SR_HPP 1
#include <RiemannSolver.hpp>

class Hll_SR : public RiemannSolver
{
public:
	Hll_SR();

	Extensive SolveRS(Primitive const& left, Primitive const& right, IdealGas const& eos, double vface) const;

	~Hll_SR();
};

#endif

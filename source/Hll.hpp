#ifndef HLL_HPP
#define HLL_HPP 1
#include <RiemannSolver.hpp>

class Hll : public RiemannSolver
{
private:
	bool iter_;
public:
	Hll(bool iter);

	Extensive SolveRS(Primitive const& left, Primitive const& right, IdealGas const& eos,double vface) const;

	~Hll();
};

#endif

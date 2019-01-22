#ifndef RS_HPP
#define RS_HPP 1
#include <Extensive.hpp>
#include <Primitive.hpp>
#include <ideal_gas.hpp>

class RiemannSolver
{
public:
	
	virtual Extensive SolveRS(Primitive const& left, Primitive const& right, IdealGas const& eos) const = 0;

	virtual ~RiemannSolver() {};
};

#endif
#ifndef HLLC_SR_HPP
#define HLLC_SR_HPP 1
#include "RiemannSolver.hpp"
class HLLC_SR : public RiemannSolver
{
public:
	HLLC_SR();
	~HLLC_SR();
	Extensive SolveRS(Primitive const& left, Primitive const& right, IdealGas const& eos, 
		double vface = 0) const;
};

#endif
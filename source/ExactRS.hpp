#ifndef EXACTRS_HPP
#define EXACTRS_HPP 1

#include "RiemannSolver.hpp"

class ExactRS : public RiemannSolver
{
private:
	const double gamma_;

public:
	ExactRS(double gama);
	~ExactRS();

	Extensive SolveRS(Primitive const& left, Primitive const& right, IdealGas const& eos) const;
};


#endif //EXACTRS_HPP
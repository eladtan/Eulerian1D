#ifndef EXTENSIVE_HPP
#define EXTENSIVE_HPP 1
class Extensive
{
public:
	double mass;
	double momentum;
	double energy;
	double entropy;
	double et;

	Extensive();
	~Extensive();

	Extensive& operator+=(const Extensive &rhs);

	Extensive& operator*=(double s);

	Extensive& operator=(Extensive const& other);
};

#endif //EXTENSIVE_HPP
#include "Extensive.hpp"

Extensive::Extensive():mass(0),momentum(0),energy(0),entropy(0),et(0)
{
}


Extensive::~Extensive()
{
}

Extensive & Extensive::operator+=(const Extensive & rhs)
{
	this->mass += rhs.mass;
	this->momentum += rhs.momentum;
	this->energy += rhs.energy;
	this->entropy += rhs.entropy;
	this->et += rhs.et;
	return *this;
}

Extensive & Extensive::operator=(Extensive const & other)
{
	this->energy = other.energy;
	this->entropy = other.entropy;
	this->et = other.et;
	this->mass = other.mass;
	this->momentum = other.momentum;
	return *this;
}

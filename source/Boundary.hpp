#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP 1

#include "Primitive.hpp"
#include "ExactRS.hpp"
#include <vector>

using namespace std;

/*class BoundarySolution
{
public:
	~BoundarySolution();

	virtual pair<RSsolution,RSsolution> GetBoundaryValues(vector<Primitive> const& cells)const =0;

	virtual pair<bool, bool> ShouldCalc()const=0;
};

class VacuumInFlow : public BoundarySolution
{
private:
	const bool calc_left_, calc_right_;
public:
	VacuumInFlow(bool calc_left, bool calc_right);

	std::pair<RSsolution, RSsolution> GetBoundaryValues(vector<Primitive> const& cells)const;

	std::pair<bool, bool> ShouldCalc()const;
};

*/
class Boundary
{
public:
	virtual ~Boundary();

	virtual vector<Primitive> GetBoundaryValues(vector<Primitive> const& cells, vector<double> const&
		edges,size_t index,double time)const=0;
};

class RigidWall : public Boundary
{
public:
	vector<Primitive> GetBoundaryValues(vector<Primitive> const& cells, vector<double> const&
		edges, size_t index,double time)const;
};

class Ratchet : public Boundary
{
public:
	vector<Primitive> GetBoundaryValues(vector<Primitive> const& cells, vector<double> const&
		edges, size_t index, double time)const;
};

class RigidWall1 : public Boundary
{
public:
	vector<Primitive> GetBoundaryValues(vector<Primitive> const& cells, vector<double> const&
		edges, size_t index, double time)const;
};


class FreeFlow : public Boundary
{
public:
	vector<Primitive> GetBoundaryValues(vector<Primitive> const& cells, vector<double> const&
		edges, size_t index, double time)const;
};

class FreeFlow2 : public Boundary
{
public:
	vector<Primitive> GetBoundaryValues(vector<Primitive> const& cells, vector<double> const&
		edges, size_t index, double time)const;
};

class Periodic : public Boundary
{
public:
	vector<Primitive> GetBoundaryValues(vector<Primitive> const& cells, vector<double> const&
		edges, size_t index, double time)const;
};

class ConstantPrimitive : public Boundary
{
private:
	Primitive outer_;
public:
	ConstantPrimitive(Primitive outer);

	vector<Primitive> GetBoundaryValues(vector<Primitive> const& cells, vector<double> const&
		edges, size_t index, double time)const;
};

class SeveralBoundary : public Boundary
{
private:
	Boundary const& left_;
	Boundary const& right_;
public:
	SeveralBoundary(Boundary const& left, Boundary const& right);

	vector<Primitive> GetBoundaryValues(vector<Primitive> const& cells, vector<double> const&
		edges, size_t index, double time)const;

};

#endif //BOUNDARY_HPP

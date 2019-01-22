#define _USE_MATH_DEFINES
#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP 1
#include <cmath>
#include <vector>

#define M_PI 3.14159265358979323846

class Geometry
{
public:
	virtual inline double GetArea(double r)const = 0;
	virtual inline double GetVolume(std::vector<double> const& edges, size_t index)const = 0;
	~Geometry() {}
};

class Planar : public Geometry
{
public:
	Planar() {}
	inline double GetArea(double r)const
	{
		return 1.0;
	}

	inline double GetVolume(std::vector<double> const& edges, size_t index)const
	{
		return edges[index + 1] - edges[index];
	}
};

class Spherical : public Geometry
{
public:
	Spherical() {}
	inline double GetArea(double r)const
	{
		return 4 * M_PI*r*r;
	}

	inline double GetVolume(std::vector<double> const& edges, size_t index)const
	{
		double r1 = edges[index + 1];
		double r0 = edges[index];
		return 4 * M_PI*0.3333333333*(r1*r1*r1 - r0*r0*r0);
	}
};

class Cylindrical : public Geometry
{
public:
	Cylindrical() {}
	inline double GetArea(double r)const
	{
		return 2 * M_PI*r;
	}

	inline double GetVolume(std::vector<double> const& edges, size_t index)const
	{
		double r1 = edges[index + 1];
		double r0 = edges[index];
		return M_PI*0.3333333333*(r1*r1 - r0*r0);
	}
};
#endif //GEOMETRY_HPP

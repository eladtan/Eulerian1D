#ifndef PRIMITIVE_HPP
#define PRIMITIVE_HPP 1

class Primitive
{
public:
	double density;
	double pressure;
	double velocity;
	double entropy;
	double energy;
	unsigned char sticker;
	double LastCool;

	Primitive();

	Primitive(Primitive const& other);

	Primitive(double Density, double Pressure, double Velocity,double Entropy,double Energy,unsigned char Sticker);

	Primitive operator-(Primitive const& other)const;

	Primitive operator+(Primitive const& other)const;

	Primitive operator*(double s);
	
	Primitive operator/(double s);

	Primitive& operator=(Primitive const& other);
};

#endif //PRIMITIVE_HPP

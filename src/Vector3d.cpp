#include "Vector3d.h"

////////////////////// CONSTRUCTORS & DESTRUCTOR ///////////////////////

/// Default constructor that sets X,Y and Z to 0.
Vector3d::Vector3d()
{
	x = 0;
	y = 0;
	z = 0;
}

/// Constructor that creates a vector with the given (X, Y, Z). 
Vector3d::Vector3d(double x, double y, double z)
{
	this -> x = x;
	this -> y = y;
	this -> z = z;
}

Vector3d::Vector3d(Vector3d const& copy_from_me)
{
	x = copy_from_me.x;
	y = copy_from_me.y;
	z = copy_from_me.z;
}

Vector3d::~Vector3d()
{
	// Nothing?
}

////////////////////// VECTOR OPERATORS ///////////////////////

/// Returns the sum of the two vectors.
Vector3d Vector3d::operator+ (const Vector3d& vect) const
{
	return Vector3d(x + vect.x, y + vect.y, z + vect.z);
}

/// Returns the difference between the two vectors.
Vector3d Vector3d::operator- (const Vector3d& vect) const
{
	return Vector3d(x - vect.x, y - vect.y, z - vect.z);
}

/// Returns the component-wise multiplication of the vectors
Vector3d Vector3d::operator* (const Vector3d& vect) const
{
	return Vector3d(x * vect.x, y * vect.y, z * vect.z);
}

/// Returns the component-wise division of the vectors
Vector3d Vector3d::operator/ (const Vector3d& vect) const
{
	return Vector3d(x / vect.x, y / vect.y, z / vect.y);
}

/// Sets the vector equal to the parameter vector. 
void Vector3d::operator= (const Vector3d& vect)
{
	x = vect.x;
	y = vect.y;
	z = vect.z;
}

////////////////////// SCALAR OPERATORS ///////////////////////

/// Returns the vector with all components multiplied by the scalar parameter
Vector3d Vector3d::operator* (const double& scalar) const
{
	return Vector3d(x * scalar, y * scalar, z * scalar); 
}

/// Returns the vector with all components divided by the scalar parameter
Vector3d Vector3d::operator/ (const double& scalar) const
{
	return Vector3d(x / scalar, y / scalar, z / scalar);
}

/// Subtracts the scalar from all components. 
void Vector3d::operator-= (const double& scalar)
{
	x = x - scalar;
	y = y - scalar;
	z = z - scalar;
}

/// Adds the scalar to all components.
void Vector3d::operator+= (const double& scalar)
{
	x = x + scalar;
	y = y + scalar;
	z = z - scalar;
}


////////////////////// PRODUCTS ///////////////////////

/// returns the DOT product between the the two vectors
/// -- sum of the component-wise multiplication
/// -- Where A is this, and B is other. 
/// -- a1*b1 + a2*b2 + a3*b3
double Vector3d::dot (const Vector3d& vect) const
{
	return (x * vect.x + y * vect.y + z * vect.z);
}

/// Returns the CROSS product between the two vectors
/// -- An orthogonal Vector3d to the two. 
/// -- Where A is this, and B is other. 
/// -- c1 = a2 * b3 - a3 * b2
/// -- c2 = a3 * b1 - a1 * b3
/// -- c3 = a1 * b2 - a2 * b1
Vector3d Vector3d::cross (const Vector3d& vect) const
{
	Vector3d result;
	result.x = y * vect.z - z * vect.y;
	result.y = z * vect.x - x * vect.z;
	result.z = x * vect.y - y * vect.x;

	return result;
}



////////////////////// LENGTHS ///////////////////////

/// Returns the square of the length of the vector. Useful when you just want to compare two vectors
/// to see which is longest, as this avoids computing the square roots. 
double Vector3d::squaredLength() const
{
	return (x * x + y * y + z * z);
}

/// returns the length of the vector
double Vector3d::length() const
{
	return sqrt(squaredLength());
}


/// Returns a vector pointing in the same direction, but with unit length (length of 1)
Vector3d Vector3d::unit() const
{
	double length = this -> length();
	return Vector3d(x / length, y / length, z / length);
}
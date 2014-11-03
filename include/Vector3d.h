//																	    //
//	Purpose:	Implement a 3D vector in the MATHEMATICAL sense, not    //
//				in the dynamic array sense (as STL vector is)		    //
// 				Adapted from my Vector2f version. (Thanks Khan Academy) //
//	Author:		Christian Holt										    //
//	Last Edit:	16/10/2014											    //
//																	    //
//////////////////////////////////////////////////////////////////////////


#pragma once											// Makes sure it is only included once per compile.

#include <math.h>

class Vector3d
{
public:
	double x;
	double y;
	double z;

	 
	// Constructors & Destructor
	Vector3d();											// Default constructor that sets X and Y to 0. 
	Vector3d(double x, double y, double z);				// Constructor that creates a vector with the given (X, Y).
	Vector3d(Vector3d const& copy_from_me);				// Copy constructor
	~Vector3d();										// Destructor. 

	// Vector Operators
	Vector3d operator+ (Vector3d vector);				// Returns the sum of the two vectors.
	Vector3d operator- (Vector3d vector);  				// Returns the difference between the two vectors.
	Vector3d operator* (Vector3d vector);				// Returns the component-wise multiplication of the vectors
	Vector3d operator/ (Vector3d vector);				// Returns the component-wise division of the vectors
	void operator= (Vector3d vector);					// Sets the vector equal to the parameter vector. 

	// Scalar operators
	Vector3d operator* (double scalar);					// Returns the vector with all components multiplied by the scalar parameter
	Vector3d operator/ (double scalar);					// Returns the vector with all components divided by the scalar parameter
	void operator-= (double scalar);					// Subtracts the scalar from all components. 
	void operator+= (double scalar);					// Adds the scalar to all components.

	// products
	double dot (Vector3d vector);						// returns the DOT product between the the two vectors
	Vector3d cross(Vector3d vector);					// Returns the CROSS product between the two vectors (an Orthogonal Vector)

	// Lengths
	double length();									// returns the length of the vector
	Vector3d unit();									// Returns a vector pointing in the same direction, but with unit length (lenght of 1)
	double squaredLength();								// Returns the square of the length of the vector. Useful when you just want to compare two vectors
														// to see which is longest, as this avoids computing the square roots. 

	// TODO: add more functionality from:
	// http://scratchapixel.com/lessons/3d-basic-lessons/lesson-4-geometry/math-operations-on-points-and-vectors/
	// http://scratchapixel.com/lessons/3d-basic-lessons/lesson-4-geometry/how-does-matrix-work-part-1/ 
};
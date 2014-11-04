#include "sphere.h"

Sphere::Sphere()
{
	// Empty Constructor
}

Sphere::Sphere(double radius, Vector3d position, Vector3d emission,
		   Vector3d color, MaterialType materialType)
{
	this -> radius = radius;
	this -> position = position;
	this -> emission = emission;
	this -> color = color;
	this -> materialType = materialType;
}

// ----- The following comment rant is mostly for my own sake ----- 
// http://en.wikipedia.org/wiki/Line–sphere_intersection
// In Vector form:

// Sphere equation => |x - c|^2 = r^2
// -- x is point on the sphere
// -- c is the center point (this.position)
// -- r is the radius (this.radius)

// Equation for a line => X = o + dL
// -- d is the distance along the line from the starting point
// -- l is the direction of the line (ray.direction)
// -- o is the origin of the line (ray.origin)
// -- x is points on the line

// Combining the two equations searching for points that 
// are present on both the line and the sphere:
// Note that * = dot product because vectors. 
// d^2(l*l) + 2d(l * (o - c)) + (o - c) * (o - c) - r^2 = 0

// We can solve for d using the quadratic formula.
// http://en.wikipedia.org/wiki/Quadratic_formula 
// X = -b +- sqroot(b^2 - 4ac) / 2a
// ad^2 + bd + c = 0
// Where:
// a = l * l
// b = 2(l * (o - c))
// c = (o - c) * (o - c) - r^2 

// since l (direction) is a unit vector, a == 1 
// final formula:
// D = -(l * (o-c)) +- sqroot((l * (o-c))^2 - |o-c|^2 +r^2)
// If the value inside sqroot is 0 we have 1 hit (just touched the sphere once)
// If the value inside sqroot is negative, no solution exists, it's a miss. 
// If the value inside sqroot is postive, we have 2 hits and our ray travels through 
// and hits the sphere "edge" twice. 
// If both solutions for D is negative, the sphere is behind the ray. 

// In code:
// Distance = -(loc) +- sqroot((loc)^2 - ocLength^2  + r2)
// Distance = -(loc) +- sqroot(det)
// @Return: 0 if miss, distance to closest hit if hit. 
// Beware -- This function relies heavily on overloaded operators from Vector3d.
double Sphere::intersect(const Ray &ray) const
{
	Vector3d oc = position - ray.origin; 		

	double distance;
	double fudge = 1e-4;
	double b = oc.dot(ray.direction);							// 1/2 b for quadratic eq. setup
	double determinant = b * b - oc.dot(oc) + radius * radius; // b^2-4ac -- a = 1 because ray is normalized. 

	if(determinant < 0)
	{
		return 0;
	}
	else
	{
		determinant = sqrt(determinant);
	}

	distance = b - determinant;

	// Return lowest positive distance, or 0 for miss. 
	if(distance > fudge)
	{
		return distance;
	}
	else
	{
		distance = b + determinant;
		if(distance > fudge)
		{
			return distance;
		}
		else
		{
			return 0;
		}
	}
}

/*double Sphere::intersect(const Ray &ray)
{
	Vector3d oc = this -> position - ray.origin; 		
	double loc = ray.direction.dot(oc); 		    
	double ocLength = oc.length();			  			
	double r2 = this -> radius * this -> radius;
	double det = (loc * loc) - (ocLength * ocLength) + r2;

	if(det < 0)
	{
		return 0;
	}
	else
	{
		det = sqrt(det);
	}

	double distance1 = -loc + det;
	double distance2 = -loc - det;
	
	if(std::abs(distance1) > std::abs(distance2))
	{
		return distance2;
	}
	else
	{
		return distance1;
	}
}*/
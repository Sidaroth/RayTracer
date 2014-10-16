#include "Ray.h"

Ray::Ray(Vector3d origin, Vector3d direction)
{
	// Beware -- this uses the overloaded assignment operator from Vector3d. 
	this.origin = origin;
	this.direction = direction;
}
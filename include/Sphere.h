// Author: Christian Holt              
// Last Edit: 04.11.14                 
// Purpose: Definition of a sphere.    
///////////////////////////////////////

#pragma once

#include <cmath>
#include "Vector3d.h"
#include "Ray.h"
#include "globals.h"

class Sphere
{
public:
	double radius;
	Vector3d position;
	Vector3d emission;
	Vector3d color;
	MaterialType materialType;

	Sphere();
	Sphere(double radius, Vector3d position, Vector3d emission,
		   Vector3d color, MaterialType materialType);

	double intersect(const Ray &ray) const;
};
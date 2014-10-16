// Purpose: Definition of a Ray for the RayTracer 

#pragma once

#include "Vector3d.h"

class Ray
{
public:
	Vector3d origin;
	Vector3d direction;

	Ray(Vector3d origin, Vector3d direction);
};
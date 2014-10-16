// Author: Christian Holt              		 
// Last Edit: 16.10.14                 		 
// Purpose: A simple Ray Tracer        		 
//									   		 
// Some code snippets by smallPt by kevin beason.
// http://kevinbeason.com/smallpt/
// https://docs.google.com/file/d/0B8g97JkuSSBwUENiWTJXeGtTOHFmSm51UC01YWtCZw/edit


#include <cstdio>
#include "Sphere.h"
#include "Vector3d.h"
#include "globals.h"

/// Used for radiance calculation. 


// Scene from smallPT. 
 Sphere spheres[] = { //Scene: radius, position, emission, color, material 
   Sphere(1e5,  Vector3d( 1e5 + 1, 40.8, 81.6),   Vector3d(), Vector3d(0.75, 0.25, 0.25), MaterialType::DIFFUSE), // Left 
   Sphere(1e5,  Vector3d(-1e5 + 99, 40.8, 81.6),  Vector3d(), Vector3d(0.25, 0.25, 0.75), MaterialType::DIFFUSE), // Right 
   Sphere(1e5,  Vector3d(50, 40.8, 1e5),      	  Vector3d(), Vector3d(0.75, 0.75, 0.75), MaterialType::DIFFUSE), // Back 
   Sphere(1e5,  Vector3d(50, 40.8, -1e5 + 170),   Vector3d(), Vector3d(),           	  MaterialType::DIFFUSE), // Front 
   Sphere(1e5,  Vector3d(50, 1e5, 81.6),     	  Vector3d(), Vector3d(0.75, 0.75, 0.75), MaterialType::DIFFUSE), // Bottomm 
   Sphere(1e5,  Vector3d(50, -1e5 + 81.6, 81.6),  Vector3d(), Vector3d(0.75, 0.75, 0.75), MaterialType::DIFFUSE), // Top 
   Sphere(16.5, Vector3d(27, 16.5, 47),        	  Vector3d(), Vector3d(1, 1, 1) * 0.999,  MaterialType::SPECULAR), // Mirror 
   Sphere(16.5, Vector3d(73, 16.5, 78),        	  Vector3d(), Vector3d(1, 1, 1) * 0.999,  MaterialType::REFRACTIVE), // Glass 
   Sphere(600,  Vector3d(50, 681.6 - 0.27, 81.6), Vector3d(12, 12, 12),  Vector3d(), 	  MaterialType::DIFFUSE)  	 // Lite 
 }; 

int main()
{
	printf("Hello World, I am a RayTracer!");

	return 0; 
}
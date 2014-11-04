// Author: Christian Holt              		 
// Last Edit: 04.02.14                 		 
// Purpose: A basic raytracer. Aspects of this may be overly commented, but this is done as a learning exercise so the comments are equally for
//          the authors benefit. Having to put what happens into words furthers our understanding.   		
// 
// Keep in mind this program does not aim for efficiency. 
//									   		 
// Some code snippets from smallPt by kevin beason, and based on concepts from the following sites and more:
// http://kevinbeason.com/smallpt/
// https://docs.google.com/file/d/0B8g97JkuSSBwUENiWTJXeGtTOHFmSm51UC01YWtCZw/edit
// http://scratchapixel.com/lessons/3d-basic-lessons/lesson-1-writing-a-simple-raytracer/ 

#include <cstdio>
#include <iostream>
#include <omp.h>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <ctime>
#include <string>
#include <fstream>
#include <random>
#include "Sphere.h"
#include "Vector3d.h"
#include "globals.h"                                                                                  

/// CONSTANTS
const double      INFINITYAPP = 1e20;                                                                         // Our approximation of infinity. 
const int         SAMPLES     = 50;                                                                           // Number of samples per subpixel
const int         WIDTH       = 1024;                                                                          // Resolution of the final image
const int         HEIGHT      = 768;
const double      FOVANGLE    = 0.5135;
const int         XSUBSAMPLES = 2;                                                                         // Number of samples to split a pixel into. 
const int         YSUBSAMPLES = 2;
const std::string FILENAME    = "image.ppm";
const long long   SEED        = 31415926535;

/// Function decl. 
bool     intersectCheck(Ray&, double&, double&);
Vector3d radiance(const Ray&, int, int E=1);
double   clamp(double);
void     writeToFile();
double   applyGamma(double);
int      convertToIntegerRange(double);


// Scene from smallPT. 
Sphere spheres[] = { // Sphere: radius, position, emission, color, material 
  Sphere(1e5,  Vector3d( 1e5 + 1, 40.8, 81.6),   Vector3d(),            Vector3d(0.75, 0.25, 0.25), MaterialType::DIFFUSE),    // Left 
  Sphere(1e5,  Vector3d(-1e5 + 99, 40.8, 81.6),  Vector3d(),            Vector3d(0.25, 0.25, 0.75), MaterialType::DIFFUSE),    // Right 
  Sphere(1e5,  Vector3d(50, 40.8, 1e5),          Vector3d(),            Vector3d(0.75, 0.75, 0.75), MaterialType::DIFFUSE),    // Back 
  Sphere(1e5,  Vector3d(50, 40.8, -1e5 + 170),   Vector3d(),            Vector3d(),                 MaterialType::DIFFUSE),    // Front 
  Sphere(1e5,  Vector3d(50, 1e5, 81.6),          Vector3d(),            Vector3d(0.75, 0.75, 0.75), MaterialType::DIFFUSE),    // Bottomm 
  Sphere(1e5,  Vector3d(50, -1e5 + 81.6, 81.6),  Vector3d(),            Vector3d(0.75, 0.75, 0.75), MaterialType::DIFFUSE),    // Top 
  Sphere(16.5, Vector3d(27, 16.5, 47),           Vector3d(),            Vector3d(1, 1, 1) * 0.999,  MaterialType::SPECULAR),   // Mirror 
  Sphere(16.5, Vector3d(73, 16.5, 78),           Vector3d(),            Vector3d(1, 1, 1) * 0.999,  MaterialType::REFRACTIVE), // Glass 
  Sphere(600,  Vector3d(50, 681.6 - 0.27, 81.6), Vector3d(12, 12, 12),  Vector3d(),                 MaterialType::DIFFUSE)     // Light 
}; 

int sphereLength = (int) (sizeof(spheres) / sizeof(Sphere));
Vector3d *image = new Vector3d[WIDTH*HEIGHT];                                            // An array holding (r,g,b) values for each pixel. 
std::mt19937 mtGenerator(SEED);
std::uniform_real_distribution<> distribution(0, 1);                                     // Uniform distribution of real numbers between 0 and 1 using a mersenne twister as generator. 


int main()
{
  /////// SETUP ////////////
  printf("\nBeginning Ray Trace.");
  std::chrono::time_point<std::chrono::system_clock> start, now, end;
  std::chrono::duration<double> elapsedSeconds;
  start = std::chrono::system_clock::now();


  Ray camera(Vector3d(50, 52, 295.6), Vector3d(0, -0.042612, -1).unit());                  // Camera position and direction. Also taken from SmallPT to match the scene. 
  Vector3d cx = Vector3d(WIDTH * FOVANGLE / HEIGHT, 0, 0);                                 // X-direction increment. (assumes a camera that is upright) - X camera direction. 
  Vector3d cy = cx.cross(camera.direction).unit() * FOVANGLE;                              // Y-direction increment. (cross product gets vector perpendicular to CX and gaze direction - the vertical up vector)

  Vector3d color;                                                                          // Color of a single sample. 

  ///////// CREATE THE IMAGE ///////////////
#pragma omp parallel for schedule(dynamic, 1) private(color)                               // Each loop iteration should run in its own thread. 
  for(int y = 0; y < HEIGHT; ++y)                                                          // Image rows. 
  {
    now = std::chrono::system_clock::now();
    elapsedSeconds = now - start;
    printf("\nRendering (%d samples): %5.2f%%. Elapsed time: %fs", SAMPLES * 4, 100.0 * y / (HEIGHT - 1), elapsedSeconds.count());    // Print progress

    for(unsigned short x = 0; x < WIDTH; ++x)                                              // Image columns. 
    {
      // XSUBSAMPLESxYSUBSAMPLES SUBSAMPLES for each pixel. Doing SAMPLES samples for each SUBSAMPLE
      for(int subY = 0, pixel = (HEIGHT - y - 1) * WIDTH + x; subY < YSUBSAMPLES; ++subY)  // Y SUBSAMPLES, 
      {
        for(int subX = 0; subX < XSUBSAMPLES; ++subX, color = Vector3d())                  // X SUBSAMPLES
        {
          for(int sample = 0; sample < SAMPLES; ++sample)                                  // Sampling. 
          {
            // TENT FILTER (r1 and r2 are random values) 
            double r1 = 2 * distribution(mtGenerator);
            double dx = r1 < 1 ? (sqrt(r1) - 1) : 1 - sqrt(2 - r1);

            double r2 = 2 * distribution(mtGenerator);
            double dy = r2 < 1 ? (sqrt(r2) - 1) : 1 - sqrt(2 - r2);

            ///////////// Determining the color ///////////////////////////
            double xDir = ((subX + 0.5 + dx) / 2 + x) / WIDTH  - 0.5;                         // Adjusting ray direction.  
            double yDir = ((subY + 0.5 + dy) / 2 + y) / HEIGHT - 0.5;
            Vector3d direction = (cx * xDir) + (cy * yDir) + camera.direction;

            color = color + radiance(Ray(camera.origin + direction * 140, direction.unit()), 0) * (1.0 / SAMPLES);    // rays are pushed forward (140) to start inside the room. 
          }

          image[pixel] = image[pixel] + Vector3d(clamp(color.x), clamp(color.y), clamp(color.z)) * 0.25;
        }
      }
    }
  }

  writeToFile();
  end = std::chrono::system_clock::now();
  elapsedSeconds = end - start;

  printf("\n\nSucessfully wrote to file (%s). Total time taken to trace: %f seconds.", FILENAME.c_str(),elapsedSeconds.count());
	
  return 0;
}

///////// FUNCTIONS /////////

/// Format: http://netpbm.sourceforge.net/doc/ppm.html#plainppm 
void writeToFile()
{
  std::ofstream file;
  file.open(FILENAME.c_str());
  file << "P3\n" << WIDTH << " " << HEIGHT << "\n255\n";

  for(int pixel = 0; pixel < WIDTH * HEIGHT; ++pixel)
  {
     int r = convertToIntegerRange(applyGamma(image[pixel].x));
     int g = convertToIntegerRange(applyGamma(image[pixel].y));
     int b = convertToIntegerRange(applyGamma(image[pixel].z));
  
     file << r << " " << g << " " << b << " ";
  }

  file.close();
}

/// Applies a gamma correction constant. 
/// A value of 2.2 is normal: http://en.wikipedia.org/wiki/Gamma_correction 
double applyGamma(double value)
{
  return pow(clamp(value), 1 / 2.2);
}

/// Takes a floating point value in the range [0.0, 1.0] and returns an integer in the range [0, 255]
int convertToIntegerRange(double value)
{
  return int(value * 255 + 0.5);
}


bool intersectCheck(Ray &ray, double &temporary, double &id)
{
  temporary = INFINITYAPP;                                       // Used to store the current closest ray hit. 
  double distance;

  for(int i = sphereLength; i >= 0; --i)                      // Loop through all spheres
  {
    distance = spheres[i].intersect(ray);

    if(distance > 0 && distance < temporary)                  // We have a ray hit.  -- Check for "temporary" to see if we intersecT with a closer sphere (i.e in front)
    {
      id = i;
      temporary = distance;
    }
  } 

  return (temporary < INFINITYAPP);
}

/// Return a value in the range [0.0, 1.0]
double clamp(double value)
{
  if(value < 0)
  {
    return 0.0;
  }
  else if(value > 1)
  {
    return 1.0;
  }
  else
  {
    return value;
  }
}

Vector3d radiance(const Ray &ray, int depth, int E)
{



  return Vector3d();
}
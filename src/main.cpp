// Author: Christian Holt              		 
// Last Edit: 04.11.14                 		 
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
#include <cassert>
#include <random>
#include "Sphere.h"
#include "Vector3d.h"
#include "globals.h"                                                                                  

/// CONSTANTS
const double      INFINITYAPP     = 1e20;                                                                         // Our approximation of infinity. 
const int         SAMPLES         = 100;                                                                           // Number of samples per subpixel
const int         WIDTH           = 512;                                                                          // Resolution of the final image
const int         HEIGHT          = 384;
const double      FOVANGLE        = 0.5135;
const int         XSUBSAMPLES     = 2;                                                                         // Number of samples to split a pixel into. 
const int         YSUBSAMPLES     = 2;
const std::string FILENAME        = "image.ppm";
const long long   SEED            = 314159265;
const int         RECURSIVE_LIMIT = 10;
const int         ROULETTE_LIMIT  = 5;
const double      M_PI            = 3.14159265358979;
const double      M_1_PI          = 1.0 / M_PI;
const int         THREADS         = 12;

/// Function decl. 
bool     intersectCheck(const Ray&, double&, int&);
Vector3d radiance(const Ray&, int, int E=1);
double   clamp(double);
void     writeToFile();
double   applyGamma(double);
int      convertToIntegerRange(double);
int      toInt(double);


// Scene from smallPT. 
Vector3d tc(0.0588, 0.361, 0.0941);
Vector3d sc = Vector3d(1,1,1)*.7;
Sphere spheres[] = {//Scene: radius, position, emission, color, material
  // center 50 40.8 62
  // floor 0
  // back  0
//  Sphere(1e5, Vector3d(50, 1e5+100, 0),  Vector3d(1,1,1)*1,Vector3d(),MaterialType::DIFFUSE), //lite
//  Sphere(1e5, Vector3d(50, -1e5, 0),  Vector3d(),Vector3d(.3,.3,.1),MaterialType::DIFFUSE), //grnd
//  Sphere(1e5, Vector3d(50, 1e5+100, 0),  Vector3d(0.761, 0.875, 1.00)*1.3,Vector3d(),MaterialType::DIFFUSE),
//  //lite
  Sphere(1e5, Vector3d(50, 1e5+130, 0),  Vector3d(1,1,1)*1.3,Vector3d(),MaterialType::DIFFUSE), //lite
  Sphere(1e2, Vector3d(50, -1e2+2, 47),  Vector3d(),Vector3d(1,1,1)*.7,MaterialType::DIFFUSE), //grnd

  Sphere(1e4, Vector3d(50, -30, 300)+Vector3d(-sin(50*M_PI/180),0,cos(50*M_PI/180))*1e4, Vector3d(), Vector3d(1,1,1)*.99,MaterialType::SPECULAR),// mirr L
  Sphere(1e4, Vector3d(50, -30, 300)+Vector3d(sin(50*M_PI/180),0,cos(50*M_PI/180))*1e4,  Vector3d(), Vector3d(1,1,1)*.99,MaterialType::SPECULAR),// mirr R
  Sphere(1e4, Vector3d(50, -30, -50)+Vector3d(-sin(30*M_PI/180),0,-cos(30*M_PI/180))*1e4,Vector3d(), Vector3d(1,1,1)*.99,MaterialType::SPECULAR),// mirr FL
  Sphere(1e4, Vector3d(50, -30, -50)+Vector3d(sin(30*M_PI/180),0,-cos(30*M_PI/180))*1e4, Vector3d(), Vector3d(1,1,1)*.99,MaterialType::SPECULAR),// mirr


  Sphere(4, Vector3d(50,6*.6,47),   Vector3d(),Vector3d(.13,.066,.033), MaterialType::DIFFUSE),//"tree"
  Sphere(16,Vector3d(50,6*2+16*.6,47),   Vector3d(), tc,  MaterialType::DIFFUSE),//"tree"
  Sphere(11,Vector3d(50,6*2+16*.6*2+11*.6,47),   Vector3d(), tc,  MaterialType::DIFFUSE),//"tree"
  Sphere(7, Vector3d(50,6*2+16*.6*2+11*.6*2+7*.6,47),   Vector3d(), tc,  MaterialType::DIFFUSE),//"tree"

  Sphere(15.5,Vector3d(50,1.8+6*2+16*.6,47),   Vector3d(), sc,  MaterialType::DIFFUSE),//"tree"
  Sphere(10.5,Vector3d(50,1.8+6*2+16*.6*2+11*.6,47),   Vector3d(), sc,  MaterialType::DIFFUSE),//"tree"
  Sphere(6.5, Vector3d(50,1.8+6*2+16*.6*2+11*.6*2+7*.6,47),   Vector3d(), sc,  MaterialType::DIFFUSE),//"tree"
};

/*Sphere spheres[] = { // Sphere: radius, position, emission, color, material 
  Sphere(1e5,  Vector3d( 1e5 + 1, 40.8, 81.6),   Vector3d(),                  Vector3d(0.75, 0.25, 0.25), MaterialType::DIFFUSE),    // Left 
  Sphere(1e5,  Vector3d(-1e5 + 99, 40.8, 81.6),  Vector3d(),                  Vector3d(0.25, 0.25, 0.75), MaterialType::DIFFUSE),    // Right 
  Sphere(1e5,  Vector3d(50, 40.8, 1e5),          Vector3d(),                  Vector3d(0.75, 0.75, 0.75), MaterialType::DIFFUSE),    // Back 
  Sphere(1e5,  Vector3d(50, 40.8, -1e5 + 170),   Vector3d(),                  Vector3d(),                 MaterialType::DIFFUSE),    // Front 
  Sphere(1e5,  Vector3d(50, 1e5, 81.6),          Vector3d(),                  Vector3d(0.75, 0.75, 0.75), MaterialType::DIFFUSE),    // Bottom 
  Sphere(1e5,  Vector3d(50, -1e5 + 81.6, 81.6),  Vector3d(),                  Vector3d(0.75, 0.75, 0.75), MaterialType::DIFFUSE),    // Top 
  Sphere(16.5, Vector3d(27, 16.5, 47),           Vector3d(),                  Vector3d(1, 1, 1) * 0.999,  MaterialType::SPECULAR),   // Mirror 
  Sphere(16.5, Vector3d(73, 16.5, 78),           Vector3d(),                  Vector3d(1, 1, 1) * 0.999,  MaterialType::REFRACTIVE), // Glass 
  Sphere(1.5,  Vector3d(50, 81.6 - 16.5, 81.6),  Vector3d(4, 4, 4) * 100,     Vector3d(),                 MaterialType::DIFFUSE)     // Light 
}; */

int sphereLength = (int) (sizeof(spheres) / sizeof(Sphere));
Vector3d *image = new Vector3d[WIDTH*HEIGHT];                                            // An array holding (r,g,b) values for each pixel. 
std::mt19937 mtGenerator(SEED);
std::uniform_real_distribution<> distribution(0, 1);                                     // Uniform distribution of real numbers between 0 and 1 using a mersenne twister as generator. 


int main()
{
  /////// SETUP ////////////
  printf("\nBeginning Ray Trace.\n");
  std::chrono::time_point<std::chrono::system_clock> start, now, end;
  std::chrono::duration<double> elapsedSeconds;
  start = std::chrono::system_clock::now();

  omp_lock_t writeLock;
  omp_set_num_threads(THREADS);
  omp_init_lock(&writeLock);

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
    omp_set_lock(&writeLock);
    printf("\rRendering (%d samples): %5.2f%%. Elapsed time: %fs", SAMPLES * 4, 100.0 * y / (HEIGHT - 1), elapsedSeconds.count());    // Print progress
    omp_unset_lock(&writeLock);

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

  printf("\n\nSucessfully wrote to file (%s). Total time taken to trace: %f seconds.\n", FILENAME.c_str(), elapsedSeconds.count());
  omp_destroy_lock(&writeLock);
	
  return 0;
}

///////// FUNCTIONS /////////

/// Radiance method needs some cleanup, should probably seperate out the different material types etc. (i.e a switch with function calls)
/// Any reference to Vector3d(); is a Vector with values (0, 0, 0) for RGB, i.e black. 
Vector3d radiance(const Ray &ray, int depth, int E)
{
  double distance;                        // Distance to intersection
  int id = 0; 
  if(!intersectCheck(ray, distance, id))
  {
    return Vector3d();                    // We missed. 
  }

  const Sphere &sphere = spheres[id];     // The sphere hit. 

  if(depth > RECURSIVE_LIMIT)                          // Recursion depth check
  {
    return Vector3d();
  }

  Vector3d intersectPoint = ray.origin + ray.direction * distance;
  Vector3d normal         = (intersectPoint - sphere.position).unit();
  Vector3d surfaceNormal  = normal.dot(ray.direction) < 0 ? normal : normal * -1; 
  Vector3d BRDFModulator  = sphere.color;

  // Use maximum reflectance for Russian Roulette. 
  double maxRefl;
  if(BRDFModulator.x > BRDFModulator.y && BRDFModulator.x > BRDFModulator.z)
  {
    maxRefl = BRDFModulator.x;
  }
  else
  {
    if(BRDFModulator.y > BRDFModulator.z)
    {
      maxRefl = BRDFModulator.y;
    }
    else
    {
      maxRefl = BRDFModulator.z;
    }
  }

  if(++depth > ROULETTE_LIMIT || !maxRefl)                                 // Only do russian roulette after a recursion depth of ROULETTE_LIMIT. 
  {
    if(distribution(mtGenerator) < maxRefl)
    {
      BRDFModulator = BRDFModulator * (1 / maxRefl);
    }
    else
    {
      return sphere.emission * E;
    }
  }

  // IDEAL DIFFUSE REFLECTION
  if(sphere.materialType == MaterialType::DIFFUSE)
  {
    double r1 = 2 * M_PI * distribution(mtGenerator);     // Angle around randomly
    double r2 = distribution(mtGenerator);                // distance from
    double r2s = sqrt(r2);                                // center - Random

    Vector3d w = surfaceNormal;                           // Vectors w, u, v are used to create an orthonormal coordinate frame
    Vector3d u;                                           // around the normal to sample the unit hemisphere. (x is perpendicular to z)

    if(fabs(w.x) > 0.1)
    {
      u = Vector3d(0, 1, 0);
    }
    else
    {
      u = Vector3d(1, 0, 0);
    }
    u = u.cross(w).unit();
    Vector3d v = w.cross(u);                              // Perpendicular to x and z. 

    Vector3d reflectanceRayDir = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).unit(); // Create a random(based on r1, r2) reflectance ray within the hemisphere. 
    
    // Sampling Lights
    Vector3d lighting;
    for(int i = 0; i < sphereLength; ++i)
    {
      const Sphere &light = spheres[i];
      if(light.emission.x <= 0 && light.emission.y <= 0 && light.emission.z <= 0)
      {
        continue;                                                 //// Skip any spheres that are not light sources. 
      }

      // Create another coordinate system sw, su, sv similar to x,y,z above.
      Vector3d sw = light.position - intersectPoint;
      Vector3d su;

      if(fabs(sw.x) > 0.1)
      {
        su = Vector3d(0, 1, 0);
      }
      else
      {
        su = Vector3d(1, 0, 0);
      }
      su = su.cross(sw).unit();
      Vector3d sv = sw.cross(su);

      double cosMax = sqrt(1 - light.radius * light.radius / (intersectPoint - light.position).dot(intersectPoint - light.position)); // Determine max angle
      double eps1 = distribution(mtGenerator);
      double eps2 = distribution(mtGenerator);            // Calculating sample direction. (This is based on an equation from Realistic Ray Tracing by Shirley et.al)
      double cosA = 1 - eps1 + eps1 * cosMax;
      double sinA = sqrt(1 - cosA * cosA);
      double phi = 2 * M_PI * eps2;

      Vector3d sampleDir = su * cos(phi) * sinA + sv * sin(phi) * sinA + sw * cosA;
      sampleDir = sampleDir.unit();

      // Shoot shadow ray
      if(intersectCheck(Ray(intersectPoint, sampleDir), distance, id) && id == i)
      {
        double omega = 2 * M_PI * (1 - cosMax);
        lighting = lighting + BRDFModulator * (light.emission * sampleDir.dot(surfaceNormal) * omega) * M_1_PI; // 1/pi for BRDF
      }

      return sphere.emission * E + lighting + BRDFModulator * (radiance(Ray(intersectPoint, reflectanceRayDir), depth, 0));
    }
  }
  else if(sphere.materialType == MaterialType::SPECULAR)  // IDEAL SPECULAR REFLECTION
  {
    // Angle of incidence == Angle of Reflection, R = D - 2(N dot D)N. 
    Vector3d recursive = radiance(Ray(intersectPoint, ray.direction - normal * 2 * normal.dot(ray.direction)), depth);
    return sphere.emission + BRDFModulator * recursive;
  }
  else                                                    // DIELECTRIC SURFACE. (Glass)
  {
    // Again, angle of incidence = angle of reflection -- IDEAL
    Ray reflectanceRay(intersectPoint, ray.direction - normal * 2 * normal.dot(ray.direction));
    bool into = normal.dot(surfaceNormal) > 0;          // entering or exiting glass

    double ior = 1.5;                                    // Index of Refraction (ior) for glass is 1.5, nnt is either 1.5 or 1/1.5 depending on entry or exit. 
    double nnt = into ? 1 / ior : ior / 1;
    double ddn = ray.direction.dot(surfaceNormal);      
    double angle = 1 - nnt * nnt * (1 - ddn * ddn);     // angle of exit from the glass. 

    // Total Internal Reflection
    if(angle < 0)                                       // if the angle is too shallow, all light is reflected. 
    {
      return sphere.emission + BRDFModulator * (radiance(reflectanceRay, depth));
    }

    // Otherwise, choose reflection or refraction with Fresnel's Term. 
    Vector3d transmissionDir = (ray.direction * nnt - normal * ((into ? 1 : -1) * (ddn * nnt + sqrt(angle)))).unit();
    double a  = ior - 1;
    double b  = ior + 1;
    double R0 = a * a / (b * b);                                     // Reflectance at normal incidence based on IOR. 
    double c  = 1 - (into ? -ddn : transmissionDir.dot(normal));     // 1 - cos(theta)
    double Re = R0 + (1 - R0) * c * c * c * c * c;                   // Fresnel Reflectance. 
    double Tr = 1 - Re;
    double P  = 0.25 + (0.5 * Re);                                     // Probability of reflecting
    double Rp = Re / P;
    double Tp = Tr / (1 - P);

    if(depth > 2)
    {
      if(distribution(mtGenerator) < P)
      {
        return sphere.emission + BRDFModulator * (radiance(reflectanceRay, depth) * Rp);                              // Reflection
      }
      else
      {
        return sphere.emission + BRDFModulator * (radiance(Ray(intersectPoint, transmissionDir), depth) * Tp);        // Refraction
      }
    }
    else
    {
      return sphere.emission + BRDFModulator * (radiance(reflectanceRay, depth) * Re + radiance(Ray(intersectPoint, transmissionDir), depth) * Tr); 
    }

  }

  return Vector3d();
}

/// Format: http://netpbm.sourceforge.net/doc/ppm.html#plainppm 
void writeToFile()
{
  std::ofstream file;
  file.open(FILENAME.c_str());
  file << "P3\n" << WIDTH << " " << HEIGHT << "\n255\n";

  for(int pixel = 0; pixel < WIDTH * HEIGHT; ++pixel)
  {
    assert(toInt(image[pixel].x) == convertToIntegerRange(applyGamma(image[pixel].x)));
    file << toInt(image[pixel].x) << " " << toInt(image[pixel].y) << " " << toInt(image[pixel].z) << " ";
   //  int r = convertToIntegerRange(applyGamma(image[pixel].x));
   //  int g = convertToIntegerRange(applyGamma(image[pixel].y));
   //  int b = convertToIntegerRange(applyGamma(image[pixel].z));
  
   //  file << r << " " << g << " " << b << " ";
  }

  file.close();
}

// Comparing with my own 2 seperate functions
int toInt(double x)
{
  return int(pow(clamp(x), 1 / 2.2) * 255 + .5);
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


bool intersectCheck(const Ray &ray, double &temporary, int &id)
{
  temporary = INFINITYAPP;                                       // Used to store the current closest ray hit. 
  double numSpheres = sizeof(spheres) / sizeof(Sphere);
  double distance;

  for(int i = int(numSpheres); i--;)                      // Loop through all spheres
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
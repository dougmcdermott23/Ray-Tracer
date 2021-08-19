/*
  CSC418 - RayTracer code - Winter 2017 - Assignment 3&4

  RayTracer header file. Structure and data definitions.
  F.J. Estrada, Dec. 9, 2010

  Goal is to make this code as simple and fast as possible...

  You DO NOT need to modify this file, but you DO NEED to
  read and understand what it does.
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<time.h>

#ifndef __RayTracer_header
#define __RayTracer_header

#define PI 3.14159265354
/* The structure below is used to hold a single RGB image */
struct image{
	void *rgbdata;
	int sx;
	int sy;
};

/* The structure below defines a point in 3D homogeneous coordinates */
struct point3D{
	double px;
	double py;
	double pz;
	double pw;
};

/* The structure below defines a ray in 3D homogeneous coordinates,
   the point corresponds to the representation r(\lambda)=p+(\lambda)d */
struct ray3D{
	struct point3D p0;	// Ray origin (at lambda=0)
	struct point3D d;		// Ray direction
	void (*rayPos)(struct ray3D *ray, double lambda, struct point3D *pos);
					// Function to return the
					// position along the ray
					// for a given lambda.
};

/*
   The structures below are used to define an object colour in terms of the
   components of the Phong illumination model. Note that we typically
   define colours together with objects, so you should not need to
   instantiate lone instances of the colour structure.

   Also, note that you can easily make your objects completely white
   (or completely monochromatic) by not being careful how the different
   components in the Phong model add up. Take a moment and think how
   you want your object to look before you set these values.
*/
struct albedosPhong{
	double ra;	// Ambient light albedo
	double rd;	// Diffuse component albedo
	double rs;	// Specular component albedo
	double rg;	// Global component albedo
};

/*
   The structure below defines an RGB colour, values are
   in [0,1]
*/
struct colourRGB{
	double R;
	double G;
	double B;
};

/*
   The structure below defines an Object within the World Coordinate Frame.
   For this ray tracer, we will use a simple linked list of objects (not
   a tree structure, so hierarchical transformations have to be fully
   specified). This is not really limiting and simplifies the data
   structure to hold the scene.

   Note that the Object3D structure is completely agnostic about the type
   of object defined, the difference is made by the intersect function.
   If you think carefully about it, in terms of raytracing objects are
   defined entirely by their surface, and the surface is defined in the
   intersect function for each object.

   Thus, to create additional objects, simply provide a suitable
   intersection function. The rest stays the same.
*/
struct object3D{
	struct albedosPhong alb;	// Object's albedos for Phong model
	struct colourRGB col;	// Object's colour in RGB
	double  T[4][4]; 	// T holds the transformation applied to this object.
	double  Tinv[4][4];      // Tinv holds the inverse transformation

        // Below we set up space for a pointer to the intersection function for this object.
        // Note that the intersection function must compute the lambda at the intersection, the
        // intersection point p, the normal at that point n, and the texture coordinates (a,b).
        // The texture coordinates are not used unless texImg!=NULL and a textureMap function
        // has been provided
	void (*intersect)(struct object3D *obj, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b);		

	// Texture mapping function. Takes normalized texture coordinates (a,b) and returns the
  	// texture colour at that point using bi-linear interpolation 
	void (*textureMap)(struct image *img, double a, double b, double *R, double *G, double *B);

        struct image *texImg;				// Pointer to structure
							// holding the texture
							// for this object
	double  alpha;		// Opacity - if less than 1 this is a semi
				// transparent object and refraction rays
				// should be implemented
	double  r_index;		// Index of refraction
        double  shinyness;	// Exponent for the Phong specular component
	int 	frontAndBack;	// Flag to indicate that both sides of the object
				// should be lit.
	int	isLightSource;	// Flag to indicate if this is an area light source
	struct object3D *next;	// Pointer to next entry in object linked list
};

/* The structure below defines a point light source */
struct pointLS{
	struct colourRGB col;	// Light source colour
	struct point3D p0;	// Light source location
	struct pointLS *next;	// Pointer to next light in the scene
};

/*
   The structure below is used to hold camera parameters. You will need
   to write code to initialize the camera position and orientation.
*/
struct view{
	struct point3D	e;	// Location of the camera center
	struct point3D	u;	// u vector
	struct point3D  v;	// v vector
	struct point3D  w;	// w vector
	double f;		// Focal length
	double wl;		// Left edge in camera coordinates
	double wt;		// Top edge in camera coordinates
	double wsize;		// Window size in distance units (not pixels!)
	double W2C[4][4];	// World2Camera conversion matrix
	double C2W[4][4];	// Camera2World conversion matrix
};

// Function definitions start here
int main(int argc, char *argv[]);									// Main raytracing function. 
void buildScene(void);											// Scene set up. Defines objects and object transformations
void rayTrace(struct ray3D *ray, int depth, struct colourRGB *col, struct object3D *Os);		// RayTracing routine
void findFirstHit(struct ray3D *ray, double *lambda, struct object3D *Os, struct object3D **obj, struct point3D *p, struct point3D *n, double *a, double *b);
void rtShade(struct object3D *obj, struct point3D *p, struct point3D *n,struct ray3D *ray, int depth, double a, double b, struct colourRGB *col);

#endif

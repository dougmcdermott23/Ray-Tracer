/*
   utils.c - F.J. Estrada, Dec. 9, 2010

   Utilities for the ray tracer. You will need to complete
   some of the functions in this file. Look for the sections
   marked "TO DO". Be sure to read the rest of the file and
   understand how the entire code works.
*/

#include "utils.h"

// A useful 4x4 identity matrix which can be used at any point to
// initialize or reset object transformations
double eye4x4[4][4]={{1.0, 0.0, 0.0, 0.0},
                    {0.0, 1.0, 0.0, 0.0},
                    {0.0, 0.0, 1.0, 0.0},
                    {0.0, 0.0, 0.0, 1.0}};

/////////////////////////////////////////////
// Primitive data structure section
/////////////////////////////////////////////
		    
void newPoint(double px, double py, double pz, struct point3D *p0) {
  p0->px = px;
  p0->py = py;
  p0->pz = pz;
  p0->pw = 1.0;
}

struct pointLS *newPLS(struct point3D *p0, double r, double g, double b)
{
 // Allocate a new point light sourse structure. Initialize the light
 // source to the specified RGB colour
 // Note that this is a point light source in that it is a single point
 // in space, if you also want a uniform direction for light over the
 // scene (a so-called directional light) you need to place the
 // light source really far away.

 struct pointLS *ls=(struct pointLS *)calloc(1,sizeof(struct pointLS));
 if (!ls) fprintf(stderr,"Out of memory allocating light source!\n");
 else
 {
  memcpy(&ls->p0,p0,sizeof(struct point3D));	// Copy light source location

  ls->col.R=r;					// Store light source colour and
  ls->col.G=g;					// intensity
  ls->col.B=b;
 }
 return(ls);
}

/////////////////////////////////////////////
// Ray and normal transforms
/////////////////////////////////////////////
inline void rayTransform(struct ray3D *ray_orig, struct ray3D *ray_transformed, struct object3D *obj)
{
 // Transforms a ray using the inverse transform for the specified object. This is so that we can
 // use the intersection test for the canonical object. Note that this has to be done carefully!
  
 //Transforming the origin for the ray to model space
 ray_transformed->p0.px = ray_orig->p0.px;
 ray_transformed->p0.py = ray_orig->p0.py;
 ray_transformed->p0.pz = ray_orig->p0.pz;
 ray_transformed->p0.pw = ray_orig->p0.pw;
 
 matVecMult(obj->Tinv, &ray_transformed->p0);
 
 //Transforming the direction vector for the ray to model space
 ray_transformed->d.px = ray_orig->d.px;
 ray_transformed->d.py = ray_orig->d.py;
 ray_transformed->d.pz = ray_orig->d.pz;
 ray_transformed->d.pw = ray_orig->d.pw;
 
 //apply inverse transform to the ray 
 matVecMult(obj->Tinv, &ray_transformed->d);
}

//Returns the transpose of the inverse matrix of an object
void Transpose (double transpose[][4], struct object3D *obj) {
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j ++) {
      transpose[i][j] = obj->Tinv[j][i];
    }
  }
  
  //set the bottom row to make sure w-component is not affected
  //we do not care about the translation component anyway 
  transpose[3][0] = 0; 
  transpose[3][1] = 0;
  transpose[3][2] = 0;
  transpose[3][3] = 1;
}

inline void normalTransform(struct point3D *n_orig, struct point3D *n_transformed, struct object3D *obj)
{
 // Computes the normal at an affinely transformed point given the original normal and the
 // object's inverse transformation. From the notes:
 // n_transformed=A^-T*n normalized.
 
 //Finding the transpose of the inverse transformation of an object
 double InverseTranspose[4][4];
 Transpose(InverseTranspose, obj);
 
 //store values of n_orig in a n_transformed
 n_transformed->px = n_orig->px;
 n_transformed->py = n_orig->py;
 n_transformed->pz = n_orig->pz;
 n_transformed->pw = 0;
 
 //calculate the normal and store the values in n_transformed
 matVecMult(InverseTranspose, n_transformed);
}

/////////////////////////////////////////////
// Object management section
/////////////////////////////////////////////
struct object3D *newPlane(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double r_index, double shiny)
{
 // Intialize a new plane with the specified parameters:
 // ra, rd, rs, rg - Albedos for the components of the Phong model
 // r, g, b, - Colour for this plane
 // alpha - Transparency, must be set to 1 unless you are doing refraction
 // r_index - Refraction index if you are doing refraction.
 // shiny - Exponent for the specular component of the Phong model
 //
 // The plane is defined by the following vertices (CCW)
 // (1,1,0), (-1,1,0), (-1,-1,0), (1,-1,0)
 // With normal vector (0,0,1) (i.e. parallel to the XY plane)

 struct object3D *plane=(struct object3D *)calloc(1,sizeof(struct object3D));

 if (!plane) fprintf(stderr,"Unable to allocate new plane, out of memory!\n");
 else
 {
  plane->alb.ra=ra;
  plane->alb.rd=rd;
  plane->alb.rs=rs;
  plane->alb.rg=rg;
  plane->col.R=r;
  plane->col.G=g;
  plane->col.B=b;
  plane->alpha=alpha;
  plane->r_index=r_index;
  plane->shinyness=shiny;
  plane->intersect=&planeIntersect;
  plane->texImg=NULL;
  memcpy(&plane->T[0][0],&eye4x4[0][0],16*sizeof(double));
  memcpy(&plane->Tinv[0][0],&eye4x4[0][0],16*sizeof(double));
  plane->textureMap=&texMap;
  plane->frontAndBack=1;
  plane->isLightSource=0;
 }
 return(plane);
}

struct object3D *newSphere(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double r_index, double shiny)
{
 // Intialize a new sphere with the specified parameters:
 // ra, rd, rs, rg - Albedos for the components of the Phong model
 // r, g, b, - Colour for this plane
 // alpha - Transparency, must be set to 1 unless you are doing refraction
 // r_index - Refraction index if you are doing refraction.
 // shiny -Exponent for the specular component of the Phong model
 //
 // This is assumed to represent a unit sphere centered at the origin.
 //

 struct object3D *sphere=(struct object3D *)calloc(1,sizeof(struct object3D));

 if (!sphere) fprintf(stderr,"Unable to allocate new sphere, out of memory!\n");
 else
 {
  sphere->alb.ra=ra;
  sphere->alb.rd=rd;
  sphere->alb.rs=rs;
  sphere->alb.rg=rg;
  sphere->col.R=r;
  sphere->col.G=g;
  sphere->col.B=b;
  sphere->alpha=alpha;
  sphere->r_index=r_index;
  sphere->shinyness=shiny;
  sphere->intersect=&sphereIntersect;
  sphere->texImg=NULL;
  memcpy(&sphere->T[0][0],&eye4x4[0][0],16*sizeof(double));
  memcpy(&sphere->Tinv[0][0],&eye4x4[0][0],16*sizeof(double));
  sphere->textureMap=&texMap;
  sphere->frontAndBack=0;
  sphere->isLightSource=0;
 }
 return(sphere);
}

struct object3D *newCylinder(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double r_index, double shiny)
{
 // Intialize a new cylinder with the specified parameters:
 // ra, rd, rs, rg - Albedos for the components of the Phong model
 // r, g, b, - Colour for this plane
 // alpha - Transparency, must be set to 1 unless you are doing refraction
 // r_index - Refraction index if you are doing refraction.
 // shiny -Exponent for the specular component of the Phong model
 //
 // This is assumed to represent a cylinder centered at the origin from y=-1 to y=1.

 struct object3D *cylinder=(struct object3D *)calloc(1,sizeof(struct object3D));

 if (!cylinder) fprintf(stderr,"Unable to allocate new sphere, out of memory!\n");
 else
 {
  cylinder->alb.ra=ra;
  cylinder->alb.rd=rd;
  cylinder->alb.rs=rs;
  cylinder->alb.rg=rg;
  cylinder->col.R=r;
  cylinder->col.G=g;
  cylinder->col.B=b;
  cylinder->alpha=alpha;
  cylinder->r_index=r_index;
  cylinder->shinyness=shiny;
  cylinder->intersect=&cylinderIntersect;
  cylinder->texImg=NULL;
  memcpy(&cylinder->T[0][0],&eye4x4[0][0],16*sizeof(double));
  memcpy(&cylinder->Tinv[0][0],&eye4x4[0][0],16*sizeof(double));
  cylinder->textureMap=&texMap;
  cylinder->frontAndBack=0;
  cylinder->isLightSource=0;
 }
 return(cylinder);
}

void planeIntersect(struct object3D *plane, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b)
{
 // Computes and returns the value of 'lambda' at the intersection
 // between the specified ray and the specified canonical plane.
  
 struct ray3D *transformedRay;			//holds the transformed ray in model space
 struct ray3D _transformedRay;			//pointer above is pointing to this object 
 struct point3D p0, d;				//origin and direction of transformed ray 
 struct point3D *planeNormal; 			//holds the normal to the plane
 struct point3D _planeNormal;			//pointer above is pointing to this object 
 
 transformedRay = &_transformedRay;
 newRay(&p0, &d, transformedRay); 
 
 //Transform the ray to model space
 rayTransform(ray, transformedRay, plane);
 transformedRay->p0.pw = 1.0; 
 transformedRay->d.pw = 0; 
 
 //determine the intersection variable
 //make sure z-component is not zero 
 if(transformedRay->d.pz != 0){
  *lambda = -(transformedRay->p0.pz) / transformedRay->d.pz;
 }
 
 //return if the object is behind the camera or parallel with the camera
 //set lambda to invalid number
 if (*lambda < 0 || transformedRay->d.pz == 0) {
   *lambda = -1;
   return;
 }
 
 //find the intersection point of the ray and the plane
 rayPosition(transformedRay, *lambda, p);
 
 //set the normal of the unit plane to the positive z axis
 planeNormal = &_planeNormal;
 newPoint(0.0, 0.0, 1.0, planeNormal); 
 planeNormal->pw = 0;
  
 //if the ray intersects with the plane, transform the vars to world space
 //checks the bounds of the plane between -1 and 1 in x and y 
 if (p->px >= -1 && p->px <= 1 && p->py >= -1 && p->py <= 1) {
   //check if there is a texture 
   //if there is set the a and b values 
   if (plane->texImg != NULL) {
     *a = (p->px + 1)/2;
     *b = (p->py + 1)/2;
   }
   //transform the normal and the point using the plane transformation 
   normalTransform(planeNormal, n, plane);
   matVecMult(plane->T, p);
   return;
 }
 
 //return a neg lambda value if the ray doesn't intersect with the unit plane
  else 
   *lambda = -1;
   
 return; 
}

//function to calculate the quadratic formula 
int QuadraticFormula (double a, double b, double c, double &t0, double &t1) {
  double discriminant = b*b - 4*a*c;		//holds the discriminant
  
  //no intersection
  if (discriminant < 0)
    return 0;	//return false 
  
  //one intersection point
  else if (discriminant == 0) {
    t0 = -0.5*b/a;
    t1 = t0;
  }
  
  //two intersection points
  else {
    t0 = (-0.5/a)*(b + sqrt(discriminant));
    t1 = (-0.5/a)*(b - sqrt(discriminant));
  }
  
  //return t0 as the closest point
  if (t0 > t1) {
    double temp = t1;
    t1 = t0;
    t0 = temp;
  }
  
  return 1;	//return true 
}

void sphereIntersect(struct object3D *sphere, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b)
{
 // Computes and returns the value of 'lambda' at the intersection
 // between the specified ray and the specified canonical sphere.
 double A, B, C, t0, t1;		//A, B, C used for quadratic formula. t0 and t1 are the intersection points found. 
 struct ray3D *transformedRay;		//holds the transformed ray 
 struct ray3D _transformedRay;		//pointer above points to this ray 
 struct point3D p0, d; 			//start point and direction of transformed ray 
 double phi, theta;			//angles used to calculate texture coordinates 
 struct point3D *normal;		//holds the normal
 struct point3D normalObj;		//pointer above points to this point 
 
 //initialize transformed ray 
 transformedRay = &_transformedRay;
 newRay(&p0, &d, transformedRay); 
 
 //Transform the ray to model space
 rayTransform(ray, transformedRay, sphere);
 
 //determine the intersection variable
 //coefficients used for the quadratic formula
 A = dot(&transformedRay->d, &transformedRay->d);
 B = 2 * dot(&transformedRay->p0, &transformedRay->d);
 C = dot(&transformedRay->p0, &transformedRay->p0) - 1;
 
 //compute the intersection points and determine the closest valid intersection point
 //if there is no valid intersection point, return with negative lambda (invalid) 
 if (QuadraticFormula(A,B,C,t0,t1) == 0) {
   *lambda = -1;
   return;
 }
 
 //select the closest valid lambda
 if (t0 <= 0) {
   t0 = t1;
   if (t0 <= 0) {
     *lambda = -1;
     return;
   }
 }
 
 //set the intersection point to thye closest intersection found above
 *lambda = t0;
 rayPosition(transformedRay, *lambda, p);

 //check if sphere has a texture 
 if (sphere->texImg != NULL) {
  theta = acos(p->pz);
  phi = atan2(p->py,p->px);
  
  //set a and b 
  *a = fmod(phi,(2*PI))/(2*PI);
  *b = (PI - theta)/PI;
 }
 
 //transform points back to world coordinates 
 normal = &normalObj;
 
 newPoint(-p->px, -p->py, -p->pz, normal);
 normal->pw = 0;
 
 //transform normal 
 normalTransform(normal, n, sphere);
 n->pw = 0; //vector 
 //transform point 
 matVecMult(sphere->T, p);
 
 return;
}

void cylinderIntersect(struct object3D *cylinder, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b){
  double A, B, C, t0, t1; 			//A, B, C used for quadratic formula coefficients. t0 and t1 are the intersection values.
  struct ray3D *transformedRay; 		//holds the transformed ray 
  struct ray3D _transformedRay; 		//pointer above points to this ray 
  struct point3D p0, d; 			//start point and direction of transformed ray 
  struct point3D *normal;
  struct point3D normalObj;
  
  //transform ray 
  transformedRay = &_transformedRay; 
  newRay(&p0, &d, transformedRay); 
  rayTransform(ray, transformedRay, cylinder); 
  
  //determine the intersection values
  //set coefficients for quadratic formula
  A = (transformedRay->d.px * transformedRay->d.px) + (transformedRay->d.pz * transformedRay->d.pz); 
  B = 2 * (transformedRay->p0.px * transformedRay->d.px) + 2 * (transformedRay->p0.pz * transformedRay->d.pz);
  C = transformedRay->p0.px * transformedRay->p0.px + transformedRay->p0.pz * transformedRay->p0.pz - 1; 
  
  //compute the intersection points and determine the closest valid intersection point
  //if there is no valid intersection point, return with negative lambda
  if (QuadraticFormula(A,B,C,t0,t1) == 0) {
    *lambda = -1;
    return;
  }
 
  normal = &normalObj;
  
  //make sure t0 holds the closest value for lambda 
  if(t0>t1){
   float temp = t0; 
   t0 = t1; 
   t1 = temp; 
  }
  
  //calculate the y values of both intersections 
  float y0 = transformedRay->p0.py + t0*transformedRay->d.py; 
  float y1 = transformedRay->p0.py + t1*transformedRay->d.py; 
  
  //check if y0 is less than -1 
  if(y0 < -1) {
      //check if y1 is less than -1 
      if(y1 < -1) {
	//both below -1, miss. Return invalid lambda.  
	*lambda = -1; 
	return; 
      }
      
      else {
	//check if hit cap of cylinder
	float cap = t0 + (t1-t0)*(y0+1)/(y0-y1);
	if(cap <= 0){
	  *lambda = -1; 
	  return; 
	}
	*lambda = cap;
	rayPosition(transformedRay, *lambda, p);
	newPoint(0, 1, 0, normal); 
	
	//check if cylinder has a texture 
	if (cylinder->texImg != NULL) {
	  *a = 0.5*p->px + 0.5;
	  *b = 0.5*p->pz + 0.5;
	}
	//transform the normal and the point 
	normalTransform(normal, n, cylinder);
	n->pw = 0; 
	matVecMult(cylinder->T, p);
	return; 
      }
    
  }
  //check if it is within the range of the cylinder 
  //hit the cylinder 
  else if (y0 >= -1 && y0 <= 1){
    if(t0 <= 0){
	*lambda = -1; 
	return; 
    }
    //set lambda, normal, and point 
    *lambda = t0;
    rayPosition(transformedRay, *lambda, p);
    newPoint(-p->px, 0, -p->pz, normal); 
    
    //check if cylinder has texture
    if (cylinder->texImg != NULL) {
      double theta = atan2(p->pz, p->px);
      
      *a = fmod(theta,(2*PI))/(2*PI);
      *b = (p->py + 1)/2;
    }
	
    normalTransform(normal, n, cylinder);
    n->pw = 0; 
    matVecMult(cylinder->T, p); 
    return; 
  }
  
  else if (y0 > 1) {
    if(y1 > 1){
      //both above 1, missed the cylinder 
	*lambda = -1; 
	return; 
    }
    else {
      //hit the other cap
      float cap = t0 + (t1-t0)*(y0-1)/(y0-y1);
      if(cap <= 0) {
	*lambda = -1; 
	return; 
      }
      *lambda = cap;
      rayPosition(transformedRay, *lambda, p);
      newPoint(0, -1, 0, normal); 
      
      //check for texture 
      if (cylinder->texImg != NULL) {
	*a = 0.5*p->px + 0.5;
	*b = 0.5*p->pz + 0.5;
      }
	
      //transform the normal and the point 
      normalTransform(normal, n, cylinder);
      n->pw = 0; //vector 
      matVecMult(cylinder->T, p);
      return; 
      
    }
    
  }
  return; 
}


void loadTexture(struct object3D *o, const char *filename)
{
 // Load a texture image from file and assign it to the
 // specified object
 if (o!=NULL)
 {
  if (o->texImg!=NULL)	// We have previously loaded a texture
  {			// for this object, need to de-allocate it
   if (o->texImg->rgbdata!=NULL) free(o->texImg->rgbdata);
   free(o->texImg);
  }
  o->texImg=readPPMimage(filename);	// Allocate new texture
 }
}

void texMap(struct image *img, double a, double b, double *R, double *G, double *B)
{
 /*
  Function to determine the colour of a textured object at
  the normalized texture coordinates (a,b).

  a and b are texture coordinates in [0 1].
  img is a pointer to the image structure holding the texture for
   a given object.

  The colour is returned in R, G, B. Uses bi-linear interpolation
  to determine texture colour.
 */

 double *rgbIm = (double *)img->rgbdata;		//pointer to the image 
 int i, j;						//coordinates to read the image
 double aPrime, bPrime;					//used for bilinear interpolation
 double Cij, Cipj, Cijp, Cipjp;				//four parts used for interpolation 
 
 //initialize i and j
 i = floor(a * img->sx);
 j = floor(b * img->sy);
 
 //check if i or j is out of bounds 
 if (i >= (img->sx) - 1)
   i-=2;
 if (j >= (img->sy) - 1)
   j-=2;
 
 //setting values for interpolation 
 aPrime = (a * img->sx) - floor(a * img->sx);
 bPrime = (b * img->sx) - floor(b * img->sy);
 
 //Red
 Cij = *(rgbIm + 3*(j*(img->sx) + i)); //C(i,j)
 Cipj = *(rgbIm + 3*(j*(img->sx) + (i+1))); //C(i+1,j)
 Cijp = *(rgbIm + 3*((j+1)*(img->sx) + i)); //C(i,j+1)
 Cipjp = *(rgbIm + 3*((j+1)*(img->sx) + (i+1))); //C(i+1,j+1)
 
 *(R) = (1 - aPrime)*(1 - bPrime)*(Cij) + aPrime*(1 - bPrime)*(Cipj) + (1 - aPrime)*bPrime*(Cijp) + aPrime*bPrime*(Cipjp);
 
 //Green
 Cij = *(rgbIm + 3*(j*(img->sx) + i) + 1); //C(i,j)
 Cipj = *(rgbIm + 3*(j*(img->sx) + (i+1)) + 1); //C(i+1,j)
 Cijp = *(rgbIm + 3*((j+1)*(img->sx) + i) + 1); //C(i,j+1)
 Cipjp = *(rgbIm + 3*((j+1)*(img->sx) + (i+1)) + 1); //C(i+1,j+1)

 *(G) = (1 - aPrime)*(1 - bPrime)*Cij + aPrime*(1 - bPrime)*Cipj + (1 - aPrime)*bPrime*Cijp + aPrime*bPrime*Cipjp;
 
 //Blue
 Cij = *(rgbIm + 3*(j*(img->sx) + i) + 2); //C(i,j)
 Cipj = *(rgbIm + 3*(j*(img->sx) + (i+1)) + 2); //C(i+1,j)
 Cijp = *(rgbIm + 3*((j+1)*(img->sx) + i) + 2); //C(i,j+1)
 Cipjp = *(rgbIm + 3*((j+1)*(img->sx) + (i+1)) + 2); //C(i+1,j+1)

 *(B) = (1 - aPrime)*(1 - bPrime)*Cij + aPrime*(1 - bPrime)*Cipj + (1 - aPrime)*bPrime*Cijp + aPrime*bPrime*Cipjp;
 
 return;
}

void insertObject(struct object3D *o, struct object3D **list)
{
 if (o==NULL) return;
 // Inserts an object into the object list.
 if (*(list)==NULL)
 {
  *(list)=o;
  (*(list))->next=NULL;
 }
 else
 {
  o->next=(*(list))->next;
  (*(list))->next=o;
 }
}

void insertPLS(struct pointLS *l, struct pointLS **list)
{
 if (l==NULL) return;
 // Inserts a light source into the list of light sources
 if (*(list)==NULL)
 {
  *(list)=l;
  (*(list))->next=NULL;
 }
 else
 {
  l->next=(*(list))->next;
  (*(list))->next=l;
 }

}

void addAreaLight(float sx, float sy, float nx, float ny, float nz,\
                  float tx, float ty, float tz, int lx, int ly,\
                  float r, float g, float b, struct object3D **o_list, struct pointLS **l_list)
{
 /*
   This function sets up and inserts a rectangular area light source
   with size (sx, sy)
   orientation given by the normal vector (nx, ny, nz)
   centered at (tx, ty, tz)
   consisting of (lx x ly) point light sources (uniformly sampled)
   and with colour (r,g,b) - which also determines intensity

   Note that the light source is visible as a uniformly colored rectangle and
   casts no shadow. If you require a lightsource to shade another, you must
   make it into a proper solid box with backing and sides of non-light-emitting
   material
 */

 struct pointLS *l;					//
 struct point3D p;					//point of area light in x direction
 struct point3D py;					//point of area light in y direction 
 struct point3D fin; 					//final value of the point light position 
 struct point3D vector1; 				//first vector that spans the area light plane 
 struct point3D *vector2; 				//second vector that spans the area light plane 
 struct point3D normal; 				//normalized version of the normal 
 
 //first vector spanning will be in x-z plane 
 vector1.px = 1; 
 vector1.py = 0; 
 vector1.pz = 1; 
 vector1.pw = 0; 
 
 //normalize to unit vector
 normalize(&vector1); 
 
 //reverse the normal 
 normal.px = -nx;
 normal.py = -ny; 
 normal.pz = -nz; 
 normal.pw = 0; 
 
 //second vector spanning the area light plane
 //will be the cross product of the normal and the first vector spanning the plane 
 vector2 = cross(&vector1, &normal); 
 vector2->pw = 0; 
 
 //normalize to unit vector
 normalize(vector2); 

 //normalize the vectors to be of size sx, sy 
 vector1.px *= sx; 
 vector1.py *= sx; 
 vector1.pz *= sx; 
 vector1.pw = 0; 
 
 vector2->px *= sy; 
 vector2->py *= sy; 
 vector2->pz *= sy; 
 vector2->pw = 0; 
 
  //place lights along these vectors
  //multiply by i/lx-1 to create a ratio as you move along the vector 
  //the same applies for the second vector 
  for(double i = (-lx/2); i < lx/2; i++) {
    p.px = vector1.px*(i/(lx-1)) + tx ;
    p.py = vector1.py*(i/(lx-1)) + ty ;
    p.pz = vector1.pz*(i/(lx-1)) + tz ;
    
      for(double j = (-ly/2); j < ly/2 ; j++) {
	py.px = vector2->px*(j/(ly-1)) + tx;
	py.py = vector2->py*(j/(ly-1)) + ty ; 
	py.pz = vector2->pz*(j/(ly-1)) + tz ;
	py.pw=1;
	
	fin.px = p.px + py.px; 
	fin.py = p.py + py.py;
	fin.pz = p.pz + py.pz;
	fin.pw = 1;
	//insert the light into the list 
	l=newPLS(&fin,r,g,b);
	insertPLS(l,l_list);
      }
  }
 //de-allocate cross product vector 
 free(vector2);
}



///////////////////////////////////
// Geometric transformation section
///////////////////////////////////

void invert(double *T, double *Tinv)
{
 // Computes the inverse of transformation matrix T.
 // the result is returned in Tinv.

 float *U, *s, *V, *rv1;
 int singFlag, i;
 float T3x3[3][3],Tinv3x3[3][3];
 double tx,ty,tz;

 // Because of the fact we're using homogeneous coordinates, we must be careful how
 // we invert the transformation matrix. What we need is the inverse of the
 // 3x3 Affine transform, and -1 * the translation component. If we just invert
 // the entire matrix, junk happens.
 // So, we need a 3x3 matrix for inversion:
 T3x3[0][0]=(float)*(T+(0*4)+0);
 T3x3[0][1]=(float)*(T+(0*4)+1);
 T3x3[0][2]=(float)*(T+(0*4)+2);
 T3x3[1][0]=(float)*(T+(1*4)+0);
 T3x3[1][1]=(float)*(T+(1*4)+1);
 T3x3[1][2]=(float)*(T+(1*4)+2);
 T3x3[2][0]=(float)*(T+(2*4)+0);
 T3x3[2][1]=(float)*(T+(2*4)+1);
 T3x3[2][2]=(float)*(T+(2*4)+2);
 // Happily, we don't need to do this often.
 // Now for the translation component:
 tx=-(*(T+(0*4)+3));
 ty=-(*(T+(1*4)+3));
 tz=-(*(T+(2*4)+3));

 // Invert the affine transform
 U=NULL;
 s=NULL;
 V=NULL;
 rv1=NULL;
 singFlag=0;

 SVD(&T3x3[0][0],3,3,&U,&s,&V,&rv1);
 if (U==NULL||s==NULL||V==NULL)
 {
  fprintf(stderr,"Error: Matrix not invertible for this object, returning identity\n");
  memcpy(Tinv,eye4x4,16*sizeof(double));
  return;
 }

 // Check for singular matrices...
 for (i=0;i<3;i++) if (*(s+i)<1e-9) singFlag=1;
 if (singFlag)
 {
  fprintf(stderr,"Error: Transformation matrix is singular, returning identity\n");
  memcpy(Tinv,eye4x4,16*sizeof(double));
  return;
 }

 // Compute and store inverse matrix
 InvertMatrix(U,s,V,3,&Tinv3x3[0][0]);

 // Now stuff the transform into Tinv
 *(Tinv+(0*4)+0)=(double)Tinv3x3[0][0];
 *(Tinv+(0*4)+1)=(double)Tinv3x3[0][1];
 *(Tinv+(0*4)+2)=(double)Tinv3x3[0][2];
 *(Tinv+(1*4)+0)=(double)Tinv3x3[1][0];
 *(Tinv+(1*4)+1)=(double)Tinv3x3[1][1];
 *(Tinv+(1*4)+2)=(double)Tinv3x3[1][2];
 *(Tinv+(2*4)+0)=(double)Tinv3x3[2][0];
 *(Tinv+(2*4)+1)=(double)Tinv3x3[2][1];
 *(Tinv+(2*4)+2)=(double)Tinv3x3[2][2];
 *(Tinv+(0*4)+3)=Tinv3x3[0][0]*tx + Tinv3x3[0][1]*ty + Tinv3x3[0][2]*tz;
 *(Tinv+(1*4)+3)=Tinv3x3[1][0]*tx + Tinv3x3[1][1]*ty + Tinv3x3[1][2]*tz;
 *(Tinv+(2*4)+3)=Tinv3x3[2][0]*tx + Tinv3x3[2][1]*ty + Tinv3x3[2][2]*tz;
 *(Tinv+(3*4)+3)=1;

 free(U);
 free(s);
 free(V);
}

void RotateX(struct object3D *o, double theta)
{
 // Multiply the current object transformation matrix T in object o
 // by a matrix that rotates the object theta *RADIANS* around the
 // X axis.

 double R[4][4];
 memset(&R[0][0],0,16*sizeof(double));

 R[0][0]=1.0;
 R[1][1]=cos(theta);
 R[1][2]=-sin(theta);
 R[2][1]=sin(theta);
 R[2][2]=cos(theta);
 R[3][3]=1.0;

 matMult(R,o->T);
}

void RotateY(struct object3D *o, double theta)
{
 // Multiply the current object transformation matrix T in object o
 // by a matrix that rotates the object theta *RADIANS* around the
 // Y axis.

 double R[4][4];
 memset(&R[0][0],0,16*sizeof(double));

 R[0][0]=cos(theta);
 R[0][2]=sin(theta);
 R[1][1]=1.0;
 R[2][0]=-sin(theta);
 R[2][2]=cos(theta);
 R[3][3]=1.0;

 matMult(R,o->T);
}

void RotateZ(struct object3D *o, double theta)
{
 // Multiply the current object transformation matrix T in object o
 // by a matrix that rotates the object theta *RADIANS* around the
 // Z axis.

 double R[4][4];
 memset(&R[0][0],0,16*sizeof(double));

 R[0][0]=cos(theta);
 R[0][1]=-sin(theta);
 R[1][0]=sin(theta);
 R[1][1]=cos(theta);
 R[2][2]=1.0;
 R[3][3]=1.0;

 matMult(R,o->T);
}

void Translate(struct object3D *o, double tx, double ty, double tz)
{
 // Multiply the current object transformation matrix T in object o
 // by a matrix that translates the object by the specified amounts.

 double tr[4][4];
 memset(&tr[0][0],0,16*sizeof(double));

 tr[0][0]=1.0;
 tr[1][1]=1.0;
 tr[2][2]=1.0;
 tr[0][3]=tx;
 tr[1][3]=ty;
 tr[2][3]=tz;
 tr[3][3]=1.0;

 matMult(tr,o->T);
}

void Scale(struct object3D *o, double sx, double sy, double sz)
{
 // Multiply the current object transformation matrix T in object o
 // by a matrix that scales the object as indicated.

 double S[4][4];
 memset(&S[0][0],0,16*sizeof(double));

 S[0][0]=sx;
 S[1][1]=sy;
 S[2][2]=sz;
 S[3][3]=1.0;

 matMult(S,o->T);
}

void printmatrix(double mat[4][4])
{
 fprintf(stderr,"Matrix contains:\n");
 fprintf(stderr,"%f %f %f %f\n",mat[0][0],mat[0][1],mat[0][2],mat[0][3]);
 fprintf(stderr,"%f %f %f %f\n",mat[1][0],mat[1][1],mat[1][2],mat[1][3]);
 fprintf(stderr,"%f %f %f %f\n",mat[2][0],mat[2][1],mat[2][2],mat[2][3]);
 fprintf(stderr,"%f %f %f %f\n",mat[3][0],mat[3][1],mat[3][2],mat[3][3]);
}

/////////////////////////////////////////
// Camera and view setup
/////////////////////////////////////////
struct view *setupView(struct point3D *e, struct point3D *g, struct point3D *up, double f, double wl, double wt, double wsize)
{
 /*
   This function sets up the camera axes and viewing direction as discussed in the
   lecture notes.
   e - Camera center
   g - Gaze direction
   up - Up vector
   fov - Fild of view in degrees
   f - focal length
 */
 struct view *c;
 struct point3D *u, *v;

 u=v=NULL;

 // Allocate space for the camera structure
 c=(struct view *)calloc(1,sizeof(struct view));
 if (c==NULL)
 {
  fprintf(stderr,"Out of memory setting up camera model!\n");
  return(NULL);
 }

 // Set up camera center and axes
 c->e.px=e->px;		// Copy camera center location, note we must make sure
 c->e.py=e->py;		// the camera center provided to this function has pw=1
 c->e.pz=e->pz;
 c->e.pw=0;

 // Set up w vector (camera's Z axis). w=-g/||g||
 c->w.px=-g->px;
 c->w.py=-g->py;
 c->w.pz=-g->pz;
 c->w.pw=0;
 normalize(&c->w);

 // Set up the horizontal direction, which must be perpenticular to w and up
 struct point3D _u;
 u = &_u;
 
 u=cross(&c->w, up);
 normalize(u);
 c->u.px=u->px;
 c->u.py=u->py;
 c->u.pz=u->pz;
 c->u.pw=0; //used to be 1

 // Set up the remaining direction, v=(u x w)  - Mind the signs
 struct point3D _v;
 v = &_v;
 
 v=cross(&c->u, &c->w);
 normalize(v);
 c->v.px=v->px;
 c->v.py=v->py;
 c->v.pz=v->pz;
 c->v.pw=0;

 // Copy focal length and window size parameters
 c->f=f;
 c->wl=wl;
 c->wt=wt;
 c->wsize=wsize;

 // Set up coordinate conversion matrices
 // Camera2World matrix (M_cw in the notes)
 // Mind the indexing convention [row][col]
 c->C2W[0][0]=c->u.px;
 c->C2W[1][0]=c->u.py;
 c->C2W[2][0]=c->u.pz;
 c->C2W[3][0]=0;

 c->C2W[0][1]=c->v.px;
 c->C2W[1][1]=c->v.py;
 c->C2W[2][1]=c->v.pz;
 c->C2W[3][1]=0;

 c->C2W[0][2]=c->w.px;
 c->C2W[1][2]=c->w.py;
 c->C2W[2][2]=c->w.pz;
 c->C2W[3][2]=0;

 c->C2W[0][3]=c->e.px;
 c->C2W[1][3]=c->e.py;
 c->C2W[2][3]=c->e.pz;
 c->C2W[3][3]=1;

 // World2Camera matrix (M_wc in the notes)
 // Mind the indexing convention [row][col]
 c->W2C[0][0]=c->u.px;
 c->W2C[1][0]=c->v.px;
 c->W2C[2][0]=c->w.px;
 c->W2C[3][0]=0;

 c->W2C[0][1]=c->u.py;
 c->W2C[1][1]=c->v.py;
 c->W2C[2][1]=c->w.py;
 c->W2C[3][1]=0;

 c->W2C[0][2]=c->u.pz;
 c->W2C[1][2]=c->v.pz;
 c->W2C[2][2]=c->w.pz;
 c->W2C[3][2]=0;

 c->W2C[0][3]=-dot(&c->u,&c->e);
 c->W2C[1][3]=-dot(&c->v,&c->e);
 c->W2C[2][3]=-dot(&c->w,&c->e);
 c->W2C[3][3]=1;

 free(u);
 free(v);
 return(c);
}

/////////////////////////////////////////
// Image I/O section
/////////////////////////////////////////
struct image *readPPMimage(const char *filename)
{
 // Reads an image from a .ppm file. A .ppm file is a very simple image representation
 // format with a text header followed by the binary RGB data at 24bits per pixel.
 // The header has the following form:
 //
 // P6
 // # One or more comment lines preceded by '#'
 // 340 200
 // 255
 //
 // The first line 'P6' is the .ppm format identifier, this is followed by one or more
 // lines with comments, typically used to inidicate which program generated the
 // .ppm file.
 // After the comments, a line with two integer values specifies the image resolution
 // as number of pixels in x and number of pixels in y.
 // The final line of the header stores the maximum value for pixels in the image,
 // usually 255.
 // After this last header line, binary data stores the RGB values for each pixel
 // in row-major order. Each pixel requires 3 bytes ordered R, G, and B.
 //
 // NOTE: Windows file handling is rather crotchetty. You may have to change the
 //       way this file is accessed if the images are being corrupted on read
 //       on Windows.
 //
 // readPPMdata converts the image colour information to floating point. This is so that
 // the texture mapping function doesn't have to do the conversion every time
 // it is asked to return the colour at a specific location.
 //

 FILE *f;
 struct image *im;
 char line[1024];
 int sizx,sizy;
 int i;
 unsigned char *tmp;
 double *fRGB;

 im=(struct image *)calloc(1,sizeof(struct image));
 if (im!=NULL)
 {
  im->rgbdata=NULL;
  f=fopen(filename,"rb+");
  if (f==NULL)
  {
   fprintf(stderr,"Unable to open file %s for reading, please check name and path\n",filename);
   free(im);
   return(NULL);
  }
  fgets(&line[0],1000,f);
  if (strcmp(&line[0],"P6\n")!=0)
  {
   fprintf(stderr,"Wrong file format, not a .ppm file or header end-of-line characters missing\n");
   free(im);
   fclose(f);
   return(NULL);
  }
  fprintf(stderr,"%s\n",line);
  // Skip over comments
  fgets(&line[0],511,f);
  while (line[0]=='#')
  {
   fprintf(stderr,"%s",line);
   fgets(&line[0],511,f);
  }
  sscanf(&line[0],"%d %d\n",&sizx,&sizy);           // Read file size
  fprintf(stderr,"nx=%d, ny=%d\n\n",sizx,sizy);
  im->sx=sizx;
  im->sy=sizy;

  fgets(&line[0],9,f);  	                // Read the remaining header line
  fprintf(stderr,"%s\n",line);
  tmp=(unsigned char *)calloc(sizx*sizy*3,sizeof(unsigned char));
  fRGB=(double *)calloc(sizx*sizy*3,sizeof(double));
  if (tmp==NULL||fRGB==NULL)
  {
   fprintf(stderr,"Out of memory allocating space for image\n");
   free(im);
   fclose(f);
   return(NULL);
  }

  fread(tmp,sizx*sizy*3*sizeof(unsigned char),1,f);
  fclose(f);

  // Conversion to floating point
  for (i=0; i<sizx*sizy*3; i++) *(fRGB+i)=((double)*(tmp+i))/255.0;
  free(tmp);
  im->rgbdata=(void *)fRGB;

  return(im);
 }

 fprintf(stderr,"Unable to allocate memory for image structure\n");
 return(NULL);
}

struct image *newImage(int size_x, int size_y)
{
 // Allocates and returns a new image with all zeros. Assumes 24 bit per pixel,
 // unsigned char array.
 struct image *im;

 im=(struct image *)calloc(1,sizeof(struct image));
 if (im!=NULL)
 {
  im->rgbdata=NULL;
  im->sx=size_x;
  im->sy=size_y;
  im->rgbdata=(void *)calloc(size_x*size_y*3,sizeof(unsigned char));
  if (im->rgbdata!=NULL) return(im);
 }
 fprintf(stderr,"Unable to allocate memory for new image\n");
 return(NULL);
}

void imageOutput(struct image *im, const char *filename)
{
 // Writes out a .ppm file from the image data contained in 'im'.
 // Note that Windows typically doesn't know how to open .ppm
 // images. Use Gimp or any other seious image processing
 // software to display .ppm images.
 // Also, note that because of Windows file format management,
 // you may have to modify this file to get image output on
 // Windows machines to work properly.
 //
 // Assumes a 24 bit per pixel image stored as unsigned chars
 //

 FILE *f;

 if (im!=NULL)
  if (im->rgbdata!=NULL)
  {
   f=fopen(filename,"wb+");
   if (f==NULL)
   {
    fprintf(stderr,"Unable to open file %s for output! No image written\n",filename);
    return;
   }
   fprintf(f,"P6\n");
   fprintf(f,"# Output from RayTracer.c\n");
   fprintf(f,"%d %d\n",im->sx,im->sy);
   fprintf(f,"255\n");
   fwrite((unsigned char *)im->rgbdata,im->sx*im->sy*3*sizeof(unsigned char),1,f);
   fclose(f);
   return;
  }
 fprintf(stderr,"imageOutput(): Specified image is empty. Nothing output\n");
}

void deleteImage(struct image *im)
{
 // De-allocates memory reserved for the image stored in 'im'
 if (im!=NULL)
 {
  if (im->rgbdata!=NULL) free(im->rgbdata);
  free(im);
 }
}

void cleanup(struct object3D *o_list, struct pointLS *l_list)
{
 // De-allocates memory reserved for the object list and the point light source
 // list. Note that *YOU* must de-allocate any memory reserved for images
 // rendered by the raytracer.
 struct object3D *p, *q;
 struct pointLS *r, *s;

 p=o_list;		// De-allocate all memory from objects in the list
 while(p!=NULL)
 {
  q=p->next;
  if (p->texImg!=NULL)
  {
   if (p->texImg->rgbdata!=NULL) free(p->texImg->rgbdata);
   free(p->texImg);
  }
  free(p);
  p=q;
 }

 r=l_list;
 while(r!=NULL)
 {
  s=r->next;
  free(r);
  r=s;
 }
}

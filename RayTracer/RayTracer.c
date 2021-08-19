/*
  CSC418 - RayTracer code - Winter 2017 - Assignment 3&4

  Written Dec. 9 2010 - Jan 20, 2011 by F. J. Estrada
  Freely distributable for adacemic purposes only.

  Uses Tom F. El-Maraghi's code for computing inverse
  matrices. You will need to compile together with
  svdDynamic.c

  You need to understand the code provided in
  this file, the corresponding header file, and the
  utils.c and utils.h files. Do not worry about
  svdDynamic.c, we need it only to compute
  inverse matrices.

  You only need to modify or add code in sections
  clearly marked "TO DO"
  
  DOUGLAS MCDERMOTT - 1001229823
  
*/

#include "utils.h"
#include <omp.h>
// A couple of global structures and data: An object list, a light list, and the
// maximum recursion depth
struct object3D *object_list;			//holds the list of objects in the scene 
struct pointLS *light_list;			//holds the list of lights in the scene 
int MAX_DEPTH;					//holds the max depth of recursion for reflection/refraction 
struct object3D *cubeMap;			//holds the environment map  
struct object3D _cubeMap;			//pointer above points to this object 

#define INF 1000000000				//infinity used for initializing lambda
#define MAXLIGHTS 16				//holds the max number of lights for the area light 
#define THREADS 8				//number of threads for multithreading 

void buildScene(void)
{
 // Sets up all objects in the scene. This involves creating each object,
 // defining the transformations needed to shape and position it as
 // desired, specifying the reflectance properties (albedos and colours)
 // and setting up textures where needed.
 // Light sources must be defined, positioned, and their colour defined.
 // All objects must be inserted in the object_list. All light sources
 // must be inserted in the light_list.
  
 struct object3D *o;		//holds the object to be inserted
 struct pointLS *l;		//holds the light to be inserted 
 struct point3D p;		//holds position 
 cubeMap = &_cubeMap;
 
 //TO SWITCH BETWEEN SCENES UNCOMMENT THE ASSIGNMENT WANTED 
 
 /*
 // **************** Simple scene for Assignment 3 ****************
 // Insert a couple of objects. A plane and two spheres
 // with some transformations.
 
 // Let's add a plane
 // Note the parameters: ra, rd, rs, rg, R, G, B, alpha, r_index, and shinyness)
 
 o=newPlane(.05,.75,.05,.05,.55,.8,.75,1,1,2);	// Note the plane is highly-reflective (rs=rg=.75) so we
						// should see some reflections if all is done properly.
						// Colour is close to cyan, and currently the plane is
						// completely opaque (alpha=1). The refraction index is
						// meaningless since alpha=1
 Scale(o,6,6,1);				// Do a few transforms...
 RotateZ(o,PI/1.20);
 RotateX(o,PI/2.25);
 Translate(o,0,-3,10);
 invert(&o->T[0][0],&o->Tinv[0][0]);		// Very important! compute
						// and store the inverse
						// transform for this object!
 insertObject(o,&object_list);			// Insert into object list
 
 // Let's add a couple spheres
 o=newSphere(.05,.95,.35,.5,1,.25,.25,1,1,50);
 Scale(o,.75,.5,1.5);
 RotateY(o,PI/2);
 Translate(o,-1.45,1.1,3.5);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
 
 o=newSphere(.05,.95,.95,.5,.75,.95,.55,1,1,50);
 Scale(o,.5,2.0,1.0);
 RotateZ(o,PI/1.5);
 Translate(o,1.75,1.25,5.0);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
 
 //**************** end of assignment 3 scene ****************
 */
 
 
 //**************** Assignment 4 scene ****************
 
 //****************  CUBE MAP  ****************
 
 loadTexture(cubeMap, "enviroMap.ppm");
 
 // ***************  FLOOR  ***************
 o=newPlane(.05, .75, .05, .05, .5, 0, 0, 1, 1, 2);	
 Scale(o,10,17,1);	
 RotateX(o,PI/2);
 Translate(o,0,-5,10);
 Translate(o, 0, 0, 1);
 loadTexture(o, "darkMarbleText.ppm"); 
 invert(&o->T[0][0],&o->Tinv[0][0]);		
 insertObject(o,&object_list);		
 
 // *************** CENTER ***************
 //center pillar
 o=newCylinder(.05,.95,.95,.05, 0.97 , 0.97 , 0.97 ,1,1,50);
 Scale(o,0.5,2.5,0.5);
 Translate(o,0,-3,7);
 Translate(o, 0, 0, 3);
 loadTexture(o, "marbleText.ppm");
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 //top part of pillar 
 o=newCylinder(.05,.75,.5, 0, 0.97 , 0.97 , 0.97 ,1,1,50);
 Scale(o,0.7,0.15,0.7);
 Translate(o, 0, -0.5, 7);
 Translate(o, 0, 0, 3);
 loadTexture(o, "marbleText.ppm");
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
 
 //bottom part of pillar 
 o=newCylinder(.05,.75,.05, 0, 0.97 , 0.97 , 0.97 ,1,1,50);
 Scale(o,0.7,0.15,0.7);
 Translate(o, 0, -4.8, 7);
 Translate(o, 0, 0, 3);
 loadTexture(o, "marbleText.ppm");
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
 
 //sphere on center pillar 
 //this sphere fully reflects the environment
 o=newSphere(.23, .27, .77, 1, 0.25, 0.25, 0.25, 1, 1, 89.6);
 RotateX(o,PI/2);
 RotateY(o,PI);
 Scale(o, 1.4, 1.4, 1.4);
 Translate(o, 0, 1.2, 7);
 Translate(o, 0, 0, 3);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
 
 // *************** BACK LEFT ***************
 //back left pillar
 o=newCylinder(.05,.95,.95,.05, 0.97 , 0.97 , 0.97 ,1,1,50);
 Scale(o,0.5,2,0.5);
 Translate(o,-3, -4, 4);
 Translate(o, 0, 0, 9);
 loadTexture(o, "marbleText.ppm");
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
 
 //top part of pillar
 o=newCylinder(.05,.75,.5, 0, 0.97 , 0.97 , 0.97 ,1,1,50);
 Scale(o,0.7,0.15,0.7);
 Translate(o, -3, -2, 4);
 Translate(o, 0, 0, 9);
 loadTexture(o, "marbleText.ppm");
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
 
 //bottom part of pillar 
 o=newCylinder(.05,.75,.05, 0, 0.97 , 0.97 , 0.97 ,1,1,50);
 Scale(o,0.7,0.15,0.7);
 Translate(o, -3, -4.8, 4);
 Translate(o, 0, 0, 9);
 loadTexture(o, "marbleText.ppm");
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
 
 //sphere on back left pillar 
 //this sphere is the water element
 o=newSphere(.05, .75, .77, 0.25, 0.25, 0.25, 0.25, 1, 1, 89.6);
 RotateX(o,PI/2);
 RotateY(o,PI);
 Scale(o, 0.75, 0.75, 0.75);
 Translate(o, -3, -1, 4);
 Translate(o, 0, 0, 9);
 loadTexture(o, "waterText.ppm");
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
 
 // *************** BACK RIGHT ***************
 //back right pillar
 o=newCylinder(.05,.95,.95,.05, 0.97 , 0.97 , 0.97 ,1,1,50);
 Scale(o,0.5,2,0.5);
 Translate(o,3,-4,4);
 Translate(o, 0, 0, 9);
 loadTexture(o, "marbleText.ppm");
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
 
 //top part of pillar 
 o=newCylinder(.05,.75,.5, 0, 0.97 , 0.97 , 0.97 ,1,1,50);
 Scale(o,0.7,0.15,0.7);
 Translate(o, 3, -2, 4);
 Translate(o, 0, 0, 9);
 loadTexture(o, "marbleText.ppm");
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
 
 //bottom part of pillar
 o=newCylinder(.05,.75,.05, 0, 0.97 , 0.97 , 0.97 ,1,1,50);
 Scale(o,0.7,0.15,0.7);
 Translate(o, 3, -4.8, 4);
 Translate(o, 0, 0, 9);
 loadTexture(o, "marbleText.ppm");
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
 
 //sphere on back right pillar 
 //this sphere is the earth element
 o=newSphere(.05, .75, .77, 0.25, 0.25, 0.25, 0.25, 1, 1, 89.6);
 RotateX(o,PI/2);
 RotateY(o,PI);
 Scale(o, 0.75, 0.75, 0.75);
 Translate(o,3,-1,4);
 Translate(o, 0, 0, 9);
 loadTexture(o, "earthText.ppm");
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
 
 // *************** MIDDLE LEFT ***************
 //middle left pillar
 o=newCylinder(.05,.95,.95,.05, 0.97 , 0.97 , 0.97 ,1,1,50);
 Scale(o,0.5,2,0.5);
 Translate(o,-6,-4,7);
 Translate(o, 0, 0, 3);
 loadTexture(o, "marbleText.ppm");
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
 
 //top part of pillar
 o=newCylinder(.05,.75,.5, 0, 0.97 , 0.97 , 0.97 ,1,1,50);
 Scale(o,0.7,0.15,0.7);
 Translate(o, -6, -2, 7);
 Translate(o, 0, 0, 3);
 loadTexture(o, "marbleText.ppm");
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
 
 //bottom part of pillar
 o=newCylinder(.05,.75,.05, 0, 0.97 , 0.97 , 0.97 ,1,1,50);
 Scale(o,0.7,0.15,0.7);
 Translate(o, -6, -4.8, 7);
 Translate(o, 0, 0, 3);
 loadTexture(o, "marbleText.ppm");
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
 
 //sphere on middle left pillar 
 //this sphere refracts 
 o=newSphere(0, 0, 0.75, 0.6, 0, 0, 0, 0, 1.05, 50);
 Scale(o, 0.75, 0.75, 0.75);
 Translate(o,-6, -1, 7);
 Translate(o, 0, 0, 3);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
 
 // *************** MIDDLE RIGHT ***************
 //middle right pillar
 o=newCylinder(.05,.95,.95,.05, 0.97 , 0.97 , 0.97 ,1,1,50);
 Scale(o,0.5,2,0.5);
 Translate(o,6,-4,7);
 Translate(o, 0, 0, 3);
 loadTexture(o, "marbleText.ppm");
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
 
 //top part of pillar
 o=newCylinder(.05,.75,.5, 0, 0.97 , 0.97 , 0.97 ,1,1,50);
 Scale(o,0.7,0.15,0.7);
 Translate(o, 6, -2, 7);
 Translate(o, 0, 0, 3);
 loadTexture(o, "marbleText.ppm");
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
 
 //bottom part of pillar
 o=newCylinder(.05,.75,.05, 0, 0.97 , 0.97 , 0.97 ,1,1,50);
 Scale(o,0.7,0.15,0.7);
 Translate(o, 6, -4.8, 7);
 Translate(o, 0, 0, 3);
 loadTexture(o, "marbleText.ppm");
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
 
 //sphere on middle right pillar 
 //this sphere refracts 
 o=newSphere(0, 0, 0.75, 0.6, 0, 0, 0, 0, 1.05, 50);
 Scale(o, 0.75, 0.75, 0.75);
 Translate(o, 6, -1, 7);
 Translate(o, 0, 0, 3);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
 
 // *************** FRONT LEFT ***************
 //front left pillar
 o=newCylinder(.05,.95,.95,.05, 0.97 , 0.97 , 0.97 ,1,1,50);
 Scale(o,0.5,2,0.5);
 Translate(o,-3,-4,1);
 Translate(o, 0, 0, 6);
 loadTexture(o, "marbleText.ppm");
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
 
 //top part of pillar
 o=newCylinder(.05,.75,.5, 0, 0.97 , 0.97 , 0.97 ,1,1,50);
 Scale(o,0.7,0.15,0.7);
 Translate(o, -3, -2, 1);
 Translate(o, 0, 0, 6);
 loadTexture(o, "marbleText.ppm");
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
 
 //bottom part of pillar
 o=newCylinder(.05,.75,.05, 0, 0.97 , 0.97 , 0.97 ,1,1,50);
 Scale(o,0.7,0.15,0.7);
 Translate(o, -3, -4.8, 1);
 Translate(o, 0, 0, 6);
 loadTexture(o, "marbleText.ppm");
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
 
 //sphere on front left pillar 
 //this sphere is the fire element 
 o=newSphere(.05, .75, .77, 0.25, 0.25, 0.25, 0.25, 1, 1, 89.6);
 RotateX(o,PI/2);
 RotateY(o,PI);
 Scale(o, 0.75, 0.75, 0.75);
 Translate(o,-3,-1,1);
 Translate(o, 0, 0, 6);
 loadTexture(o, "fireText.ppm");
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
 
 // *************** FRONT RIGHT ***************
 //front right pillar
 o=newCylinder(.05,.75,.5,.0, 0.97 , 0.97 , 0.97 ,1,1,50);
 Scale(o,0.5,2,0.5);
 Translate(o,3,-4,1);
 Translate(o, 0, 0, 6);
 loadTexture(o, "marbleText.ppm");
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
 
 //top part of pillar
 o=newCylinder(.05,.75,.5, 0, 0.97 , 0.97 , 0.97 ,1,1,50);
 Scale(o,0.7,0.15,0.7);
 Translate(o, 3, -2, 1);
 Translate(o, 0, 0, 6);
 loadTexture(o, "marbleText.ppm");
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
 
 //bottom part of pillar
 o=newCylinder(.05,.75,.05, 0, 0.97 , 0.97 , 0.97 ,1,1,50);
 Scale(o,0.7,0.15,0.7);
 Translate(o, 3, -4.8, 1);
 Translate(o, 0, 0, 6);
 loadTexture(o, "marbleText.ppm");
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
 
 //sphere on front right pillar 
 //this sphere is the wind element
 o=newSphere(.05, .75, .77, 0.25, 0.25, 0.25, 0.25, 1, 1, 89.6);
 RotateX(o,PI/2);
 RotateY(o,PI);
 Scale(o, 0.75, 0.75, 0.75);
 Translate(o,3,-1,1);
 Translate(o, 0, 0, 6);
 loadTexture(o, "windText.ppm");
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 // **************** end of assignment 4 scene ****************
 
 
 //call function to add area lights 
 addAreaLight(4, 4, 3, 15.5, -5.5, 3, 15.5, -5.5, 16, 16, .95, .95, .95, &object_list, &light_list); 
 
}

void LoadCubeMapTexture(struct image *img, int index, double a, double b, double *R, double *G, double *B) {
  
 double *rgbIm = (double *)img->rgbdata;		//holds the image
 int i, j;						//used to retrieve texture 
 int ioff;						//holds the offset needed to retrieve the correct image in env map 
 
 //calculate the offset for the given index
 ioff = (index * img->sx)/6;
 i = floor(a * (img->sx / 6)) + ioff;
 j = floor(b * img->sy);
 
 //get the correct pixel colour
 *(R) = *(rgbIm + 3*(j*(img->sx) + i));
 *(G) = *(rgbIm + 3*(j*(img->sx) + i) + 1);
 *(B) = *(rgbIm + 3*(j*(img->sx) + i) + 2);
  
 return;
}

void CubeMap(struct point3D dir, int *index, double *u, double *v) {
  double x = dir.px;		//x value of the direction of ray 
  double y = dir.py;		//y value of the direction of ray 
  double z = dir.pz;		//z value of the direction of ray 
  
  double absX = fabs(x);	//absolute value of x-value
  double absY = fabs(y);	//absolute value of y-value
  double absZ = fabs(z);	//absolute value of z-value
  
  double maxAxis, uc, vc;	//maxAxis holds the largest component of ray. uc, vc are the axis used to traverse the texture 
  
  //find the wall that the ray hits 
  
  // X is the highest component 
  if (absX >= absY && absX >= absZ) {
    maxAxis = absX;
    // POSITIVE X
    if (x > 0) {
      *index = 1;
      uc = -z;
      vc = -y;
    }
    // NEGATIVE X
    else {
      *index = 0;
      uc = z;
      vc = -y;
    }
   }
  
  //Y is the highest component 
  if (absY >= absX && absY >= absZ) {
    maxAxis = absY;
    // POSITIVE Y
    if (y > 0) {
      *index = 2;
      uc = -x;
      vc = -z;
    }
    // NEGATIVE Y
    else {
      *index = 3;
      uc = x;
      vc = -z;
    }
  }
  
  //Z is the highest component 
  if (absZ >= absX && absZ >= absY) {
    maxAxis = absZ;
    // POSITIVE Z
    if (z > 0) {
      *index = 5;
      uc = x;
      vc = -y;
    }
    // NEGATIVE Z
    else {
      *index = 4;
      uc = -x;
      vc = -y;
    }
  }

  //calculate the u-v texture coordinates and puts it within [0,1] 
  *u = 0.5f * (uc / maxAxis + 1.0f);
  *v = 0.5f * (vc / maxAxis + 1.0f);
  
  return;
}

void PhongIlluminationModel (struct object3D *obj, struct point3D *p, struct point3D *n, struct ray3D *ray, struct pointLS *lightSource, double &ambient, double &diffuse, double &specular) {

  struct point3D *normal;						//normalizes the normal
  struct point3D *viewpoint;						//holds the normalized ray 
  struct point3D *incidentLight;					//holds the normalized incident 
  struct point3D *reflectedLight;					//holds the normalized reflected 
  struct point3D _normal, _viewpoint, _reflectedLight, _incidentLight;	//pointers above point to these points 
  
  //normalize normal
  normal = &_normal;
  newPoint(n->px, n->py, n->pz, normal);
  normalize(normal); 
  normal->pw = 0; 
   
  //normalize viewpoint
  viewpoint = &_viewpoint;
  newPoint(-ray->d.px, -ray->d.py, -ray->d.pz, viewpoint);
  normalize(viewpoint); 
  viewpoint->pw = 0; 
  
  //find the incident light vector
  incidentLight = &_incidentLight;
  newPoint((p->px - lightSource->p0.px), (p->py - lightSource->p0.py), (p->pz - lightSource->p0.pz), incidentLight);
  incidentLight->pw = 0;
  
  //normalize incident light
  normalize(incidentLight); 
  
  //find the reflected light vector
  reflectedLight = &_reflectedLight;
  newPoint((2*dot(incidentLight, normal)*normal->px) - incidentLight->px, (2*dot(incidentLight, normal)*normal->py) - incidentLight->py, (2*dot(incidentLight, normal)*normal->pz) - incidentLight->pz, reflectedLight);
  reflectedLight->pw = 0;
  
  //normalize reflected light vector
  normalize(reflectedLight); 
  
  //ambient
  ambient = obj->alb.ra;
  
  //diffuse
  diffuse = 0;
  double diffuseReflection = dot(normal, incidentLight);
  
  //max of (0, N dot L)
  if (diffuseReflection > 0) {
    diffuse = obj->alb.rd * diffuseReflection;
  }
  
  //specular
  specular = 0;
  double dotProd = dot(viewpoint, reflectedLight); 
  double specularReflection = 0; 
  
  if(dotProd < 0) {
    specularReflection = pow(dotProd, obj->shinyness);
  }
  
  //max of (0, (V dot R))^alpha
  if (specularReflection > 0) {
    specular = obj->alb.rs * specularReflection;
  }
  
  return;
}

void rtShade(struct object3D *obj, struct point3D *p, struct point3D *n, struct ray3D *ray, int depth, double a, double b, struct colourRGB *col)
{
 // This function implements the shading model as described in lecture. It takes
 // - A pointer to the first object intersected by the ray (to get the colour properties)
 // - The coordinates of the intersection point (in world coordinates)
 // - The normal at the point
 // - The ray (needed to determine the reflection direction to use for the global component, as well as for
 //   the Phong specular component)
 // - The current racursion depth
 // - The (a,b) texture coordinates (meaningless unless texture is enabled)
 //
 // Returns:
 // - The colour for this ray (using the col pointer)
 //
  
 struct colourRGB tmp_col;	// Accumulator for colour components
 struct colourRGB tex_col;	//texture colour 
 double R,G,B;			// Colour for the object in R G and B

 // This will hold the colour as we process all the components of
 // the Phong illumination model
 tmp_col.R=0;
 tmp_col.G=0;
 tmp_col.B=0;
 
 double ambient;					//ambient component 
 double diffuse;					//diffuse component 
 double specular;					//specular component 
 double tempA, tempB;					//holds temp a and b values so we don't overwrite texture coords 
 struct point3D _tempP, _tempN, _movedP, _viewpoint;	//pointers below will point to these points 
 struct ray3D _shadowRay;				//pointer below will point to this ray 
 struct pointLS *lightPtr; 				//pointer to the light to loop through all the lights 
 
 //initialize light pointer 
 lightPtr = light_list; 
 
 struct ray3D *shadowRay;				//holds the ray to be sent to calculate shadows 
 struct point3D *viewpoint;				//holds the viewpoint vector 
 double *lambdaShadow; 					//holds the lambda returned by the shadow ray 
 double value; 						//lambdaShadow points to this variable 
 //initialize lambda
 lambdaShadow = &value;
 
 struct object3D *tempobj;				//temp object to send with shadow ray 
 struct object3D *objHit; 				//object to use with shadow ray 
 struct point3D *tempP; 				//used to send to find first hit 
 struct point3D *tempN; 				//used to send to find first hit 
 struct point3D *movedP;				//point slightly above intersection point 
 
 //temporary variables used to send to findFirstHit 
 tempP = &_tempP;
 newPoint(p->px, p->py, p->pz, tempP); 
 tempP->pw = 1; 
 
 tempN = &_tempN;
 newPoint (n->px, n->py, n->pz, tempN);
 tempN->pw = 0; 
 
 //move point slightly above to eliminate self intersections 
 movedP = &_movedP;
 newPoint(p->px + abs(n->px*0.1), p->py + abs(n->py*0.1), p->pz + abs(n->pz*0.1), movedP); 
 movedP->pw = 1; 
   
 //initialize ray/point to be used later 
 viewpoint = &_viewpoint;
 shadowRay = &_shadowRay;
  
 //check if an object was hit - not null
 if(obj != NULL) {
   //loop through all light sources
   while(lightPtr != NULL){
     //create vector from the point to the light 
     newPoint((lightPtr->p0.px - p->px), (lightPtr->p0.py - p->py), (lightPtr->p0.pz - p->pz), viewpoint);
     viewpoint->pw = 0; 
     
     //ray to the light to check for shadows 
     newRay(movedP, viewpoint, shadowRay); 
     *lambdaShadow = INF; 
     //shoot ray from this point to find if it is in shadow
     findFirstHit(shadowRay, lambdaShadow, obj, &tempobj, tempP, tempN, &tempA, &tempB);
     //check if shadowRay intersected anything 
     //if it didn't, set the color to the phongIllumination 
     //else keep it's color 0 (black)
     if(*lambdaShadow < 0){
       PhongIlluminationModel(obj, p, n, ray, lightPtr, ambient, diffuse, specular);
       
       //if the object has a texture attached to it
       if (obj->texImg != NULL) {
	 texMap(obj->texImg, a, b, &tex_col.R, &tex_col.G, &tex_col.B);
	 tmp_col.R += ((ambient + diffuse)*tex_col.R*lightPtr->col.R + specular) /(MAXLIGHTS*MAXLIGHTS); 
	 tmp_col.G += ((ambient + diffuse)*tex_col.G*lightPtr->col.G + specular) /(MAXLIGHTS*MAXLIGHTS);
	 tmp_col.B += ((ambient + diffuse)*tex_col.B*lightPtr->col.B + specular) /(MAXLIGHTS*MAXLIGHTS);
       }
       //else use the objects RGB value
       else {
	tmp_col.R += ((ambient + diffuse)*obj->col.R*lightPtr->col.R + specular) /(MAXLIGHTS*MAXLIGHTS); 
	tmp_col.G += ((ambient + diffuse)*obj->col.G*lightPtr->col.G + specular) /(MAXLIGHTS*MAXLIGHTS);
	tmp_col.B += ((ambient + diffuse)*obj->col.B*lightPtr->col.B + specular) /(MAXLIGHTS*MAXLIGHTS);
       }
     }
     
     //move to next light
     lightPtr = lightPtr->next; 
   }
 }
 //if it did not intersect any object sample from the environment 
 else if (cubeMap->texImg != NULL) {
   int index;
   point3D direction;
   newPoint(ray->d.px, ray->d.py, ray->d.pz, &direction);
   
   CubeMap(direction, &index, &a, &b);
   LoadCubeMapTexture(cubeMap->texImg, index, a, b, &tex_col.R, &tex_col.G, &tex_col.B);
   tmp_col.R += tex_col.R; 
   tmp_col.G += tex_col.G;
   tmp_col.B += tex_col.B;
 }
 
 //update 'col' with the final colour computed
 col->R += tmp_col.R*255; 
 col->G += tmp_col.G*255; 
 col->B += tmp_col.B*255; 
 
 //cap the colours to 255 
   if (col->R > 255) {
     col->R = 255;
   }
   if (col->G > 255) {
     col->G = 255;
   }
   if (col->B > 255) {
     col->B = 255;
   }
   
 return;
}

void findFirstHit(struct ray3D *ray, double *lambda, struct object3D *Os, struct object3D **obj, struct point3D *p, struct point3D *n, double *a, double *b)
{
 // Find the closest intersection between the ray and any objects in the scene.
 // It returns:
 //   - The lambda at the intersection (or < 0 if no intersection)
 //   - The pointer to the object at the intersection (so we can evaluate the colour in the shading function)
 //   - The location of the intersection point (in p)
 //   - The normal at the intersection point (in n)
 //
 // Os is the 'source' object for the ray we are processing, can be NULL, and is used to ensure we don't 
 // return a self-intersection due to numerical errors for recursive raytrace calls.
 
 double *currentLambda;				//holds the currentLambda of the object being checked
 double value; 					//current lambda points to this value 
 //initialize current lambda 
 currentLambda = &value;
 
 double *currentA, *currentB;			//temp a and b coords to not overwrite the texture coordinates
 double _a, _b;					//pointers above will point to these values
 //initialize values 
 currentA = &_a;
 currentB = &_b;
 
 struct point3D _currentP, _currentN;		//pointers below will point to these values 
 struct point3D *currentP; 			//temp p to send to find first hit 
 struct point3D *currentN;  			//temp n to send to find first hit 
 
 //variables used to hold the value returned by the intersection function
 currentP = &_currentP;
 newPoint(0,0,0,currentP);
 currentN = &_currentN; 
 newPoint(0,0,0,currentN);
 
 //pointer used to loop through all objects
 struct object3D *objectPtr;
 objectPtr = object_list;
  
 //go through every object in the object list and determine the closest intersection point of all objects
 while (objectPtr != NULL) {
   //make sure we don't check the object that we came from - Os
   if(objectPtr->T != Os->T) {
      (objectPtr->intersect)(objectPtr, ray, currentLambda, currentP, currentN, currentA, currentB);
      //if the currentLambda is valid and is lower than the one currently help update 
      if (*currentLambda > 0 && *currentLambda < *lambda) {
	*lambda = *currentLambda;
	*obj = objectPtr;
	p->px = currentP->px; 
	p->py = currentP->py; 
	p->pz = currentP->pz; 
	p->pw = 1; 
	
	n->px = currentN->px;
	n->py = currentN->py; 
	n->pz = currentN->pz; 
	n->pw = 0;
	
	*a = *currentA;
	*b = *currentB;
      }
   }
  
    //move to the next object in the list
      objectPtr = objectPtr->next;
 }

 //if it didn't intersect anything return invalid lambda 
 if (*lambda == INF){
   *lambda = -1;
 }
 
 
}

void rayTrace(struct ray3D *ray, int depth, struct colourRGB *col, struct object3D *Os)
{
 // Ray-Tracing function. It finds the closest intersection between
 // the ray and any scene objects, calls the shading function to
 // determine the colour at this intersection, and returns the
 // colour.
 //
 // Os is needed for recursive calls to ensure that findFirstHit will
 // not simply return a self-intersection due to numerical
 // errors. For the top level call, Os should be NULL. And thereafter
 // it will correspond to the object from which the recursive
 // ray originates.
 
 double lambda;		// Lambda at intersection
 double a,b;		// Texture coordinates
 
 struct colourRGB tmp_col;
 struct object3D *obj;											// Pointer to object at intersection
 struct point3D p;											// Intersection point
 struct point3D n;											// Normal at intersection
 struct colourRGB I;											// Colour returned by shading function
 struct ray3D *reflectedRay;										// Ray used for reflection
 struct ray3D *refractedRay; 										// Ray used for refraction
 struct point3D *normalizedRay;										// Incident ray normalized 
 struct point3D *normal;										// Normalized n used for reflecting/refracting 
 struct point3D *reflectedDir;										// Direction of the reflected ray 
 struct point3D *reflectedP; 										// Start point of the reflected ray 
 struct point3D *refractedDir; 										// Direction of the refracted ray 
 struct point3D *refractedP; 										// Start point of the refracted ray 
 struct point3D _normalizedRay, _normal, _reflectedDir, _reflectedP, _refractedDir, _refractedP;	// Variables used for the pointers above 
 struct ray3D _reflectedRay, _refractedRay;								// Rays used for the pointers above
 double index1, index2, temp; 										// Index holds the index of refraction. Temp is used for swapping. 
 int reflect = 1; 											// Boolean that holds when to reflect, initialized to true. 
 
 //initialize variables 
 lambda = INF;
 obj = NULL;
 
 //local illumination
 findFirstHit(ray, &lambda, Os, &obj, &p, &n, &a, &b);
 rtShade(obj, &p, &n, ray, depth, a, b, col);
 
 //if object is not null
 if(obj != NULL) {
   //check for transparent material and send refractive ray 
    if(obj->alpha != 1) {
	index1 = 1; 
	index2 = obj->r_index; 
      
	//normalize incident ray 
	normalizedRay = &_normalizedRay;
	newPoint(ray->d.px, ray->d.py, ray->d.pz, normalizedRay); 
	normalize(normalizedRay); 
	normalizedRay->pw = 0; 
	
	//normalize the normal 
	normal = &_normal;
	newPoint(n.px, n.py, n.pz, normal);
	normalize(normal); 
	normal->pw = 0; 
	
	//calculate the dot product to check 
	float dotProd = dot(normalizedRay, normal); 
	
	//cap the values between -1 and 1
	if(dotProd < -1)
	  dotProd = -1; 
	if(dotProd > 1)
	  dotProd = 1; 
	
	//if the dotProd is neg, make it pos
	if(dotProd < 0) 
	  dotProd = -dotProd; 
	//otherwise we are going from inside to outside, flip the indexes and flip the direction of the normal 
	else {
	  temp = index1; 
	  index1 = index2; 
	  index2 = temp; 
	  normal->px = -normal->px; 
	  normal->py = -normal->py; 
	  normal->pz = -normal->pz; 
	}
	
	//calculate the ratio of the indexes 
	float indexRatio = index1/index2; 
	//holds the value of the square root term in refraction equation 
	double dirConst; 
	//if it is negative, this is total internal reflection
	//only reflect do not refract
	//set the dirConst to an invalid value 
	if( (1 - ( (1 - pow(dotProd, 2))*pow(indexRatio, 2) )) < 0 ){
	    dirConst = -1;  
	    reflect = 1; 
	}
	//otherwise set the constant normally 
	else {
	  dirConst = sqrt(1 - ( (1 - pow(dotProd, 2))*pow(indexRatio, 2) ) ); 
	}
	
	//start point of the refracted ray is slightly above the intersection point 
	refractedP = &_refractedP;
	newPoint(p.px + n.px*0.001, p.py + n.py*0.001, p.pz + n.pz*0.001, refractedP); 
	refractedP->pw = 1; 
	
	//calculate the direction of the refracted ray based on snell's law 
	refractedDir = &_refractedDir; 
	newPoint((normalizedRay->px - dotProd*normal->px)*indexRatio, (normalizedRay->py - dotProd*normal->py)*indexRatio, (normalizedRay->pz - dotProd*normal->pz)*indexRatio, refractedDir);  
	refractedDir->px -= normal->px * (dirConst);  
	refractedDir->py -= normal->py * (dirConst);  
	refractedDir->pz -= normal->pz * (dirConst);  
	refractedDir->pw = 0; 
	
	//set the refracted ray 
	refractedRay = &_refractedRay; 
	newRay(refractedP, refractedDir, refractedRay); 
	
	//initialize colour to 0
	tmp_col.R = 0;
	tmp_col.G = 0;
	tmp_col.B = 0;
	
	//make sure dirConst is valid, if it is send refracted ray, do not reflect 
	if(dirConst >= 0) {
	  rayTrace(refractedRay, depth+1, &tmp_col, obj);
	  reflect = 0; 
	}
	
	//add the colour returned to the colour of the pixel 
	col->R += tmp_col.R; 
	col->G += tmp_col.G; 
	col->B += tmp_col.B; 
    
	//cap the colour to 255 
	if (col->R > 255)
	  col->R = 255;
	if (col->G > 255)
	  col->G = 255;
	if (col->B > 255)
	  col->B = 255;
    }
   
   //check depth and send reflected ray if reflect is still 1 (no refraction) 
    if(depth < MAX_DEPTH && reflect == 1) {
	//create a reflected ray about the normal at the point of intersection 
	//normalize the incident ray 
	normalizedRay = &_normalizedRay;
	newPoint(ray->d.px, ray->d.py, ray->d.pz, normalizedRay); 
	normalize(normalizedRay); 
	normalizedRay->pw = 0; 
	
	//normalize the normal 
	normal = &_normal;
	newPoint(n.px, n.py, n.pz, normal);
	normalize(normal); 
	normal->pw = 0; 

	//calculate the reflected direction
	reflectedDir = &_reflectedDir;
	newPoint(-((2*dot(normalizedRay, normal)*normal->px) - normalizedRay->px), -((2*dot(normalizedRay, normal)*normal->py) - normalizedRay->py), -((2*dot(normalizedRay, normal)*normal->pz) - normalizedRay->pz), reflectedDir);
	reflectedDir->pw = 0;
	
	//the start point is the intersection point
	//add a little piece of the normal to prevent self intersection 
	reflectedP = &_reflectedP;
	newPoint(p.px + n.px*0.001, p.py + n.py*0.001, p.pz + n.pz*0.001, reflectedP); 
	reflectedP->pw = 1; 
	
	//set the reflected ray 
	reflectedRay = &_reflectedRay;
	newRay(reflectedP, reflectedDir, reflectedRay);
	
	//initialize colour to zero
	tmp_col.R = 0;
	tmp_col.G = 0;
	tmp_col.B = 0;
	
	//send ray recursively, adding to the depth 
	rayTrace(reflectedRay, depth+1, &tmp_col, obj);
	
	//multiply the colour by the rg component of the object 
	tmp_col.R *= obj->alb.rg; 
	tmp_col.G *= obj->alb.rg; 
	tmp_col.B *= obj->alb.rg;
	
	//add to the colour of the pixel 
	col->R += tmp_col.R; 
	col->G += tmp_col.G; 
	col->B += tmp_col.B; 
    
	//cap the colour at 255
	if (col->R > 255)
	  col->R = 255;
	if (col->G > 255)
	  col->G = 255;
	if (col->B > 255)
	  col->B = 255;
    } 
 }
}


int main(int argc, char *argv[])
{
  
  //random value used for anti-aliasing 
 srand48(time(NULL));
  
 // Main function for the raytracer. Parses input parameters,
 // sets up the initial blank image, and calls the functions
 // that set up the scene and do the raytracing.
 struct image *im;		// Will hold the raytraced image
 struct view *cam;		// Camera and view for this scene
 int sx;			// Size of the raytraced image
 int antialiasing;		// Flag to determine whether antialiaing is enabled or disabled
 char output_name[1024];	// Name of the output file for the raytraced .ppm image
 struct point3D e;		// Camera view parameters 'e', 'g', and 'up'
 struct point3D g;
 struct point3D up;
 double du, dv;			// Increase along u and v directions for pixel coordinates
 struct point3D pc,d;		// Point structures to keep the coordinates of a pixel and
				// the direction or a ray
 struct colourRGB background;   // Background colour
 //int i,j;			// Counters for pixel coordinates
 unsigned char *rgbIm;		// Pointer to fill the image 
 
 if (argc<5)
 {
  fprintf(stderr,"RayTracer: Can not parse input parameters\n");
  fprintf(stderr,"USAGE: RayTracer size rec_depth antialias output_name\n");
  fprintf(stderr,"   size = Image size (both along x and y)\n");
  fprintf(stderr,"   rec_depth = Recursion depth\n");
  fprintf(stderr,"   antialias = A single digit, 0 disables antialiasing. Anything else enables antialiasing\n");
  fprintf(stderr,"   output_name = Name of the output file, e.g. MyRender.ppm\n");
  exit(0);
 }
 sx=atoi(argv[1]);
 MAX_DEPTH=atoi(argv[2]);
 if (atoi(argv[3])==0) antialiasing=0; else antialiasing=1;
 strcpy(&output_name[0],argv[4]);

 fprintf(stderr,"Rendering image at %d x %d\n",sx,sx);
 fprintf(stderr,"Recursion depth = %d\n",MAX_DEPTH);
 if (!antialiasing) fprintf(stderr,"Antialising is off\n");
 else fprintf(stderr,"Antialising is on\n");
 fprintf(stderr,"Output file name: %s\n",output_name);
  
 object_list=NULL;
 light_list=NULL;

 // Allocate memory for the new image
 im=newImage(sx, sx);
 if (!im)
 {
  fprintf(stderr,"Unable to allocate memory for raytraced image\n");
  exit(0);
 }
 else rgbIm=(unsigned char *)im->rgbdata;

 
 buildScene();		// Create a scene. This defines all the
			// objects in the world of the raytracer

 // Camera center is at (0,0,-1)
 e.px=0;
 e.py=0;
 e.pz=-3;
 e.pw=1;

 // To define the gaze vector, we choose a point 'pc' in the scene that
 // the camera is looking at, and do the vector subtraction pc-e.
 // Here we set up the camera to be looking at the origin, so g=(0,0,0)-(0,0,-1)
 g.px=0;
 g.py=0;
 g.pz=1;
 g.pw=1;

 // Define the 'up' vector to be the Y axis
 up.px=0;
 up.py=1;
 up.pz=0;
 up.pw=1;

 // Set up view with given the above vectors, a 4x4 window,
 // and a focal length of -1 (why? where is the image plane?)
 // Note that the top-left corner of the window is at (-2, 2)
 // in camera coordinates.
 cam=setupView(&e, &g, &up, -3, -2, 2, 4);

 if (cam==NULL)
 {
  fprintf(stderr,"Unable to set up the view and camera parameters. Our of memory!\n");
  cleanup(object_list,light_list);
  deleteImage(im);
  exit(0);
 }

 // Set up background colour here
 background.R=0;
 background.G=0;
 background.B=0;
 
 // Do the raytracing
 du=cam->wsize/(sx-1);		// du and dv. In the notes in terms of wl and wr, wt and wb,
 dv=-cam->wsize/(sx-1);		// here we use wl, wt, and wsize. du=dv since the image is
				// and dv is negative since y increases downward in pixel
				// coordinates and upward in camera coordinates.

 fprintf(stderr,"View parameters:\n");
 fprintf(stderr,"Left=%f, Top=%f, Width=%f, f=%f\n",cam->wl,cam->wt,cam->wsize,cam->f);
 fprintf(stderr,"Camera to world conversion matrix (make sure it makes sense!):\n");
 printmatrix(cam->C2W);
 fprintf(stderr,"World to camera conversion matrix\n");
 printmatrix(cam->W2C);
 fprintf(stderr,"\n");

 fprintf(stderr,"Rendering rows: ");
 
 #pragma omp parallel for schedule(dynamic) num_threads(THREADS)
 for (int j=0; j<sx; j++)		// For each of the pixels in the image
 {
  fprintf(stderr,"%d/%d, ",j,sx);
  for (int i=0; i<sx; i++)
  {
    struct colourRGB col;		// Return colour for raytraced pixels
    struct ray3D *ray;			// Structure to keep the ray from e to a pixel
    struct point3D imagePlane; 		// Holds image plane values
    struct point3D direction; 		// Holds direction of ray 
    struct point3D origin;		// Holds origin of ray 
    
    origin.px = 0;
    origin.py = 0; 
    origin.pz = 0; 
    origin.pw = 1;
     
    //check if antialiasing is on 
    if (!antialiasing) {
      //pixels
      imagePlane.px = cam->wl + ((cam->wl + cam->wsize) - cam->wl)*(i + 0.5)/sx; 
      imagePlane.py =  cam->wt + ((cam->wt - cam->wsize) - cam->wt)*(j + 0.5)/sx;
      imagePlane.pz = -3; 
      
      //ray direction
      direction.px = imagePlane.px - origin.px; 
      direction.py = imagePlane.py - origin.py; 
      direction.pz = imagePlane.pz - origin.pz; 
      direction.pw = 0; //vector
      
      //convert to world-space
      matVecMult(cam->C2W, &direction);
      matVecMult(cam->C2W, &origin);
      
      //initialize ray to be sent 
      struct ray3D _ray;
      ray = &_ray;
      newRay(&origin, &direction, ray);
      //initialize colour of pixel to 0 
      col.R = 0;
      col.B = 0; 
      col.G = 0; 
      
      //send ray, source is NULL for first call 
      rayTrace(ray, 0, &col, NULL);
   
      //set the colour of the pixel to the one returned 
      *(rgbIm + 3*(j*sx + i)) = col.R; 
      *(rgbIm + 3*(j*sx + i) + 1) = col.G; 
      *(rgbIm + 3*(j*sx + i) + 2) = col.B;
      
    }
    //antialiasing is turned on 
    else {
      double supersamples = 20;		//set the number of rays to be sent 
      struct colourRGB tempCol;		//holds the colour of the supersamples 
      
      //initialize colours to 0
      col.R = 0;
      col.G = 0;
      col.B = 0;
      
      //loop through all supersamples 
      for (int k = 0; k < supersamples; k++) {
	//initialize colour and origin of this sample 
	tempCol.R = 0;
	tempCol.G = 0;
	tempCol.B = 0;
	
	origin.px = 0;
	origin.py = 0; 
	origin.pz = 0; 
	origin.pw = 1; 		//point
	
	//set image plane using a random number added to i and j instead of just 0.5 
	//this is the jitter method 
	imagePlane.px = cam->wl + ((cam->wl + cam->wsize) - cam->wl)*(i + drand48())/(sx);
	imagePlane.py =  cam->wt + ((cam->wt - cam->wsize) - cam->wt)*(j + drand48())/(sx);
	imagePlane.pz = -3;
	
	//ray direction
	direction.px = imagePlane.px - origin.px; 
	direction.py = imagePlane.py - origin.py; 
	direction.pz = imagePlane.pz - origin.pz; 
	direction.pw = 0; 	//vector
	
	//convert to world-space
	matVecMult(cam->C2W, &direction);
	matVecMult(cam->C2W, &origin);

	//initialize ray with found origin and direction 
	struct ray3D _ray;
	ray = &_ray;
	newRay(&origin, &direction, ray);
	
	//send ray, source is NULL for first call 
	rayTrace(ray, 0, &tempCol, NULL);
	
	//add the tempCol of this sample to the total colour 
	col.R += tempCol.R;
	col.G += tempCol.G;
	col.B += tempCol.B;
      }
      
      //divide by number of samples to obtain an average 
      col.R /= supersamples;
      col.G /= supersamples;
      col.B /= supersamples;
      
      //set the colour of the pixel 
      *(rgbIm + 3*(j*sx + i)) = col.R; 
      *(rgbIm + 3*(j*sx + i) + 1) = col.G; 
      *(rgbIm + 3*(j*sx + i) + 2) = col.B;
    }
    
  } // end for i
 } // end for j
 
 fprintf(stderr,"\nDone!\n");

 // Output rendered image
 imageOutput(im,output_name);

 // Exit section. Clean up and return.
 cleanup(object_list,light_list);		// Object and light lists
 deleteImage(im);				// Rendered image
 free(cam);					// camera view
 exit(0);
}

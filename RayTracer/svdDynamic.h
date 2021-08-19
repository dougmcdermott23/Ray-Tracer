/*
 * svdDynamic.h
 * Copyright (c) 2000
 * Thomas F. El-Maraghi
 *
 * Singular value decomposition (SVD) routines.
 *
 */

// Modified for stand-alone use, FEG, Jul 18, 2006

#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#ifndef __SVD_dynamic
#define __SVD_dynamic

#define MAX_SVD_ITERATIONS 100
#define max(A,B) ((A)<(B)?(B):(A))

#define signof(A,B)    (((B)>=0)? (fabs(A)) : (-fabs(A)))

int SVDHelper( const int m, const int n,
	       float *U, float *w, float *V,
	       float *rv1 );
static float SVD_PYTHAG( const float a, const float b );
int SVD( const float *A, const int m, const int n,
	 float **U, float **w, float **V, float **rv1 );
void SortSV( int *svPerm, float *w, const int n );
int SolveLinearSystem( const float *A, const float *b,
		       const int m, const int n,
		       float **x, float **w );
void InvertMatrix( const float *U, const float *w, const float *V,
		   const int n, float *I );
#endif


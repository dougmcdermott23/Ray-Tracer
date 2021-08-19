/*
 * svdDynamic.c
 * Copyright (c) 2000
 * Thomas F. El-Maraghi
 *
 * Singular value decomposition (SVD) routines.
 *
 */

// Modified for stand-alone use, FEG, Jul 18, 2006

#include "svdDynamic.h"

#define signof(A,B)    (((B)>=0)? (fabs(A)) : (-fabs(A)))

/*
 * Returns the Singular Value Decomposition U*w*V^T of the matrix A
 * (passed to the function in U).
 *
 * U orthogonal mxn matrix (mxn, m >= n)
 * w n-vector of signular values (nx1)
 * V orthogonal nxn matrix (nxn)
 * rv1  superdiagonal of singular value matrix
 *
 * This subroutine is a translation of the algol procedure svd, 
 * num. math. 14, 403-420(1970) by golub and reinsch. 
 * handbook for auto. comp., vol ii-linear algebra, 134-151(1971). 
 * See http://www.netlib.org/  for Eispack svd.f Fortran77 version 
 */
int SVDHelper( const int m, const int n,
	       float *U, float *w, float *V, 
	       float *rv1 )
{
  int flag, i, its, j, jj, k, l = 0, nm = 0;
  float c, f, h, s, x, y, z;
  float anorm = (float)0;
  float g = (float)0;
  float scale = (float)0;
  float tst;

  if( m < n )
  {
//    FATAL( "SVD: You must augment A with extra zero rows" );
   fprintf(stderr,"SVD: Thou shouldst augment A with extra zero rows, silly!\n");
   exit(0);
  }
  
  for( i = 0; i < n; i++ ) {
    l = i + 1;
    rv1[i] = (float)(scale * g);
    g = s = scale = (float)0;
    if( i < m ) {
      for( k = i; k < m; k++ ) 
	scale += fabs( U[k*n+i] );
      if( scale ) {
	for( k = i; k < m; k++ ) {
	  U[k*n + i] /= scale;
	  s += (float)(U[k*n+i] * U[k*n+i]);
	}
	f = U[i*n + i];
	g = (float)-signof(sqrt(s),f);
	h = (float)(f * g - s);
	U[i*n+i] = (float)(f - g);
	if( i != n - 1 ) {
	  for( j = l; j < n; j++ ) {
	    for( s = (float)0, k = i; k < m; k++ ) 
	      s += (float)(U[k*n+i] * U[k*n+j]);
	    f = (float)(s / h);
	    for( k = i; k < m; k++ )
	      U[k*n+j] += (float)(f * U[k*n+i]);
	  }
	}
	for( k = i; k < m; k++ ) 
	  U[k*n+i] *= scale;
      }
    }
    w[i] = (float)(scale * g);
    g = s = scale = (float)0;
    if( i < m && i != n - 1) {
      for( k = l; k < n; k++ )
	scale += fabs(U[i*n+k]);
      if( scale ) {
	for( k = l; k < n; k++ ) {
	  U[i*n+k] /= scale;
	  s += (float)(U[i*n+k] * U[i*n+k]);
	}
	f = U[i*n+l];
	g = (float)-signof(sqrt(s),f);
	h = (float)(f * g - s);
	U[i*n+l] = f - g;
	for( k = l; k < n; k++ ) 
	  rv1[k] = (float)(U[i*n+k] / h);
	if( i != m - 1) {
	  for( j = l; j < m; j++ ) {
	    for( s = (float)0, k = l; k < n; k++ )
	      s += (float)(U[j*n+k] * U[i*n+k]);
	    for( k = l; k < n; k++ )
	      U[j*n+k] += (float)(s * rv1[k]);
	  }
	}
	for( k = l; k < n; k++ )
	  U[i*n+k] *= scale;
      }
    }
    anorm = max( anorm, (fabs(w[i]) + fabs(rv1[i])) );
  }
  for( i = n - 1; i >= 0; i-- ) {
    if( i < n - 1 ) {
      if( g ) {
	for( j = l; j < n; j++ )
	  V[j*n+i] = (float)((U[i*n+j] / U[i*n+l]) / g);
	for( j = l; j < n; j++ ) {
	  for( s = (float)0, k = l; k < n; k++ )
	    s += (float)(U[i*n+k] * V[k*n+j]);
	  for( k = l; k < n; k++ )
	    V[k*n+j] += (float)(s * V[k*n+i]);
	}
      }
      for( j = l; j < n; j++ )
	V[i*n+j] = V[j*n+i] = (float)0;
    }
    V[i*n+i] = (float)1;
    g = rv1[i];
    l = i;
  }
  for( i = n - 1; i >= 0; i-- ) {
    l = i + 1;
    g = w[i];
    if( i < n - 1 )
      for( j = l; j < n; j++ )
	U[i*n+j] = (float)0;
    if( g ) {
      g = (float)((float)1 / g);
      if( i != n - 1 ) {
	for( j = l; j < n; j++ ) {
	  for( s = (float)0, k = l; k < m; k++ )
	    s += (float)(U[k*n+i] * U[k*n+j]);
	  f = (float)((s / U[i*n+i]) * g);
	  for( k = i; k < m; k++ )
	    U[k*n+j] += (float)(f * U[k*n+i]);
	}
      }
      for( j = i; j < m; j++ )
	U[j*n+i] *= g;
    } else {
      for( j = i; j < m; j++ )
	U[j*n+i] = (float)0;
    }
    ++U[i*n+i];
  }

  for( k = n - 1; k >= 0; k-- ) {
    for( its = 1; its <= MAX_SVD_ITERATIONS; its++ ) {
      flag = 1;
      for( l = k; l >= 0; l-- ) {
	nm = l - 1;
	tst = fabs(rv1[l]) + anorm;
	if( tst == anorm ) {
	  flag = 0;
	  break;
	}
	tst = fabs(w[nm]) + anorm;
	if( tst == anorm )
	  break;
      }
      /* Found a zero diagonal element w[nm] */
      if( flag ) {
	c = (float)0;
	s = (float)1;
	for( i = l; i <= k; i++ ) {
	  f = (float)(s * rv1[i]);
	  rv1[i] = (float)(c * rv1[i]);
	  tst = fabs(f) + anorm;
	  if( tst == anorm ) 
	    break;
	  else {
	    g = w[i];
	    h = SVD_PYTHAG(f,g);
	    w[i] = h;
	    h = (float)((float)1 / h);
	    c = (float)(g * h);
	    s = (float)(-f * h);
	    for( j = 0; j < m; j++ ) {
	      y = U[j*n+nm];
	      z = U[j*n+i];
	      U[j*n+nm] = (float)(y * c + z * s);
	      U[j*n+i] = (float)(z * c - y * s);
	    }
	  }
	}
      }
      z = w[k];
      if( l == k ) {
     	if( z < (float)0 ) {
	  w[k] = (float)-z;
	  for( j = 0; j < n; j++ )
	    V[j*n+k] = (float)(-V[j*n+k]);
	}
	break;
      }
      if( its >= MAX_SVD_ITERATIONS ) {
	return( -1 );
      }
      
      x = w[l];
      nm = k - 1;
      y = w[nm];
      g = rv1[nm];
      h = rv1[k];
      f = (float)(0.5 * (((g + z) / h) * ((g - z) / y) + y / h - h / y));
      g = SVD_PYTHAG(f,(float)1);
      f = (float)(x - (z / x) * z + (h / x) * (y / (f + signof(g,f)) - h));
      c = s = (float)1;
      for( j = l; j <= nm; j++ ) {
	i = j + 1;
	g = rv1[i];
	y = w[i];
	h = (float)(s * g);
	g = (float)(c * g);
	z = SVD_PYTHAG(f,h);
	rv1[j] = z;
	c = (float)(f / z);
	s = (float)(h / z);
	f = (float)(x * c + g * s);
	g = (float)(g * c - x * s);
	h = (float)(y * s);
	y = (float)(y * c);
	for( jj = 0; jj < n; jj++ ) {
	  x = V[jj*n+j];
	  z = V[jj*n+i];
	  V[jj*n+j] = (float)(x * c + z * s);
	  V[jj*n+i] = (float)(z * c - x * s);
	}
	z = SVD_PYTHAG(f,h);
	w[j] = z;
	if( z != (float)0 ) {
	  c = (float)(f / z);
	  s = (float)(h / z);
	} 
	f = (float)((c * g) + (s * y));
	x = (float)((c * y) - (s * g));
	for( jj = 0; jj < m; jj++ ) {
	  y = U[jj*n+j];
	  z = U[jj*n+i];
	  U[jj*n+j] = (float)(y * c + z * s);
	  U[jj*n+i] = (float)(z * c - y * s);
	}
      }
      rv1[l] = (float)0;
      rv1[k] = f;
      w[k] = x;
    }
  }
  return( 0 );
}

static float SVD_PYTHAG( const float a, const float b )
{
  float at = fabs(a);
  float bt = fabs(b);
  float ct;
  if( at > bt ) {
    ct = bt / at;
    return (float)( at * sqrt( 1.0 + ct * ct ) );
  }
  if( bt <= 0.0 )
    return 0.0;
  ct = at / bt;
  return (float)( bt * sqrt( 1.0 + ct * ct ) );
}


/*
 * Returns the Singular Value Decomposition A = U*w*V^T of the matrix A
 *
 * A input mxn matrix
 * U orthogonal mxn matrix (mxn)
 * w n-vector of signular values (nx1)
 * V orthogonal nxn matrix (nxn)
 * rv1  superdiagonal of singular value matrix
 */
int SVD( const float *A, const int m, const int n,
	 float **U, float **w, float **V, float **rv1 )
{
  float *tmp;
  int r, c, k, svdReturn;
  int N = m * n;
  int ownRv1 = 0;

  if( m >= n ) {

    /* allocate memory if necessary */
    if( *U == NULL ) *U = (float *)calloc(N,sizeof(float));
    if( *w == NULL ) *w = (float *)calloc(n,sizeof(float));
    if( *V == NULL ) *V = (float *)calloc(n*n,sizeof(float));

    if( *rv1 == NULL ) {
      *rv1 = (float *)calloc(n,sizeof(float));
      ownRv1 = 1;
    }

    /* copy A to U if necessary */
    if( *U != A )
      for( k = 0; k < N; k++ )
	(*U)[k] = A[k];

    svdReturn = SVDHelper( m, n, *U, *w, *V, *rv1 );
    
  } else {
    
    /* allocate memory if necessary */
    if( *U == NULL ) *U = (float *)calloc(N,sizeof(float));
    if( *w == NULL ) *w = (float *)calloc(m,sizeof(float));
    if( *V == NULL ) *V = (float *)calloc(m*m,sizeof(float));

    if( *rv1 == NULL ) {
      *rv1 = (float *)calloc(m,sizeof(float));
      ownRv1 = 1;
    }

    /* transpose A and store in U */
    tmp = *U;
    for( c = 0; c < n; c++ )
      for( r = 0; r < m; r++, tmp++ ) 
	*tmp = A[r*n + c];

    svdReturn = SVDHelper( n, m, *U, *w, *V, *rv1 );
    
    /* swap U and V */
    tmp = *U;
    *U = *V;
    *V = tmp;
  }   

  if( ownRv1 ) free(*rv1);

  return svdReturn;
}

/*
 * Sort singular values into decreasing order,
 * return permutation array svPerm *... ith sorted
 * singular value is then w[svPerm[i]] 
 */
void SortSV( int *svPerm, float *w, const int n )
{
  int i, j, iTmp;
  float tmp;

  for( i = 0; i < n; i++ )
    svPerm[i] = i;

  for( i = 0; i < n; i++ ) {
    /* Find max in remaining set i..numSamp */
    tmp = w[svPerm[i]];
    iTmp = i;
    for( j = i + 1; j < n; j++ )
      if( w[svPerm[j]] > tmp ) {
	tmp = w[svPerm[j]];
	iTmp = j;
      }
    /* Switch */
    j = svPerm[i];
    svPerm[i] = svPerm[iTmp];
    svPerm[iTmp] = j;
  }
}


/*
 * Solve the linear system Ax=b
 */
int SolveLinearSystem( const float *A, const float *b,
		       const int m, const int n,
		       float **x, float **w )
{ 
  float *U, *V, *s;
  int i, j, svdReturn;

  if( *x == NULL ) *x = (float *)calloc(n,sizeof(float));

  svdReturn = SVD( A, m, n, &U, w, &V, &s );

  for( i = 0; i < n; i++ ) {      
    s[i] = 0.0;
    for( j = 0; j < n; j++ )          
      s[i] += U[j*n+i] * b[j-1];                                   
    s[i] /= (*w)[i];                        
  }                                                    
                                   
  for( i = 0; i < n; i++ ) {        
    (*x)[i-1] = 0.0;      
    for( j = 0; j < n; j++ )                                            
      (*x)[i-1] += V[i*n+j] * s[j];
  }
 
  free(U);
  free(V);
  free(s);

  return( svdReturn );                  
}
      

/*
 * Given svd decomposition U w V^T of a matrix, compute
 * the inverse (or pseudoinverse) I = V w^-1 U^T 
 */
void InvertMatrix( const float *U, const float *w, const float *V, 
		   const int n, float *I )
{ 
  int i, j, k;
  float *scr;

  scr = (float *)calloc(n,sizeof(float));
  
  for( k = 0; k < n; k++ ) {
    /* Compute scr = kth column of (w^-1 U^T)  */
    for( i = 0; i < n; i++ )
      scr[i] = (float)(U[k*n+i] / w[i]);
      
    /* Compute kth col of Ainv = V * scr */ 
    for( i = 0; i < n; i++ ) {
      I[i*n+k] = (float)0;      
      for( j = 0; j < n; j++ )
	I[i*n+k] += (float)(V[i*n+j] * scr[j]);
    }
  }
  
  free(scr);
}

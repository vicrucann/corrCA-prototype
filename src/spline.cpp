/* Spline interpolation.
    Copyright (C) 2007 Lionel Moisan <Lionel.Moisan@parisdescartes.fr>
    Copyright (C) 2010 Pascal Monasse <monasse@imagine.enpc.fr>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
// WARNING: 
// This file implements an algorithm possibly linked to the patent
//
// David Lowe  "Method and apparatus for identifying scale invariant 
// features in an image and use of same for locating an object in an 
// image",  U.S. Patent 6,711,293.
//
// This file is made available for the exclusive aim of serving as
// scientific tool to verify of the soundness and
// completeness of the algorithm description. Compilation,
// execution and redistribution of this file may violate exclusive
// patents rights in certain countries.
// The situation being different for every country and changing
// over time, it is your responsibility to determine which patent
// rights restrictions apply to you before you compile, use,
// modify, or redistribute this file. A patent lawyer is qualified
// to make this determination.
// If and only if they don't conflict with any patent terms, you
// can benefit from the following license terms attached to this
// file.
//
// This program is provided for scientific and educational only:
// you can use and/or modify it for these purposes, but you are
// not allowed to redistribute this work or derivative works in
// source or executable form. A license must be obtained from the
// patent right holders for any other use.

#include "spline.h"
#include "library.h"
#include <string.h>
//#include "libIO/nan.h"
#include <cmath>
#include <cfloat>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>

static double initcausal(double* c, int step, int n, double z)
{
    double zk,z2k,iz,sum;

    zk = z; iz = 1./z;
    z2k = pow(z,(double)n-1.);
    sum = c[0] + z2k * c[step*(n-1)];
    z2k = z2k*z2k*iz;
    for(int k = 1; k <= n-2; k++) {
        sum += (zk+z2k)*c[step*k];
        zk *= z;
        z2k *= iz;
    }
    return (sum/(1.-zk*zk));
}

static double initanticausal(double* c, int step, int n, double z)
{
    return (z/(z*z-1.)) * (z*c[step*(n-2)]+c[step*(n-1)]);
}

double initcausal(double *c,int n,double z)
{
  double zk,z2k,iz,sum;
  int k;

  zk = z; iz = 1./z;
  z2k = pow(z,(double)n-1.);
  sum = c[0] + z2k * c[n-1];
  z2k = z2k*z2k*iz;
  for (k=1;k<=n-2;k++) {
    sum += (zk+z2k)*c[k];
    zk *= z;
    z2k *= iz;
  }
  return (sum/(1.-zk*zk));
}

double initanticausal(double *c,int n,double z)
{
  return((z/(z*z-1.))*(z*c[n-2]+c[n-1]));
}

static void invspline1D(double* c, int step, int size, double* z, int npoles)
{
    /* normalization */
    double lambda=1.0;
    for(int k = npoles-1; k >= 0; k--)
        lambda *= (1.-z[k])*(1.-1./z[k]);
    for(int n = size-1; n >= 0; n--)
        c[step*n] *= static_cast<double>(lambda);

    for(int k=0 ; k < npoles; k++) { // Loop on poles
        /* forward recursion */
        c[0] = static_cast<double>(initcausal(c, step, size, z[k]));
        for(int n=1; n < size; n++)
            c[step*n] += static_cast<double>(z[k]*c[step*(n-1)]);
        /* backward recursion */
        c[step*(size-1)] = static_cast<double>(initanticausal(c, step, size, z[k]));
        for(int n=size-2; n >= 0; n--)
            c[step*n] = static_cast<double>(z[k]*(c[step*(n+1)]-c[step*n]));
    }
}

void invspline1D(double *c,int size,double *z,int npoles)
{
  double lambda;
  int n,k;

  /* normalization */
  for (k=npoles,lambda=1.;k--;) lambda *= (1.-z[k])*(1.-1./z[k]);
  for (n=size;n--;) c[n] *= lambda;

  /*----- Loop on poles -----*/
  for (k=0;k<npoles;k++) {

    /* forward recursion */
    c[0] = initcausal(c,size,z[k]);
    for (n=1;n<size;n++) 
      c[n] += z[k]*c[n-1];

    /* backwards recursion */
    c[size-1] = initanticausal(c,size,z[k]);
    for (n=size-1;n--;) 
      c[n] = z[k]*(c[n+1]-c[n]);
    
  }
}

/// Put in array \a z the poles of the spline of given \a order.
static bool fill_poles(double* z, int order)
{
    switch(order) {
    case 1:
        break;
        // case 2: z[0]=-0.17157288;  /* sqrt(8)-3 */
        // break;
    case 3: z[0]=-0.26794919;  /* sqrt(3)-2 */ 
        break;
        // case 4: z[0]=-0.361341; z[1]=-0.0137254;
        // break;
    case 5: z[0]=-0.430575; z[1]=-0.0430963;
        break;
        // case 6: z[0]=-0.488295; z[1]=-0.0816793; z[2]=-0.00141415;
        // break;
    case 7: z[0]=-0.53528; z[1]=-0.122555; z[2]=-0.00914869;
        break;
        // case 8: z[0]=-0.574687; z[1]=-0.163035; z[2]=-0.0236323;
        // z[3]=-0.000153821;
        // break;
    case 9: z[0]=-0.607997; z[1]=-0.201751; z[2]=-0.0432226; z[3]=-0.00212131;
        break;
        // case 10: z[0]=-0.636551; z[1]=-0.238183; z[2]=-0.065727;
        // z[3]=-0.00752819; z[4]=-0.0000169828;
        // break;
    case 11: z[0]=-0.661266; z[1]=-0.27218; z[2]=-0.0897596; z[3]=-0.0166696; 
        z[4]=-0.000510558;
        break;
    default:
        return false;
    }
    return true;
}

/// Prepare image for cardinal spline interpolation.
bool prepare_spline(image_double& im, int order)
{
    if(order < 3)
        return true;

    // Init poles of associated z-filter
    double z[5];
    if(! fill_poles(z, order))
        return false;
    int npoles = order/2;

	for(int y = 0; y < im->ysize; y++) // Filter on lines
		invspline1D(im->data+y*im->xsize, 1, im->xsize, z, npoles);
	for(int x = 0; x < im->xsize; x++) // Filter on columns
		invspline1D(im->data+x, 1*im->xsize, im->ysize, z, npoles);
    return true;
}

/* c[] = values of interpolation function at ...,t-2,t-1,t,t+1,... */

/* coefficients for cubic interpolant (Keys' function) */
void keys(float *c,float t,float a)
{
  float t2,at;

  t2 = t*t;
  at = a*t;
  c[0] = a*t2*(1-t);
  c[1] = (2*a+3 - (a+2)*t)*t2 - at;
  c[2] = ((a+2)*t - a-3)*t2 + 1;
  c[3] = a*(t-2)*t2 + at;
}

/* c[] = values of interpolation function at ...,t-2,t-1,t,t+1,... */

/* coefficients for cubic interpolant (Keys' function) */
static void keys(double* c, double t, double a)
{
    double t2 = t*t;
    double at = a*t;
    c[0] = a*t2*(1.0f-t);
    c[1] = (2.0f*a+3.0f - (a+2.0f)*t)*t2 - at;
    c[2] = ((a+2.0f)*t - a-3.0f)*t2 + 1.0f;
    c[3] = a*(t-2.0f)*t2 + at;
}

/* coefficients for cubic spline */
void spline3(float *c,float t)
{
  float tmp;

  tmp = 1-t;
  c[0] = 0.1666666666f*t*t*t;
  c[1] = 0.6666666666f-0.5f*tmp*tmp*(1+t);
  c[2] = 0.6666666666f-0.5f*t*t*(2-t);
  c[3] = 0.1666666666f*tmp*tmp*tmp;
}

/* coefficients for cubic spline */
static void spline3(double* c, double t)
{
    double tmp = 1.0f-t;
    c[0] = 0.1666666666f *t*t*t;
    c[1] = 0.6666666666f -0.5f*tmp*tmp*(1.0f+t);
    c[2] = 0.6666666666f -0.5f*t*t*(2.0f-t);
    c[3] = 0.1666666666f *tmp*tmp*tmp;
}

/* pre-computation for spline of order >3 */
void init_splinen(float *a,int n)
{
  int k;

  a[0] = 1.;
  for (k=2;k<=n;k++) a[0]/=(float)k;
  for (k=1;k<=n+1;k++)
    a[k] = - a[k-1] *(float)(n+2-k)/(float)k;
}

/* pre-computation for spline of order >3 */
static void init_splinen(double* a, int n)
{
    a[0] = 1.;
    for(int k=2; k <= n; k++)
        a[0]/=(double)k;
    for(int k=1; k <= n+1; k++)
        a[k] = - a[k-1] *(n+2-k)/k;
}

/* fast integral power function */
static float ipow(float x, int n)
{
    float res;
    for(res = 1.; n; n>>=1) {
        if(n&1)
            res *= x;
        x *= x;
    }
    return res;
}

/* coefficients for spline of order >3 */
void splinen(float *c,float t,float *a,int n)
{
  int i,k;
  float xn;
  
  memset((void *)c,0,(n+1)*sizeof(float));
  for (k=0;k<=n+1;k++) { 
    xn = ipow(t+(float)k,n);
    for (i=k;i<=n;i++) 
      c[i] += a[i-k]*xn;
  }
}

/* coefficients for spline of order >3 */
static void splinen(double* c, double t, double* a, int n)
{
    memset((void*)c, 0, (n+1)*sizeof(double));
    for(int k=0; k <= n+1; k++) { 
        double xn = ipow(t+(double)k, n);
        for(int i=k; i <= n; i++) 
            c[i] += a[i-k]*xn;
    }
}

/// Spline interpolation of given \a order of image \a im at point (x,y).
/// \a out must be an array of size the number of components.
/// Supported orders: 0(replication), 1(bilinear), -3(Keys's bicubic), 3, 5, 7,
/// 9, 11.
/// \a paramKeys is Keys's parameter, only used for order -3.
/// Success means a valid order and pixel in image.
bool interpolate_spline( image_double& im, int order, double x, double y, double& out, double paramKeys)
{
	double  cx[12],cy[12];

	/* CHECK ORDER */
	if(order != 0 && order != 1 && order != -3 &&
		order != 3 && order != 5 && order != 7 && order != 9 && order != 11)
		return false;

	double ak[13];
	if(order > 3)
		init_splinen(ak, order);

	bool bInside = false;
	/* INTERPOLATION */
	if(order == 0) { /* zero order interpolation (pixel replication) */
		int xi = (int)floor((double)x);
		int yi = (int)floor((double)y);
		bInside = valid_image_double(im, xi, yi);//im.valid(xi, yi);
		if(bInside) {
			double p = im->data[xi+yi*im->xsize]; //im.pixel(xi, yi);
			out = p;
		}
	} else { /* higher order interpolations */
		bInside = (x>=0.0f && x<=(double)im->xsize && y>=0.0f && y<=(double)im->ysize);
		if(bInside) {
			x -= 0.5f; y -= 0.5f;
			int xi = (x<0)? -1: (int)x;
			int yi = (y<0)? -1: (int)y;
			double ux = x - (double)xi;
			double uy = y - (double)yi;
			switch(order)  {
			case 1: /* first order interpolation (bilinear) */
				cx[0] = ux; cx[1] = 1.0f-ux;
				cy[0] = uy; cy[1] = 1.0f-uy;
				break;
			case -3: /* third order interpolation (bicubic Keys' function) */
				keys(cx, ux, paramKeys);
				keys(cy, uy, paramKeys);
				break;
			case 3: /* spline of order 3 */
				spline3(cx, ux);
				spline3(cy, uy);
				break;
			default: /* spline of order >3 */
				splinen(cx, ux, ak, order);
				splinen(cy, uy, ak, order);
				break;
			}
			int n2 = (order==-3)? 2: (order+1)/2;
			int n1 = 1-n2;
			/* this test saves computation time */
			if(valid_image_double(im, xi+n1, yi+n1) && valid_image_double(im, xi+n2,yi+n2)) {
				out = 0.0f;
				for(int dy = n1; dy <= n2; dy++) {
					int adrs = (xi+n1) + (yi+dy) * im->xsize;
					for(int dx = n1; dx <= n2; dx++) {
						double f = im->data[adrs];
						out += cy[n2-dy]*cx[n2-dx] * f;
						adrs++;
					}
				}

			} else
				out = 0.0f;
			for(int dy = n1; dy <= n2; dy++)
				for(int dx = n1; dx <= n2; dx++) {
					double v = 0.0f; // the image is not infinite, there is no data outside
					out += cy[n2-dy]*cx[n2-dx]*v;
				}

		}
	}
	return bInside;
}

/*------------------------------ MAIN MODULE ------------------------------*/

void finvspline(float *in,int order,float *out, int width, int height)
{
  double *c,*d,z[5];
  int npoles,nx,ny,x,y;
 
  ny = height; nx = width;

  /* initialize poles of associated z-filter */
  switch (order) 
    {
    case 2: z[0]=-0.17157288;  /* sqrt(8)-3 */
      break;

    case 3: z[0]=-0.26794919;  /* sqrt(3)-2 */ 
      break;

    case 4: z[0]=-0.361341; z[1]=-0.0137254;
      break;

    case 5: z[0]=-0.430575; z[1]=-0.0430963;
      break;
      
    case 6: z[0]=-0.488295; z[1]=-0.0816793; z[2]=-0.00141415;
      break;

    case 7: z[0]=-0.53528; z[1]=-0.122555; z[2]=-0.00914869;
      break;
      
    case 8: z[0]=-0.574687; z[1]=-0.163035; z[2]=-0.0236323; z[3]=-0.000153821;
      break;

    case 9: z[0]=-0.607997; z[1]=-0.201751; z[2]=-0.0432226; z[3]=-0.00212131;
      break;
      
    case 10: z[0]=-0.636551; z[1]=-0.238183; z[2]=-0.065727; z[3]=-0.00752819;
      z[4]=-0.0000169828;
      break;
      
    case 11: z[0]=-0.661266; z[1]=-0.27218; z[2]=-0.0897596; z[3]=-0.0166696; 
      z[4]=-0.000510558;
      break;
      
     default:
      printf("finvspline: order should be in 2..11.\n");
      exit(-1);
    }

  npoles = order/2;

  /* initialize double array containing image */
  c = (double *)malloc(nx*ny*sizeof(double));
  d = (double *)malloc(nx*ny*sizeof(double));
  for (x=nx*ny;x--;) 
    c[x] = (double)in[x];

  /* apply filter on lines */
  for (y=0;y<ny;y++) 
    invspline1D(c+y*nx,nx,z,npoles);

  /* transpose */
  for (x=0;x<nx;x++)
    for (y=0;y<ny;y++) 
      d[x*ny+y] = c[y*nx+x];
      
  /* apply filter on columns */
  for (x=0;x<nx;x++) 
    invspline1D(d+x*ny,ny,z,npoles);

  /* transpose directy into image */
  for (x=0;x<nx;x++)
    for (y=0;y<ny;y++) 
      out[y*nx+x] = (float)(d[x*ny+y]);

  /* free array */
  free(d);
  free(c);
}

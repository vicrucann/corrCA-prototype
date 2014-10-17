/* LMTacheC derived class from LMA base class
 *
Copyright (C) 2014 Victoria Rudakova <vicrucann@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CENTERS_H
#define CENTERS_H

#include "LM.h"
#include "pgm_io.h"

using namespace libNumerics;

template <typename T>
T model_luminance_alternative(int x, int y, int wi, int he, vector<T> P, bool clr)
{
	T err = 0, lum = 0;
	T dist = -512;
	T l1=P[0], l2=P[1], theta=P[2], tu=P[3], tv=P[4], r1=P[5], p=P[6], r2=P[7], vh=P[8], vb=P[9], c=P[10];
	if (x >= 1 && x <= wi-2 && y >= 1 && y <= he-2)
	{
		T coordX = l1 * (x-tu)*cos(theta) - l1 * (y-tv)*sin(theta);
		T coordY = l2 * (x-tu)*sin(theta) + l2 * (y-tv)*cos(theta);
		T d = coordX*coordX + coordY*coordY;
		if (d < 2.0*2.0)
			dist = std::sqrt(d);
	}

	T d = dist;
	if (d == -512) // does not belong to circle - return the 'not belongin color'.
	{
		lum = 255;
		vh = 255;
		vb = 0;
		if (!clr)
			lum = vb + lum;
		else 
		lum = vh - lum;	
		return lum;
	}
	if (!clr)
		d = -d+c;
	else
		d = d-c;
	r1 = 0;//abs(r1);
	r2 = 0; //abs(r2);
	T y20 = -1.0;
	T x21 = y20/p;
	T x20 = x21;
	if (d >= x20)
		lum = 0;
	else
	{
		T y10 = 1.0;
		T x11 = y10/p;
		T x10 = x11;
		if (d <= x10)
			lum = 1;
		else
		{
			lum = p*d;
			lum = 0.5*lum+0.5;
		}
	}
	
	lum = lum * (P[8] - P[9]);
	if (!clr)
		lum = P[9] + lum;
	else 
		lum = P[8] - lum;
	return lum;
}

template <typename T>
void jacobian_alternative(matrix<T>& J, const vector<T>& P, int wi, int he, bool clr, int xbegin, int xend, int ybegin, int yend)
{
	T l1=P[0], l2=P[1], theta=P[2], tu=P[3], tv=P[4], r1=P[5], p=P[6], r2=P[7], vh=P[8], vb=P[9], c=P[10];
	int idx = 0;
	for (int v = ybegin; v <= yend; v++) {
		for (int u = xbegin; u <= xend; u++) {

			for (int k = 0; k < P.size(); k++) // initialize
				J(idx, k) = 0;

			T x = 0, y = 0, d = -512;
			if (u >= 1 && u <= wi-2 && v >= 1 && v <= he-2) {
				x = l1 * (u-tu)*cos(theta) - l1 * (v-tv)*sin(theta);
				y = l2 * (u-tu)*sin(theta) + l2 * (v-tv)*cos(theta);
				T dist = x*x + y*y;
				if (dist <= 2.0*2.0)
					d = sqrt(dist);
			}

			if (d != -512) {
				if (!clr) d = -d+c;
				else d = d-c;
				r1 = 0;//abs(r1);
				r2 = 0;//abs(r2);

				T y20 = -1.0;
				T x21 = y20/p;
				T x20 = x21;

				T y10 = 1.0;
				T x11 = y10/p;
				T x10 = x11;

				T s2=0;
				if(d>=x20)
					s2 = 0.0;
				else if(d<=x10)
					s2 = 1.0;
				else //if (d<=x21)
				{
					s2=p*d;
					s2 = 0.5*s2+0.5;}
				if (d < x20 && d > x10) {
					T dXdL1 = (u-tu)*cos(theta) - (v-tv)*sin(theta); // "-" - for black taches, because d = -d+1!
					T dYdL2 = (u-tu)*sin(theta) + (v-tv)*cos(theta);
					T dDdL1 = -x*dXdL1 / sqrt(x*x + y*y);
					T dDdL2 = -y*dYdL2 / sqrt(x*x + y*y);
					T dXdTheta = -l1*(sin(theta)*(u-tu) + cos(theta)*(v-tv));
					T dYdTheta = l2*(cos(theta)*(u-tu) - sin(theta)*(v-tv));
					T dDdTheta = -(x*dXdTheta + y*dYdTheta) / sqrt(x*x + y*y);
					T dXdTu = -l1*cos(theta);
					T dYdTu = -l2*sin(theta);
					T dDdTu = -(x*dXdTu + y*dYdTu) / sqrt(x*x + y*y);
					T dXdTv = l1*sin(theta);
					T dYdTv = -l2*cos(theta);
					T dDdTv = -(x*dXdTv + y*dYdTv) / sqrt(x*x + y*y);
					T dDdC = 1;
					J(idx, 0) = p*dDdL1;
					J(idx, 1) = p*dDdL2;
					J(idx, 2) = p*dDdTheta;
					J(idx, 3) = p*dDdTu;
					J(idx, 4) = p*dDdTv;
					//J(idx, 5) = 0;
					J(idx, 6) = d;
					//J(idx, 7) = 0;
					J(idx, 10) = p*dDdC;
					for (int k = 0; k < P.size(); k++)
						J(idx,k) *= 0.5*(vh-vb);
				}
				J(idx, 8) = s2;
				J(idx, 9) = 1-s2;
			}
			idx++;
		}
	}
	//J *= 0.5*(vh-vb);
}

template <typename T>
T model_luminance(int x, int y, int wi, int he, vector<T> P, bool clr)
{
	T err = 0, lum = 0;
	T dist = -512;
	T l1=P[0], l2=P[1], theta=P[2], tu=P[3], tv=P[4], r1=P[5], p=P[6], r2=P[7], vh=P[8], vb=P[9], c=P[10];
	if (x >= 1 && x <= wi-2 && y >= 1 && y <= he-2)
	{
		T coordX = l1 * ((x-tu)*cos(theta) - (y-tv)*sin(theta));
		T coordY = l2 * ((x-tu)*sin(theta) + (y-tv)*cos(theta));
		T d = coordX*coordX + coordY*coordY;
		if (d < 2.0*2.0)
			dist = std::sqrt(d);
	}

	T d = dist;
	if (d == -512) // does not belong to circle - return the 'not belongin color'.
	{
		lum = 255;
		vh = 255;
		vb = 0;
		if (!clr)
			lum = vb + lum;
		else 
		lum = vh - lum;	
		return lum;
	}
	if (!clr)
		d = -d+c;
	else
		d = d-c;
	r1 = abs(r1);
	r2 = abs(r2);
	T val = 1.0 / sqrt(p*p+1.0);
	T y20 = -(1.0-r2);
	T x21 = (y20-r2*val)/p;
	T x20 = x21-r2*p*val;
	if (d >= x20)
		lum = 0;
	else
	{
		T y10 = 1.0-r1;
		T x11 = (y10+r1*val)/p;
		T x10 = x11+r1*p*val;
		if (d <= x10)
			lum = 1;
		else
		{
			if (d <= x11)
				lum = y10+std::sqrt(r1*r1-(d-x10)*(d-x10));
			else if (d <= x21)
				lum = p*d;
			else
				lum = (y20-std::sqrt(r2*r2-(d-x20)*(d-x20)));
			lum = 0.5*lum+0.5;
		}
	}
	
	lum = lum * (P[8] - P[9]);
	if (!clr)
		lum = P[9] + lum;
	else 
		lum = P[8] - lum;
	return lum;
}

template <typename T>
void jacobian(matrix<T>& J, const vector<T>& P, int wi, int he, bool clr, int xbegin, int xend, int ybegin, int yend)
{
	T l1=P[0], l2=P[1], theta=P[2], tu=P[3], tv=P[4], r1=P[5], p=P[6], r2=P[7], vh=P[8], vb=P[9], c=P[10];
	int idx = 0;
	for (int v = ybegin; v <= yend; v++) {
		for (int u = xbegin; u <= xend; u++) {

			for (int k = 0; k < P.size(); k++) // initialize
				J(idx, k) = 0;

			T x = 0, y = 0, d = -512;
			if (u >= 1 && u <= wi-2 && v >= 1 && v <= he-2) {
				x = l1*(u-tu)*cos(theta) - l1*(v-tv)*sin(theta);
				y = l2*(u-tu)*sin(theta) + l2*(v-tv)*cos(theta);
				T dist = x*x + y*y;
				if (dist <= 2.0*2.0)
					d = sqrt(dist);
			}

			if (d != -512) {
				if (!clr) d = -d+c;
				else d = d-c;
				r1 = abs(r1);
				r2 = abs(r2);
				T val = 1.0 / std::sqrt(p*p+1.0);

				T y20 = -(1.0-r2);
				T x21 = (y20-r2*val)/p;
				T x20 = x21-r2*p*val;

				T y10 = 1.0-r1;
				T x11 = (y10+r1*val)/p;
				T x10 = x11+r1*p*val;

				T s2=0;
				if(d>=x20)
					s2 = 0.0;
				else if(d<=x10)
					s2 = 1.0;
				else if (d<=x11) {
					s2=y10+sqrt(r1*r1-(d-x10)*(d-x10));
					s2 = 0.5*s2+0.5;}
				else if (d<=x21){
					s2=p*d;
					s2 = 0.5*s2+0.5;}
				else{
					s2=(y20-sqrt(r2*r2-(d-x20)*(d-x20)));
					s2 = 0.5*s2+0.5;}

				if (d < x20 && d > x10) {
					T dXdL1 = (u-tu)*cos(theta) - (v-tv)*sin(theta);
					T dYdL2 = (u-tu)*sin(theta) + (v-tv)*cos(theta);
					T dDdL1 = -x*dXdL1 / sqrt(x*x + y*y);
					T dDdL2 = -y*dYdL2 / sqrt(x*x + y*y);
					T dXdTheta = -l1*(sin(theta)*(u-tu) + cos(theta)*(v-tv));
					T dYdTheta = l2*(cos(theta)*(u-tu) - sin(theta)*(v-tv));
					T dDdTheta = -(x*dXdTheta + y*dYdTheta) / sqrt(x*x + y*y);
					T dXdTu = -l1*cos(theta);
					T dYdTu = -l2*sin(theta);
					T dDdTu = -(x*dXdTu + y*dYdTu) / sqrt(x*x + y*y);
					T dXdTv = l1*sin(theta);
					T dYdTv = -l2*cos(theta);
					T dDdTv = -(x*dXdTv + y*dYdTv) / sqrt(x*x + y*y);
					T dDdC = 1;
					if (d <= x21 && d >= x11) {
						J(idx, 0) = p*dDdL1;
						J(idx, 1) = p*dDdL2;
						J(idx, 2) = p*dDdTheta;
						J(idx, 3) = p*dDdTu;
						J(idx, 4) = p*dDdTv;
						//J(idx, 5) = 0;
						J(idx, 6) = d;
						//J(idx, 7) = 0;
						J(idx, 10) = p*dDdC;
					}
					else if (d < x11) {
						J(idx, 0) = -dDdL1 * (d-x10) / sqrt(r1*r1-(d-x10)*(d-x10));
						J(idx, 1) = -dDdL2 * (d-x10) / sqrt(r1*r1-(d-x10)*(d-x10));
						J(idx, 2) = -dDdTheta * (d-x10) / sqrt(r1*r1-(d-x10)*(d-x10));
						J(idx, 3) = -dDdTu * (d-x10) / sqrt(r1*r1-(d-x10)*(d-x10));
						J(idx, 4) = -dDdTv * (d-x10) / sqrt(r1*r1-(d-x10)*(d-x10));
						T dX11dR1 = 1/p * (-1+1/sqrt(1+p*p));
						T dX10dR1 = dX11dR1 + p/sqrt(1+p*p);
						J(idx, 5) = -1 + 1/(sqrt(r1*r1 - (d-x10)*(d-x10))) * (r1 + (d-x10)*dX10dR1);
						T dX10dP = -1/(p*p) * (-r1+1+r1*sqrt(1+p*p)) + r1/sqrt(1+p*p);
						J(idx, 6) = 1/(sqrt(r1*r1-(d-x10)*(d-x10))) * (d-x10) * dX10dP;
						//J(idx, 7) = 0;
						J(idx,10) = -dDdC * (d-x10) / sqrt(r1*r1-(d-x10)*(d-x10));
					}
					else {
						J(idx, 0) = dDdL1 * (d-x20) / sqrt(r2*r2-(d-x20)*(d-x20));
						J(idx, 1) = dDdL2 * (d-x20) / sqrt(r2*r2-(d-x20)*(d-x20));
						J(idx, 2) = dDdTheta * (d-x20) / sqrt(r2*r2-(d-x20)*(d-x20));
						J(idx, 3) = dDdTu * (d-x20) / sqrt(r2*r2-(d-x20)*(d-x20));
						J(idx, 4) = dDdTv * (d-x20) / sqrt(r2*r2-(d-x20)*(d-x20));
						//J(idx, 5) = 0;
						T dX20dP = 1/(p*p) * (-r2+1+r2*sqrt(1+p*p)) - r2/sqrt(1+p*p);
						J(idx, 6) = -1/(sqrt(r2*r2-(d-x20)*(d-x20))) * (d-x20) * dX20dP;
						T dX21dR2 = 1/p*(1-1/sqrt(1+p*p));
						T dX20dR2 = dX21dR2 - p*1/sqrt(1+p*p);
						J(idx, 7) = 1 - 1/(sqrt(r2*r2 - (d-x20)*(d-x20))) * (r2 + (d-x20) * dX20dR2);
						J(idx, 10) = dDdC * (d-x20) / sqrt(r2*r2 - (d-x20)*(d-x20));
					}
					for (int k = 0; k < P.size(); k++)
						J(idx,k) *= 0.5*(vh-vb);
				}
				J(idx, 8) = s2;
				J(idx, 9) = 1-s2;
			}
			idx++;
		}
	}
	//J *= 0.5*(vh-vb);
}

template <typename T>
class LMTacheC : public MinLM<T> {

public:
	LMTacheC(image_double Im, T cx, T cy, T delta, bool tache_color, int img_width, int img_height) {
		clr = tache_color; 
		wi = img_width;
		he = img_height;
		xbegin = cx+0.5-delta;
		xend = cx+0.5+delta;
		ybegin = cy+0.5-delta;
		yend = cy+0.5+delta;
		im = Im;
	}

private:
	bool clr;
	int wi, he;
	int xbegin, xend, ybegin, yend;
	image_double im;
	
public:
	virtual void modelData(const vector<T>& P, vector<T>& ymodel) const
	{
		int idx = 0;
		for (int i = ybegin; i <= yend; i++) {
			for (int j = xbegin; j <= xend; j++) {
				ymodel[idx] = model_luminance_alternative(j, i, wi, he, P, clr); 
				idx++;
			} }	
	}
	
	virtual void modelJacobian(const vector<T>& P, matrix<T>& J) const
	{
		jacobian_alternative(J, P, wi, he, clr, xbegin, xend, ybegin, yend);
	}
};


#endif

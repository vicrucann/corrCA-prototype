/* Chromatic aberration correction prototype.
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

#include "centers.h"
#include "abberation.h"
#include "distortion.h"
#include "pgm_io.h"
#include "spline.h"
#include "correction.h"

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

image_double average_image(image_double img) {
	int w = img->xsize; 
	int h = img->ysize;
	image_double img_avg = new_image_double_copy(img);
	for (int v = 1; v < h-1; v++) {
		for (int u = 1; u < w-1; u++) {
			double pix = img->data[u-1+(v-1)*w] + img->data[u+(v-1)*w] + img->data[u+1+(v-1)*w] + 
				img->data[u-1+v*w] + img->data[u+v*w] + img->data[u+1+v*w] + 
				img->data[u-1+(v+1)*w] + img->data[u+(v+1)*w] + img->data[u+1+(v+1)*w] ;
			img_avg->data[u+v*w] = pix/9; 
		}
	}
	return img_avg;
}

template <typename T>
int initial_tache(image_double I, vector<T>& h, T& rayon, bool color, T x, T y) {
	int COL_IMA = I->xsize; 
	int LIG_IMA = I->ysize;
	int j = x;
	int i = y;
	int d = 2*rayon;
	if (2*d+1 > LIG_IMA) 
		d=(LIG_IMA-1)/2;
	if(2*d+1 > COL_IMA) 
		d=(COL_IMA-1)/2;
	if(i<d)
		i=d+1;
	if(i>LIG_IMA-1-d) 
		i=LIG_IMA-2-d;
	if(j<d)
		j=d+1;
	if(j>COL_IMA-1-d) 
		j=COL_IMA-2-d;
	int val_haut=0;
	int val_bas=255;
	for (int k = -d; k <= d; k++) {
		for (int l = -d;  l <= d; l++) {
			T lum = I->data[j+l+(i+k)*COL_IMA];
			if (lum > val_haut)
				val_haut = lum;
			else if (lum<val_bas)
				val_bas=lum;
		} }
	T seuil = 0;
	if (!color)
		seuil = val_bas + (val_haut - val_bas)/3 * 2;
	else if (color)
		seuil = val_bas + (val_haut - val_bas)/3;
	
	matrix<T> tab = matrix<T>::zeros(2*d+1, 2*d+1);
	int label = 1;
	
	for (int k = -d+1; k <= d; k++){
		for(int l = -d+1; l <= d-1; l++){
			T lum= I->data[j+l+(i+k)*COL_IMA];
			if (lum < seuil) {
				int imin=l;
				int imax=l+1;          
				while( (I->data[j+imax+(i+k)*COL_IMA] <= seuil) && (imax <= d-1) ){
					imax++; }
				int vallab=0;
				for(int m = imin; m <= imax; m++){
					if (tab(k+d-1,m+d) != 0){
						vallab = tab(k+d-1,m+d);
					}
				}
				if (vallab == 0){
					vallab=label; 
					label++;
				}
				for(int m = imin; m <= imax; m++){
					tab(k+d,m+d)=vallab;
				}
				l=imax;
			}
		}
	}
	matrix<T> bary = matrix<T>::zeros(label, 4);
	for(int k = -d; k <= d; k++){
		for(int l = -d; l <= d; l++){
			if(tab(k+d,l+d)!=0){
				T lum= I->data[j+l+(i+k)*COL_IMA];
				bary(tab(k+d,l+d),0)+=(255-lum)*(j+l);
				bary(tab(k+d,l+d),1)+=(255-lum)*(i+k);
				bary(tab(k+d,l+d),2)+=(255-lum);
				bary(tab(k+d,l+d),3)++; 
			}
		}
	}
	int distmin=100;
	int labelmin=0;
	for(int k = 1; k < label; k++){
		T dist = std::sqrt( (bary(k,0)/bary(k,2)-x) * (bary(k,0)/bary(k,2)-x)+
			(bary(k,1)/bary(k,2)-y) * (bary(k,1)/bary(k,2)-y));
		if(dist < distmin && bary(k,3) > 25 ){ /* 25 =  surface min*/
			distmin=dist;
			labelmin=k;
		}      
	}
	if(labelmin == 0) {
		printf("pb tache trop petite (<=25 pixels)\n");
		return 1;
	}
	x=bary(labelmin,0)/bary(labelmin,2);
	y=bary(labelmin,1)/bary(labelmin,2);
	T sx2 = 0, sy2 = 0, sxy = 0, ss = 0;
	for(int k = -d; k <= d; k++){
		for(int l = -d; l <= d; l++){
			if(tab(k+d,l+d) == labelmin){
				sx2+=(j+l-x)*(j+l-x);
				sy2+=(i+k-y)*(i+k-y);
				sxy+=(i+k-y)*(j+l-x);
				ss++;
			}
		}
	}
	T lambda1 =  ((sx2+sy2)/ss + std::sqrt(( (sx2+sy2)*(sx2+sy2)+4*(sxy*sxy-sx2*sy2)))/ss)/2.0; 
	T lambda2 =  ((sx2*sy2-sxy*sxy)/(ss*ss))/lambda1;
	rayon=std::sqrt(lambda1)*2;
	h[0] = 1.0/rayon; 	                        /* lambda1 		*/
	h[1] = (std::sqrt(lambda1/lambda2))/rayon; 	/* lambda2		*/
	h[2] = std::atan2(sx2/ss-lambda1, -sxy/ss);    	/* alpha  	 	*/
	h[3] = x;        	/* tu     		*/
	h[4] = y;        	/* tv      		*/
	h[5] = 0.25;       	/* rayon cercle 1   	*/
	h[6] = -2.0;      	/* pente 	      	*/
	h[7] = 0.25;       	/* rayon cercle 2   	*/
	h[8] = val_haut; 	/* val_haut   		*/
	h[9] = val_bas; 	/* val_bas   		*/
	h[10] = 1.0; 	/* position step  	*/
	return 0;
}

template <typename T>
vector<T> trgtDataCalc(image_double img_avg, T cx, T cy, T delta) {
	int xbegin = cx+0.5-delta;
	int xend = cx+0.5+delta;
	int ybegin = cy+0.5-delta;
	int yend = cy+0.5+delta;
	int nerr = (yend-ybegin+1)*(xend-xbegin+1);
	vector<T> trgData = vector<T>::zeros(nerr);
	int wi = img_avg->xsize;
	int he = img_avg->ysize;
	int idx = 0;
	for (int v = ybegin; v <= yend; v++) {
		for (int u = xbegin; u <= xend; u++) {
			if (u >= 1 && u <= wi-2 && v >= 1 && v <= he-2)
				trgData[idx] = img_avg->data[u+v*img_avg->xsize]; 
			else
				trgData[idx] = 255; // for 'tache noire'
			idx++;
		}
	}
	return trgData;
}

template <typename T>
T centerLMA(image_double sub_img, bool clr, T& centerX, T& centerY)
{
	image_double img_avg = average_image(sub_img);
	int w = sub_img->xsize; 
	int h = sub_img->ysize;
	T cx = w/2, cy = h/2, radi = 0.4*w;
	vector<T> P(11);
	initial_tache(sub_img, P, radi, clr, cx, cy);
	vector<T> trgData = trgtDataCalc<T>(img_avg, P[3], P[4], radi*2);
	LMTacheC<T> ellipseLMA(img_avg, P[3], P[4], radi*2, clr, w, h);
	T rmse = ellipseLMA.minimize(P, trgData, 0.001);
	free_image_double(img_avg);
    //T lambda1 = P[0]; T lambda2 = P[1]; T theta = P[2];
	centerX = P[3];
	centerY = P[4];
	return rmse;
}

template <typename T>
image_double takeSubImg(image_double IMG, T cx, T cy, T radi, int& x0, int& y0)
{
	int size = 2.5 * radi;
	int x1 = cx - 0.5*size, x2 = cx + 0.5*size; 
	int y1 = cy - 0.5*size, y2 = cy + 0.5*size;
	x0 = x1;
	y0 = y1;
	if (y2-y1 != x2-x1)
		y2 = y1+x2-x1;
	image_double img = new_image_double(x2-x1, y2-y1);
	for (int i = 0; i < img->xsize; i++) {
		for (int j = 0; j < img->ysize; j++) 
			if (x1+i >= 0 && y1+j >= 0 && x1+i<IMG->xsize && y1+j<IMG->ysize)
				img->data[i+j*img->xsize] = IMG->data[x1+i+(y1+j)*IMG->xsize];
			else
				img->data[i+j*img->xsize] = 255;
	}
	return img;
}

template <typename T>
bool loadKeypts(const char* fname, std::vector<CCStats>& cc_green, std::vector<CCStats>& cc_red)
{
	std::ifstream f(fname);
    while( f.good() ) {
        std::string str;
        std::getline(f, str);
        if( f.good() ) {
            std::istringstream s(str);
            CCStats red, green;
			s >> green.centerX >> green.centerY >> red.centerX >> red.centerY;
            if(!s.fail() )
			{
                cc_green.push_back(green);
				cc_red.push_back(red);
			}
        }
    }
	
	return true;
}

template <typename T>
image_double image_rotate_left(image_double img)
{
	image_double res = new_image_double_ini(img->ysize, img->xsize, 0);
	int idx = 0;
	for (int x = img->xsize-1; x >= 0; x--)
	{
		for (int y = 0; y < img->ysize; y++)	
		{
			res->data[idx] = img->data[x+y*img->xsize];
			idx++;
		}
	}
	return res;
}

template <typename T>
image_double image_rotate_right(image_double img)
{
	image_double res = new_image_double_ini(img->ysize, img->xsize, 0);
	int idx = 0;
	for (int x = 0; x < img->xsize; x++)
	{
		for (int y = img->ysize-1; y >= 0; y--)	
		{
			res->data[idx] = img->data[x+y*img->xsize];
			idx++;
		}
	}
	return res;
}

template <typename T>
void binarization(image_double& imgbiR, image_double& imgbiG, image_double& imgbiB,
	image_double& imgR, image_double& imgG, image_double& imgB,
	T threR, T threG, T threB)
{
	int wiRB = imgR->xsize;
	int heRB = imgR->ysize;
	int wiG = imgG->xsize;
	int heG = imgG->ysize;
	T red, blue, green;
	for (int i = 0; i < wiRB; i++) {
		for (int j = 0; j < heRB; j++) {
			red = imgR->data[i+j*wiRB];
			if (red <= threR) imgbiR->data[i+j*wiRB] = 0;

			blue = imgB->data[i+j*wiRB];
			if (blue <= threB) imgbiB->data[i+j*wiRB] = 0;

			
			if (wiG == wiRB && heG == heRB) {
				green = imgG->data[i+j*wiRB];
				if (green <= threG) imgbiG->data[i+j*wiRB] = 0;	}
			else {
				green = imgG->data[i*2+1+j*2*imgbiG->xsize];
				if (green <= threG) imgbiG->data[i*2+1+j*2*imgbiG->xsize] = 0;
				green = imgG->data[i*2+j*2*imgbiG->xsize];
				if (green <= threG) imgbiG->data[i*2+j*2*imgbiG->xsize] = 0;
				green = imgG->data[i*2+(j*2+1)*imgbiG->xsize];
				if (green <= threG) imgbiG->data[i*2+(j*2+1)*imgbiG->xsize] = 0;
				green = imgG->data[i*2+1+(j*2+1)*imgbiG->xsize];
				if (green <= threG) imgbiG->data[i*2+1+(j*2+1)*imgbiG->xsize] = 0; }
		}
	}
}

template <typename T>
void img_extremas(image_double& img, T& min, T& max) {
	min = 255;
	max = 0;
	for (int i = 4; i < img->xsize-4; i++) {
		for (int j = 4; j < img->ysize-4; j++) {
			T clr = img->data[j*img->xsize + i];
			if (clr > max) max = clr;
			if (clr < min) min = clr; }
	}
}


template <typename T>
void circle_redefine(image_double& imgR, image_double& imgG, image_double& imgB, 
	vector<T>& xR, vector<T>& yR, vector<T>& rR, 
	vector<T>& xGr, vector<T>& yGr, vector<T>& rG, 
	vector<T>& xB, vector<T>& yB, vector<T>& rB,
	vector<T>& xGb, vector<T>& yGb,
    int scale, bool clr, bool green = true)
{
	printf("\nLMA center redefinition for the channels... \n");
	int ntaches = xR.size();
	for (int i = 0; i < ntaches; i++) {
		int x0R, y0R, x0G, y0G, x0B, y0B;
		image_double sub_imgR = takeSubImg(imgR, xR[i], yR[i], rR[i], x0R, y0R);
		image_double sub_imgB = takeSubImg(imgB, xB[i], yB[i], rB[i], x0B, y0B);

		T cxR=0, cyR=0, cxB=0, cyB=0;
		centerLMA<T>(sub_imgR, clr, cxR, cyR);
		centerLMA<T>(sub_imgB, clr, cxB, cyB);

		free_image_double(sub_imgR);
		free_image_double(sub_imgB);
	
		xR[i] = scale*(x0R + cxR); 
		yR[i] = scale*(y0R + cyR);
		xB[i] = scale*(x0B + cxB); 
		yB[i] = scale*(y0B + cyB);

		if (green) {
			image_double sub_imgG = takeSubImg(imgG, xGr[i], yGr[i], rG[i], x0G, y0G);
			T cxG=0, cyG=0;
			centerLMA<T>(sub_imgG, clr, cxG, cyG);
			free_image_double(sub_imgG);
			xGr[i] = x0G + cxG; 
			yGr[i] = y0G + cyG;
			xGb[i] = x0G + cxG; 
			yGb[i] = y0G + cyG; }

		double percent = ((double)i / (double)ntaches)*100;
		if (!(i % (int)(0.2*ntaches))) printf("%i%c", (int)percent+1, '%');
		else if (!(i % (int)(0.04*ntaches))) printf(".");
	}
}

template <typename T>
void keypnts_circle(image_double& imgR, image_double& imgG, image_double& imgB,
	vector<T>& xR, vector<T>& yR, vector<T>& rR,
	vector<T>& xGr, vector<T>& yGr, vector<T>& rG, 
	vector<T>& xB, vector<T>& yB, vector<T>& rB,
	vector<T>& xGb, vector<T>& yGb,
	T scale, bool clr)
{
	int wiRB = imgR->xsize, heRB = imgR->ysize;
	int wiG = wiRB*scale, heG = heRB*scale;

	T maxR = 0, maxG = 0, maxB = 0;
	T minR = 255, minG = 255, minB = 255;
	img_extremas(imgR, minR, maxR);
	img_extremas(imgG, minG, maxG);
	img_extremas(imgB, minB, maxB);
	T threR = 0.5*(maxR-minR);
	T threG = 0.54*(maxG-minG);
	T threB = 0.55*(maxB-minB);

	image_double imgbiR = new_image_double_ini(wiRB, heRB, 255);
	image_double imgbiG = new_image_double_ini(wiG, heG, 255);
	image_double imgbiB = new_image_double_ini(wiRB, heRB, 255);
	binarization(imgbiR, imgbiG, imgbiB, imgR, imgG, imgB, threR, threG, threB);
	//write_pgm_image_double(imgbiB, "rawdata/b.pgm");
	
	printf("finding connected components... \n");
	std::vector<CCStats> ccstatsR, ccstatsG, ccstatsB;
	CC(ccstatsR, imgbiR, 'R'); printf("number = %i ", ccstatsR.size());
	CC(ccstatsG, imgbiG, 'G'); printf("number = %i ", ccstatsG.size());
	CC(ccstatsB, imgbiB, 'B'); printf("number = %i ", ccstatsB.size());
	
	printf("\nnumber of connected components per channel R=%i G=%i B=%i\n", ccstatsR.size(), ccstatsG.size(), ccstatsB.size());
	assert(ccstatsR.size() == ccstatsG.size() && ccstatsG.size() == ccstatsB.size());
	printf("centers initialization for channels is done \n");

	printf("\nMatching the centers... ");
	int ntaches = ccstatsG.size();
	xR = xR.ones(ntaches); yR = yR.ones(ntaches); rR = rR.ones(ntaches);
	xB = xR.ones(ntaches); yB = yR.ones(ntaches); rB = rB.ones(ntaches);
	xGr = xGr.ones(ntaches); yGr = yGr.ones(ntaches); rG = rG.ones(ntaches);
	xGb = xGb.ones(ntaches); yGb = yGb.ones(ntaches); 
	xR = -1; xB = -1;
	yR = -1; yB = -1;
	for (int i = 0; i < ntaches; i++) {
		T xg = ccstatsG[i].centerX;
		T yg = ccstatsG[i].centerY;
		int idxR = findMatch(xg, yg, ccstatsR, scale);
		int idxB = findMatch(xg, yg, ccstatsB, scale);
		xGr[i] = xg; yGr[i] = yg; 
		xGb[i] = xg; yGb[i] = yg; 
		rG[i] = 0.5*(ccstatsG[i].radius1+ccstatsG[i].radius2);
		if (idxR != -1) { 
			xR[i] = ccstatsR[idxR].centerX; 
			yR[i] = ccstatsR[idxR].centerY; 
			rR[i] = 0.5*(ccstatsR[idxR].radius1+ccstatsR[idxR].radius2); }
		if (idxB != -1) { 
			xB[i] = ccstatsB[idxB].centerX; 
			yB[i] = ccstatsB[idxB].centerY; 
			rB[i] = 0.5*(ccstatsB[idxB].radius1+ccstatsB[idxB].radius2); }
	}
	printf("done.\n");
	circle_redefine(imgR, imgG, imgB, xR, yR, rR, xGr, yGr, rG, xB, yB, rB, xGb, yGb, scale, clr);
	free_image_double(imgbiR); 
	free_image_double(imgbiG); 
	free_image_double(imgbiB);
}

template <typename T>
void raw2rgb(image_double& img_bayer, image_double& imgR, image_double& imgG, image_double& imgB)
{
	printf("raw data processing... ");
	int wiRB = imgR->xsize;
	int heRB = imgR->ysize;
	T red, blue, green;
	for (int i = 1; i < wiRB-1; i++) {
		for (int j = 1; j < heRB-1; j++) {
			red = img_bayer->data[i*2+j*2*img_bayer->xsize];
			//if (red < minR) minR = red;
			//if (red > maxR) maxR = red;
			imgR->data[i+j*wiRB] = red;

			blue = img_bayer->data[i*2+1+(j*2+1)*img_bayer->xsize];
			//if (blue < minB) minB = blue;
			//if (blue > maxB) maxB = blue;
			imgB->data[i+j*wiRB] = blue;
			
			// //// if scale is 1:
			//green = 0.5*(img_bayer->data[i*2+1+j*2*wi] + img_bayer->data[i*2+(j*2+1)*wi]);
			//if (green < minG) minG = green;
			//if (green > maxG) maxG = green;
			//imgG->data[i+j*wiRB] = green;
			green = img_bayer->data[i*2+1+j*2*img_bayer->xsize];
			imgG->data[i*2+1+j*2*imgG->xsize] = green;
			//if (green < minG) minG = green;
			//if (green > maxG) maxG = green;
			green = img_bayer->data[i*2+(j*2+1)*img_bayer->xsize];
			//if (green < minG) minG = green;
			//if (green > maxG) maxG = green;
			imgG->data[i*2+(j*2+1)*imgG->xsize] = green;
			imgG->data[i*2+j*2*imgG->xsize] = 0.25* (img_bayer->data[(i*2+1)+j*2*img_bayer->xsize] + img_bayer->data[i*2+(j*2+1)*img_bayer->xsize] + 
				img_bayer->data[(i*2-1)+j*2*img_bayer->xsize] + img_bayer->data[i*2+(j*2-1)*img_bayer->xsize]);
			imgG->data[i*2+1+(j*2+1)*imgG->xsize] = 0.25* (img_bayer->data[(i*2+1)+j*2*img_bayer->xsize] + img_bayer->data[i*2+(j*2+1)*img_bayer->xsize] + 
				img_bayer->data[(i*2+2)+(j*2+1)*img_bayer->xsize] + img_bayer->data[(i*2+1)+(j*2+2)*img_bayer->xsize]);
		}
	}
	printf("done\n");
}

template <typename T>
void print_RMSE(vector<T>& xR, vector<T>& yR, vector<T>& xGr, vector<T>& yGr, vector<T>& xB, vector<T>& yB, vector<T>& xGb, vector<T>& yGb)
{
	printf("Data stats calculation... \n");
	int ntachesr = xR.size();
	int ntachesb = xB.size();
	T RMSE_red_dist = 0, RMSE_blue_dist = 0;
	T mean_red_dist = 0, mean_blue_dist = 0, stddev_red_dist = 0, stddev_blue_dist = 0, dispR_dist = 0, dispB_dist = 0;
	for (int i = 0; i < ntachesr; i++)
	{
		T dispR = std::sqrt((xR[i]-xGr[i])*(xR[i]-xGr[i]) + (yR[i]-yGr[i])*(yR[i]-yGr[i]));
		if (dispR_dist < dispR) dispR_dist = dispR;
		RMSE_red_dist += dispR*dispR;
		mean_red_dist += dispR;
	}
	for (int i = 0; i < ntachesb; i++)
	{
		T dispB = std::sqrt((xB[i]-xGb[i])*(xB[i]-xGb[i]) + (yB[i]-yGb[i])*(yB[i]-yGb[i]));
		if (dispB_dist < dispB) dispB_dist = dispB;
		RMSE_blue_dist += dispB*dispB;
		mean_blue_dist += dispB;
	}
	mean_red_dist /= ntachesr;
	mean_blue_dist /= ntachesb;
	for (int i = 0; i < ntachesr; i++)
		stddev_red_dist += (std::sqrt((xR[i]-xGr[i])*(xR[i]-xGr[i]) + (yR[i]-yGr[i])*(yR[i]-yGr[i])) - mean_red_dist) * 
			(std::sqrt((xR[i]-xGr[i])*(xR[i]-xGr[i]) + (yR[i]-yGr[i])*(yR[i]-yGr[i])) - mean_red_dist);
	for (int i = 0; i < ntachesb; i++)
		stddev_blue_dist += (std::sqrt((xB[i]-xGb[i])*(xB[i]-xGb[i]) + (yB[i]-yGb[i])*(yB[i]-yGb[i])) - mean_blue_dist) * 
			(std::sqrt((xB[i]-xGb[i])*(xB[i]-xGb[i]) + (yB[i]-yGb[i])*(yB[i]-yGb[i])) - mean_blue_dist);
	printf(" done.\n");

	printf("RMSE R-G & B-G are:	%f	%f \n", std::sqrt(RMSE_red_dist/ntachesr), std::sqrt(RMSE_blue_dist/ntachesb));
	printf("mean R-G & B-G are:	%f	%f \n", mean_red_dist, mean_blue_dist);
	printf("stdd R-G & B-G are:	%f	%f \n", std::sqrt(stddev_red_dist/ntachesr), std::sqrt(stddev_blue_dist/ntachesb));
	printf("maxd R-G & B-G are:	%f	%f \n", dispR_dist, dispB_dist);	
}

template <typename T>
void keypnts2file(const char* fnameXYdist, 
	vector<T>& xR, vector<T>& yR, vector<T>& xGr, vector<T>& yGr, vector<T>& xB, vector<T>& yB, vector<T>& xGb, vector<T>& yGb)
{
	printf("Saving corrected keypoints to file... ");
	FILE *pfile_dist;
	pfile_dist = fopen(fnameXYdist, "wt");
	if(pfile_dist == NULL) {
		printf("cannot open file %s.\n", fnameXYdist);
		exit(1);
	}
	int lenR = xR.size();
	int lenB = xB.size();
	int lenmax = std::max(lenR, lenB);
	int lenmin = std::min(lenR, lenB);
	for (int i = 0; i < lenmax; i++)
	{
		if (i < lenmin)
			fprintf(pfile_dist, "%f %f %f %f %f %f %f %f\n", xR[i], yR[i], xGr[i], yGr[i], xGb[i], yGb[i], xB[i], yB[i]);
		else
		{
			if (lenR < lenB)
				fprintf(pfile_dist, "%f %f %f %f %f %f %f %f\n", 0, 0, 0, 0, xGb[i], yGb[i], xB[i], yB[i]);
			else
				fprintf(pfile_dist, "%f %f %f %f %f %f %f %f\n", xR[i], yR[i], xGr[i], yGr[i], 0, 0, 0, 0);
		}
	}
	fclose(pfile_dist);
	printf("done.\n");
}

template <typename T>
void get_polynom(vector<T>& xF, vector<T>& yF, vector<T>& xGf, vector<T>& yGf,
	vector<T>& paramsXF, vector<T>& paramsYF, int degX, int degY, T xp, T yp)
{
	printf("Obtaining correction polynomials... ");
	//vector<T> polyB = getParamsCorrection(xB, yB, xGr, yGr, degX, degY, xp, yp); 
	vector<T> polyF = getParamsCorrection(xF, yF, xGf, yGf, degX, degY, xp, yp);
	int sizex = (degX + 1) * (degX + 2) / 2;
	int sizey = (degY + 1) * (degY + 2) / 2;
	paramsXF = polyF.copyRef(0, sizex-1);
	paramsYF = polyF.copyRef(sizex, sizex+sizey-1);
	//paramsXB = polyB.copyRef(0, sizex-1);
	//paramsYB = polyB.copyRef(sizex, sizex+sizey-1);
	printf("done.\n");
}

template <typename T>
void correct_channel(image_double& imgF, image_double& imgFz, vector<T>& paramsXF, vector<T>& paramsYF, 
	int spline_order, int degX, int degY, T xp, T yp, int wiG, int heG, int scale)
{
	printf("\nCorrected channel is being calculated... \n");
	prepare_spline(imgF, spline_order);
	for (int i = 0; i < wiG; i++) {
		for (int j = 0; j < heG; j++) {
			T p1=0, p2=0;
			undistortPixel(p1, p2, paramsXF, paramsYF, i, j, xp, yp, degX, degY);
			T clr = interpolate_image_double(imgF, spline_order, p1/scale+0.5, p2/scale+0.5); // +0.5 to compensate -0.5 in interpolation function
			if (clr < 0) clr = 0; 
			if (clr > 255) clr = 255;
			imgFz->data[i+j*imgFz->xsize] = clr;
		}
		double percent = ((double)i / (double)wiG)*100;
		if (!(i % (int)(0.2*wiG))) printf("%i%c", (int)percent+1, '%');
		else if (!(i % (int)(0.04*wiG))) printf(".");
	}
	printf("done.\n");
}

template <typename T>
void circuit(int argc, char ** argv, bool clr, bool test = false)
{
	char* fnameRGB = argv[1];
	char* fnameR = argv[2]; 
	char* fnameG = argv[3];
	char* fnameB = argv[4];
	char* fnameXYdist = argv[5]; 
	char* fnameXYcorr = argv[6];
	int nImgs = (argc-7)/4;

	T scale = 2;
	image_double img_bayer = read_pgm_image_double(fnameRGB);
	// // Uncomment below to replace the line above, make sure you know which direction to rotate (left or right)
	//image_double img_bayer2 = read_pgm_image_double(fnameRGB);
	//image_double img_bayer = image_rotate_right<T>(img_bayer2); free_image_double(img_bayer2);
	int wi = img_bayer->xsize, he = img_bayer->ysize;
	int wiRB = wi/2, heRB = he/2;
	int wiG = wiRB*scale, heG = heRB*scale;
	image_double imgR = new_image_double_ini(wiRB, heRB, 255);
	image_double imgG = new_image_double_ini(wiG, heG, 255);
	image_double imgB = new_image_double_ini(wiRB, heRB, 255);
	raw2rgb<T>(img_bayer, imgR, imgG, imgB);
	vector<T> xR, yR, xGr, yGr, xB, yB, xGb, yGb, rR, rG, rB;
	keypnts_circle<T>(imgR, imgG, imgB, xR, yR, rR, xGr, yGr, rG, xB, yB, rB, xGb, yGb, scale, clr);
	//keypnts_sift<T>(imgR, imgG, imgB, xR, yR, xGr, yGr, xB, yB, xGb, yGb, scale, clr);
	keypnts2file(fnameXYdist, xR, yR, xGr, yGr, xB, yB, xGb, yGb);
	print_RMSE(xR, yR, xGr, yGr, xB, yB, xGb, yGb);
	
	vector<T> paramsXR, paramsYR, paramsXB, paramsYB;
	int degX = 11, degY = 11;
	T xp = (T)imgG->xsize/2+0.2, yp = (T)imgG->ysize/2+0.2;
	get_polynom<T>(xR, yR, xGr, yGr, paramsXR, paramsYR, degX, degY, xp, yp);
	get_polynom<T>(xB, yB, xGb, yGb, paramsXB, paramsYB, degX, degY, xp, yp);	
	
	int spline_order = 3;
	image_double imgRz = new_image_double_ini(wiG, heG, 255);
	image_double imgBz = new_image_double_ini(wiG, heG, 255);
	correct_channel<T>(imgR, imgRz, paramsXR, paramsYR, spline_order, degX, degY, xp, yp, wiG, heG, scale);
	correct_channel<T>(imgB, imgBz, paramsXB, paramsYB, spline_order, degX, degY, xp, yp, wiG, heG, scale);
	
	printf("\nSaving images to file... \n");
	write_pgm_image_double(imgRz, fnameR);
	write_pgm_image_double(imgG, fnameG);
	write_pgm_image_double(imgBz, fnameB);
	
    bool green_proc = false;
    vector<T> rB_scale = rB*scale;
    vector<T> rR_scale = rR*scale;
    circle_redefine<T>(imgRz, imgG, imgBz, xR, yR, rR_scale,
                       xGr, yGr, rG, xB, yB, rB_scale, xGb, yGb, 1, clr, green_proc);
	//keypnts_sift<T>(imgRz, imgG, imgBz, xR, yR,  xGr, yGr, xB, yB, xGb, yGb, 1, clr);
	keypnts2file(fnameXYcorr, xR, yR, xGr, yGr, xB, yB, xGb, yGb);
	print_RMSE(xR, yR, xGr, yGr, xB, yB, xGb, yGb);

	printf("\nCorrecting blue and red channels for other input images... \n");
	for (int i = 0; i < nImgs; i++)
	{
		image_double imgn_bayer = read_pgm_image_double(argv[7+i*4+0]);
		image_double imgnR = new_image_double_ini(wiRB, heRB, 255);
		image_double imgnG = new_image_double_ini(wiG, heG, 255);
		image_double imgnB = new_image_double_ini(wiRB, heRB, 255);
		//separate the channels
		raw2rgb<T>(imgn_bayer, imgnR, imgnG, imgnB);
		// measure test image RMSE if necessary
		vector<T> xnR, ynR, xnGr, ynGr, xnB, ynB, xnGb, ynGb, rnR, rnG, rnB;
		if (test) {
			//keypnts_sift<T>(imgnR, imgnG, imgnB, xnR, ynR,  xnGr, ynGr, xnB, ynB, xnGb, ynGb, scale, clr);
			keypnts_circle<T>(imgnR, imgnG, imgnB, xnR, ynR, rnR, xnGr, ynGr, rnG, xnB, ynB, rnB, xnGb, ynGb, scale, clr);
			print_RMSE(xnR, ynR, xnGr, ynGr, xnB, ynB, xnGb, ynGb);
		}
		// perform the correction
		image_double imgnRz = new_image_double_ini(wiG, heG, 255);
		image_double imgnBz = new_image_double_ini(wiG, heG, 255);
		correct_channel<T>(imgnR, imgnRz, paramsXR, paramsYR, spline_order, degX, degY, xp, yp, wiG, heG, scale);
		correct_channel<T>(imgnB, imgnBz, paramsXB, paramsYB, spline_order, degX, degY, xp, yp, wiG, heG, scale);
		// save corrected images to files	
		printf("\Saving images to file... \n");
		write_pgm_image_double(imgnRz, argv[7+i*4+1]); 
		write_pgm_image_double(imgnG, argv[7+i*4+2]);
		write_pgm_image_double(imgnBz, argv[7+i*4+3]);
		// if its test image, measure RMSE
		if (test) {
			//keypnts_sift<T>(imgnRz, imgnG, imgnBz, xnR, ynR,  xnGr, ynGr, xnB, ynB, xnGb, ynGb, 1, clr);
            bool green_proc = false;
            vector<T> rnR_scale = rnR*scale;
            vector<T> rnB_scale = rnB*scale;
            circle_redefine<T>(imgnRz, imgnG, imgnBz, xnR, ynR,
                               rnR_scale, xnGr, ynGr, rnG, xnB, ynB, rnB_scale, xnGb, ynGb, 1, clr, green_proc);
			print_RMSE(xnR, ynR, xnGr, ynGr, xnB, ynB, xnGb, ynGb);
		}
		// free memory
		free_image_double(imgnR); free_image_double(imgnG); free_image_double(imgnB);
		free_image_double(imgnRz); free_image_double(imgnBz);
	}
	// free memory
	free_image_double(img_bayer);
	free_image_double(imgR); free_image_double(imgG); free_image_double(imgB);
	free_image_double(imgRz); free_image_double(imgBz);
}

template <typename T>
void printMono(FILE* pfile, T mono, int degX, int degY) {
	if (degX != 0 && degY != 0) {
		if (mono > 0)	fprintf(pfile, "+ %.16g * x^%i * y^%i\n", mono, degX, degY);
		else			fprintf(pfile, "- %.16g * x^%i * y^%i\n", -1*mono, degX, degY); }
	else if (degX == 0 && degY == 0) {
		if (mono > 0)	fprintf(pfile, "+ %.16g\n", mono);
		else			fprintf(pfile, "- %.16g\n", -1*mono); }
	else if (degX == 0) {
		if (mono > 0)	fprintf(pfile, "+ %.16g * y^%i\n", mono, degY);
		else			fprintf(pfile, "- %.16g * y^%i\n", -1*mono, degY); }
	else {
		if (mono > 0)	fprintf(pfile, "+ %.16g * x^%i\n", mono, degX);
		else			fprintf(pfile, "- %.16g * x^%i\n", -1*mono, degX); }
}

template <typename T>
void save_poly(char* fname, vector<T>& paramsX, vector<T>& paramsY, const int degX, const int degY) {
	int sizex = (degX + 1) * (degX + 2) / 2;
	int sizey = (degY + 1) * (degY + 2) / 2;
	//vector<T> paramsX = poly_params.copyRef(0, sizex-1);
	//vector<T> paramsY = poly_params.copyRef(sizex, sizex+sizey-1);
	FILE *pfile;
	pfile = fopen(fname, "wt");
	if(pfile == NULL)
		printf("cannot open file %s.\n", fname);
	int idx = 0;
	fprintf(pfile, "# polyX(x,y): \n");
	for (int i = degX; i >= 0; i--) {
		for (int j = 0; j <= i; j++) {
			printMono(pfile, paramsX[idx], i-j, j);
			idx++; 	} }
	idx = 0;
	fprintf(pfile, "# polyY(x,y): \n");
	for (int i = degY; i >= 0; i--) {
		for (int j = 0; j <= i; j++) {
			printMono(pfile, paramsY[idx], i-j, j);
			idx++; } }
	fclose(pfile);
}

template <typename T>
void polyEstimation(int argc, char ** argv, bool clr) {
	printf("\Polynomial estimation... \n");
	char* fnameRGB = argv[1];
	char* fnamePolyR = argv[2]; 
	char* fnamePolyB = argv[3]; 
	T scale = 2;
	image_double img_bayer = read_pgm_image_double(fnameRGB);
	//image_double img_bayer2 = read_pgm_image_double(fnameRGB);
	//image_double img_bayer = image_rotate_right<T>(img_bayer2); free_image_double(img_bayer2);
	int wi = img_bayer->xsize, he = img_bayer->ysize;
	int wiRB = wi/2, heRB = he/2;
	int wiG = wiRB*scale, heG = heRB*scale;
	image_double imgR = new_image_double_ini(wiRB, heRB, 255);
	image_double imgG = new_image_double_ini(wiG, heG, 255);
	image_double imgB = new_image_double_ini(wiRB, heRB, 255);
	
	raw2rgb<T>(img_bayer, imgR, imgG, imgB);
	vector<T> xR, yR, xGr, yGr, xB, yB, xGb, yGb, rR, rG, rB;
	keypnts_circle<T>(imgR, imgG, imgB, xR, yR, rR, xGr, yGr, rG, xB, yB, rB, xGb, yGb, scale, clr);
	print_RMSE(xR, yR, xGr, yGr, xB, yB, xGb, yGb);
	vector<T> paramsXR, paramsYR, paramsXB, paramsYB;
	int degX = 11, degY = 11;
	T xp = (T)imgG->xsize/2+0.2, yp = (T)imgG->ysize/2+0.2;
	get_polynom<T>(xR, yR, xGr, yGr, paramsXR, paramsYR, degX, degY, xp, yp);
	get_polynom<T>(xB, yB, xGb, yGb, paramsXB, paramsYB, degX, degY, xp, yp);

	save_poly(fnamePolyR, paramsXR, paramsYR, degX, degY);
	save_poly(fnamePolyB, paramsXB, paramsYB, degX, degY);

	free_image_double(img_bayer);
	free_image_double(imgR); free_image_double(imgG); free_image_double(imgB);
}

/* Pop out one character from char array */
char popchar( char *c, int idx, int size) {
	char res = c[idx];
	for (int i = idx; i < size; i++)
		if (i > idx) c[i-1] = c[i];
	return res;
}

/* Difines the degee values for X and Y polinomials from file */
void read_degree(FILE *pfile, int& degX, int& degY) {
	degX = 0; degY = 0;
	char buffer[250], sign[2];
	double coef = 0;
	const char* grid = "#", *star = "*" , *plus = "+", *minu = "-";
	bool flagX = false, flagY = false;
	while (!feof(pfile)) {
		/* expects each monimial to have a form of: "+/- coef * x^deg1 * y^deg2 " */
		/* deg1 and/or deg2 can be zeros, leaving just a coefficent value */
		fscanf(pfile, "%s", sign);
		/* if sign "#" is met - it's a comment, we may skip it. */
		if (strcmp(sign, grid) != 0) {
			char mono1[5], mono2[5];
			fscanf(pfile, "%lf", &coef);
			long currPos = ftell (pfile);
			fscanf(pfile, "%s", sign);
			if (strcmp(sign, star) != 0) {
				if (!feof(pfile)) {
					fseek(pfile, currPos, SEEK_SET);
					assert(strcmp(sign, plus) == 0 || strcmp(sign, minu) == 0 || strcmp(sign, grid) == 0); 	} }
			else {
				fscanf(pfile, "%s", mono1);
				currPos = ftell (pfile);
				fscanf(pfile, "%s", sign);
				if (strcmp(sign, star) != 0) {
					if (!feof(pfile)) {
						fseek(pfile, currPos, SEEK_SET);
						assert(strcmp(sign, plus) == 0 || strcmp(sign, minu) == 0 || strcmp(sign, grid) == 0); 	}
				}
				else fscanf(pfile, "%s", mono2); }
			/* pop-out "x^" and "y^" so that to calculate degrees. */	
			popchar(mono1, 0, 5); popchar(mono1, 0, 5);
			popchar(mono2, 0, 5); popchar(mono2, 0, 5);
			int tmpdeg = std::max(atoi(mono1), atoi(mono2));
			/* there are two polynomials in file: degree of X is defined by flagX, second - flagY. */
			/* assumed first poly is for X direction, second is for Y. */
			if (flagX) { if (tmpdeg > degX) degX = tmpdeg; }
			else { if (tmpdeg > degY) degY = tmpdeg; }		
		}
		else {
			fscanf(pfile, "%s", buffer);
			if (degX == 0 && degY == 0) flagX = true;
			else {flagX = false; flagY = true;}
		}
	}
	fseek(pfile, 0, SEEK_SET);
}

/* Finds an index of coefTerm[] vector based on given x and y degrees. */
int coefIdx(int degree, int x, int y) {
	int a1 = degree + 1;
	int n = std::abs(x+y-degree)+1;
	int an = x+y+1;
	int Sn = n * (a1 + an);
	Sn /= 2;
	return Sn-an+y; }

/* Reads the poly coefficients and insert them into vector of coefficients - coefTerm[]. */
/* Also returns the degrees for each polynomial. */
template <typename T>
vector<T> read_poly(char* fname, int& degX, int& degY) {
	FILE *pfile;
	pfile = fopen(fname, "r");
	if(pfile == NULL) printf("unable to open file %s.\n", fname);
	/* get the degrees for each polynomial - need it to know the size of vectors paramsX and paramsY. */
	read_degree(pfile, degX, degY);
	int sizex = (degX + 1) * (degX + 2) / 2;
	int sizey = (degY + 1) * (degY + 2) / 2;
	vector<T> paramsX(sizex), paramsY(sizey);
	char buffer[500], sign[2];
	const char* grid = "#", *star = "*", *plus = "+", *minu = "-";
	bool flagX = false, flagY = false;
	while (!feof(pfile)) {
		/* reading is done by the same manner as in "read_degree(pfile, degX, degY);" */
		fscanf(pfile, "%s", sign);
		if (strcmp(sign, grid) != 0 && (strcmp(sign, plus) == 0 || strcmp(sign, minu) == 0 ) ) {
			char mono1[5], mono2[5];
			double coef = 0;
			fscanf(pfile, "%lf", &coef);
			long currPos = ftell (pfile);
			/* save the coef sign. */
			if (strcmp(sign, minu) == 0) coef *= -1;
			fscanf(pfile, "%s", sign);
			if (strcmp(sign, star) != 0) {
				if (!feof(pfile)) {
					fseek(pfile, currPos, SEEK_SET);
					assert(strcmp(sign, plus) == 0 || strcmp(sign, minu) == 0 || strcmp(sign, grid) == 0); 	} 	}
			else {
				fscanf(pfile, "%s", mono1);
				currPos = ftell (pfile);
				fscanf(pfile, "%s", sign);
				if (strcmp(sign, star) != 0) {
					if (!feof(pfile)) {
						fseek(pfile, currPos, SEEK_SET);
						assert(strcmp(sign, plus) == 0 || strcmp(sign, minu) == 0 || strcmp(sign, grid) == 0); 	} 	}
				else fscanf(pfile, "%s", mono2); }
			char x = popchar(mono1, 0, 5); popchar(mono1, 0, 5);
			char y = popchar(mono2, 0, 5); popchar(mono2, 0, 5);
			/* see which degree belongs to which variable. */
			int tmpdegX = 0, tmpdegY = 0;
			if ( x == 'x')		tmpdegX = atoi(mono1);
			else if ( x == 'y')	tmpdegY = atoi(mono1);
			if ( y == 'y') 	tmpdegY = atoi(mono2);
			/* save the coef to accoring paramsX/Y; fist poly in file belongs to paramsX, second - paramsY. */
			if (flagX && !flagY) { 
				int idx = coefIdx(degX, tmpdegX, tmpdegY);
				paramsX[idx] = coef; }
			else {
				int idx = coefIdx(degY, tmpdegX, tmpdegY);
				paramsY[idx] = coef; }		
		}
		else {
			fscanf(pfile, "%s", buffer);
			if (!flagX) flagX = true;
			else flagY = true; }
	}
	fclose(pfile);
	vector<T> poly_params(sizex+sizey);
	/* copy the paramsX and paramsY to one vector. */
	/* later we will be able to separate them since the degrees for each poly are known. */
	for (int k = 0; k < sizex; k++) poly_params[k] = paramsX[k];
	for (int k = 0; k < sizey; k++) poly_params[k+sizex] = paramsY[k];
	return poly_params;
}

template <typename T>
void aberCorrection(int argc, char ** argv, bool clr)
{
	printf("\Aberration correction... \n");
	char* fnameRGB = argv[1];
	char* fnamePolyR = argv[2]; 
	char* fnamePolyB = argv[3]; 
	char* fnameR = argv[4]; 
	char* fnameG = argv[5];
	char* fnameB = argv[6];
	
	int degX = 11, degY = 11;
	int sizex = (degX + 1) * (degX + 2) / 2;
	int sizey = (degY + 1) * (degY + 2) / 2;

	T scale = 2;
	image_double imgn_bayer = read_pgm_image_double(fnameRGB);
	int wi = imgn_bayer->xsize, he = imgn_bayer->ysize;
	int wiRB = wi/2, heRB = he/2;
	int wiG = wiRB*scale, heG = heRB*scale;
	image_double imgnR = new_image_double_ini(wiRB, heRB, 255);
	image_double imgnG = new_image_double_ini(wiG, heG, 255);
	image_double imgnB = new_image_double_ini(wiRB, heRB, 255);
	raw2rgb<T>(imgn_bayer, imgnR, imgnG, imgnB);

	vector<T> paramsR = read_poly<T>(fnamePolyR, degX, degY);
	vector<T> paramsB = read_poly<T>(fnamePolyB, degX, degY);
	vector<T> paramsXR = paramsR.copyRef(0, sizex-1);
	vector<T> paramsYR = paramsR.copyRef(sizex, sizex+sizey-1);
	vector<T> paramsXB = paramsB.copyRef(0, sizex-1);
	vector<T> paramsYB = paramsB.copyRef(sizex, sizex+sizey-1);

	int spline_order = 3;
	image_double imgnRz = new_image_double_ini(wiG, heG, 255);
	image_double imgnBz = new_image_double_ini(wiG, heG, 255);
	T xp = (T)imgnG->xsize/2+0.2, yp = (T)imgnG->ysize/2+0.2;
	correct_channel<T>(imgnR, imgnRz, paramsXR, paramsYR, spline_order, degX, degY, xp, yp, wiG, heG, scale);
	correct_channel<T>(imgnB, imgnBz, paramsXB, paramsYB, spline_order, degX, degY, xp, yp, wiG, heG, scale);
	printf("\Saving images to file... \n");
	write_pgm_image_double(imgnRz, fnameR); 
	write_pgm_image_double(imgnG, fnameG);
	write_pgm_image_double(imgnBz, fnameB);

	free_image_double(imgn_bayer);
	free_image_double(imgnR); free_image_double(imgnG); free_image_double(imgnB);
	free_image_double(imgnRz); free_image_double(imgnBz);
}

int main(int argc, char ** argv)
{
	bool clr = false; // deals with black circles on white background
	bool test = false; // true if the image to correct is a test image to measure the correction RMSE
    printf("Program usage: \n");
    printf("Polynomial estimation:\n");
    printf("chromaberrat fname_raw_calib.pgm fname_poly_red.txt fname_poly_blue.txt \n\n");
    printf("CA correction using estimated polynomial:\n");
    printf("chromaberrat fname_raw.pgm fname_poly_red.txt fname_poly_blue.txt fname_raw_red_corr.pgm fname_raw_green_corr.pgm fname_raw_blue_corr.pgm \n\n");
    printf("Running all circuit (polynomial estimation - image correction):\n");
    printf("chromaberrat fname_raw_calib.pgm fname_raw_calib_red_corr.pgm fname_raw_calib_green_corr.pgm fname_raw_calib_blue_corr.pgm fname_raw_calib_keyp_dist.txt fname_raw_calib_keyp_corr.txt [fname_img_n.pgm fname_img_n_red_corr.pgm fname_img_n_green_corr.pgm fname_img_n_blue_corr.pgm, ...]\n\n");

	if (argc > 7)  // runs all circuit, change settings inside
        circuit<double>(argc, argv, clr, test);

	if (argc == 4) // only estimates and saves polynomial
		polyEstimation<double>(argc, argv, clr);

	if (argc == 7) // reads image and poly, corrects input and saves three corrected channels separately
		aberCorrection<double>(argc, argv, clr);

	return 0; 	
}

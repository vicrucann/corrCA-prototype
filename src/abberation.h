#ifndef ABBERATION_H
#define ABBERATION_H

#include <stack>
#include "image.h"

#define PI 3.14159

class Pixel {
public: int x, y;
		Pixel(int X, int Y) {x = X, y = Y;}
};

class CCStats {
public:
	int nPoints;
	double centerX, centerY, radius1, radius2;
	int perimeter; // need to calculate compactness measure of shape to eliminate noise
};

int white_neighbors(Pixel& p, image_double img)
{
	int nn = 0;
	double color;
	Pixel p1(p.x-1, p.y); 
	if (p1.x >= 0 && p1.x < img->xsize && p1.y >= 0 && p1.y < img->ysize) {
		if (img->data[p1.x+p1.y*img->xsize] == 255)
			nn++; }
	else nn++;

	Pixel p2(p.x+1, p.y);
	if (p2.x >= 0 && p2.x < img->xsize && p2.y >= 0 && p2.y < img->ysize) {
		if (img->data[p2.x+p2.y*img->xsize] == 255)
			nn++; }
	else nn++;

	Pixel p3(p.x, p.y-1);
	if (p3.x >= 0 && p3.x < img->xsize && p3.y >= 0 && p3.y < img->ysize) {
		if (img->data[p3.x+p3.y*img->xsize] == 255)
			nn++; }
	else nn++;

	Pixel p4(p.x, p.y+1);
	if (p4.x >= 0 && p4.x < img->xsize && p4.y >= 0 && p4.y < img->ysize) {
		if (img->data[p4.x+p4.y*img->xsize] == 255)
			nn++; }
	else nn++;
	
	if (nn > 0) nn = 1;
	return nn;
}

void extract_CCStats(std::vector<Pixel>& cc, CCStats& stats, image_double img)
{
	stats.nPoints = cc.size();
	stats.perimeter = 0;
	double meanX = 0, meanY = 0;
	double minX = img->xsize, minY = img->ysize, maxX = 0, maxY = 0;
	for (int i = 0; i < cc.size(); i++) {
		meanX += cc[i].x;
		meanY += cc[i].y;
		if (cc[i].x > maxX) maxX = cc[i].x;
		if (cc[i].x < minX) minX = cc[i].x;
		if (cc[i].y > maxY) maxY = cc[i].y;
		if (cc[i].y < minY) minY = cc[i].y;
		stats.perimeter += white_neighbors(Pixel(cc[i].x, cc[i].y), img);
	}
	stats.centerX = meanX / cc.size();
	stats.centerY = meanY / cc.size();
	stats.radius1 = 0.5 * (maxX - minX);
	stats.radius2 = 0.5 * (maxY - minY);
}

int extract_cc_(Pixel& p, std::vector<Pixel>& cc, image_double& img)
{
	std::stack<Pixel> s;
	double* color;
	if (p.x >= 0 && p.x < img->xsize && p.y >= 0 && p.y < img->ysize)
		color = &img->data[p.x+p.y*img->xsize];
	else 
		return 0;
	while (*color == 0 || !s.empty())
	{
		if (*color == 0)
		{
			cc.push_back(p);
			*color = 255;
			s.push(Pixel(p.x+1, p.y));
			s.push(Pixel(p.x-1, p.y));
			s.push(Pixel(p.x, p.y+1));
			s.push(Pixel(p.x, p.y-1));
		}
		p = s.top();
		s.pop();
		if (p.x >= 0 && p.x < img->xsize && p.y >= 0 && p.y < img->ysize)
			color = &img->data[p.x+p.y*img->xsize];
		else 
			*color = 255;
	}
	return cc.size();
}

bool extract_cc(Pixel& p, std::vector<Pixel>& cc, image_double& img)
{
	bool result = false;
	std::stack<Pixel> s;
	double* color;
	if (p.x >= 0 && p.x < img->xsize && p.y >= 0 && p.y < img->ysize)
		color = &img->data[p.x+p.y*img->xsize];
	else
		return result;
	if (*color == 0) {
		result = true;
		cc.push_back(p);
		*color = 255;
		s.push(Pixel(p.x+1, p.y));
		s.push(Pixel(p.x-1, p.y));
		s.push(Pixel(p.x, p.y+1));
		s.push(Pixel(p.x, p.y-1));
		while (!s.empty()) {
			extract_cc(s.top(), cc, img);
			s.pop(); }
	}
	return result;
}

void CC(std::vector<CCStats>& ccstats, image_double& imgbi, char channel)
{
	image_double img_copy = new_image_double_copy(imgbi);
	printf("\nchannel %c: ", channel);
	double meansize = 0;
	for (int i = 0; i < imgbi->xsize; i++) {
		for (int j = 0; j < imgbi->ysize; j++) {
			std::vector<Pixel> ccC;
			CCStats stats;
			int npix = extract_cc_(Pixel(i, j), ccC, imgbi);
			if (npix > 180) {			
				extract_CCStats(ccC, stats, img_copy);
				double compactness = 4*PI*stats.nPoints / (stats.perimeter*stats.perimeter);
				if (std::min(stats.radius1, stats.radius2) > 8 && compactness < 1.3 && compactness > 0.7)
				{ 
					ccstats.push_back(stats);
					meansize += stats.nPoints;
				}
			}
		}
		double percent = ((double)i / (double)imgbi->xsize)*100;
		if (!(i % (int)(0.2*imgbi->xsize))) printf("%i%c", (int)percent+1, '%');
		else if (!(i % (int)(0.04*imgbi->xsize))) printf(".");
	}

	// retrieve min_size and max_size to build a size histogram
	int max_val = 0;
	int rad_thre = 7;
	for (int i = 0; i < ccstats.size(); i++) {
		if (ccstats[i].nPoints > max_val)
			max_val = ccstats[i].nPoints;
	}

	std::vector<int> hist(max_val);
	std::vector<std::stack<int> > hist_stack(max_val); // to keep indeces of all the circles for given size
	//run through all the sizes and build frequency histogram 
	for (int i = 0; i < ccstats.size(); i++) {
		int val = ccstats[i].nPoints-1;
		//double compactness = 4*PI*ccstats[i].nPoints / (ccstats[i].perimeter*ccstats[i].perimeter);
		//double r1 = ccstats[i].radius1;
		//double r2 = ccstats[i].radius2;
		//if (std::min(r1, r2) >= rad_thre && compactness < 1.4 && compactness > 0.6)
		//{
		hist[val]++;
		hist_stack[val].push(i);
		//}
	}
	int pick = 0, commonsize = 0;
	for (int i = 0; i < hist.size(); i++)
	{
		if (pick < hist[i])
		{
			pick = hist[i];
			commonsize = i;
		}
	}
	
	// collect the inliers in positive direction from commonsize idx
	int zerosgap = 1700;
	int count = 0;
	int flag = zerosgap;
	std::vector<int> inliers(ccstats.size());
	while (flag != 0 && commonsize+count < hist.size()) {
		int hist_idx = commonsize + count;
		int onesizecircles = hist[hist_idx];
		if (onesizecircles == 0)
			flag--;
		else  {
			while (!hist_stack[hist_idx].empty()) {
				inliers[hist_stack[hist_idx].top()] = 1;
				//inliers.push(hist_stack[hist_idx].top());
				hist_stack[hist_idx].pop();
			}
		}
		count++;
	}
	// collect the inliers in negative direction from commonsize idx
	count = -1;
	flag = zerosgap;
	while (flag != 0 && commonsize+count >= 0 ) {
		int hist_idx = commonsize + count;
		int onesizecircles = hist[hist_idx];
		if (onesizecircles == 0)
			flag--;
		else {
			while (!hist_stack[hist_idx].empty()) {
				inliers[hist_stack[hist_idx].top()] = 1;
				hist_stack[hist_idx].pop();
			}
		}
		count--;
	}
			
	//meansize /= ccstats.size();
	int idx = 0;
	while (idx < ccstats.size()) {
		//double p = ccstats[idx].perimeter;
		//double r1 = ccstats[idx].radius1;
		//double r2 = ccstats[idx].radius2;
		//double q = std::max(r1,r2) / std::min(r1,r2);
		//double area = PI*r1*r2;
		//double area_ = ccstats[idx].nPoints;
		//double compactness = 4*PI*ccstats[idx].nPoints / (ccstats[idx].perimeter*ccstats[idx].perimeter);
		//double length = 2*PI* ccstats[idx].radius;
		//if (length >= ccstats[idx].perimeter-3 && length <= ccstats[idx].perimeter+3)
		//if (compactness > 1.01 && compactness < 0.99 || std::min(r1,r2) <= 8)
		//if (ccstats[idx].nPoints < 0.7*meansize || ccstats[idx].nPoints > 2.5*meansize) 
		//if ((ccstats[idx].nPoints > commonsize + 150 || ccstats[idx].nPoints < commonsize - 150) && (compactness > 1.4 || compactness < 0.6) || std::min(r1, r2) < rad_thre)
		if (inliers[idx] == 0)
		{
			ccstats.erase(ccstats.begin() + idx);
			inliers.erase(inliers.begin() + idx);
		}
		else 
			idx++;
	}
	free_image_double(img_copy);
}

// finds corresponding circle center match among two channels
template <typename T>
int findMatch(T xg, T yg, std::vector<CCStats>& trgstats, T scale = 1)
{
	int minidx = -1;
	T mindist = 5; // look for center in radius of 10 pixels
	for (int i = 0; i < trgstats.size(); i++)
	{
		T eucdist = std::sqrt((scale*trgstats[i].centerX - xg)*(scale*trgstats[i].centerX - xg) + (scale*trgstats[i].centerY - yg)*(scale*trgstats[i].centerY - yg));
		if (eucdist <= mindist)
		{
			mindist = eucdist;
			minidx = i;
		}
	}
	return minidx;
}

#endif
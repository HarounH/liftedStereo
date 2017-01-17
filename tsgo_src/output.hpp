#ifndef _OUTPUT_HPP_
#define _OUTPUT_HPP_


// Haroun: Systems includes... I/O stuff and PNG include. 
#include <fstream>
#include <png++/png.hpp>
#include <stdlib.h>
#include <cmath>
void outputGrayDisparities(float *data, int width, int height, const char* filename) {
	png::image< png::ga_pixel > img(width, height);
	float mi,ma;
	mi = ma = data[0];
	float val;
	for (int y = 0; y < height; ++y)
	{
		for (int x = 0; x < width; ++x)
		{
			val = data[x + y*width];
			if(val<mi) mi = val;
			if(val>ma) ma = val;
		}
	}
	float scalefactor = 4; 255.0/(ma-mi);
	std::cout << "scalefactor=" << scalefactor << std::endl;
	std::cout << "Saved image, with min=" << mi << " & max=" << ma << std::endl;

	for (int y = 0; y < height; ++y)
	{
		for (int x = 0; x < width; ++x)
		{
			// std::cout << "writing val=" << val << " rounded to " << round(val) << endl;
			if (data[x+y*width] < 0) {
				img[y][x] = png::ga_pixel(0);
			} else {
				img[y][x] = png::ga_pixel(1 + round(scalefactor*(data[x + y*width])));
			}
		}
	}
	img.write(filename);
}

void outputGrayDisparities(float *data, int width, int height, const char* filename, float scalefactor=2.0) {
	png::image< png::ga_pixel > img(width, height);
	float mi,ma;
	mi = ma = scalefactor*data[0];
	float val;

	for (int y = 0; y < height; ++y)
	{
		for (int x = 0; x < width; ++x)
		{
			val = scalefactor*data[x + y*width];
			if(val<mi) mi = val;
			if(val>ma) ma = val;
			// std::cout << "writing val=" << val << " rounded to " << round(val) << endl;
			if (data[x+y*width] < 0) {
				img[y][x] = png::ga_pixel(0);
			} else {
				img[y][x] = png::ga_pixel(1 + round(scalefactor*(data[x + y*width])));
			}
		}
	}
	img.write(filename);
	std::cout << "scalefactor=" << scalefactor << std::endl;
	std::cout << "Saved image, with min=" << mi << " & max=" << ma << std::endl;
}

void saveSymmetries(int nx, int ny, std::vector<std::vector<int> >& syms, std::string symsimagename) {
	using namespace std;
	png::image<png::rgb_pixel> output(nx,ny);
	symsimagename = "./sym/" + symsimagename;
	auto variableIndex = [&nx](const int x, const int y) {return x + nx*y;};
	auto vx = [&nx](const int vid) { return vid%nx; };
	auto vy = [&nx](const int vid) { return (int)(vid/nx); };

	int nGroupsDisplayed=1;

	for(int grp=0; grp<syms.size(); ++grp) {
		if (syms[grp].size()>1) {
			for(int i=0; i<syms[grp].size(); ++i) {
				int vid=syms[grp][i];
				int r = (nGroupsDisplayed%3)*(nGroupsDisplayed);
				int g = ((nGroupsDisplayed+1)%3)*(nGroupsDisplayed);
				int b = ((nGroupsDisplayed+2)%3)*(nGroupsDisplayed);
				output[vy(vid)][vx(vid)] = png::rgb_pixel(r,g,b);
			}
			nGroupsDisplayed += 1;
		}
	}
	cout << "\tnGroupsDisplayed=" << nGroupsDisplayed << endl;
	output.write(symsimagename);
}


void copyVecLabels2Rez(std::vector<int>& veclabels, int* rez) {
	for(int i=0; i<veclabels.size();++i) {
		rez[i] = veclabels[i];
	}
}

void copyVecLabels2Rez(std::vector<int>& veclabels, std::vector<int>& invsyms, int* rez) {
	for(int i=0; i<invsyms.size();++i) {
		rez[i] = veclabels[ invsyms[i] ];
	}
}

#endif
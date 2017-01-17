// Systems includes... I/O stuff and PNG include. 
#include <iostream>
#include <fstream>
#include <png++/png.hpp>
#include <stdlib.h>
#include <typeinfo>


void saveSymmetries(int nx, int ny, std::vector<std::vector<size_t> >& syms, std::string symsimagename) {
	using namespace std;
	png::image<png::rgb_pixel> output(nx,ny);
	
	auto variableIndex = [&nx](const size_t x, const size_t y) {return x + nx*y;};
	auto vx = [&nx](const size_t vid) { return vid%nx; };
	auto vy = [&nx](const size_t vid) { return (int)(vid/nx); };

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
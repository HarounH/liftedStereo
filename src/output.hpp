#ifndef _OUTPUT_HPP_
#define _OUTPUT_HPP_

#include <iostream>
#include <fstream>
#include <png++/png.hpp>
#include <stdlib.h>
#include <cmath>
#include <algorithm>

void labels2grayimage(std::vector<size_t>& labels, int nx, int ny, double scale, std::string filename) {
	// functions to help with indexing.
	auto variableIndex = [&nx](const size_t x, const size_t y) {return x + nx*y;};
	auto vx = [&nx](const size_t vid) { return vid%nx; };
	auto vy = [&nx](const size_t vid) { return (int)(vid/nx); };
	
	png::image< png::ga_pixel > img(nx, ny);

	for(int vid=0; vid< labels.size(); ++vid) {
		img[vy(vid)][vx(vid)] = png::ga_pixel( (int)(labels[vid]*scale) );
	}
	img.write(filename);
}

#endif
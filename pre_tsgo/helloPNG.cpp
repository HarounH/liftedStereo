/**
@author Haroun
This file is an hello world for PNG++

IMPORTANT: compilcation : g++ this.cpp -lpng
*/

#include <png++/png.hpp> // libpng++-dev
#include <iostream>

using namespace std;

int main(int argc, char** argv) {
	png::image< png::rgb_pixel > image("../../data/stereo2006/Aloe/view1.png");
	png::image< png::ga_pixel > newimg(image.get_width(), image.get_height());
	size_t nx = image.get_width();
	size_t ny = image.get_height();
	for(size_t x=0; x<nx; ++x) {
		for(size_t y=0; y<ny; ++y) {
			newimg[y][x] = png::ga_pixel((int)(255*(x+y)/(nx+ny)));
		}
	}
	newimg.write("output.png");
	return 0;
}
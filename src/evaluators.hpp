#ifndef _EVALUATORS_HPP_
#define _EVALUATORS_HPP_

// Systems includes... I/O stuff and PNG include. 
#include <iostream>
#include <fstream>
#include <png++/png.hpp>
#include <stdlib.h>


double badPixelScore(
	std::string outputFileName,
	std::string goldFileName,
	int nDisparities,
	double threshold=2
	) {

	png::image< png::ga_pixel > output(outputFileName);
	png::image< png::ga_pixel > gold(goldFileName);

	int nx = output.get_width();
	int ny = output.get_height();
	double dout, dgold, diff, total;
	for(int y = 0; y<ny; ++y) {
		for(int x = 0; x<nx; ++x) {
			dout = double(nDisparities)*(output[y][x].value)/255.0;
			dgold = double(nDisparities)*(gold[y][x].value)/255.0;
			diff = dout-dgold;
			diff = diff>0?diff:-diff;
			total += diff>threshold?1:0;
		}
	}
	return total/(nx*ny);
}

double badPixelScore(
	std::vector<size_t>& labels,
	png::image< png::ga_pixel >& gold,
	double goldScale,
	double threshold=2
	) {

	int nx = gold.get_width();
	int ny = gold.get_height();
	auto variableIndex = [&nx](const size_t x, const size_t y) {return x + nx*y;};
	
	double dout, dgold, diff, total;
	for(int y = 0; y<ny; ++y) {
		for(int x = 0; x<nx; ++x) {
			dout = labels[variableIndex(x,y)];
			dgold = (gold[y][x].value)/goldScale;
			diff = dout-dgold;
			diff = diff>0?diff:-diff;
			total += diff>threshold?1:0;
		}
	}
	return total/(nx*ny);
}
#endif
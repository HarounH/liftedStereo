#ifndef _EVALUATOR_HPP_
#define _EVALUATOR_HPP_

// Systems includes... I/O stuff and PNG include. 
#include <iostream>
#include <fstream>
#include <png++/png.hpp>
#include <stdlib.h>
#include <cmath> // for sqrt()

double rmsScore(
	// arguments
		png::image< png::ga_pixel >& output,
		png::image< png::ga_pixel >& gold,
		size_t nDisparities,
		double threshold
	) { // function begins
	size_t nx = output.get_width();
	size_t ny = output.get_height();
	
	int total = 0;
	double dOut, dGold, diff;
	for(int x=0; x<nx; ++x) {
		for(int y=0; y<ny; ++y) {
			dOut = nDisparities*((output[y][x].value)/255.0);
			dGold = nDisparities*((gold[y][x].value)/255.0);
			diff = dOut-dGold;
			total += diff*diff;
		}
	}
	return (1.0 - (sqrt(total/(nx*ny))));

}

// incorrect function. avoid at all costs.
double badPixelScore(
	// arguments
		png::image< png::ga_pixel >& output,
		png::image< png::ga_pixel >& gold,
		size_t nDisparities,
		double threshold
	) { // function begins
	size_t nx = output.get_width();
	size_t ny = output.get_height();
	
	int total = 0;
	double dOut, dGold, diff;
	for(int x=0; x<nx; ++x) {
		for(int y=0; y<ny; ++y) {
			dOut = nDisparities*((output[y][x].value)/255.0);
			dGold = nDisparities*((gold[y][x].value)/255.0);
			diff = dOut-dGold;
			diff = diff>0.0?diff:-diff;
			total += diff>threshold?1:0;
		}
	}
	return (1.0 - (double(total)/(nx*ny)));
}

double calculatePixelError(
	// arguments
		png::image< png::ga_pixel >& output,
		png::image< png::ga_pixel >& gold,
		size_t nOutputDisparities,
		size_t goldScale,
		double threshold // the threshold beyond which a pixel is counted as bad
	) {
		size_t nx = output.get_width();
		size_t ny = output.get_height();
		
		
		double dOut, dGold, diff, total=0.0;

		// iterate over pixels.
		for(int x=0; x<nx; ++x) {
			for(int y=0; y<ny; ++y) {
				dOut = nOutputDisparities*((output[y][x].value)/255.0);
				dGold = goldScale*((gold[y][x].value));
				diff = dOut-dGold;
				diff = diff>0.0?diff:-diff;
				total += diff>threshold?1:0;
			}
		}
		return total/(nx*ny);
}

#endif
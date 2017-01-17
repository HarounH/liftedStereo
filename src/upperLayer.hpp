#ifndef _UPPER_LAYER_HPP_
#define _UPPER_LAYER_HPP_

// Systems includes... I/O stuff and PNG include. 
#include <iostream>
#include <fstream>
#include <png++/png.hpp>
#include <stdlib.h>

#include "Image.h"
#include "Image.cpp"
#include "permutohedral.h"

// Implementes bilateral filtering. 
// Following the TSGO naming scheme as religiously as possible.
void bilateral_filter(
	// arguments
		png::image< png::rgb_pixel >& img,
		std::vector<std::vector<std::vector<double> > >& cost,
		double sigmax, double sigmaf
	) {
	

	int ny = img.get_height();
	int nx = img.get_width();
	int nLabels = cost[0][0].size();
	
	double invSigmax = 1.0/sigmax;
	double invSigmaf = 1.0/sigmaf;

	Image ref(1, nx, ny, 5); // 5 channels.
	Image cost_image(1, nx, ny, nLabels); // one channel per disparity label.
	
	for(int y=0; y<ny; ++y) {
		for(int x=0; x<nx; ++x) {
			// set up ref.
			ref(x,y)[0] = x*invSigmax;
			ref(x,y)[1] = y*invSigmax;
			ref(x,y)[2] = img[y][x].red * invSigmaf;
			ref(x,y)[3] = img[y][x].green * invSigmaf;
			ref(x,y)[4] = img[y][x].blue * invSigmaf;
			
			// set up cost.
			for(int l=0; l<nLabels; ++l) {
				cost_image(x,y)[l] = cost[y][x][l];
			}
		}
	}

	// invoke filter.
	Image filtered_costs = PermutohedralLattice::filter(cost_image, ref);

	// translate back.
	double ncost;
	for(int y=0; y<ny; ++y) {
		for(int x=0; x<nx; ++x) {
			for(int l=0; l<nLabels; ++l) {
				ncost = filtered_costs(x,y)[l];
				// cout << "Filter(" << cost[y][x][l] << ")=" << ncost << endl;
				cost[y][x][l] = ncost;
			}
		}
	}
}
#endif
#ifndef _GAUSSIAN_HPP_
#define _GAUSSIAN_HPP_

#include <algorithm>
#include <cmath>


void gaussianFilterCosts(
	// arguments
		png::image< png::rgb_pixel >& img, // reference image
		std::vector<std::vector<std::vector<double> > >& costs,
		int radius,
		double sqsigmax,
		double sqsigmaf
	) { // function begins
	int ny = costs.size();
	int nx = costs[0].size();
	int nLabels = costs[0][0].size();

	// Things we need for a gaussian filter.
	std::vector<double> gaussx(radius+1,0.0);
	for(int i=0; i<gaussx.size(); ++i) {
		gaussx[i] = exp(-0.5*(i*i)/sqsigmax);
	}
	std::vector<double> gaussf(255*3, 0.0);
	for(int i=0; i<gaussf.size(); ++i) {
		gaussf[i] = exp(-0.5*(i*i)/sqsigmaf);
	}

	// the filtration itself.
	std::vector<std::vector<std::vector<double> > > ocosts(costs); // make a copy.
	int ylower, yupper, xlower, xupper;
	double w, tw;
	for(int y1=0; y1<ny; ++y1) {
		for(int x1=0; x1<nx; ++x1) {
			for(int l=0; l<nLabels; ++l) {
				costs[y1][x1][l] = 0.0;
				tw = 0.0;
				// Walk over neighbourhood.
				ylower = std::max(y1-radius,0);
				yupper = std::min(ny,y1+radius+1);
				xlower = std::max(x1-radius,0);
				xupper = std::min(nx,x1+radius+1);				
				for(int y2= ylower; y2<yupper; ++y2) {
					for(int x2=xlower; x2<xupper; ++x2) {
						w = gaussx[abs(x2-x1)]*gaussx[abs(y2-y1)]*gaussf[absolute_difference_color(img[y1][x1], img[y2][x2])];
						tw += w;
						costs[y1][x1][l] += w*ocosts[y2][x2][l];
					}
				}
				costs[y1][x1][l] /= tw;
			}
		}
	}
}




#endif
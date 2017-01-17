#ifndef UNARY_POTENTIALS
#define UNARY_POTENTIALS

// Systems includes... I/O stuff and PNG include. 
#include <iostream>
#include <fstream>
#include <png++/png.hpp>
#include <stdlib.h>
#include <cmath>
#include <algorithm>

// OpenGM includes
#include <opengm/opengm.hxx>
#include <opengm/graphicalmodel/graphicalmodel.hxx>
#include <opengm/operations/minimizer.hxx>
#include <opengm/operations/adder.hxx>
#include <opengm/inference/graphcut.hxx>
#include <opengm/inference/alphaexpansion.hxx>
#include <opengm/inference/alphabetaswap.hxx>
#include <opengm/inference/alphaexpansionfusion.hxx>
#include <opengm/inference/auxiliary/minstcutboost.hxx>
#include <opengm/graphicalmodel/space/simplediscretespace.hxx>
#include <opengm/functions/potts.hxx>
#include <opengm/functions/pottsn.hxx>
#include <opengm/functions/truncated_absolute_difference.hxx>
#include <opengm/inference/messagepassing/messagepassing.hxx>
#include <opengm/inference/gibbs.hxx>

#include <typeinfo>

using namespace std;
#define N_LABELS 15
#define UNARY_POTENTIAL_RADIUS 7
#define RANK_SCALE 1000.0

size_t nLabels = N_LABELS; // Modify this.
size_t unaryPotentialRadius = UNARY_POTENTIAL_RADIUS; // modify this.
double inv_255square = 1.0/(255.0*255.0);
// quick and dirty function to get difference between two pixels.
inline double diff_square(png::rgb_pixel const& l, png::rgb_pixel const& r) {
	return (l.red - r.red)*(l.red - r.red) + (l.green - r.green)*(l.green - r.green) + (l.blue - r.blue)*(l.blue - r.blue);
}

// quick and dirty helper function to calculate intensity of a pixel.
inline double unnormalized_squared_intensity(png::rgb_pixel const& p) {
	return (p.red*p.red) + (p.green*p.green) + (p.blue*p.blue);
}
inline double unnormalized_intensity(png::rgb_pixel const& p) {
	return sqrt((p.red*p.red) + (p.green*p.green) + (p.blue*p.blue));
}

// computes a unary potential function, using the left right images and UNARY_POTENTIAL_RADIUS
double unaryPotentialMeanAbsoluteDifference(png::image< png::rgb_pixel >& left, png::image< png::rgb_pixel >& right, size_t x, size_t y, size_t d) {
	// TODO try different variants of this.
	double potential = 0.0;
	double redDiff = 0.0;
	double greenDiff = 0.0;
	double blueDiff = 0.0;
	size_t nx = left.get_width();
	size_t ny = right.get_height();
	size_t count=1;
	for(size_t xiter = x - unaryPotentialRadius + 1; xiter<= x+unaryPotentialRadius-1; ++xiter) {
		for(size_t yiter = y - unaryPotentialRadius + 1; yiter <= y+unaryPotentialRadius-1; ++yiter) {
			if ((yiter>0) && (yiter<ny) && (xiter>0) && (xiter<nx) && ((xiter-d)>0) && ((xiter-d)<nx)) {
				auto lpixel = left[yiter][xiter];
				auto rpixel = right[yiter][xiter-d];
				redDiff += abs(lpixel.red - rpixel.red);
				greenDiff += abs(lpixel.green - rpixel.green);
				blueDiff += abs(lpixel.blue - rpixel.blue);
				count++;
			}
		}
	}
	potential += redDiff*redDiff + greenDiff*greenDiff + blueDiff*blueDiff;
	return potential/(count*count);
}



double unaryPotentialMeanFilter(png::image< png::rgb_pixel >& left, png::image< png::rgb_pixel >& right, size_t& x, size_t& y, size_t& d) {
	// constants needed
	size_t nDim = 3;
	size_t nx = left.get_width();
	size_t ny = left.get_height();

	std::vector<double> leftMean(nDim, 0.0);
	std::vector<double> rightMean(nDim, 0.0);
	std::vector<size_t> nPixels(2, 0);
	// iterate over window.
	for(size_t xiter = x - (unaryPotentialRadius-1); xiter <= x + (unaryPotentialRadius-1); ++xiter) {
		for(size_t yiter = y - (unaryPotentialRadius-1); yiter <= y + (unaryPotentialRadius-1); ++yiter) {
			// get left mean
			if ((xiter>0) && (xiter<nx) && (yiter>0) && (yiter<ny)) {
				nPixels[0] += 1;
				leftMean[0] += left[yiter][xiter].red;
				leftMean[1] += left[yiter][xiter].green;
				leftMean[2] += left[yiter][xiter].blue;
			}

			// get right mean
			if (((xiter-d)>0) && ((xiter-d)<nx) && (yiter>0) && (yiter<ny)) {
				nPixels[0] += 1;
				rightMean[0] += right[yiter][xiter-d].red;
				rightMean[1] += right[yiter][xiter-d].green;
				rightMean[2] += right[yiter][xiter-d].blue;
			}
		}
	}
	double ans = 0.0;
	for(size_t dim=0; dim<nDim; ++dim) {
		leftMean[dim] /= nPixels[0];
		rightMean[dim] /= nPixels[1];
		ans += (leftMean[dim]-rightMean[dim])*(leftMean[dim]-rightMean[dim]);
	}
	return ans/(nDim*nDim);
}

// Rank/BT as of the cost evaluation paper.
double unaryPotentialIntensityRankMeanAbsoluteDifference(png::image< png::rgb_pixel >& left, png::image< png::rgb_pixel >& right, size_t& x, size_t& y, size_t& d) {
	// constants needed
	size_t nDim = 3;
	size_t nx = left.get_width();
	size_t ny = left.get_height();


	double lIntensity = unnormalized_squared_intensity(left[y][x]);
	double rIntensity;
	if ((x+d)<nx) {
		rIntensity = unnormalized_squared_intensity(right[y][x-d]);		
	} else {
		return 0.0; // FIXME what do i return here?
	}

	int lgreater = 0;
	int rgreater = 0;
	int windowSize = 0;
	double newLIntensity, newRIntensity;
	// iterate over window.
	for(size_t xiter = x - (unaryPotentialRadius-1); xiter <= x + (unaryPotentialRadius-1); ++xiter) {
		for(size_t yiter = y - (unaryPotentialRadius-1); yiter <= y + (unaryPotentialRadius-1); ++yiter) {
			if ((xiter>0) && (xiter<nx) && (yiter>0) && (yiter<ny)
				&& ((xiter-d)>0) && ((xiter-d)<nx) && (yiter>0) && (yiter<ny) ) {

				windowSize++;
				newLIntensity = unnormalized_squared_intensity(left[yiter][xiter]);
				newRIntensity = unnormalized_squared_intensity(right[yiter][xiter-d]);

				if (newLIntensity>lIntensity) {
					lgreater++;
				}
				if (newRIntensity>rIntensity) {
					rgreater++;
				}
			}
		}
	}
	if (windowSize==0) return 0.0;
	return (RANK_SCALE*(1.0 + abs(double(rgreater-lgreater)))/windowSize);
}


double unaryPotentialSumIntensityDifference(png::image< png::rgb_pixel >& left, png::image< png::rgb_pixel >& right, size_t& x, size_t& y, size_t& d) {
	// constants needed
	size_t nDim = 3;
	size_t nx = left.get_width();
	size_t ny = left.get_height();

	double ans = 0.0;
	int windowSize=0;
	// iterate over window.
	for(size_t xiter = x - (unaryPotentialRadius-1); xiter <= x + (unaryPotentialRadius-1); ++xiter) {
		for(size_t yiter = y - (unaryPotentialRadius-1); yiter <= y + (unaryPotentialRadius-1); ++yiter) {
			if ((xiter>0) && (xiter<nx) && (yiter>0) && (yiter<ny)
				&& ((xiter-d)>0) && ((xiter-d)<nx) && (yiter>0) && (yiter<ny) ) {

				windowSize++;
				
				ans += abs(unnormalized_intensity(left[yiter][xiter]) - unnormalized_intensity(right[yiter][xiter-d]));
			}
		}
	}
	if (windowSize==0) return 0.0;
	return ans/(windowSize*inv_255square);
}

// pixel dissimalirty invariant to image sampling.
double unaryPotentialBirchfieldNoWindow(png::image< png::rgb_pixel >& left, png::image< png::rgb_pixel >& right, size_t& x, size_t& y, size_t& d) {
	double ilm, il, ilp, irm, ir, irp;

	int nx = left.get_width();

	il = unnormalized_intensity(left[y][x]);
	ir = unnormalized_intensity(right[y][x-d]);
	// left stuff.
	if (x>0) {
		ilm = unnormalized_intensity(left[y][x-1]);
	} else {
		ilm = il;
	}
	if (x<(nx-1)) {
		ilp = unnormalized_intensity(left[y][x+1]);
	} else {
		ilp = il;
	}
	// right stuff.
	if (x>d) {
		irm = unnormalized_intensity(right[y][x-d-1]);
	} else {
		irm = ir;
	}
	if (x<(d+nx-1)) {
		irp = unnormalized_intensity(right[y][x-d+1]);
	} else {
		irp = ir;
	}

	// calculating the value itself.
	double ilmax, ilmin, irmax, irmin, dl, dr;
	ilmax = max(il , max(0.5*(il + ilp) , 0.5*(il + ilm)) );
	ilmin = max(il , max(0.5*(il + ilp) , 0.5*(il + ilm)) );
	irmax = max(ir , max(0.5*(ir + irp) , 0.5*(ir + irm)) );
	irmin = max(ir , max(0.5*(ir + irp) , 0.5*(ir + irm)) );
	dl = max(0.0, max(il - irmax, irmin - il));
	dr = max(0.0, max(ir - ilmax, ilmin - ir));
	return min(dl, dr);
}






double unaryPotentialBirchwoodTomiskiWindowed(png::image< png::rgb_pixel >& left, png::image< png::rgb_pixel >& right, size_t& x, size_t& y, size_t& d) {
	// constants needed
	size_t nDim = 3;
	size_t nx = left.get_width();
	size_t ny = left.get_height();

	double ans = 0.0;
	int windowSize=0;

	// iterate over window.
	for(size_t xiter = x - (unaryPotentialRadius-1); xiter <= x + (unaryPotentialRadius-1); ++xiter) {
		for(size_t yiter = y - (unaryPotentialRadius-1); yiter <= y + (unaryPotentialRadius-1); ++yiter) {
			if ((xiter>0) && (xiter<nx) && (yiter>0) && (yiter<ny)
				&& ((xiter-d)>0) && ((xiter-d)<nx) && (yiter>0) && (yiter<ny) ) {

				windowSize++;
				ans += unaryPotentialBirchfieldNoWindow(left, right, xiter, yiter, d);
			}
		}
	}
	if (windowSize==0) return 0.0;
	return ans/(windowSize);
}
#endif
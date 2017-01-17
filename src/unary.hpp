#ifndef _UNARY_HPP_
#define _UNARY_HPP_

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
const double one_by_sqrt2 = 1.0/sqrt(2.0);

// ERF constants
const double erf_a1 =  0.254829592;
const double erf_a2 = -0.284496736;
const double erf_a3 =  1.421413741;
const double erf_a4 = -1.453152027;
const double erf_a5 =  1.061405429;
const double erf_p  =  0.3275911;

double erf(double z) { // computes the gauss error function of z.
	double sign=1.0;
	if (z<0) sign = -1.0;
	z = fabs(z);
	// Famous formula, from stackoverflow for ERF.
	double t = 1.0/(1.0 + erf_p*z);
	double y = 1.0 - ((((((erf_a5*t + erf_a4)*t) + erf_a3)*t + erf_a2)*t + erf_a1)*t*exp(-z*z));
	// cout << sign*y << endl;
	return sign*y;
}


inline double square_intensity(png::rgb_pixel const& p) {
	return (p.red*p.red) + (p.green*p.green) + (p.blue*p.blue);
}
inline double sum_of_colors(png::rgb_pixel const& p) {
	return p.red + p.green + p.blue;
}
inline int absolute_difference_color(png::rgb_pixel& p, png::rgb_pixel& q) {
	return abs(p.red - q.red) + abs(p.green - q.green) + abs(p.blue - q.blue);
}

// function to calculate minmean, as per tsgo.
double minmean(std::vector<std::vector<std::vector<double> > >& cost) {
	int ny = cost.size();
	int nx = cost[0].size();
	int nLabels = cost[0][0].size();
	double mean = 0.0;
	double min;
	for(int y=0; y<ny; ++y) {
		for(int x=0; x<nx; ++x) {
			min = cost[y][x][0];
			for(int l=0; l<nLabels; ++l) {
				min = (cost[y][x][l]<min)?cost[y][x][l]:min;
			}
			mean += min;
		}
	}
	return mean/(nx)/(ny);
}

// function to calculate minvar, as per tsgo.
double minvar(std::vector<std::vector<std::vector<double> > >& cost) {
	int ny = cost.size();
	int nx = cost[0].size();
	int nLabels = cost[0][0].size();
	
	double min, ans;
	ans = 0.0;
	for(int y=0; y<ny; ++y) {
		for(int x=0; x<nx; ++x) {
			min = cost[y][x][0];
			for(int l=0; l<nLabels; ++l) {
				min = (cost[y][x][l]<min)?cost[y][x][l]:min;
			}

			for(int l=0; l<nLabels; ++l) {
				ans += ( cost[y][x][l] - min );
			}
		}
	}
	return ans/(nLabels)/(nx)/(ny);
}

void doERFTransform(std::vector<std::vector<std::vector<double> > >& cost, bool makeInt = false) {
	double fancyv, theta;
	int ny = cost.size();
	int nx = cost[0].size();
	int nLabels = cost[0][0].size();
	theta = minmean(cost);
	fancyv = minvar(cost);
	fancyv = 0.00095*(fancyv-theta)*(fancyv-theta);
	if(fancyv>2.5) fancyv*=0.91; // what even. EDP.cpp does it, so will i.
	cout << "theta: " << theta << endl;
	cout << "fancyv:" << fancyv << endl;
	
	if(makeInt) { // Keeping it here is ugly, but avoids x*y*l checks.
		for(int y=0; y<ny; ++y) {
			for(int x=0; x<nx; ++x) {
				for(int l=0; l<nLabels; ++l) {
					// cout << "erf(" << fancyv*(cost[y][x][l]- theta)/theta << ")=";
					cost[y][x][l] = 0.5*(1.0 + erf(fancyv*(cost[y][x][l]- theta)/theta));
					cost[y][x][l] = (cost[y][x][l]>0.5)?1:0;
				}
			}
		}
		return;		
	} else {
		for(int y=0; y<ny; ++y) {
			for(int x=0; x<nx; ++x) {
				for(int l=0; l<nLabels; ++l) {
					// cout << "erf(" << fancyv*(cost[y][x][l]- theta)/theta << ")=";
					cost[y][x][l] = 0.5*(1.0 + erf(fancyv*(cost[y][x][l]- theta)/theta));
				}
			}
		}
		return;		
	}
}

/* this function replicates Get_mul_gr in TSGO:EDP.cpp */
double getGradientMultiplier(
	// arguments
		const png::image< png::rgb_pixel >& img,
		const std::vector<std::vector<double> > & grad_xr,
		const std::vector<std::vector<double> > & grad_xg,
		const std::vector<std::vector<double> > & grad_xb,
		const std::vector<std::vector<double> > & grad_yr,
		const std::vector<std::vector<double> > & grad_yg,
		const std::vector<std::vector<double> > & grad_yb
	) {
	int nx = img.get_width();
	int ny = img.get_height();
	// deviation of img
	double devImage = 0.0;
	std::vector<double> vecMeans(3,0.0); // Means for color.
	for(int y=0; y< ny; ++y) {
		for(int x=0; x<nx; ++x) {
			vecMeans[0] += img[y][x].red;
			vecMeans[1] += img[y][x].green;
			vecMeans[2] += img[y][x].blue;
		}
	}
	for(int i=0; i<vecMeans.size(); ++i) {
		vecMeans[i] /= ny;
		vecMeans[i] /= nx;
	}
	for(int y=0; y< ny; ++y) {
		for(int x=0; x<nx; ++x) {
			devImage += (vecMeans[0] - img[y][x].red)*(vecMeans[0] - img[y][x].red);
			devImage += (vecMeans[1] - img[y][x].green)*(vecMeans[1] - img[y][x].green);
			devImage += (vecMeans[2] - img[y][x].blue)*(vecMeans[2] - img[y][x].blue);
		}
	}
	devImage = sqrt(devImage/ny/nx);

	// "deviation" of gradients
	double devGrad = 0.0;
	for(int y=0; y< ny; ++y) {
		for(int x=0; x<nx; ++x) {
			devGrad += (grad_xr[y][x]*grad_xr[y][x]);
			devGrad += (grad_xg[y][x]*grad_xg[y][x]);
			devGrad += (grad_xb[y][x]*grad_xb[y][x]);
			devGrad += (grad_yr[y][x]*grad_yr[y][x]);
			devGrad += (grad_yg[y][x]*grad_yg[y][x]);
			devGrad += (grad_yb[y][x]*grad_yb[y][x]);
		}
	}
	// cout << "devImage=" << devImage << " devGrad=" << devGrad << endl;
	devGrad = sqrt(devGrad/ny/nx);
	return (devImage/devGrad);
}
// a function that computes all gradients for each pixel.
// It returns the multiplier.
double getAllGradients(
	// arguments
		png::image< png::rgb_pixel >& img,
		std::vector<std::vector<double> > & grad_xr,
		std::vector<std::vector<double> > & grad_xg,
		std::vector<std::vector<double> > & grad_xb,
		std::vector<std::vector<double> > & grad_yr,
		std::vector<std::vector<double> > & grad_yg,
		std::vector<std::vector<double> > & grad_yb,
		double multiplier = 0.0
	) { // function begins
	int nx = img.get_width();
	int ny = img.get_height();
	int xf,yf,xb,yb; // the positions using which to calculate gradients.
	double dx, dy, dx1, dy1;
	double temp;
	double bias;

	for(int y=0; y<ny; ++y) {
		for (int x = 0; x < nx; ++x) {
			xf = (x+1 < nx)?x+1 : nx-1;
			yf = (y+1 < ny)?y+1 : ny-1;
			xb = (x-1>=0)?(x-1):0;
			yb = (y-1>=0)?(y-1):0;

			// need to do stuff for each color.
			// red
			dx = img[y][xf].red-img[y][xb].red;
			dy = img[yf][x].red - img[yb][x].red;
			dx1 = 0.5*(img[yb][xf].red - img[yf][xb].red);
			dy1 = 0.5*(img[yf][xf].red - img[yb][xb].red);
			grad_xr[y][x] = 4*(dx + (dx1+dy1));
			grad_yr[y][x] = 4*(dy + (-dx1+dy1));

			// blue
			dx = img[y][xf].blue-img[y][xb].blue;
			dy = img[yf][x].blue - img[yb][x].blue;
			dx1 = 0.5*(img[yb][xf].blue - img[yf][xb].blue);
			dy1 = 0.5*(img[yf][xf].blue - img[yb][xb].blue);
			grad_xb[y][x] = 4*(dx + (dx1+dy1));
			grad_yb[y][x] = 4*(dy + (-dx1+dy1));
			
			// green
			dx = img[y][xf].green-img[y][xb].green;
			dy = img[yf][x].green - img[yb][x].green;
			dx1 = 0.5*(img[yb][xf].green - img[yf][xb].green);
			dy1 = 0.5*(img[yf][xf].green - img[yb][xb].green);
			grad_xg[y][x] = 4*(dx + (dx1+dy1));
			grad_yg[y][x] = 4*(dy + (-dx1+dy1));			
		}
	} // End of equivalent of GrBuf().

	// sqrt everything.
	for(int y=0; y<ny; ++y) {
		for(int x=0; x<nx; ++x) {
			grad_xr[y][x] = (grad_xr[y][x]>0)?(sqrt(grad_xr[y][x])):(-sqrt(-grad_xr[y][x]));
			grad_yr[y][x] = (grad_yr[y][x]>0)?(sqrt(grad_yr[y][x])):(-sqrt(-grad_yr[y][x]));
			grad_xb[y][x] = (grad_xb[y][x]>0)?(sqrt(grad_xb[y][x])):(-sqrt(-grad_xb[y][x]));
			grad_yb[y][x] = (grad_yb[y][x]>0)?(sqrt(grad_yb[y][x])):(-sqrt(-grad_yb[y][x]));
			grad_xg[y][x] = (grad_xg[y][x]>0)?(sqrt(grad_xg[y][x])):(-sqrt(-grad_xg[y][x]));
			grad_yg[y][x] = (grad_yg[y][x]>0)?(sqrt(grad_yg[y][x])):(-sqrt(-grad_yg[y][x]));
		}
	}

	if(multiplier==0.0) {
		multiplier = getGradientMultiplier(img, 
											grad_xr,
											grad_yr,
											grad_xb,
											grad_yb,
											grad_xg,
											grad_yg);		
	}
	// cout << "gradient multiplier=" << multiplier << endl;
	bias = 128;

	// threshold everything
	for(int y=0; y<ny; ++y) {
		for(int x=0; x<nx; ++x) {
			grad_xr[y][x] *= multiplier;
			grad_xr[y][x] += bias;
			grad_xr[y][x] = (grad_xr[y][x]<0)?(0):(grad_xr[y][x]>255?255:grad_xr[y][x]);
			grad_yr[y][x] *= multiplier;
			grad_yr[y][x] += bias;
			grad_yr[y][x] = (grad_yr[y][x]<0)?(0):(grad_yr[y][x]>255?255:grad_yr[y][x]);
			grad_xb[y][x] *= multiplier;
			grad_xb[y][x] += bias;
			grad_xb[y][x] = (grad_xb[y][x]<0)?(0):(grad_xb[y][x]>255?255:grad_xb[y][x]);
			grad_yb[y][x] *= multiplier;
			grad_yb[y][x] += bias;
			grad_yb[y][x] = (grad_yb[y][x]<0)?(0):(grad_yb[y][x]>255?255:grad_yb[y][x]);
			grad_xg[y][x] *= multiplier;
			grad_xg[y][x] += bias;
			grad_xg[y][x] = (grad_xg[y][x]<0)?(0):(grad_xg[y][x]>255?255:grad_xg[y][x]);
			grad_yg[y][x] *= multiplier;
			grad_yg[y][x] += bias;
			grad_yg[y][x] = (grad_yg[y][x]<0)?(0):(grad_yg[y][x]>255?255:grad_yg[y][x]);
		
			// cout << grad_xr[y][x] << " "
			// 	 <<	grad_yr[y][x] << " "
			// 	 <<	grad_xb[y][x] << " "
			// 	 <<	grad_yb[y][x] << " "
			// 	 <<	grad_xg[y][x] << " "
			// 	 <<	grad_yg[y][x] << endl;
		}
	}

	return multiplier;
}


/*
	Function to get dissimilarity between 
	limg[y][x] and rimg[y][x+l]

	The function computes
	min( (birchfield) ,tau_l)

	We use intensity instead of color itself.
*/
double getTSGOIntensityCost_faulty(
	// arguments
		png::image< png::rgb_pixel >& left,
		png::image< png::rgb_pixel >& right,
		int x,
		int y,
		int l, // displacement
		double taul=90.0
	) {
	double ilm, il, ilp, irm, ir, irp;
	int nx = left.get_width();
	il = sum_of_colors(left[y][x]); // sqrt(square_intensity(left[y][x]));
	if ((x-l)>=0 && (x-l)<nx) {
		ir = sum_of_colors(right[y][x-l]); // sqrt(square_intensity(right[y][x-l]));
	} else {
		return taul;
	}
	// left stuff.
	if (x>0) {
		ilm = sum_of_colors(left[y][x-1]); // sqrt(square_intensity(left[y][x-1]));
	} else {
		ilm = il;
	}
	if (x<(nx-1)) {
		ilp = sum_of_colors(left[y][x+1]); // sqrt(square_intensity(left[y][x+1]));
	} else {
		ilp = il;
	}
	// right stuff.
	if (((l>=0) && (x>l)) || ((l<0) && (x<(l+nx-1)))) {
		irm = sum_of_colors(right[y][x-l-1]); // sqrt(square_intensity(right[y][x-l-1]));
	} else {
		irm = ir;
	}
	if (((l<0) && (x>l)) || ((l>=0) && (x<(l+nx-1)))) {
		irp = sum_of_colors(right[y][x-l+1]); // sqrt(square_intensity(right[y][x-l+1]));
	} else {
		irp = ir;
	}

	// calculating the value itself.
	double ilmax, ilmin, irmax, irmin, dl, dr;
	ilmax = std::max(il , std::max(0.5*(il + ilp) , 0.5*(il + ilm)) );
	ilmin = std::max(il , std::max(0.5*(il + ilp) , 0.5*(il + ilm)) );
	irmax = std::max(ir , std::max(0.5*(ir + irp) , 0.5*(ir + irm)) );
	irmin = std::max(ir , std::max(0.5*(ir + irp) , 0.5*(ir + irm)) );
	dl = std::max(0.0, std::max(il - irmax, irmin - il));
	dr = std::max(0.0, std::max(ir - ilmax, ilmin - ir));
	return std::min(/*std::min(dl, dr)*/dl, taul);
}

/*
	Computes 
		min( taul, 
			min_l-0.5<=d<=l+0.5 (\sum_(c \in {r,g,b}) |left[y][x].c - right[y][x-d].c| )
			)
		That is similar to Birchfield, but its not exactly birchfield.
*/
double getTSGOIntensityCost(
	// arguments
		png::image< png::rgb_pixel >& left,
		png::image< png::rgb_pixel >& right,
		int x,
		int y,
		int l, // displacement
		double taul=90.0,
		bool birchfield = false
	) {
	double val = 0.0; // return value, essentially.
	int nx = left.get_width();
	int ny = left.get_height();
	int xr; // x-l, essentially.
	int xp, xrp; // x+1, xr+1 essentially.
	// some preparation for indices.
	xr = x-l;
	if(xr<0) xr=0;
	if(xr>=nx) xr = nx-1;

	xp = x+1;
	if(xp<0) xp=0;
	if(xp>=nx) xp = nx-1;

	xrp = xr+1;
	if(xrp<0) xrp=0;
	if(xrp>=nx) xrp = nx-1;


	if(!birchfield) { // go for a far more simplistic energy.
		val += abs(left[y][x].red - right[y][xr].red);
		val += abs(left[y][x].green - right[y][xr].green);
		val += abs(left[y][x].blue - right[y][xr].blue);
		return std::min(val, taul);
	}

	// val should become
	// 		min_l-0.5<=d<=l+0.5 (\sum_(c \in {r,g,b}) |left[y][x].c - right[y][x-d].c|
	std::vector<double> sumdfb(4,0.0); // As EDP.cpp calls it.

	// red
	sumdfb[0] += abs(left[y][x].red + left[y][xp].red - right[y][x].red - right[y][xrp].red)/2.0;
	sumdfb[1] += abs(left[y][x].red - right[y][xr].red);
	sumdfb[2] += abs((left[y][x].red + left[y][xp].red)/2.0 - right[y][xr].red);
	sumdfb[3] += abs(left[y][x].red - (right[y][xr].red - right[y][xrp].red)/2.0);
	// green
	sumdfb[0] += abs(left[y][x].green + left[y][xp].green - right[y][x].green - right[y][xrp].green)/2.0;
	sumdfb[1] += abs(left[y][x].green - right[y][xr].green);
	sumdfb[2] += abs((left[y][x].green + left[y][xp].green)/2.0 - right[y][xr].green);
	sumdfb[3] += abs(left[y][x].red - (right[y][xr].red - right[y][xrp].red)/2.0);
	// blue
	sumdfb[0] += abs(left[y][x].blue + left[y][xp].blue - right[y][x].blue - right[y][xrp].blue)/2.0;
	sumdfb[1] += abs(left[y][x].blue - right[y][xr].blue);
	sumdfb[2] += abs((left[y][x].blue + left[y][xp].blue)/2.0 - right[y][xr].blue);
	sumdfb[3] += abs(left[y][x].red - (right[y][xr].red - right[y][xrp].red)/2.0);
	val = std::min( std::min(sumdfb[2], sumdfb[3]), std::min(sumdfb[0], sumdfb[1]) );
	return std::min(taul, val);
}



double getTSGOGradientCost(
	// arguments
		std::vector<std::vector<std::vector<double> > >& grads_xr,
		std::vector<std::vector<std::vector<double> > >& grads_xg,
		std::vector<std::vector<std::vector<double> > >& grads_xb,
		std::vector<std::vector<std::vector<double> > >& grads_yr,
		std::vector<std::vector<std::vector<double> > >& grads_yg,
		std::vector<std::vector<std::vector<double> > >& grads_yb,
		bool lr, // or rl ... represents which direction to get the cost for.
		int x, int y,
		int l,
		double taug=180.0,
		bool birchfield = false
	) {
	int fromid, toid;
	if (lr) {
		fromid=0;
		toid=1;
	} else {
		fromid = 1;
		toid = 0;
	}
	int nx = grads_xr[0][0].size();
	if(!birchfield) {
		if ((x-l) > 0 && ((x-l) < nx)) {
			return std::min(
					fabs( grads_xr[fromid][y][x] - grads_xr[toid][y][x-l] ) +
					fabs( grads_xg[fromid][y][x] - grads_xg[toid][y][x-l] ) +
					fabs( grads_xb[fromid][y][x] - grads_xb[toid][y][x-l] ) +
					fabs( grads_yr[fromid][y][x] - grads_yr[toid][y][x-l] ) +
					fabs( grads_yg[fromid][y][x] - grads_yg[toid][y][x-l] ) +
					fabs( grads_yb[fromid][y][x] - grads_yb[toid][y][x-l] ),
					taug);
		} else {
			return taug;
		}
	} else { // Birchfield.
		int xr; // x-l, essentially.
		int xp, xrp; // x+1, xr+1 essentially.
		// some preparation for indices.
		xr = x-l;
		if(xr<0) xr=0;
		if(xr>=nx) xr = nx-1;

		xp = x+1;
		if(xp<0) xp=0;
		if(xp>=nx) xp = nx-1;

		xrp = xr+1;
		if(xrp<0) xrp=0;
		if(xrp>=nx) xrp = nx-1;

		double val = 0.0;

		std::vector<double> sumdfb(4,0.0); // As EDP.cpp calls it.

		// gotta add it up over 6 gradients.
		sumdfb[0] += abs(grads_xr[fromid][y][x] + grads_xr[fromid][y][xp] - grads_xr[toid][y][x] - grads_xr[toid][y][xrp])/2.0;
		sumdfb[1] += abs(grads_xr[fromid][y][x] - grads_xr[toid][y][xr]);
		sumdfb[2] += abs((grads_xr[fromid][y][x] + grads_xr[fromid][y][xp])/2.0 - grads_xr[toid][y][xr]);
		sumdfb[3] += abs(grads_xr[fromid][y][x] - (grads_xr[toid][y][xr] - grads_xr[toid][y][xrp])/2.0);

		sumdfb[0] += abs(grads_xg[fromid][y][x] + grads_xg[fromid][y][xp] - grads_xg[toid][y][x] - grads_xg[toid][y][xrp])/2.0;
		sumdfb[1] += abs(grads_xg[fromid][y][x] - grads_xg[toid][y][xr]);
		sumdfb[2] += abs((grads_xg[fromid][y][x] + grads_xg[fromid][y][xp])/2.0 - grads_xg[toid][y][xr]);
		sumdfb[3] += abs(grads_xg[fromid][y][x] - (grads_xg[toid][y][xr] - grads_xg[toid][y][xrp])/2.0);

		sumdfb[0] += abs(grads_xb[fromid][y][x] + grads_xb[fromid][y][xp] - grads_xb[toid][y][x] - grads_xb[toid][y][xrp])/2.0;
		sumdfb[1] += abs(grads_xb[fromid][y][x] - grads_xb[toid][y][xr]);
		sumdfb[2] += abs((grads_xb[fromid][y][x] + grads_xb[fromid][y][xp])/2.0 - grads_xb[toid][y][xr]);
		sumdfb[3] += abs(grads_xb[fromid][y][x] - (grads_xb[toid][y][xr] - grads_xb[toid][y][xrp])/2.0);

		sumdfb[0] += abs(grads_yr[fromid][y][x] + grads_yr[fromid][y][xp] - grads_yr[toid][y][x] - grads_yr[toid][y][xrp])/2.0;
		sumdfb[1] += abs(grads_yr[fromid][y][x] - grads_yr[toid][y][xr]);
		sumdfb[2] += abs((grads_yr[fromid][y][x] + grads_yr[fromid][y][xp])/2.0 - grads_yr[toid][y][xr]);
		sumdfb[3] += abs(grads_yr[fromid][y][x] - (grads_yr[toid][y][xr] - grads_yr[toid][y][xrp])/2.0);

		sumdfb[0] += abs(grads_yg[fromid][y][x] + grads_yg[fromid][y][xp] - grads_yg[toid][y][x] - grads_yg[toid][y][xrp])/2.0;
		sumdfb[1] += abs(grads_yg[fromid][y][x] - grads_yg[toid][y][xr]);
		sumdfb[2] += abs((grads_yg[fromid][y][x] + grads_yg[fromid][y][xp])/2.0 - grads_yg[toid][y][xr]);
		sumdfb[3] += abs(grads_yg[fromid][y][x] - (grads_yg[toid][y][xr] - grads_yg[toid][y][xrp])/2.0);

		sumdfb[0] += abs(grads_yb[fromid][y][x] + grads_yb[fromid][y][xp] - grads_yb[toid][y][x] - grads_yb[toid][y][xrp])/2.0;
		sumdfb[1] += abs(grads_yb[fromid][y][x] - grads_yb[toid][y][xr]);
		sumdfb[2] += abs((grads_yb[fromid][y][x] + grads_yb[fromid][y][xp])/2.0 - grads_yb[toid][y][xr]);
		sumdfb[3] += abs(grads_yb[fromid][y][x] - (grads_yb[toid][y][xr] - grads_yb[toid][y][xrp])/2.0);
		

		val = std::min( std::min(sumdfb[2], sumdfb[3]), std::min(sumdfb[0], sumdfb[1]) );
		return std::min(taug, val);
	}

}


void getTSGODissimilarityMeasures(
	// arguments
		png::image< png::rgb_pixel >& limg,
		png::image< png::rgb_pixel >& rimg,
		std::vector<std::vector<std::vector<double> > >& lcost,
		std::vector<std::vector<std::vector<double> > >& rcost,
		int nLabels
	) {
	// ASSERT: dimensions of lcost, rcost are nx*ny*nLabels
	int nx = limg.get_width();
	int ny = limg.get_height();
	
	std::vector<std::vector<std::vector<double> > > grads_xr(2, std::vector<std::vector<double> >(ny, std::vector<double>(nx, 0.0)));
	std::vector<std::vector<std::vector<double> > > grads_xg(2, std::vector<std::vector<double> >(ny, std::vector<double>(nx, 0.0)));
	std::vector<std::vector<std::vector<double> > > grads_xb(2, std::vector<std::vector<double> >(ny, std::vector<double>(nx, 0.0)));
	std::vector<std::vector<std::vector<double> > > grads_yr(2, std::vector<std::vector<double> >(ny, std::vector<double>(nx, 0.0)));
	std::vector<std::vector<std::vector<double> > > grads_yg(2, std::vector<std::vector<double> >(ny, std::vector<double>(nx, 0.0)));
	std::vector<std::vector<std::vector<double> > > grads_yb(2, std::vector<std::vector<double> >(ny, std::vector<double>(nx, 0.0)));
	double multiplier = getAllGradients(limg,
			grads_xr[0],
			grads_xg[0],
			grads_xb[0],
			grads_yr[0],
			grads_yg[0],
			grads_yb[0]
		);
	getAllGradients(rimg,
			grads_xr[1],
			grads_xg[1],
			grads_xb[1],
			grads_yr[1],
			grads_yg[1],
			grads_yb[1], multiplier
		);


	// unary costs
	std::vector<std::vector<std::vector<double> > > lucost(lcost);
	std::vector<std::vector<std::vector<double> > > rucost(rcost);

	// gradient costs
	std::vector<std::vector<std::vector<double> > > lgcost(lcost);
	std::vector<std::vector<std::vector<double> > > rgcost(rcost);
	
	double laplha, ralpha;
	double lumin, rumin, lgmin, rgmin;
	double lu, ru, lg, rg;

	for(int l=0; l<nLabels; ++l) {
		for(int y=0; y<ny; ++y) {
			for(int x=0; x<nx; ++x) {
				// left
					// intensity
					lucost[y][x][l] = getTSGOIntensityCost(limg, rimg, x, y, l);
					// gradients
					lgcost[y][x][l] = getTSGOGradientCost(
												grads_xr,
												grads_xg,
												grads_xb,
												grads_yr,
												grads_yg,
												grads_yb,
												true, 
												x,y,
												l);
				// right
					// intensity
					rucost[y][x][l] = getTSGOIntensityCost(rimg, limg, x, y, -l);
					// gradients
					rgcost[y][x][l] = getTSGOGradientCost(
												grads_xr,
												grads_xg,
												grads_xb,
												grads_yr,
												grads_yg,
												grads_yb,
												false, 
												x,y,
												-l);				
			}
		}
	}

	// need to calculate alphas, and for that, we need to calculate the minimas. as per the TSGO paper
	for(int y=0; y<ny; ++y) {
		for(int x=0; x<nx; ++x) {
			// get all the minimums.
			lumin = lucost[y][x][0];
			rumin = rucost[y][x][0];
			lgmin = lgcost[y][x][0];
			rgmin = rgcost[y][x][0];
			for(int l=0; l<nLabels; ++l) {
				lumin = (lucost[y][x][l]<lumin)?(lucost[y][x][l]):lumin;
				rumin = (rucost[y][x][l]<rumin)?(rucost[y][x][l]):rumin;
				lgmin = (lgcost[y][x][l]<lgmin)?(lgcost[y][x][l]):lgmin;
				rgmin = (rgcost[y][x][l]<rgmin)?(rgcost[y][x][l]):rgmin;
			}

			// add the energies appropriately.
			lu += lumin;
			ru += rumin;
			lg += lgmin;
			rg += rgmin;
		}
	}
	laplha = 3.5*lu/lg;
	ralpha = 3.5*ru/rg;
	for(int y=0; y<ny; ++y) {
		for(int x=0; x<nx; ++x) {
			for(int l=0; l<nLabels; ++l) {
				// FIXME 
				lcost[y][x][l] = lucost[y][x][l] + laplha*lgcost[y][x][l];
				// cout << "costs: " << lcost[y][x][l] << endl;
				rcost[y][x][l] = rucost[y][x][l] + ralpha*rgcost[y][x][l];
			}
		}
	}
}


void getLowestCostAssignment(
	std::vector<std::vector<std::vector<double> > >& cost,
	std::vector<size_t>& labels
	){
	int ny = cost.size();
	int nx = cost[0].size();
	int nl = cost[0][0].size();

	auto variableIndex = [&nx](const size_t x, const size_t y) {return x + nx*y;};
	auto vx = [&nx](const size_t vid) { return vid%nx; };
	auto vy = [&nx](const size_t vid) { return (int)(vid/nx); };

	for(int y = 0; y<ny; ++y) {
		for(int x=0; x<nx; ++x) {
			for(int l=0; l<nl; ++l) {
				if ((cost[y][x][l] < cost[y][x][labels[variableIndex(x,y)]]))
					labels[variableIndex(x,y)] = l;
			}
		}
	}
}


void addUnaryPotentials(
	// arguments
		std::vector<std::vector<std::vector<double> > >& costs,
		opengm::GraphicalModel<double, opengm::Adder, OPENGM_TYPELIST_3(opengm::ExplicitFunction<double> , opengm::TruncatedAbsoluteDifferenceFunction<double> , opengm::PottsNFunction<double> ) , opengm::SimpleDiscreteSpace<size_t, size_t> >& gm
	) {
	typedef opengm::GraphicalModel<double, opengm::Adder, OPENGM_TYPELIST_3(opengm::ExplicitFunction<double> , opengm::TruncatedAbsoluteDifferenceFunction<double> , opengm::PottsNFunction<double> ) , opengm::SimpleDiscreteSpace<size_t, size_t> > Model;

	size_t ny = costs.size();
	size_t nx = costs[0].size();
	size_t nLabels = costs[0][0].size();

		auto variableIndex = [&nx](const size_t x, const size_t y) {return x + nx*y;};
		auto vx = [&nx](const size_t vid) { return vid%nx; };
		auto vy = [&nx](const size_t vid) { return (int)(vid/nx); };

	const size_t unaryShape[] = {nLabels};	
	for(size_t y=0; y<ny; ++y) {
		for(size_t x=0; x<nx; ++x) {
			opengm::ExplicitFunction<double> f(unaryShape, unaryShape+1);
			for(size_t l=0; l<nLabels; ++l) {
				f(l) = costs[y][x][l];
			}
			Model::FunctionIdentifier fid = gm.addFunction(f);
			size_t variables[] = {variableIndex(x,y)};
			gm.addFactor(fid, variables, variables+1);
		}
	}
}
#endif
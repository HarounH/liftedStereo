#ifndef _PAIRWISE_HPP_
#define _PAIRWISE_HPP_

// Systems includes... I/O stuff and PNG include. 
#include <iostream>
#include <fstream>
#include <png++/png.hpp>
#include <stdlib.h>

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

inline double absolute_color_difference(png::rgb_pixel& p, png::rgb_pixel& q) {
	return abs(p.red - q.red) + abs(p.green - q.green) + abs(p.blue - q.blue);
}

void addPairwisePotentials(
	// arguments
		png::image< png::rgb_pixel >& limg,
		png::image< png::rgb_pixel >& rimg,
		opengm::GraphicalModel<double, opengm::Adder, OPENGM_TYPELIST_3(opengm::ExplicitFunction<double> , opengm::TruncatedAbsoluteDifferenceFunction<double> , opengm::PottsNFunction<double> ) , opengm::SimpleDiscreteSpace<size_t, size_t> >& gm
	) {
	typedef opengm::GraphicalModel<double, opengm::Adder, OPENGM_TYPELIST_3(opengm::ExplicitFunction<double> , opengm::TruncatedAbsoluteDifferenceFunction<double> , opengm::PottsNFunction<double> ) , opengm::SimpleDiscreteSpace<size_t, size_t> > Model;
	size_t nx = limg.get_width();
	size_t ny = limg.get_height();
	size_t nLabels = gm.numberOfLabels(0);
		auto variableIndex = [&nx](const size_t x, const size_t y) {return x + nx*y;};
		auto vx = [&nx](const size_t vid) { return vid%nx; };
		auto vy = [&nx](const size_t vid) { return (int)(vid/nx); };

	double truncationThreshold = 2.0;
	double lambda1=3.5*10; // 2.5*0.125;
	double lambda2=0.6; // 1.25*0.125; 
	double lambda3=0.2; // 1.0*0.125;
	double u1 = 7.0;
	double u2 = 15.0;
	double deltaf = 0.0;

	opengm::TruncatedAbsoluteDifferenceFunction<double> f1(nLabels, nLabels, truncationThreshold, lambda1/truncationThreshold);
	opengm::TruncatedAbsoluteDifferenceFunction<double> f2(nLabels, nLabels, truncationThreshold, lambda2/truncationThreshold);
	opengm::TruncatedAbsoluteDifferenceFunction<double> f3(nLabels, nLabels, truncationThreshold, lambda3/truncationThreshold);
	Model::FunctionIdentifier fid1 = gm.addFunction(f1);
	Model::FunctionIdentifier fid2 = gm.addFunction(f2);
	Model::FunctionIdentifier fid3 = gm.addFunction(f3);

	int c1,c2,c3;
	c1 = c2 = c3 = 0;

	for(size_t y=0; y<ny; ++y) {
		for(size_t x=0; x<nx; ++x) {
			// Two neighbours - x+1,y and x,y+1
			if(x+1 < nx) {
				deltaf = absolute_color_difference(limg[y][x], limg[y][x+1]);
				size_t vars[] = { variableIndex(x,y) , variableIndex(x+1,y) };
				if (deltaf<u1) {
					gm.addFactor(fid1, vars, vars+2);
					c1++;
				} else if (deltaf < u2) {
					gm.addFactor(fid2, vars, vars+2);
					c2++;
				} else {
					gm.addFactor(fid3, vars, vars+2);
					c3++;
				}
			}

			if(y+1 < ny) {
				deltaf = absolute_color_difference(limg[y][x], limg[y+1][x]);
				size_t vars[] = { variableIndex(x,y) , variableIndex(x,y+1) };
				if (deltaf<u1) {
					gm.addFactor(fid1, vars, vars+2);
					c1++;
				} else if (deltaf < u2) {
					gm.addFactor(fid2, vars, vars+2);
					c2++;
				} else {
					gm.addFactor(fid3, vars, vars+2);
					c3++;
				}
			}
		}
	}

	std::cout << "Factor count: \n";
	std::cout << "\tType1: " << c1 << endl;
	std::cout << "\tType2: " << c2 << endl;
	std::cout << "\tType3: " << c3 << endl;
}
#endif
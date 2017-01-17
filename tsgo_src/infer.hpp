#ifndef _INFER_HPP_
#define _INFER_HPP_

// Systems includes... I/O stuff and PNG include. 
#include <iostream>
#include <fstream>
#include <png++/png.hpp>
#include <stdlib.h>
#include <cmath>
#include <ctime>

#include "image.h" // Contains uchar
#include "output.hpp"
// #include "symmetries.hpp"

// OpenGM includes
#include <opengm/opengm.hxx>
#include <opengm/graphicalmodel/graphicalmodel.hxx>
#include <opengm/operations/minimizer.hxx>
#include <opengm/operations/adder.hxx>
#include <opengm/inference/graphcut.hxx>
#include <opengm/inference/alphaexpansion.hxx>
#include <opengm/inference/alphabetaswap.hxx>
#include <opengm/inference/alphaexpansionfusion.hxx>
#include <opengm/inference/trws/trws_trws.hxx>
#include <opengm/inference/auxiliary/minstcutboost.hxx>
#include <opengm/inference/external/trws.hxx>
#include <opengm/graphicalmodel/space/simplediscretespace.hxx>
#include <opengm/functions/potts.hxx>
#include <opengm/functions/pottsn.hxx>
#include <opengm/functions/truncated_absolute_difference.hxx>
#include <opengm/inference/messagepassing/messagepassing.hxx>
#include <opengm/inference/gibbs.hxx>
#include <typeinfo>


void copyDisparities(int* disparities, std::vector<int>& vecDisparities, int nx, int ny) {
	// x + y*nx
	for(int i=0; i<vecDisparities.size(); ++i) {
		disparities[i] = vecDisparities[i];
	}
}

void infer(
	// arguments
		int maxIter, double* mcost, int* disparities, int lr,
		unsigned char * inImage, int ny, int nx, int nLabels,
		bool measureError = false
	) {
	clock_t start, end;
	using namespace std;
	// Step 1 : Reconstruct image.
	png::image< png::rgb_pixel > img(nx, ny);
	start = clock();
	ucharArray2pngImage(img, inImage, nx, ny);
	end = clock();
	std::cout << "uchar2png took:" << double(end-start)/CLOCKS_PER_SEC << std::endl;
	// Step 2 : Create graphical model
	typedef opengm::SimpleDiscreteSpace<int, int> Space;
	Space space(nx * ny, nLabels); // we're doing pixel labelling afterall.
		// functions to aid in indexing.
		auto variableIndex = [&nx](const int x, const int y) {return x + nx*y;};
		auto vx = [&nx](const int vid) { return vid%nx; };
		auto vy = [&nx](const int vid) { return (int)(vid/nx); };

	typedef opengm::GraphicalModel<double, opengm::Adder, OPENGM_TYPELIST_3(opengm::ExplicitFunction<double> , opengm::TruncatedAbsoluteDifferenceFunction<double> , opengm::PottsNFunction<double> ) , Space> Model;
	Model gm(space);

	start = clock();
	// Step 3 : Add potentials to graphical model
	addUnaryPotentials(gm, mcost, nx, ny, nLabels);
	addPairwisePotentials(gm, img, nx, ny, nLabels);
	end = clock();
	cout << "Adding potentials took:" << double(end-start)/CLOCKS_PER_SEC << endl;
	// Step 4 : Create inferencer
	typedef opengm::AlphaExpansionFusion<Model, opengm::Minimizer> AEFInferType;
	AEFInferType::Parameter params(maxIter); // we're gonna step through this.
	AEFInferType::TimingVisitorType visitor;	
	AEFInferType ae_gm(gm, params);

	// Step 5 : Infer
	int nv=nx*ny;
	start = clock();
	std::vector<int> vecDisparities(nx*ny, 0);
	ae_gm.infer(visitor);
	ae_gm.arg(vecDisparities);

	// saveVecDisparities(vecDisparities, nx, ny, 255.0/nLabels);
	
	// Step 6 : Decode inference into disparities.
	copyDisparities(disparities, vecDisparities, nx, ny);
	return;
}

#include "symmetries.hpp"
#include "symmetries2.hpp"
void symmetries_infer( // TODO
	// arguments
		int maxIter, double* mcost, int* disparities, int lr,
		unsigned char * inImage, int ny, int nx, int nLabels,
		int nRanksConsidered=5, int nCBPIterations=1
	) {
	clock_t start, end;	
	using namespace std;
	int nv = nx*ny;
	// Step 1 : Reconstruct image.
	png::image< png::rgb_pixel > img(nx, ny);
	start = clock();
	ucharArray2pngImage(img, inImage, nx, ny);
	end = clock();
	std::cout << "uchar2png took:" << double(end-start)/CLOCKS_PER_SEC << std::endl;
		
	// Step 2 : Create graphical model
	typedef opengm::SimpleDiscreteSpace<int, int> Space;
	Space space(nx * ny, nLabels); // we're doing pixel labelling afterall.
		// functions to aid in indexing.
		auto variableIndex = [&nx](const int x, const int y) {return x + nx*y;};
		auto vx = [&nx](const int vid) { return vid%nx; };
		auto vy = [&nx](const int vid) { return (int)(vid/nx); };

	typedef opengm::GraphicalModel<double, opengm::Adder, OPENGM_TYPELIST_3(opengm::ExplicitFunction<double> , opengm::TruncatedAbsoluteDifferenceFunction<double> , opengm::PottsNFunction<double> ) , Space> Model;
	Model gm(space);

	start = clock();
	// Step 3 : Add potentials to graphical model
	addUnaryPotentials(gm, mcost, nx, ny, nLabels);
	addPairwisePotentials(gm, img, nx, ny, nLabels);
	end = clock();
	cout << "Adding potentials took:" << double(-(start-end))/CLOCKS_PER_SEC << endl;
	
	// Step 3.5 : Reduce graph.
	Model *rgm; // We'll make this point to a reduced GM.
	std::vector<std::vector<int> > syms;
	std::vector<int> invsyms(nx*ny, -1);
	
	processSymmetries(gm, rgm, syms, invsyms, nRanksConsidered, nCBPIterations);
	
	time_t t = time(0);
	struct tm now;
	now = *localtime(&t);
	char fname[80];
	strftime(fname, sizeof(fname),"%X", &now);
	saveSymmetries(nx, ny, syms, strcat(fname,".png"));
	
	// Step 4 : Create inferencer
	typedef opengm::AlphaExpansionFusion<Model, opengm::Minimizer> AEFInferType;
	AEFInferType::Parameter params(maxIter); // we're gonna step through this.
	AEFInferType::TimingVisitorType visitor;	
	AEFInferType ae_rgm(*rgm, params);

	std::cout << "rgm:\n\tnVar=" << rgm->numberOfVariables() << "\n\tnFactors" << rgm->numberOfFactors() << std::endl;

	// Step 5 : Infer
	std::vector<int> reducedDisparities(rgm->numberOfVariables(), 0);
	ae_rgm.infer(visitor);
	std::vector<int> vecDisparities(nx*ny, 0);
	ae_rgm.arg(reducedDisparities);
	for(int i=0; i<vecDisparities.size(); ++i) {
		vecDisparities[i] = reducedDisparities[ invsyms[i] ];
	}
	// Step 6 : Decode inference into disparities.
	copyDisparities(disparities, vecDisparities, nx, ny);

	if(rgm) delete rgm;
	
	return;
}


#endif
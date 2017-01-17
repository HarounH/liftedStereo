// This file runs alphaexpansionfusion using QPBO 
// Reports the change in accuracy as iterations change.

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

// Symmetries include
#include "useSym.hpp" // contains reduceGraphUsingSymmetries and also necessary inference stuff.
#include "mapSym.hpp" // contains compress factor graph... which is what we want.
#include "evaluator.hpp"
#include "outputSymmetries.hpp"
#include "unaryPotentials.hpp"
#include "outputgroups.hpp"

using namespace std;

/*
	The main function does everything. It's a very hacky code, sorry.
*/
int main(int argc, char** argv) {
	// STEP 1 : Get parameters

	// params - files
	string leftImageFilename = argv[1];
	string rightImageFilename= argv[2];
	string noSymsOutputFile = argv[3];
	string symsOutputFile = argv[4];
	string symsObtained = symsOutputFile + "_syms.png";
	string goldImageFilename = argv[5];

	// params - problem parameters, 
	int nLabels = atoi(argv[6]);
	double scalingFactor = 255.0/nLabels; // least appropriate variable name, I think.
	double goldScale = atof(argv[7]);
	double badPixelThreshold = atof(argv[8]);

	// params - energy parameters
	int unaryPotentialRadius = atoi(argv[9]);
	double truncationThreshold = atof(argv[10]);
	double truncationPenalty = atof(argv[11]);

	// params - inference parameters
	int nIterations = atoi(argv[12]);

	// params - symmetry parameters
	int nRanksUsed = atoi(argv[13]);
	int nIterationsOfCBP = atoi(argv[14]);

	// variables
	clock_t begin, end; // use these for timing.

	// STEP 2 : Print all the parameters
	cout << "Parameters are:\n";
	cout << "lImage: " << leftImageFilename << "\nrImage: " << rightImageFilename << "\n";
	cout << "gold: " << goldImageFilename << "\n";
	cout << "predicting using " << nLabels << " labels, scale=" << scalingFactor << "\n"; 
	cout << "\tgoldScale=" << goldScale << "\n";
	cout << "badPixelThreshold=" << badPixelThreshold << "\n";
	cout << "Potentials using:\n";
	cout << "\tunaryPotentialRadius: " << unaryPotentialRadius << "\n";
	cout << "\ttruncationThreshold: " << truncationThreshold << "\n";
	cout << "\ttruncationPenalty: " << truncationPenalty << "\n";
	cout << "Alphaexpansionfusion using " << nIterations << "iterations\n";
	cout << "Symmetries using:\n";
	cout << "nRanksUsed= " << nRanksUsed << "\n";
	cout << "nIterationsOfCBP= " << nIterationsOfCBP << "\n";

	// STEP 3 : load images
	png::image< png::rgb_pixel > left(leftImageFilename);
	png::image< png::rgb_pixel > right(rightImageFilename);
	png::image<png::ga_pixel> gold(goldImageFilename);
	size_t nx = left.get_width();
	size_t ny = left.get_height();
	
	
	// STEP 4 : create nosymmetries graphical model 'gm'
	typedef opengm::SimpleDiscreteSpace<size_t, size_t> Space;
	Space space(nx * ny, nLabels); // we're doing pixel labelling afterall.
		// (x,y) -> flat array index
		auto variableIndex = [&nx](const size_t x, const size_t y) {return x + nx*y;};
		auto vx = [&nx](const size_t vid) { return vid%nx; };
		auto vy = [&nx](const size_t vid) { return (int)(vid/nx); };
	typedef opengm::GraphicalModel<double, opengm::Adder, OPENGM_TYPELIST_3(opengm::ExplicitFunction<double> , opengm::TruncatedAbsoluteDifferenceFunction<double> , opengm::PottsNFunction<double> ) , Space> Model;
	Model gm(space); // No potential functions yet.
	
	// unary potentials first
	const int unaryShape[] = {nLabels};
	for(size_t x=0; x<nx; ++x) {
		for(size_t y=0; y<ny; ++y) {
			opengm::ExplicitFunction<double> f(unaryShape, unaryShape+1);
			for(size_t d=0; d<nLabels; ++d) { // d is the disparity that we assign to it.
				// TODO: Refine this metric.
				if (x+d<nx) {
					f(d) = unaryPotentialBirchfieldNoWindow(left, right, x, y, d);
					// cout << "f(d)=" << f(d) << endl;
					// f(d) = diff_square(left[y][x], right[y][x+1+d]); // IMPORTANT: LABEL=0 -> 1 pixel away.
				} else {
					f(d) = nLabels*nLabels*255*255; // FIXME : What do I do with this edge case? :(
					break; // break out of the loop on 'd'.
				}
			}
			Model::FunctionIdentifier fid = gm.addFunction(f);
			size_t variables[] = {variableIndex(x,y)};
			gm.addFactor(fid, variables, variables+1);
		}
	}

	// pairwise potentials then
	opengm::TruncatedAbsoluteDifferenceFunction<double> f(nLabels, nLabels, truncationThreshold, truncationPenalty);
	Model::FunctionIdentifier fid = gm.addFunction(f);
	for(size_t x=0; x<nx; ++x) {
		for(size_t y=0; y<ny; ++y) {
			if (x+1 < nx) { // add y,x+1 as neighbourdouble lambda = LAMBDA; 
				size_t variables[] = {variableIndex(x,y), variableIndex(x+1, y)};
				sort(variables, variables+2);
				gm.addFactor(fid, variables, variables+2);
			}
			if (y+1 < ny) { // add y+1,x as neighbour
				size_t variables[] = {variableIndex(x,y), variableIndex(x, y+1)};
				sort(variables, variables+2);
				gm.addFactor(fid, variables, variables+2);
			}
		}
	}
	
	// no higher order potentials.


	// STEP 5 : reduce graphical model to create 'rgm'
	begin = clock_t();
	std::vector< std::vector<size_t> > syms = getVarGroupingsUsingMapSyms
												<double, opengm::Adder, OPENGM_TYPELIST_3(opengm::ExplicitFunction<double> , opengm::TruncatedAbsoluteDifferenceFunction<double> , opengm::PottsNFunction<double> ) , Space>
													(gm, nIterationsOfCBP, nRanksUsed);
	std::vector<size_t> invsyms(gm.numberOfVariables());
	getVar2GroupID(syms, invsyms);
	end = clock_t();
	cout << "Took " << double(end-begin)/CLOCKS_PER_SEC << " s to find symmetries" << endl;	
	
	cout << "Saving symmetries to " << symsObtained << "\n";
	saveSymmetries(nx, ny, syms, symsObtained);

	begin = clock();
	Model rgm = reduceGraphUsingSymmetriesEfficienctly<double, opengm::Adder, OPENGM_TYPELIST_3(opengm::ExplicitFunction<double> , opengm::TruncatedAbsoluteDifferenceFunction<double> , opengm::PottsNFunction<double> ) , Space>(
		gm, syms, invsyms, truncationThreshold, truncationPenalty, 0.0);
	end = clock();
	
	cout << "Reducing graph using symmetries took " << double(end-begin)/CLOCKS_PER_SEC << " s" << endl;
	
	cout << "From gm with:" << endl;
	cout << "\tnVar/Grp = " << gm.numberOfVariables() << endl;
	cout << "\tnFactors = " << gm.numberOfFactors() << endl;
	cout << "Using rgm with:" << endl;
	cout << "\tnVar/Grp = " << rgm.numberOfVariables() << endl;
	cout << "\tnFactors = " << rgm.numberOfFactors() << endl;


	// STEP 6 : create inferencer.
	typedef opengm::MinSTCutBoost<size_t, double, opengm::KOLMOGOROV> MinStCutType;
	typedef opengm::GraphCut<Model, opengm::Minimizer, MinStCutType> MinGraphCut;
	typedef opengm::AlphaExpansion<Model, MinGraphCut> AlphaExpansionTypeDef;
	typedef opengm::AlphaExpansionFusion<Model, opengm::Minimizer> AEFInferType;
	AEFInferType::Parameter params(1); // we're gonna step through this.

	


	
	// STEP 7 : loop for iterations.
		// no symmetries
		AEFInferType ae_gm(gm, params);
		png::image<png::ga_pixel> noSymsOutput(nx,ny);
		std::vector<size_t> gm_labels(gm.numberOfVariables());
		for(int iter=0; iter<nIterations; ++iter) {
			AEFInferType::TimingVisitorType visitor;
			if(iter!=0) { // we have output from previous run :)
				ae_gm.setStartingPoint(gm_labels.begin());
			}
			ae_gm.infer(visitor); // that will print some stuff out.
			ae_gm.arg(gm_labels);
			// copy labels into noSymsOutput.
			for(size_t vid=0; vid<nx*ny; ++vid) {
				noSymsOutput[vy(vid)][vx(vid)] = png::ga_pixel( (int)((gm_labels[vid])*scalingFactor) );
			}
			cout << "NoSyms iteration " << iter << " error=" << calculatePixelError(noSymsOutput, gold, nLabels, goldScale, badPixelThreshold);
		}

		// with symmetries.
		AEFInferType ae_rgm(rgm, params);
		png::image<png::ga_pixel> symsOutput(nx,ny);
		std::vector<size_t> rgm_labels(rgm.numberOfVariables());
		for(int iter=0; iter<nIterations; ++iter) {
			AEFInferType::TimingVisitorType visitor;
			if(iter!=0) { // we have output from previous run :)
				ae_rgm.setStartingPoint(rgm_labels.begin());
			}
			ae_rgm.infer(visitor); // that will print some stuff out.
			ae_rgm.arg(rgm_labels);
			// copy labels into noSymsOutput.
			for(size_t vid=0; vid<nx*ny; ++vid) {
				symsOutput[vy(vid)][vx(vid)] = png::ga_pixel( (int)((rgm_labels[invsyms[vid]])*scalingFactor) );
			}
			cout << "syms iteration " << iter << " error=" << calculatePixelError(symsOutput, gold, nLabels, goldScale, badPixelThreshold);
		}


	// STEP 8 : save outputs
	// nosymmetries image
	cout << "saving no symmetries output to " << noSymsOutputFile << "\n";
	noSymsOutput.write(noSymsOutputFile);
	// symmetries used image
	cout << "saving symmetries output to " << symsOutputFile << "\n";
	symsOutput.write(symsOutputFile);
	
	cout << "completed." << endl;
	return 0;
}
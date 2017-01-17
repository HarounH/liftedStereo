/**
	@author Haroun
	IMPORTANT: symmetries are used in stereoInfer2.cpp
	This file takes in two PNG files (left, right), and outputs a stereo matching

	Requirements:
		Needs OpenGM 2
		Needs png++

	Compile:
		g++ <filename> -std=c++11 -lpng
			Need c++11 because I'm using 'auto'... being the lazy ass that I am. 
	Execute:
		singleInstance: <./a.out> single <path-to-left> <path-to-right> <output-file=output.png> <nlabels=10 lambda=1000 unaryPotentialRadius=3>
			all files are assumed to be PNGs
		multipleInstances: <./a.out> multi <path to file contaning list of instances> <nlabels=10 lambda=1000 unaryPotentialRadius=3>
			file format is as follows:
				N (number of instances)
				<path-to-left> <path-to-right> <output-file=output.png>
				<path-to-left> <path-to-right> <output-file=output.png>
				...
				<path-to-left> <path-to-right> <output-file=output.png>

	Evaluation of results:
		A seperate script is run to evaluate accuracies.
*/

// Systems includes... I/O stuff and PNG include. 
#include <iostream>
#include <fstream>
#include <png++/png.hpp>
#include <stdlib.h>

// OpenGM includes
#include <iostream>
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
#include "binVariables.hpp" // contains compress factor graph... which is what we want.

#include "evaluator.hpp"

using namespace std;

#define N_LABELS 15
#define LAMBDA 100
#define UNARY_POTENTIAL_RADIUS 7
#define PN_CLIQUE_DIAM 2
#define N_STEPS 100

size_t nLabels = N_LABELS; // Modify this.
double lambda = LAMBDA; // Modify this.
size_t unaryPotentialRadius = UNARY_POTENTIAL_RADIUS; // modify this.
size_t nSteps = N_STEPS;
double pnLambda = LAMBDA;

int nBinsPerLabel = 16;
int nMaxIterations = 4;

// quick and dirty function to get difference between two pixels.
inline double diff_square(png::rgb_pixel const& l, png::rgb_pixel const& r) {
	return (l.red - r.red)*(l.red - r.red) + (l.green - r.green)*(l.green - r.green) + (l.blue - r.blue)*(l.blue - r.blue);
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

double unaryPotentialMeanFilter(png::image< png::rgb_pixel >& left, png::image< png::rgb_pixel >& right, size_t x, size_t y, size_t d) {
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

void singleInstance(string leftFile, string rightFile, string outputFile, string groundTruth) {
	// read the left and right files using PNG++
	// cout << "running with " << endl << "nlabels=" << nLabels << " and lambda=" << lambda << " and upr=" << unaryPotentialRadius << " and nSteps=" << nSteps <<endl;
	
	png::image< png::rgb_pixel > left(leftFile);
	png::image< png::rgb_pixel > right(rightFile);
	
	size_t nx = left.get_width();
	size_t ny = left.get_height();
	cout << "Image size=" << nx << "," << ny << endl;
	// contruct a graphical model... use the energies appropriately.
	typedef opengm::SimpleDiscreteSpace<size_t, size_t> Space;
	Space space(nx * ny, nLabels); // we're doing pixel labelling afterall.
		// (x,y) -> flat array index
		auto variableIndex = [&nx](const size_t x, const size_t y) {return x + nx*y;};
		auto vx = [&nx](const size_t vid) { return vid%nx; };
		auto vy = [&nx](const size_t vid) { return (int)(vid/nx); };
	// TODO: Worry about functions in the model.
	typedef opengm::GraphicalModel<double, opengm::Adder, OPENGM_TYPELIST_3(opengm::ExplicitFunction<double> , opengm::TruncatedAbsoluteDifferenceFunction<double> , opengm::PottsNFunction<double> ) , Space> Model;
	Model gm(space);
	
	// start with a simple Potts based segmentation model.

	// Unary energy
	const size_t unaryShape[] = {nLabels};
	for(size_t x=0; x<nx; ++x) {
		for(size_t y=0; y<ny; ++y) {
			opengm::ExplicitFunction<double> f(unaryShape, unaryShape+1);
			for(size_t d=0; d<nLabels; ++d) { // d is the disparity that we assign to it.
				// TODO: Refine this metric.
				if (x+d<nx) {
					f(d) = unaryPotentialMeanAbsoluteDifference(left, right, x, y, d);
					// cout << "f(d)=" << f(d) << endl;
					// f(d) = diff_square(left[y][x], right[y][x+1+d]); // IMPORTANT: LABEL=0 -> 1 pixel away.
				} else {
					f(d) = unaryPotentialRadius*lambda*nLabels*d; // FIXME : What do I do with this edge case? :(
					break; // break out of the loop on 'd'.
				}
			}
			Model::FunctionIdentifier fid = gm.addFunction(f);
			size_t variables[] = {variableIndex(x,y)};
			gm.addFactor(fid, variables, variables+1);
		}
	}

	// Pott's energy.
	// opengm::PottsFunction<double> f(nLabels, nLabels, 0.0, lambda);
	// Model::FunctionIdentifier fid = gm.addFunction(f);
	
	// truncated absolute difference... because why not.
	opengm::TruncatedAbsoluteDifferenceFunction<double> f(nLabels, nLabels, nLabels/2, lambda);
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


	// higher order potential functions
	size_t pnShape[PN_CLIQUE_DIAM*PN_CLIQUE_DIAM];
	for(size_t i=0; i<PN_CLIQUE_DIAM*PN_CLIQUE_DIAM; ++i) {
		pnShape[i] = nLabels;
	}
	opengm::PottsNFunction<double> pn(pnShape, pnShape+PN_CLIQUE_DIAM, 0.0, pnLambda); // Need to decide shape etc.
	fid = gm.addFunction(pn);
	bool temp = true;
	for(size_t x=0; x<nx; ++x) {
		for(size_t y=0; y<ny; ++y) {
			// add the function on stuff.
			if(((x+PN_CLIQUE_DIAM)<nx) && ((y+PN_CLIQUE_DIAM)<ny)) {
				std::vector<size_t> vars;
				for(size_t xiter=x; xiter<(x+PN_CLIQUE_DIAM); ++xiter) {
					for(size_t yiter=y; yiter<(y+PN_CLIQUE_DIAM); ++yiter) {
						vars.push_back(variableIndex(xiter, yiter));
					}
				}
				if (temp) {
					// cout << "pottsn_size=" << vars.size();
					temp = false;
				}
				sort(vars.begin(), vars.end());
				gm.addFactor(fid, vars.begin(), vars.end());
			}
		}
	}
	// cout << "gm.nOF=" << gm.numberOfFactors() << endl;

	// reduce the graphical model.
	clock_t begin = clock();
	
	// need to create syms, invsyms
	std::vector< std::vector<size_t> > syms = getVarGroupingsUsingBins
												<double, opengm::Adder, OPENGM_TYPELIST_3(opengm::ExplicitFunction<double> , opengm::TruncatedAbsoluteDifferenceFunction<double> , opengm::PottsNFunction<double> ) , Space>
													(gm, nMaxIterations, nBinsPerLabel);

	std::vector<size_t> invsyms(gm.numberOfVariables());
	getVar2GroupID(syms, invsyms);

	

	clock_t end = clock();
	cout << "Symmetry finding took " << double(end-begin)/CLOCKS_PER_SEC << "s" << endl;	

	int nGroups = syms.size();
	std::cout << "nGroups = " << nGroups << std::endl;
	
	begin = clock();
	
	Model rgm = reduceGraphUsingSymmetries<double, opengm::Adder, OPENGM_TYPELIST_3(opengm::ExplicitFunction<double> , opengm::TruncatedAbsoluteDifferenceFunction<double> , opengm::PottsNFunction<double> ) , Space>(
		gm, syms, invsyms);
	
	end = clock();
	cout << "Reducing graph using symmetries took " << double(end-begin)/CLOCKS_PER_SEC << " s" << endl;

	cout << "From gm with:" << endl;
	cout << "\tnVar/Grp = " << gm.numberOfVariables() << endl;
	cout << "\tnFactors = " << gm.numberOfFactors() << endl;
	cout << "Using rgm with:" << endl;
	cout << "\tnVar/Grp = " << rgm.numberOfVariables() << endl;
	cout << "\tnFactors = " << rgm.numberOfFactors() << endl;

	
	
	// call the inference function.
	// typedef opengm::MinSTCutKolmogorov<size_t, double> MinStCutType;
	typedef opengm::MinSTCutBoost<size_t, double, opengm::KOLMOGOROV> MinStCutType;
	typedef opengm::GraphCut<Model, opengm::Minimizer, MinStCutType> MinGraphCut;
	typedef opengm::AlphaExpansion<Model, MinGraphCut> AlphaExpansionTypeDef;
	typedef opengm::AlphaExpansionFusion<Model, opengm::Minimizer> AEFInferType;
	// typedef opengm::AlphaExpansionFusion<Model, opengm::Minimizer> AlphaExpansionTypeDef; 
	AEFInferType::Parameter params(nSteps);

	AEFInferType ae(rgm, params);

	AEFInferType::TimingVisitorType visitor;
	ae.infer(visitor);

	vector<size_t> labels(rgm.numberOfVariables());
	ae.arg(labels);
	cout << "inference complete." << endl;
	// convert the output to png++ format.
	cout << "saving MAP estimate" << endl;
	
	png::image<png::ga_pixel> output(nx,ny);
	double scalingFactor = 255.0/nLabels; // scaling 0...nLabels to 0...255
	for(size_t vid=0; vid<nx*ny; ++vid) {
		output[vy(vid)][vx(vid)] = png::ga_pixel( (int)((labels[invsyms[vid]])*scalingFactor) );
	}
	// save the output
	cout << "writing nx*ny=" << nx*ny << " pixels to" << outputFile << endl;
	output.write(outputFile);

	// evaluate bad pixels.
	png::image<png::ga_pixel> truth(groundTruth);
	cout << "bad pixel score=" << (double)badPixelScore(output, truth, nLabels, 1.0) << endl;
}

void multipleInstances(string instancefile) {
	// Repeatedly call singleInstance, really
	std::vector< std::vector<std::string> > instances;
	std::ifstream file(instancefile);
	std::string line;

	while(std::getline(file, line)) {
		std::stringstream sstream(line);
		std::string left,right,output, truth;

		std::getline(sstream, left, ' ');
		std::getline(sstream, right, ' ');
		std::getline(sstream, output, ' ');
		std::getline(sstream, truth, ' ');
		instances.push_back({left, right, output, truth});
	}

	for(int id=0; id< instances.size(); ++id) {
		singleInstance(instances[id][0],instances[id][1],instances[id][2],instances[id][3]);
	}
}

int main(int argc, char** argv) {
	if (argc>1) { // need to do better, but ah well.
		if (strcmp("single",argv[1]) == 0) { // singleInstance
			std::string golden;
			if (argc>5) {
				// assume that N_LABELS and LAMBDA have been given in that order
				nLabels = atoi(argv[5]);
				lambda = atof(argv[6]);
				pnLambda = lambda;
				unaryPotentialRadius=atoi(argv[7]);
				nSteps = atoi(argv[8]);
				golden = argv[9];
				
				nBinsPerLabel = atoi(argv[10]);
				nMaxIterations= atoi(argv[11]);
			}
			cout << "running as follows:" << endl;
			cout << "\tnLabels=" << nLabels << "\n\tlambda=" << lambda << "\n\tunaryPotentialRadius=" << unaryPotentialRadius << "\n\tnSteps=" << nSteps;
			cout << "\n\tnBinsPerLabel=" << nBinsPerLabel << "\n\tnMaxIterations of CBP=" << nMaxIterations << endl;
			singleInstance(argv[2], argv[3], argv[4], golden);
		} else if (strcmp("multi",argv[1]) == 0) { // multipleInstances
			if (argc>3) {
				// assume that N_LABELS and LAMBDA have been given in that order
				nLabels = atoi(argv[3]);
				lambda = atof(argv[4]);
				unaryPotentialRadius=atoi(argv[5]);			
				nSteps = atoi(argv[6]);
			}
			multipleInstances(argv[2]);
		} else if (strcmp("eval", argv[1])) {
			
		}
	} else {
		cout << "Please read usage instructions!" << endl;
		cout << "\nUSAGE" << endl;
		cout << "singleInstance: <./a.out> single <path-to-left> <path-to-right> <output-file=output.png> <nlabels=10 lambda=1000 unaryPotentialRadius=3> <truth-file>" << endl;
		cout << "multipleInstances: <./a.out> multi <path to file contaning list of instances> <nlabels=10 lambda=1000 unaryPotentialRadius=3>" << endl << endl;
		return 1;
	}
	cout << "Exiting now" << argv[1] << endl;

	return 0;
}
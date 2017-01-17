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

// no symmetries here :)
#include "evaluator.hpp"
#include "unaryPotentials.hpp"

using namespace std;

#define LAMBDA 100
#define PN_CLIQUE_DIAM 2
#define N_STEPS 100

double lambda = LAMBDA; // Modify this.
size_t nSteps = N_STEPS;
double pnLambda = LAMBDA;
double truncationThreshold = nLabels/2;
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
				// cout << x << " " << y << " " << d << endl;
				if (x+d<nx) {
					f(d) = unaryPotentialBirchwoodTomiskiWindowed(left, right, x, y, d);
				} else {
					if ((d-1)>0) {// FIXME : What do I do with this edge case? :(
						f(d) = f(d-1);
					} else {
						f(d) = 0.5;
					}
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
	opengm::TruncatedAbsoluteDifferenceFunction<double> f(nLabels, nLabels, truncationThreshold, lambda);
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
				sort(vars.begin(), vars.end());
				gm.addFactor(fid, vars.begin(), vars.end());
			}
		}
	}
	cout << "gm.nOF=" << gm.numberOfFactors() << endl;
	// call the inference function.
	// typedef opengm::MinSTCutKolmogorov<size_t, double> MinStCutType;
	typedef opengm::MinSTCutBoost<size_t, double, opengm::KOLMOGOROV> MinStCutType;
	typedef opengm::GraphCut<Model, opengm::Minimizer, MinStCutType> MinGraphCut;
	typedef opengm::AlphaExpansion<Model, MinGraphCut> AlphaExpansionTypeDef;
	typedef opengm::AlphaExpansionFusion<Model, opengm::Minimizer> AEFInferType;
	// typedef opengm::AlphaExpansionFusion<Model, opengm::Minimizer> AlphaExpansionTypeDef; 
	AEFInferType::Parameter params(nSteps);

	AEFInferType ae(gm, params);

	AEFInferType::TimingVisitorType visitor;
	ae.infer(visitor);

	vector<size_t> labels(nx*ny);
	ae.arg(labels);
	cout << "inference complete." << endl;
	// convert the output to png++ format.
	cout << "saving MAP estimate" << endl;
	
	png::image<png::ga_pixel> output(nx,ny);
	double scalingFactor = 255.0/nLabels; // scaling 0...nLabels to 0...255
	for(size_t vid=0; vid<nx*ny; ++vid) {
		output[vy(vid)][vx(vid)] = png::ga_pixel( (int)((labels[vid])*scalingFactor) );
	}	
	// save the output
	output.write(outputFile);

	// evaluate bad pixels.
	// png::image<png::ga_pixel> truth(groundTruth);
	// cout << "bad pixel score=" << (double)badPixelScore(output, truth, nLabels, 1.0) << endl;
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
				truncationThreshold=atof(argv[6]);
				lambda = atof(argv[7]);
				pnLambda = atof(argv[8]);
				unaryPotentialRadius=atoi(argv[9]);
				nSteps = atoi(argv[10]);
				golden = argv[11];
			}
			cout << "running as follows:" << endl;
			cout << "\tnLabels=" << nLabels << "\n\ttruncationThreshold" << truncationThreshold <<"\n\tlambda=" << lambda << "\n\tunaryPotentialRadius=" << unaryPotentialRadius << "\n\tnSteps=" << nSteps << endl;
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
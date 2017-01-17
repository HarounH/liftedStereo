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

// user defined includes...
#include "preprocessImages.hpp"

#include "unary.hpp"

#include "imageStatistics.hpp"

#include "deprecated_upperLayer.hpp"

#include "upperLayer.hpp"

#include "pairwise.hpp"

#include "gaussian.hpp"
// #include "getSymmetries.hpp"
// #include "useSymmetries.hpp"
#include "output.hpp"

#include "postProcess.hpp"

#include "evaluators.hpp"

using namespace std;
clock_t start, stop;
int main(int argc, char const *argv[]) {
	// Step 1 : Parameters

	// params - file parameters
	string leftImageFilename = argv[1];
	string rightImageFilename= argv[2];
	string noSymsOutputFile = argv[3];
	string symsOutputFile = argv[4];
	string symsObtained = symsOutputFile + "_syms.png";
	string goldImageFilename = argv[5];

	// params - problem parameters
	int nLabels = atoi(argv[6]);
	double scalingFactor = 255.0/nLabels; // least appropriate variable name, I think.
	double goldScale = atof(argv[7]);
	double badPixelThreshold = atof(argv[8]);

	// params - energy
	double truncationThreshold = atof(argv[9]);
	double truncationPenalty = atof(argv[10]);

	// params - inference parameters
	int nIterations = atoi(argv[11]);

	// params - symmetry parameters
	int nRanksUsed = atoi(argv[12]);
	int nIterationsOfCBP = atoi(argv[13]);


	// Step 2 : Load images
	png::image< png::rgb_pixel > limg(leftImageFilename);
	png::image< png::rgb_pixel > rimg(rightImageFilename);
	int nx = limg.get_width();
	int ny = limg.get_height();
	normalize(limg, rimg);
	cout << "loaded images " << endl;

	// Step 3 : Upper layer of TSGO
	std::vector<std::vector<std::vector<double> > > lcost(ny, std::vector<std::vector<double> >(nx, std::vector<double>(nLabels,0.0)));
	std::vector<std::vector<std::vector<double> > > rcost(ny, std::vector<std::vector<double> >(nx, std::vector<double>(nLabels,0.0)));
	start = clock();
	getTSGODissimilarityMeasures(limg, rimg, lcost, rcost, nLabels);
	doERFTransform(lcost, false);
	doERFTransform(rcost, false);
	int gauss_radius = 4;
	stop = clock();
	cout << "got original dissimalirity measures in " << (double(stop - start)/CLOCKS_PER_SEC) << " s" << endl;
	// the sqsigmas are set by the paper...so...
	start = clock();
	bilateral_filter(limg, lcost, 14.0, 1.55);
	bilateral_filter(rimg, rcost, 14.0, 1.55);
	stop = clock();
	cout << "bilateral filter took " << (double(stop - start)/CLOCKS_PER_SEC) << " s" << endl;

	// Step 4 : create graphical model of lower layer
	typedef opengm::SimpleDiscreteSpace<size_t, size_t> Space;
	Space space(nx * ny, nLabels); // we're doing pixel labelling afterall.
		// (x,y) -> flat array index
		auto variableIndex = [&nx](const size_t x, const size_t y) {return x + nx*y;};
		auto vx = [&nx](const size_t vid) { return vid%nx; };
		auto vy = [&nx](const size_t vid) { return (int)(vid/nx); };
	typedef opengm::GraphicalModel<double, opengm::Adder, OPENGM_TYPELIST_3(opengm::ExplicitFunction<double> , opengm::TruncatedAbsoluteDifferenceFunction<double> , opengm::PottsNFunction<double> ) , Space> Model;
	Model lgm(space); // No potential functions yet.
	Model rgm(space); // No potential functions yet.

	addUnaryPotentials(lcost, lgm);
	addUnaryPotentials(rcost, rgm); // Defined in unary.hpp
	cout << "added unary potentials " << endl;
	addPairwisePotentials(limg, rimg, lgm);
	addPairwisePotentials(rimg, limg, rgm); // Defined in pairwise.hpp
	cout << "added pairwise potentials " << endl;
	// Step 5 : Reduce or don't reduce graphical model.

	// Step 6 : Infer
	typedef opengm::AlphaExpansionFusion<Model, opengm::Minimizer> AEFInferType;
	AEFInferType::Parameter params(nIterations); // we're gonna step through this.
	AEFInferType::TimingVisitorType lvisitor;
	AEFInferType::TimingVisitorType rvisitor;	

	std::vector<size_t> llabels(nx*ny,0);
	AEFInferType ae_lgm(lgm, params);
	cout << "Starting left image inference " << endl;
	ae_lgm.infer(lvisitor);	
	ae_lgm.arg(llabels);


	std::vector<size_t> rlabels(nx*ny,0);
	AEFInferType ae_rgm(rgm, params);
	cout << "Starting right image inference " << endl;
	ae_rgm.infer(rvisitor);	
	ae_rgm.arg(rlabels);
	
	// Step 7 : Post processing.
	start = clock();
	double sqsigmafl = getSquaredSigmaOfLabels(limg, llabels);
	double sqsigmafr = getSquaredSigmaOfLabels(rimg, rlabels);
	cout << "sqsigmafl: " << sqsigmafl << endl;
	cout << "sqsigmafr: " << sqsigmafr << endl;
	
	cout << "Attempting post processing. " << endl;
	std::vector<bool> luncertain(llabels.size(), false);
	std::vector<bool> runcertain(rlabels.size(), false);

	left_to_right_cross_checking(llabels, rlabels, luncertain, nx, ny, 1);
	// gaussian_median_filter(limg, llabels, luncertain, nLabels, gauss_radius, (gauss_radius*0.25*gauss_radius), /*255*9**/sqsigmafl/nLabels/nLabels);
	// gaussian_median_filter(rimg, rlabels, runcertain, nLabels, gauss_radius, (gauss_radius*0.25*gauss_radius), /*255*9**/sqsigmafr/nLabels/nLabels);
	// left_to_right_cross_checking(llabels, rlabels, luncertain, nx, ny, 1);
	// gaussian_median_filter(limg, llabels, luncertain, nLabels, gauss_radius, (gauss_radius*0.25*gauss_radius), /*255*9**/sqsigmafl/nLabels/nLabels);
	// left_to_right_cross_checking(rlabels, llabels, runcertain, nx, ny, -1);
	// left_to_right_cross_checking(llabels, rlabels, luncertain, nx, ny, 1);
	stop = clock();
	cout << "post processing took " << (double(double(stop - start)/CLOCKS_PER_SEC)) << " s" << endl;
	// Step 8 : Evaluate & Save output
	cout << "nx: " << nx << " ny: " << ny << endl;
	png::image< png::ga_pixel > gold(goldImageFilename);
	cout << "left image error" <<  badPixelThreshold << ": " << badPixelScore(llabels, gold, goldScale, badPixelThreshold) << endl;
	labels2grayimage(llabels, nx, ny, scalingFactor, noSymsOutputFile+"_l.png");
	labels2grayimage(rlabels, nx, ny, scalingFactor, noSymsOutputFile+"_r.png");
	return 0;
}
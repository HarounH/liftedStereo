#ifndef _MODEL_CREATION_HPP_
#define _MODEL_CREATION_HPP_
// Haroun: Systems includes... I/O stuff and PNG include. 
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

#define LAMBDA1 2.5*0.125
#define LAMBDA2 1.25*0.125
#define LAMBDA3 1*0.125
const double u1 = 7.0;
const double u2 = 15.0;

typedef opengm::SimpleDiscreteSpace<int, int> Space;
typedef opengm::GraphicalModel<double, opengm::Adder, OPENGM_TYPELIST_3(opengm::ExplicitFunction<double> , opengm::TruncatedAbsoluteDifferenceFunction<double> , opengm::PottsNFunction<double> ) , Space> Model;

void addUnaryPotentials(
	// arguments
		opengm::GraphicalModel<double, opengm::Adder, OPENGM_TYPELIST_3(opengm::ExplicitFunction<double> , opengm::TruncatedAbsoluteDifferenceFunction<double> , opengm::PottsNFunction<double> ) , opengm::SimpleDiscreteSpace< int, int > >& gm,
		double * mcost,
		int nx, int ny, int nLabels
	) { // function begins

	const int unaryShape[] = {nLabels};	
	int nPixels = nx*ny;
	for(int y=0; y<ny; ++y) {
		for(int x=0; x<nx; ++x) {
			opengm::ExplicitFunction<double> f(unaryShape, unaryShape+1);
			for(int l=0; l<nLabels; ++l) {
				f(l) = mcost[y*nx + x  + l*nPixels];
				if (f(l) <0) {
					std::cout << "Super weird. negative unary\n";
				}
			}
			Model::FunctionIdentifier fid = gm.addFunction(f);
			int variables[] = {x + y*nx};
			gm.addFactor(fid, variables, variables+1);
		}
	}
}


void addUnaryPotentials(
	// arguments
		opengm::GraphicalModel<double, opengm::Adder, OPENGM_TYPELIST_3(opengm::ExplicitFunction<double> , opengm::TruncatedAbsoluteDifferenceFunction<double> , opengm::PottsNFunction<double> ) , opengm::SimpleDiscreteSpace< int, int > >& gm,
		double * mcost,
		int nx, int ny, int nLabels,
		std::vector< std::vector<int> >& syms, std::vector<int>& invsyms
	) { // function begins

	const int unaryShape[] = {nLabels};	
	int nPixels = nx*ny;
	int x,y;
	for(int v=0; v<syms.size(); ++v) {
		opengm::ExplicitFunction<double> f(unaryShape, unaryShape+1);
		for(int i=0; i<syms[v].size(); ++i) {
			y = syms[v][i]/nx;
			x = syms[v][i]%nx;
			for(int l=0; l<nLabels; ++l) {
				f(l) += mcost[y*nx + x  + l*nPixels];
			}
		}
		Model::FunctionIdentifier fid = gm.addFunction(f);
		int variables[] = {v};
		gm.addFactor(fid, variables, variables+1);
	}
}

inline double absolute_color_difference(png::rgb_pixel& p, png::rgb_pixel& q) {
	return 0.0 + abs(p.red - q.red) + abs(p.green - q.green) + abs(p.blue - q.blue);
}


void addPairwisePotentials(
	// arguments
		opengm::GraphicalModel<double, opengm::Adder, OPENGM_TYPELIST_3(opengm::ExplicitFunction<double> , opengm::TruncatedAbsoluteDifferenceFunction<double> , opengm::PottsNFunction<double> ) , opengm::SimpleDiscreteSpace< int, int > >& gm,
		png::image< png::rgb_pixel >& img,
		int nx, int ny, int nLabels
	) {
	auto variableIndex = [&nx](const int x, const int y) {return x + nx*y;};
	
	// Create functions.
	const int shape[] = {nLabels, nLabels};
	double lambda1, lambda2, lambda3;
	lambda1 = LAMBDA1;
	lambda2 = LAMBDA2;
	lambda3 = LAMBDA3; // hash defined at top of the file.
	double deltaf;
	int c1,c2,c3;
	c1=c2=c3=0;
	opengm::ExplicitFunction<double> f1(shape, shape+2);
	opengm::ExplicitFunction<double> f2(shape, shape+2);
	opengm::ExplicitFunction<double> f3(shape, shape+2);
	int labels[] = {0,0};
	for(int l1=0; l1<nLabels; ++l1) {
		for(int l2=0; l2<nLabels; ++l2) {
			labels[0] = l1;
			labels[1] = l2;
			if(abs(l1-l2)==0) {
				f1.template operator()<int*>(labels) = 0;
				f2.template operator()<int*>(labels) = 0;
				f3.template operator()<int*>(labels) = 0;
			} else if ( abs(l1-l2)==1 ) {
				f1.template operator()<int*>(labels) = lambda1*1.0/6;
				f2.template operator()<int*>(labels) = lambda2*1.0/6;
				f3.template operator()<int*>(labels) = lambda3*1.0/6;
			} else {
				f1.template operator()<int*>(labels) = lambda1*1;
				f2.template operator()<int*>(labels) = lambda2*1;
				f3.template operator()<int*>(labels) = lambda3*1;
			}
		}
	}
	Model::FunctionIdentifier fid1 = gm.addFunction(f1);
	Model::FunctionIdentifier fid2 = gm.addFunction(f2);
	Model::FunctionIdentifier fid3 = gm.addFunction(f3);

	for(int y=0; y<ny; ++y) {
		for(int x=0; x<nx; ++x) {
			// down neighbour.
			if(x+1 < nx) {
				deltaf = absolute_color_difference(img[y][x], img[y][x+1]);
				int vars[] = { variableIndex(x,y) , variableIndex(x+1,y) };
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
				deltaf = absolute_color_difference(img[y][x], img[y+1][x]);
				int vars[] = { variableIndex(x,y) , variableIndex(x,y+1) };
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
	// Print out function count.
	std::cout << "Number of functions of different types:\n\tType1 (<7)=" << c1 << "\n\tType2 (<15)=" << c2 << "\n\tType3 (ow)=" << c3 << std::endl;
}

void addPairwisePotentials(
	// arguments
		opengm::GraphicalModel<double, opengm::Adder, OPENGM_TYPELIST_3(opengm::ExplicitFunction<double> , opengm::TruncatedAbsoluteDifferenceFunction<double> , opengm::PottsNFunction<double> ) , opengm::SimpleDiscreteSpace< int, int > >& gm,
		png::image< png::rgb_pixel >& img,
		int nx, int ny, int nLabels,
		std::vector< std::vector<int> >& syms, std::vector<int>& invsyms
	) {
	auto variableIndex = [&nx](const int x, const int y) {return x + nx*y;};
	
	// Create functions.
	const int shape[] = {nLabels, nLabels};

	double lambda1, lambda2, lambda3;
	lambda1 = LAMBDA1;
	lambda2 = LAMBDA2;
	lambda3 = LAMBDA3; // hash defined at top of the file.
	int c;
	double deltaf;
	int labels[] = {0,0};

	// Key is a vector of size 2. value is a vector of size 3.
	std::map< std::vector<int> , std::vector<int> > pair2weights;
	std::map< std::vector<int> , std::vector<int> >::iterator findVars;
	int g1, g2;
	for(int y=0; y<ny; ++y) {
		for(int x=0; x<nx; ++x) {
			g1 = invsyms[variableIndex(x,y)];
			// Rightward potential
			if( x<(nx-1) ) {
				g2 = invsyms[variableIndex(x+1, y)];
				if(g1!=g2) {
					std::vector<int> vars = {g1,g2};
					std::sort(vars.begin(), vars.end());
					findVars = pair2weights.find(vars);
					if(findVars == pair2weights.end()) {
						findVars = pair2weights.insert(std::pair< std::vector<int>, std::vector<int> >(vars, std::vector<int>(3,0))).first;
					}
					deltaf = absolute_color_difference(img[y][x], img[y][x+1]);
					if (deltaf<u1) {
						(findVars->second)[0] = 1 + findVars->second[0];
					} else if (deltaf < u2) {
						(findVars->second)[1] = 1 + findVars->second[1];
					} else {
						(findVars->second)[2] = 1 + findVars->second[2];
					}
				}
			}
			// Downward potential
			if( y<(ny-1) ) {
				g2 = invsyms[variableIndex(x,y+1)];
				if(g1!=g2) {
					std::vector<int> vars = {g1,g2};
					std::sort(vars.begin(), vars.end());
					findVars = pair2weights.find(vars);
					if(findVars == pair2weights.end()) {
						findVars = pair2weights.insert(std::pair< std::vector<int>, std::vector<int> >(vars, std::vector<int>(3,0))).first;
					}
					deltaf = absolute_color_difference(img[y][x], img[y][x+1]);
					if (deltaf<u1) {
						(findVars->second)[0] = 1 + findVars->second[0];
					} else if (deltaf < u2) {
						(findVars->second)[1] = 1 + findVars->second[1];
					} else {
						(findVars->second)[2] = 1 + findVars->second[2];
					}
				}		
			}
		}
	}
	// Now to iterate over pair2weights
	typedef std::map< std::vector<int> , std::vector<int> >::iterator map_iterator_type;
	std::vector<int> vars(2);
	std::map< std::vector<int> , Model::FunctionIdentifier > weight2fid;	
	double lambda;
	double one_sixth_of_lambda = lambda/6.0;
	Model::FunctionIdentifier fid;	
	const int pairwiseShape[] = {nLabels, nLabels};
	for(map_iterator_type iter = pair2weights.begin(); iter!=pair2weights.end(); ++iter) {
		vars[0] = (iter->first)[0];
		vars[1] = (iter->first)[1];
		if (weight2fid.find(iter->second) != weight2fid.end()) {
			fid = weight2fid[iter->second];
		} else {
			opengm::ExplicitFunction<double> f(pairwiseShape, pairwiseShape+2);			
			int labels[] = {0,0};
			lambda = LAMBDA1*iter->second[0] + LAMBDA2*iter->second[1] + LAMBDA3*iter->second[2];
			one_sixth_of_lambda = lambda/6.0;
			for(int l1=0; l1<nLabels; ++l1) {
				for(int l2=0; l2<nLabels; ++l2) {
					labels[0] = l1;
					labels[1] = l2;
					if(abs(l1-l2)>1) {
						f.template operator()<int*>(labels) = lambda;
					} else if (abs(l1-l2)==1) {
						f.template operator()<int*>(labels) = one_sixth_of_lambda;
					} else {
						f.template operator()<int*>(labels) = 0;
					}
				}
			}
			fid = gm.addFunction(f);
			weight2fid[iter->second] = fid;			
		}
		gm.addFactor(fid, vars.begin(), vars.end());
	}
	c = weight2fid.size();
	std::cout << "Number of pairwise potentials is=" << c << std::endl;
}


// This function isn't used.
void getPairwiseTypes(
	// arguments
		png::image< png::rgb_pixel >& img,
		std::vector<std::vector<int> >& factorType2right,
		std::vector<std::vector<int> >& factorType2down
	) {

	int nx = img.get_width();
	int ny = img.get_height();

	// ASSERT: factorType2right and factorType2down are of the correct size
	double deltaf;

	for(int y=0; y<ny; ++y) {
		for(int x=0; x<nx; ++x) {
			// Figure out the factor types.
			if(x+1 < nx) {
				deltaf = absolute_color_difference(img[y][x], img[y][x+1]);
				if (deltaf<u1) {
					factorType2right[y][x] = 1;
				} else if (deltaf < u2) {
					factorType2right[y][x] = 2;
				} else {
					factorType2right[y][x] = 3;
				}
			}

			if(y+1 < ny) {
				deltaf = absolute_color_difference(img[y][x], img[y+1][x]);
				if (deltaf<u1) {
					factorType2down[y][x] = 1;
				} else if (deltaf < u2) {
					factorType2down[y][x] = 2;
				} else {
					factorType2down[y][x] = 3;
				}
			}
		}
	}
}

Model constructGraphicalModel(png::image< png::rgb_pixel >& img, int nLabels, double* mcost) {
	int nx=img.get_width();
	int ny=img.get_height();
	Space space(nx*ny, nLabels);
	Model gm(space);
	addUnaryPotentials(gm, mcost, nx, ny, nLabels); // overloaded
	addPairwisePotentials(gm, img, nx, ny, nLabels); // overloaded
	return gm; // Use RVO please :<
}

Model constructGraphicalModel(
		png::image< png::rgb_pixel >& img, int nLabels, double* mcost, 
		std::vector< std::vector<int> >& syms, std::vector<int>& invsyms
	) {
	clock_t start, end;
	int nx=img.get_width();
	int ny=img.get_height();	
	int nVariables = syms.size();
	Space space(nVariables, nLabels);
	Model gm(space);
	// Don't have to worry about pairwise potentials within the same group because they're identical.
	std::cout << "starting potential addition" << std::endl;
	start=clock();
	addUnaryPotentials(gm, mcost, nx, ny, nLabels, syms, invsyms); // overloaded
	end = clock();
	std::cout << "done adding unary potentials " << double(end-start)/CLOCKS_PER_SEC << "s" << std::endl;
	start=clock();
	addPairwisePotentials(gm, img, nx, ny, nLabels, syms, invsyms); // overloaded
	end = clock();
	std::cout << "done adding pairwise potentials " << double(end-start)/CLOCKS_PER_SEC << "s" << std::endl;
	return gm;
}

#endif
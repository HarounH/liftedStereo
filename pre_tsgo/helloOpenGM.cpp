/**
	@author Haroun
	It's been a while, dear C++.
*/

#include <iostream>
#include <opengm/opengm.hxx>
#include <opengm/graphicalmodel/graphicalmodel.hxx>
#include <opengm/operations/minimizer.hxx>
#include <opengm/operations/adder.hxx>
#include <opengm/inference/icm.hxx>
#include <opengm/graphicalmodel/space/simplediscretespace.hxx>
#include <opengm/functions/potts.hxx>
#include <opengm/functions/pottsn.hxx>
#include <opengm/inference/messagepassing/messagepassing.hxx>
#include <opengm/inference/gibbs.hxx>
#include <typeinfo>

#include "useSym.hpp" // contains reduceGraphUsingSymmetries
#include "compressFactorGraph.hpp" // why nut.
using namespace std;

//This is the example code from openGM's html documentation.
int hello1() {

	// construct a graphical model with 
   // - 5 variables with {5, 5, 2, 2, 10} labels
   // - addition as the operation (template parameter Adder)
   typedef opengm::GraphicalModel<float, opengm::Adder> Model;  
   size_t numbersOfLabels[] = {5, 5, 2, 2, 10};
   Model gm(opengm::DiscreteSpace<>(numbersOfLabels, numbersOfLabels + 5));

   // add 1st order functions and factors to the model
   typedef opengm::ExplicitFunction<float> ExplicitFunction;
   typedef Model::FunctionIdentifier FunctionIdentifier;
   for(size_t variable = 0; variable < gm.numberOfVariables(); ++variable) {
      // construct 1st order function
      const size_t shape[] = {gm.numberOfLabels(variable)};
      ExplicitFunction f(shape, shape + 1);
      for(size_t state = 0; state < gm.numberOfLabels(variable); ++state) {
         f(state) = float(rand()) / RAND_MAX; // random toy data
      }
      // add function
      FunctionIdentifier id = gm.addFunction(f);
      // add factor
      size_t variableIndex[] = {variable};
      gm.addFactor(id, variableIndex, variableIndex + 1);
   }

   // add 3rd order functions and factors to the model
   for(size_t variable1 = 0; variable1 < gm.numberOfVariables(); ++variable1)
   for(size_t variable2 = variable1 + 1; variable2 < gm.numberOfVariables(); ++variable2)
   for(size_t variable3 = variable2 + 1; variable3 < gm.numberOfVariables(); ++variable3) {
      const size_t shape[] = {
         gm.numberOfLabels(variable1),
         gm.numberOfLabels(variable2),
         gm.numberOfLabels(variable3)
      };
      // construct 3rd order function
      ExplicitFunction f(shape, shape + 3);
      for(size_t state1 = 0; state1 < gm.numberOfLabels(variable1); ++state1)
      for(size_t state2 = 0; state2 < gm.numberOfLabels(variable2); ++state2)
      for(size_t state3 = 0; state3 < gm.numberOfLabels(variable3); ++state3) {         
         f(state1, state2, state3) = float(rand()) / RAND_MAX; // random toy data
      }
      FunctionIdentifier id = gm.addFunction(f);
      // sequences of variable indices need to be (and in this case are) sorted
      size_t variableIndexSequence[] = {variable1, variable2, variable3};
      gm.addFactor(id, variableIndexSequence, variableIndexSequence + 3);
   }

   // set up the optimizer (ICM)
   typedef opengm::ICM<Model, opengm::Minimizer> IcmType;
   typedef IcmType::VerboseVisitorType VerboseVisitorType;
   IcmType icm(gm);

   // obtain the (approximate) argmin
   VerboseVisitorType verboseVisitor;
   icm.infer(verboseVisitor);

   // output the (approximate) argmin
   vector<size_t> argmin;
   icm.arg(argmin);
   for(size_t variable = 0; variable < gm.numberOfVariables(); ++variable) {
      cout << "x" << variable << "=" << argmin[variable] << "\n";
   }
}


	
void hello2() {
	using namespace opengm;
	const size_t nx = 30; // width of the grid
	const size_t ny = 30; // height of the grid
	const size_t numberOfLabels = 5;

	// (x,y) -> flat array index
	auto variableIndex = [](const size_t x, const size_t y) {return x + nx*y;};
	
	double lambda = 0.1; // coupling strength of the Potts model
	// construct a label space with
	// - nx * ny variables 
	// - each having numberOfLabels many labels
	typedef SimpleDiscreteSpace<size_t, size_t> Space;
	Space space(nx * ny, numberOfLabels);

	// construct a graphical model with 
	// - addition as the operation (template parameter Adder)
	// - support for Potts functions (template parameter PottsFunction<double>)
	typedef GraphicalModel<double, Adder, OPENGM_TYPELIST_3(ExplicitFunction<double> , PottsFunction<double> , opengm::PottsNFunction<double> ) , Space> Model;
	Model gm(space);

	// for each node (x, y) in the grid, i.e. for each variable 
	// variableIndex(x, y) of the model, add one 1st order functions 
	// and one 1st order factor
	for(size_t y = 0; y < ny; ++y) 
	for(size_t x = 0; x < nx; ++x) {
	  // function
	  const size_t shape[] = {numberOfLabels};
	  ExplicitFunction<double> f(shape, shape + 1);
	  for(size_t s = 0; s < numberOfLabels; ++s) {
	     f(s) = (1.0 - lambda) * rand() / RAND_MAX;
	  }
	  Model::FunctionIdentifier fid = gm.addFunction(f);

	  // factor
	  size_t variableIndices[] = {variableIndex(x, y)};
	  gm.addFactor(fid, variableIndices, variableIndices + 1);
	}

	// add one (!) 2nd order Potts function
	PottsFunction<double> f(numberOfLabels, numberOfLabels, 0.0, lambda);
	Model::FunctionIdentifier fid = gm.addFunction(f);

	// for each pair of nodes (x1, y1), (x2, y2) which are adjacent on the grid,
	// add one factor that connects the corresponding variable indices and 
	// refers to the Potts function
	for(size_t y = 0; y < ny; ++y) 
	for(size_t x = 0; x < nx; ++x) {
	  if(x + 1 < nx) { // (x, y) -- (x + 1, y)
	     size_t variableIndices[] = {variableIndex(x, y), variableIndex(x + 1, y)};
	     sort(variableIndices, variableIndices + 2);
	     gm.addFactor(fid, variableIndices, variableIndices + 2);
	  }
	  if(y + 1 < ny) { // (x, y) -- (x, y + 1)
	     size_t variableIndices[] = {variableIndex(x, y), variableIndex(x, y + 1)};
	     sort(variableIndices, variableIndices + 2);
	     gm.addFactor(fid, variableIndices, variableIndices + 2);
	  }
	}    

	// set up the optimizer (loopy belief propagation)
	typedef BeliefPropagationUpdateRules<Model, opengm::Minimizer> UpdateRules;
	typedef MessagePassing<Model, opengm::Minimizer, UpdateRules, opengm::MaxDistance> BeliefPropagation;
	const size_t maxNumberOfIterations = 40;
	const double convergenceBound = 1e-7;
	const double damping = 0.5;
	BeliefPropagation::Parameter parameter(maxNumberOfIterations, convergenceBound, damping);
	BeliefPropagation bp(gm, parameter);

	// optimize (approximately)
	BeliefPropagation::VerboseVisitorType visitor;
	bp.infer(visitor);

	// obtain the (approximate) argmin
	vector<size_t> labeling(nx * ny);
	bp.arg(labeling);

	// output the (approximate) argmin
	size_t variableIdx = 0;
	for(size_t y = 0; y < ny; ++y) {
	  for(size_t x = 0; x < nx; ++x) {
	     cout << labeling[variableIdx] << ' ';
	     ++variableIdx;
	  }
	  cout << endl;
	}
}

// Trying to iterate over graph elements.
void hello3() {
	using namespace opengm;
	const size_t nx = 30; // width of the grid
	const size_t ny = 30; // height of the grid
	const size_t numberOfLabels = 5;

	// (x,y) -> flat array index
	auto variableIndex = [](const size_t x, const size_t y) {return x + nx*y;};
	
	double lambda = 0.1; // coupling strength of the Potts model
	// construct a label space with
	// - nx * ny variables 
	// - each having numberOfLabels many labels
	typedef SimpleDiscreteSpace<size_t, size_t> Space;
	Space space(nx * ny, numberOfLabels);

	// construct a graphical model with 
	// - addition as the operation (template parameter Adder)
	// - support for Potts functions (template parameter PottsFunction<double>)
	typedef GraphicalModel<double, Adder, OPENGM_TYPELIST_3(ExplicitFunction<double> , PottsFunction<double> , opengm::PottsNFunction<double>) , Space> Model;
	Model gm(space);

	cout << gm.space().numberOfVariables() << endl; return;
	// for each node (x, y) in the grid, i.e. for each variable 
	// variableIndex(x, y) of the model, add one 1st order functions 
	// and one 1st order factor
	for(size_t y = 0; y < ny; ++y) 
	for(size_t x = 0; x < nx; ++x) {
	  // function
	  const size_t shape[] = {numberOfLabels};
	  ExplicitFunction<double> f(shape, shape + 1);
	  for(size_t s = 0; s < numberOfLabels; ++s) {
	     f(s) = (1.0 - lambda) * rand() / RAND_MAX;
	  }
	  Model::FunctionIdentifier fid = gm.addFunction(f);

	  // factor
	  size_t variableIndices[] = {variableIndex(x, y)};
	  gm.addFactor(fid, variableIndices, variableIndices + 1);
	}

	// add one (!) 2nd order Potts function
	PottsFunction<double> f(numberOfLabels, numberOfLabels, 0.0, lambda);
	Model::FunctionIdentifier fid = gm.addFunction(f);



	// for each pair of nodes (x1, y1), (x2, y2) which are adjacent on the grid,
	// add one factor that connects the corresponding variable indices and 
	// refers to the Potts function
	for(size_t y = 0; y < ny; ++y) 
	for(size_t x = 0; x < nx; ++x) {
	  if(x + 1 < nx) { // (x, y) -- (x + 1, y)
	     size_t variableIndices[] = {variableIndex(x, y), variableIndex(x + 1, y)};
	     sort(variableIndices, variableIndices + 2);
	     gm.addFactor(fid, variableIndices, variableIndices + 2);
	  }
	  if(y + 1 < ny) { // (x, y) -- (x, y + 1)
	     size_t variableIndices[] = {variableIndex(x, y), variableIndex(x, y + 1)};
	     sort(variableIndices, variableIndices + 2);
	     gm.addFactor(fid, variableIndices, variableIndices + 2);
	  }
	}

	// the following is just messing around code by H.
	cout << "Graph created." << endl;
	cout << "nVariables = " << gm.numberOfVariables() << endl;
	cout << "nFactors = " << gm.numberOfFactors() << endl;
	int v = 3;
	cout << "nFactors of variable " << v << " = " << gm.numberOfFactors(v) << endl;
	cout << "factor order=" << gm.factorOrder() << endl;
	int k = 4;
	auto factorK = gm[gm.factorOfVariable(v,k-1)];

	const size_t shape[] = {numberOfLabels};
	
	ExplicitFunction<double> f2(shape, shape + 1);
	for(size_t s = 0; s < numberOfLabels; ++s) {
		f2(s) = (1.0 - lambda) * rand() / RAND_MAX;
	}
	
	ExplicitFunction<double> f3(shape, shape + 1);
	for(size_t s = 0; s < numberOfLabels; ++s) {
		f3(s) = (1.0 - lambda) * rand() / RAND_MAX;
	}

	ExplicitFunction<double> f4(shape, shape+1);
	for(size_t s = 0; s < numberOfLabels; ++s) {
		f4(s) = f2(s) + f3(s);
	}

	cout << f2.asString() << endl;
	

	// how do i deal with factorK?
}

void hello4() {
	typedef opengm::GraphicalModel<float, opengm::Adder> Model;  
	size_t numbersOfLabels[] = {2,2};
	Model gm(opengm::DiscreteSpace<>(numbersOfLabels, numbersOfLabels + 2));
	
	typedef opengm::ExplicitFunction<float> ExplicitFunction;
	typedef Model::FunctionIdentifier FunctionIdentifier;
	for(size_t variable = 0; variable < gm.numberOfVariables(); ++variable) {
		// construct 1st order function
		const size_t shape[] = {gm.numberOfLabels(variable)};
		ExplicitFunction f(shape, shape + 1);
		for(size_t state = 0; state < gm.numberOfLabels(variable); ++state) {
			f(state) = float(rand()) / RAND_MAX; // random toy data
			cout << f(state) << " ";
		}
		cout << endl;
		// add function
		FunctionIdentifier id = gm.addFunction(f);
		// add factor
		size_t variableIndex[] = {variable};
		gm.addFactor(id, variableIndex, variableIndex + 1);
	}
	for(size_t fid=0; fid< gm.numberOfFactors(); ++fid) {
		for(size_t s=0; s<numbersOfLabels[fid]; ++s) {
			cout << gm[fid](s) << " ";
		}
		cout << endl;
	}
}

void hello5() {
	typedef opengm::GraphicalModel<double, opengm::Adder, OPENGM_TYPELIST_3(opengm::ExplicitFunction<double> , opengm::PottsFunction<double> , opengm::PottsNFunction<double> ) > Model;
	size_t numbersOfLabels[] = {2,2,2};
	Model gm(opengm::DiscreteSpace<>(numbersOfLabels, numbersOfLabels + 3));
	size_t nLabels = 2;
	double lambda = 10.0;

	opengm::PottsFunction<double> f(nLabels, nLabels, 0.0, lambda);
	Model::FunctionIdentifier fid = gm.addFunction(f);
	size_t vars[] = {0,1};
	gm.addFactor(fid, vars, vars+2);

	cout << "nFactors=" << gm.numberOfFactors() << endl;
	// opengm::PottsFunction<double>::
	for(size_t s1=0; s1<nLabels; ++s1) {
		for(size_t s2=0; s2<nLabels; ++s2) {
			size_t vals[] = {s1,s2};
			cout << gm[0].operator()<size_t*>(vals) << " ";
		}
		cout << endl;
	}
}

void hello6() {
	typedef opengm::SimpleDiscreteSpace<size_t, size_t> Space;
	Space space(6, 2); // we're doing pixel labelling afterall.
	typedef opengm::GraphicalModel<double, opengm::Adder, OPENGM_TYPELIST_3(opengm::ExplicitFunction<double> , opengm::PottsFunction<double> , opengm::PottsNFunction<double>) , Space> Model;
	Model gm(space);

	// make 1,2,3 symmetric.
	std::vector< std::vector<size_t> > groups;
	std::vector<size_t> gp0 = {0,1,2};
	std::vector<size_t> gp1 = {3,4,5};

	groups.push_back(gp0);
	groups.push_back(gp1);
	
	std::vector<size_t> invsyms = {0,0,0,1,1,1};

	// Model rgm = reduceGraphUsingSymmetries<double, opengm::Adder, OPENGM_TYPELIST_3(opengm::ExplicitFunction<double> , opengm::PottsFunction<double> ) , Space>(gm, groups, invsyms);

	// step1 : verify variables in rgm.
	// cout << "rgm.nvar=" << rgm.numberOfVariables() << endl;


	// step2 : add some unary factors to gm, see if they are reflected in rgm
	size_t nLabels = 2;
	const size_t unaryShape[] = {nLabels};
	for(size_t v=0; v<gm.numberOfVariables(); ++v) {
		opengm::ExplicitFunction<double> f(unaryShape, unaryShape+1);
		for (size_t l=0; l<nLabels; ++l) {
			f(l) = (float(rand()) / RAND_MAX);
		}
		Model::FunctionIdentifier fid = gm.addFunction(f);
		size_t vars[] = {v};
		gm.addFactor(fid, vars, vars+1);
	}


	// step3: add some pairwise factors.
	opengm::PottsFunction<double> f(nLabels, nLabels, 10.0, 100.0);
	Model::FunctionIdentifier fid = gm.addFunction(f);
	size_t g1int[] = {0,1};
	size_t g2int[] = {3,4};
	size_t g3int[] = {2,5};
	gm.addFactor(fid, g1int, g1int+2);
	gm.addFactor(fid, g2int, g2int+2);
	gm.addFactor(fid, g3int, g3int+2);



	Model rgm = reduceGraphUsingSymmetries<double, opengm::Adder, OPENGM_TYPELIST_3(opengm::ExplicitFunction<double> , opengm::PottsFunction<double> , opengm::PottsNFunction<double>) , Space>(gm, groups, invsyms);
	cout << "gm.nfac = " << gm.numberOfFactors() << endl;
	cout << "rgm.nfac= " << rgm.numberOfFactors() << endl;

	// print gm factors
	cout << "gm factors:" << endl;
	for(size_t fid=0; fid<gm.numberOfFactors(); ++fid) {
		if(gm[fid].dimension()==1) {
			cout << "var=" << gm.variableOfFactor(fid, 0) << " "; 
			for(size_t l=0; l<nLabels; ++l) {
				size_t labels[] = {l};
				cout << gm[fid].template operator()<size_t*>(labels) << " ";
			}
		} else {
			cout << "var=" << gm.variableOfFactor(fid,0) << "," << gm.variableOfFactor(fid,1) << "\n";
			for(size_t l1=0; l1<nLabels; ++l1) {
				cout << "\t";
				for(size_t l2=0; l2<nLabels; ++l2) {
					size_t ls[] = {l1,l2};
					cout << gm[fid].template operator()<size_t*>(ls) << " ";
				}
				cout << endl;
			}
		}
		cout << endl;
	}

	// print rgm factors
	cout << "rgm factors" << endl;
	for(size_t fid=0; fid<rgm.numberOfFactors(); ++fid) {
		if(rgm[fid].dimension()==1) {
			cout << "var=" << rgm.variableOfFactor(fid, 0) << " "; 
			for(size_t l=0; l<nLabels; ++l) {
				size_t labels[] = {l};
				cout << rgm[fid].template operator()<size_t*>(labels) << " ";
			}
		} else {
			cout << "var=" << rgm.variableOfFactor(fid,0) << "," << rgm.variableOfFactor(fid,1) << "\n";
			for(size_t l1=0; l1<nLabels; ++l1) {
				cout << "\t";
				for(size_t l2=0; l2<nLabels; ++l2) {
					size_t ls[] = {l1,l2};
					cout << rgm[fid].template operator()<size_t*>(ls) << " ";
				}
				cout << endl;
			}
		}
		cout << endl;
	}
}

/*
	function to test compressFactorGraph.hpp
*/
void hello7() {
	// construct normal graph.
	typedef opengm::SimpleDiscreteSpace<size_t, size_t> Space;
	Space space(6, 3); // we're doing pixel labelling afterall.
	typedef opengm::GraphicalModel<double, opengm::Adder, OPENGM_TYPELIST_3(opengm::ExplicitFunction<double> , opengm::PottsFunction<double> , opengm::PottsNFunction<double>) , Space> Model;
	Model gm(space);
	size_t nLabels = 3;
	size_t nVars = 6;

	// unaries
	size_t unaryShape[] = {nLabels};
	for(size_t vid=0; vid<nVars; ++vid) {
		if(vid%2==1) {
			// copy the unary for vid-1 ... because why not.
			continue;
		} else { // generate a random unary.
			opengm::ExplicitFunction<double> f(unaryShape, unaryShape+1);
			cout << "vars={" << vid << "," << vid+1 << "}" << endl;
			for(size_t l=0; l<nLabels; ++l) {
				f(l) = (float(rand()) / RAND_MAX);
				cout << f(l) << " ";
			}
			cout << endl;
			Model::FunctionIdentifier fid1 = gm.addFunction(f);
			Model::FunctionIdentifier fid2 = gm.addFunction(f);
			size_t vars[] = {vid};
			gm.addFactor(fid1, vars, vars+1);
			vars[0] = vid+1;
			gm.addFactor(fid2, vars, vars+1);
		}
	}

	// pairwise.
	opengm::PottsFunction<double> fpotts(nLabels, nLabels, 0.0, 0.5);
	Model::FunctionIdentifier fid = gm.addFunction(fpotts);
	// add fid over 6 different pairs.
	size_t vars[] = {2,3};
	gm.addFactor(fid, vars, vars+2);
	vars[1] = 4; gm.addFactor(fid, vars, vars+2);
	vars[1] = 5; gm.addFactor(fid, vars, vars+2);
	vars[0] = 3; gm.addFactor(fid, vars, vars+2);
	vars[0] = 4; gm.addFactor(fid, vars, vars+2);
	// vars[0] = 4 ;vars[1] = 5; gm.addFactor(fid, vars, vars+2);
	
	// pass it to compressFactorGraph
	std::vector< std::vector<size_t> > varGroups = getVarGroupings<double, opengm::Adder, OPENGM_TYPELIST_3(opengm::ExplicitFunction<double> , opengm::PottsFunction<double> , opengm::PottsNFunction<double>) , Space>(gm, 1);
	
	// read symmetries
	for(size_t gid=0; gid< varGroups.size(); ++gid) {
		// group gid needs to be printed, along with its factor.
		cout << "group " << gid << " ";
		for(size_t vidid=0; vidid<varGroups[gid].size(); ++vidid) {
			cout << varGroups[gid][vidid] << " ";
		}
		cout << endl;
	}
}


// trying out higher order potts functions.
void hello8() {
	typedef opengm::GraphicalModel<double, opengm::Adder, OPENGM_TYPELIST_3(opengm::ExplicitFunction<double> , opengm::PottsFunction<double> , opengm::PottsNFunction<double> ) , opengm::SimpleDiscreteSpace<size_t, size_t> > Model;
	size_t numbersOfLabels[] = {2,2,2,2,2};
	size_t nLabels = 2;
	Model gm(opengm::SimpleDiscreteSpace<size_t, size_t>(6, nLabels));
	double lambda = 10.0;

	
	
	size_t shape[] = {2,2,2};
	opengm::PottsNFunction<double> f(shape, shape+3, 0.0, lambda);
	Model::FunctionIdentifier fid = gm.addFunction(f);
	size_t vars[] = {0,1,2};
	std::sort(vars, vars+3);
	gm.addFactor(fid, vars, vars+3);
	cout << "nFactors=" << gm.numberOfFactors() << endl;
	// for(size_t s1=0; s1<nLabels; ++s1) {
	// 	for(size_t s2=0; s2<nLabels; ++s2) {
	// 		for(size_t s3=0; s3<nLabels; ++s3) {
	// 			size_t s[] = {s1,s2,s3};				
	// 			cout << s1 << "," << s2 << "," << s3 << ":" << gm[0].template operator()<size_t*>(s) << endl;
	// 		}
	// 	}
	// }
	vars[0] = 3;
	vars[1] = 4;
	vars[2] = 5;
	gm.addFactor(fid, vars, vars+3);
	
	// reduce the graph.
	std::vector< std::vector<size_t> > syms;
	syms.push_back(std::vector<size_t>()); syms.push_back(std::vector<size_t>()); syms.push_back(std::vector<size_t>());
	syms[0].push_back(0); syms[0].push_back(5);
	syms[1].push_back(1); syms[1].push_back(4);
	syms[2].push_back(2); syms[2].push_back(3);
	std::vector<size_t> invsyms(6);
	invsyms[0] = 0;
	invsyms[1] = 1;
	invsyms[2] = 2;
	invsyms[3] = 1;
	invsyms[4] = 1;
	invsyms[5] = 1;

	Model rgm = reduceGraphUsingSymmetries<double, opengm::Adder, OPENGM_TYPELIST_3(opengm::ExplicitFunction<double> , opengm::PottsFunction<double> , opengm::PottsNFunction<double> ) , opengm::SimpleDiscreteSpace<size_t, size_t> >
					(gm, syms, invsyms);

	cout << "nf= " << rgm.numberOfFactors() << endl;
	for(int fid=0; fid< rgm.numberOfFactors(); ++fid) {
		int dim = rgm[fid].dimension();
		cout << "factor " << fid << " of dim=" << dim << endl;
		
		std::vector<size_t> labels(dim,0);
		// TODO ... implement this loop nicely.
		rgm[fid].template operator()< std::vector<size_t>::iterator >(labels.begin());
	}
}

int main(int argc, char** argv) {
	// hello3();
	hello8();
	return 0;
}
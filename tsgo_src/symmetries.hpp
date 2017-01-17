#ifndef _SYMMETRIES_HPP_
#define _SYMMETRIES_HPP_
#include <boost/unordered_map.hpp>

struct sort_pred {
    bool operator()(const std::pair<double,int> &left, const std::pair<double,int> &right) {
        return left.first < right.first;
    }
};

int mostLikelyStatesSymmetries(
	// arguments
		opengm::GraphicalModel<double, opengm::Adder, OPENGM_TYPELIST_3(opengm::ExplicitFunction<double> , opengm::TruncatedAbsoluteDifferenceFunction<double> , opengm::PottsNFunction<double> ) , opengm::SimpleDiscreteSpace<int, int> >& gm,
		std::vector<int>& factorColors,
		int nRanksConsidered		
	) {
	clock_t start = clock();
	int nVariables = gm.numberOfVariables();
	int nFactors = gm.numberOfFactors();
	int nLabels = gm.numberOfLabels(0);
	boost::unordered_map< std::vector<int> , int > rank2color;
	boost::unordered_map< std::vector<int> , int >::const_iterator findRanks;
	// Colors 0,1,2 are given to the three pairwise potential functions.
	int nColorsUsed = 3;
	for(int fid=0; fid<nFactors; ++fid) {
		if( gm[fid].dimension()==1 ) {
			// get the nRanksConsidered most likely states.
			std::vector< std::pair<double, int> > potentialAndIdx(nLabels);
			int ls[] = {0};
			for(int i=0; i<nLabels; ++i) {
				ls[0] = i;
				potentialAndIdx[i] = std::pair<double, int>(gm[fid].template operator()<int*>(ls) , i);
			}
			std::sort(potentialAndIdx.begin(), potentialAndIdx.end(),sort_pred());
			std::vector<int> ranks;
			for(int i=0; i<nRanksConsidered; ++i) {
				ranks.push_back(potentialAndIdx[i].second);	
			}
			// assign it a color.
			findRanks = rank2color.find(ranks);
			if ( findRanks != rank2color.end() ) {
				factorColors[fid] = findRanks->second;
			} else {
				rank2color[ranks] = nColorsUsed;
				factorColors[fid] = nColorsUsed;
				nColorsUsed++;
			}
		} else if (gm[fid].dimension()==2) { // Give it a color based on 
			int labels[] = {0,2};
			double val = gm[fid].template operator()<int*>(labels);
			if (val == LAMBDA1) {
				factorColors[fid] = 0;
			} else if (val == LAMBDA2) {
				factorColors[fid] = 1;
			} else if (val == LAMBDA3) {
				factorColors[fid] = 2;
			} else {
				std::cout << "Houston, we have a problem here" << std::endl;
				std::cout << "Stumbled upong a pairwise potential which isn't of any of the types given" << std::endl;
				exit(1);
			}
		}
	}
	clock_t end = clock();
	std::cout << "Got most likely initialization using " << nColorsUsed << " colors" << std::endl;
	std::cout << "\tTime taken= " << double(end-start)/CLOCKS_PER_SEC << std::endl;
	return nColorsUsed;
}

int assignColorBySignature(
	// arguments
		std::vector< std::vector<int> >& signature , // IN
		std::vector<int>& color // OUT
	) {
	int nColorsOut = 0;
	int helperVarColor;	
	std::map< std::vector<int> , int > sig2color; // takes far less memory.
	for(int idx=0; idx<signature.size(); ++idx) {
		if (signature[idx].size()==1) { // disconnected part of factor graph.
			continue;
		}
		sort(signature[idx].begin()+1, signature[idx].end());
		if (sig2color.count(signature[idx]) != 0) {
			color[idx] = sig2color[signature[idx]];
		} else {
			color[idx] = nColorsOut;
			sig2color[signature[idx]] = nColorsOut;
			nColorsOut+=1;
		}
	}
	return nColorsOut;
}
void constructGroupsFromColors(std::vector<int>& colors, std::vector< std::vector<int> >& group) {
	std::vector<int> colorToGroup(colors.size(), -1);
	int nColorsSeen = 0;
	group.clear();
	for(int vid=0; vid<colors.size(); ++vid) {
		int vcolor = colors[vid];
		if(vcolor==-1) continue; // invalid object.
		int gid = colorToGroup[vcolor];
		if (gid==-1) {
			// this is the first time this color is being seen.
			colorToGroup[vcolor] = nColorsSeen;
			gid = nColorsSeen;
			group.push_back( std::vector<int>() );
			group[gid].push_back(vid);
			nColorsSeen++;
		} else {
			group[gid].push_back(vid);
		}
	}
}
bool detectGroupsChanged(std::vector< std::vector<int> >& newGrouping, std::vector< std::vector<int> >& oldGrouping) {
	if (newGrouping.size()!=oldGrouping.size()) { // trivial exit.
		return true;
	}
	int nGroups = newGrouping.size(); // ASSERT: newGrouping.size() == oldGrouping.size()
	int nVarsInGroup;
	for(int gid=0; gid<nGroups; ++gid) {
		if (newGrouping[gid].size()!=oldGrouping[gid].size()) {
			return true;
		}
		nVarsInGroup = newGrouping[gid].size();
		for(int vidid=0; vidid< nVarsInGroup; ++vidid) {
			if (newGrouping[gid][vidid] != oldGrouping[gid][vidid]) {
				return true;
			}
		}
	}
	return false;
}

void getSymmetries(
	// arguments
		opengm::GraphicalModel<double, opengm::Adder, OPENGM_TYPELIST_3(opengm::ExplicitFunction<double> , opengm::TruncatedAbsoluteDifferenceFunction<double> , opengm::PottsNFunction<double> ) , opengm::SimpleDiscreteSpace<int, int> >& gm,
		std::vector<std::vector<int> >& syms,
		std::vector<int>& invsyms,
		int nRanksConsidered,
		int nCBPIterations
	) { // function begins.
	int nVariables = gm.numberOfVariables();
	int nLabels = gm.numberOfLabels(0);
	int nFactors= gm.numberOfFactors();
	std::vector<std::vector<std::vector<int> > > factorGroups(2);
	std::vector<std::vector<std::vector<int> > > variableGroups(2);
	std::vector<int> factorColors(nFactors, 0);
	std::vector<int> variableColors(nVariables, 0);
	std::vector<std::vector<int> > factorSignatures(nFactors);
	std::vector<std::vector<int> > variableSignatures(nVariables);

	// Initiatelize signatures.
	for(int i=0; i<nFactors; ++i) {
		factorSignatures[i] = std::vector<int>(1 + gm.numberOfVariables(i));
		factorSignatures[i][0] = 0; // default color.
	}
	for(int i=0; i<nVariables; ++i) {
		variableSignatures[i] = std::vector<int>(1 + gm.numberOfFactors(i));
		variableSignatures[i][0] = 0; // 0th color.
	}

	int nFactorColors = mostLikelyStatesSymmetries(gm, factorColors, nRanksConsidered);
	for(int i=0; i<factorColors.size(); ++i) factorSignatures[i][0] = factorColors[i];
	
	int newid = 0;
	constructGroupsFromColors(factorColors, factorGroups[newid]);
	std::cout << "nFactorColors:" << nFactorColors << " factorGroups.size(): " << factorGroups[newid].size() << std::endl;
	// Send 1st iteration of messages from factors to variables.
	for(int vid=0; vid<nVariables; ++vid) {
		int max_fidid = gm.numberOfFactors(vid);
		for(int fidid=0; fidid<max_fidid; ++fidid) {
			int fid = gm.factorOfVariable(vid, fidid);
			variableSignatures[vid][fidid+1] = 0; // factorColors[fid];
			if ( gm[fid].dimension()==1 ) { // Its a unary!
				variableColors[vid] = factorColors[fid];
				variableSignatures[vid][fidid+1] = factorColors[fid];
				variableSignatures[vid][0] = factorColors[fid];
			}
		}
	}

	int nVariableColors = assignColorBySignature(variableSignatures, variableColors);
	constructGroupsFromColors(variableColors, variableGroups[newid]);
	int cbp_iter=0;
	while(cbp_iter<nCBPIterations) {
		newid = (newid+1)%2;
		std::cout << "Counting BP iteration #" << cbp_iter << ", " << variableGroups[ (newid+1) % 2].size() << " groups now" << std::endl;

		// Messages from variables to factors
		for(int fid=0; fid<nFactors; ++fid) {
			factorSignatures[fid][0] = factorColors[fid];
			int vid = -1;
			for(int vidid=0; vidid<gm.numberOfVariables(fid); ++vidid) {
				vid = gm.variableOfFactor(fid, vidid);
				factorSignatures[fid][1+vidid] = variableColors[vid]; // +1 because of the 
			}
		}

		assignColorBySignature(factorSignatures, factorColors);// assignColors by signatures
		constructGroupsFromColors(factorColors, factorGroups[newid]); // construct groups

		// Messages from factors to variables
		for(int vid=0; vid<nVariables; ++vid) {
			variableSignatures[vid][0] = variableColors[vid];
			int fid = -1;
			for(int fidid=0; fidid<gm.numberOfFactors(vid); ++fidid) {
				fid = gm.factorOfVariable(vid, fidid);
				variableSignatures[vid][1+fidid] = factorColors[fid];
			}
		}
		assignColorBySignature(variableSignatures, variableColors); // assignColors by signatures
		constructGroupsFromColors(variableColors, variableGroups[newid]);// construct groups

		if (!detectGroupsChanged(variableGroups[0], variableGroups[1])) {
			break;
		}
		cbp_iter++;
	}
	std::cout << "Completed CBP. " << variableGroups[newid].size() << " groups now" << std::endl;
	for(int i=0; i<variableGroups[newid].size(); ++i) {
		syms.push_back(std::vector<int>());
		for(int j=0; j<variableGroups[newid][i].size(); ++j) {
			syms[i].push_back(variableGroups[newid][i][j]);
			invsyms[ variableGroups[newid][i][j] ] = i;
		}
	}
	return;
}

void reduceMRF(
	// arguments
		opengm::GraphicalModel<double, opengm::Adder, OPENGM_TYPELIST_3(opengm::ExplicitFunction<double> , opengm::TruncatedAbsoluteDifferenceFunction<double> , opengm::PottsNFunction<double> ) , opengm::SimpleDiscreteSpace<int, int> >& gm,
		std::vector<std::vector<int> >& syms,
		std::vector<int>& invsyms,
		opengm::GraphicalModel<double, opengm::Adder, OPENGM_TYPELIST_3(opengm::ExplicitFunction<double> , opengm::TruncatedAbsoluteDifferenceFunction<double> , opengm::PottsNFunction<double> ) , opengm::SimpleDiscreteSpace<int, int> >*& rgm
	) { // function begins
	int nLabels = gm.numberOfLabels(0);
	int nFactors = gm.numberOfFactors();
	
	rgm = new Model(opengm::SimpleDiscreteSpace<int, int>(syms.size(), nLabels));
	
	// declare unaries.
	int unaryShape[] = {nLabels};
	std::vector< opengm::ExplicitFunction<double> > unaries(syms.size(), opengm::ExplicitFunction<double>(unaryShape, unaryShape+1));
	// declare pairwise.

	// Key is a vector of size 2. value is a vector of size 3.
	std::map< std::vector<int> , std::vector<int> > pair2weights;
	
	int dim;
	int vid1, vid2, gid1, gid2, potentialtype;
	for(int fid=0; fid<nFactors; ++fid) {
		dim = gm[fid].dimension();
		if(dim==1) { // definitely a unary function.
			int labels[] = {0};
			int gid = invsyms[ gm.variableOfFactor(fid,0) ];
			for(int l=0; l<nLabels; ++l) {
				labels[0] = l;
				unaries[gid](l) += gm[fid].template operator()<int*>(labels);
			}
		} else if (dim==2){ // decide based on the variables involved.
			vid1 = gm.variableOfFactor(fid, 0);
			vid2 = gm.variableOfFactor(fid, 1);
			gid1 = invsyms[vid1];
			gid2 = invsyms[vid2];
			if (gid1==gid2) { // Becomes a unary potential. Specifically, it doesn't matter.
				continue;
			} else { // Add it as a pairwise factor.
				std::vector<int> vars = {gid1, gid2};
				std::sort(vars.begin(), vars.end());
				if( pair2weights.find(vars)==pair2weights.end() ) {
					pair2weights[vars] = std::vector<int>(3,0);
				}
				int labels[] = {0,2};
				double val = gm[fid].template operator()<int*>(labels);
				if(val == LAMBDA1) {
					potentialtype = 0;
				} else if (val == LAMBDA2) {
					potentialtype = 1;
				} else if (val == LAMBDA3){
					potentialtype = 2;
				} else {
					std::cout << "unrecognized pairwise potential\n";
					exit(1);
				}
				pair2weights[vars][potentialtype] += 1;
			}
		} else {
			std::cout << "unrecognized higher order potential\n";
			exit(1);
		}
	}
	std::cout <<  "processed weights for reduceMRF. Now going to add it to rgm\n";
	// Now to add said symmetries to model.

	// unaries
	for(int gid=0; gid<unaries.size(); ++gid) {
		Model::FunctionIdentifier fid = rgm->addFunction( unaries[gid] );
		int vars[] = {gid};
		rgm->addFactor(fid, vars, vars+1);
	}

	// pairwise
	typedef std::map< std::vector<int> , std::vector<int> >::iterator pairIterator;
	std::map< std::vector<int> , Model::FunctionIdentifier > weight2fid;
	Model::FunctionIdentifier fid;
	const int pairwiseShape[] = {nLabels, nLabels};
	for(pairIterator iter = pair2weights.begin(); iter != pair2weights.end(); ++iter) {
		auto key = iter->first;
		auto val = iter->second;
		if ( weight2fid.find(val)!=weight2fid.end() ) {
			fid = weight2fid[val];
		} else {
			// Create an explicit function and add it.
			opengm::ExplicitFunction<double> f(pairwiseShape, pairwiseShape+2);
			int labels[] = {0,0};
			double lambda = val[0]*LAMBDA1 + val[1]*LAMBDA2 + val[2]*LAMBDA3;
			double one_sixth_of_lambda = lambda/6.0;
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
			fid = rgm->addFunction(f);
			weight2fid[val] = fid;
		}
		rgm->addFactor(fid, key.begin(), key.end());
	}
}

void processSymmetries(
	// arguments
		opengm::GraphicalModel<double, opengm::Adder, OPENGM_TYPELIST_3(opengm::ExplicitFunction<double> , opengm::TruncatedAbsoluteDifferenceFunction<double> , opengm::PottsNFunction<double> ) , opengm::SimpleDiscreteSpace<int, int> >& gm,
		opengm::GraphicalModel<double, opengm::Adder, OPENGM_TYPELIST_3(opengm::ExplicitFunction<double> , opengm::TruncatedAbsoluteDifferenceFunction<double> , opengm::PottsNFunction<double> ) , opengm::SimpleDiscreteSpace<int, int> >*& rgm,
		std::vector<std::vector<int> >& syms,
		std::vector<int>& invsyms,
		int nRanksConsidered = 5,
		int nCBPIterations = 1
	) {
	// Get symmetries
	std::cout << "Getting symmetries by \n\tLooking at " << nRanksConsidered << "\n\tIterations of CBP=" << nCBPIterations << std::endl;
	clock_t start = clock();
	getSymmetries(gm, syms, invsyms, nRanksConsidered, nCBPIterations);
	clock_t end = clock();
	
	std::cout << "Got symmetries in " << (double(end-start)/CLOCKS_PER_SEC) << "s" << std::endl;
	
	start = clock();
	// Create a reduced graphical model
	reduceMRF(gm, syms, invsyms, rgm);
	end = clock();
	std::cout << "Reduced MRF in " << (double(end-start)/CLOCKS_PER_SEC) << "s" << std::endl;
	std::cout << "Original MRF has " << gm.numberOfVariables() << " variables and " << gm.numberOfFactors() << " factors\n";	
	std::cout << "Reduced MRF has " << rgm->numberOfVariables() << " variables and " << rgm->numberOfFactors() << " factors\n";
}

#endif
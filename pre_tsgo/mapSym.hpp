#ifndef MAP_SYM_HPP
#define MAP_SYM_HPP


#include <iostream>
#include <opengm/opengm.hxx>
#include <opengm/graphicalmodel/graphicalmodel.hxx>

#include <opengm/operations/minimizer.hxx>
#include <opengm/operations/adder.hxx>

#include <opengm/graphicalmodel/space/simplediscretespace.hxx>
#include <opengm/functions/potts.hxx>

#include <vector>
#include <map>

#ifndef _GET_SYM_HPP_
#define _GET_SYM_HPP_
/**
	This file haphazardly implements binning of variables in graphical model given the graphical model.
*/

#include <iostream>
#include <opengm/opengm.hxx>
#include <opengm/graphicalmodel/graphicalmodel.hxx>

#include <opengm/operations/minimizer.hxx>
#include <opengm/operations/adder.hxx>

#include <opengm/graphicalmodel/space/simplediscretespace.hxx>
#include <opengm/functions/potts.hxx>

#include <vector>
#include <utility>
#include <map>

struct sort_pred {
    bool operator()(const std::pair<double,int> &left, const std::pair<double,int> &right) {
        return left.first < right.first;
    }
};

/*

*/
template<class T, class OPERATOR, class FUNCTION_TYPE_LIST, class SPACE>
size_t mapSymmetries(
	// arguments
		const opengm::GraphicalModel<T, OPERATOR, FUNCTION_TYPE_LIST, SPACE>& gm, 
		int nRanksConsidered, // this argument should get more special as time goes on.
		std::vector<int>& colors
	) { // function begins TODO

	// two variables are symmetric if their unaries have similar MAP thingies.
	int nVar = gm.numberOfVariables();
	int nLabels= gm.numberOfLabels(0);
	std::vector<int> var2unary(nVar,-1);
	std::map< std::vector<int> , int > rank2color;
	size_t nColorsUsed = 0;
	std::cout << "nRanksConsidered=" << nRanksConsidered << std::endl;
	int nFactorsOfVid = 0;
	for(int vid=0; vid<nVar; ++vid) {
		nFactorsOfVid = gm.numberOfFactors(vid);
		for(int fidid=0; fidid<nFactorsOfVid; ++fidid) {
			int fid = gm.factorOfVariable(vid, fidid);
			if(gm[fid].dimension()==1) {
				var2unary[vid] = fid;

				std::vector< std::pair<double, int> > potentialAndIdx(nLabels);
				size_t ls[] = {0};
				for(int i=0; i<nLabels; ++i) {
					ls[0] = i;
					potentialAndIdx[i] = std::pair<double, int>(gm[fid].template operator()<size_t*>(ls) , i);
				}

				std::sort(potentialAndIdx.begin(),
						potentialAndIdx.end(),
						sort_pred()
						);

				std::vector< int > ranks;

				
				for(int i=0; i<nRanksConsidered; ++i) {
					ranks.push_back(potentialAndIdx[i].second);
				}

				if ( rank2color.count(ranks) != 0 ) {
					colors[fid] = rank2color[ranks];
				} else {
					rank2color[ranks] = nColorsUsed;
					colors[fid] = nColorsUsed;
					nColorsUsed++;
				}
				break;
			}
		}
	}
	for(int i=0; i<colors.size();++i) {
		if (colors[i]<0) {
			colors[i]=nColorsUsed;
		}
	}
	return (nColorsUsed+1);
}

/*
	IMPORTANT 
		groups are assigned such that the smallest variables are assigned a group first.
		signatures are built such that sig[0] is object's old color, 
			which was used to obtain object's present color

		sig[1...] are the colors of neighbours in factor graph

		for each node in factor graph, 
			we must have:
				signature - facSignature, varSignature
				color - facColor, varColor
				group - facGroupings, varGroupings ... contain oldGroups and newGroups. indicated using iterationParity.
					groupings[(iterationParity+1)%2] is the latest one.
					groupings[iterationParity] is the one being constructed.
*/
					

/* 
	helper function for compressFactorGraph. 
	@params
		IN signature: the signature of each factor or variable
		IN nColorsIn: the number of colors present in 'sig[1...]'
		OUT colors: the colors assigned to each object for which a signature is given

	@returns
		OUT nColorsOut: the number of colors used to color the objects given.

	Since clusters and colors are the same for me, this function is all I need.
*/
size_t assignColorBySignature(
	// arguments
		std::vector< std::vector<int> >& signature , // IN
		size_t nColorsIn, // IN
		std::vector<int>& color // OUT
	) {
	size_t nColorsOut = 0;
	size_t helperVarColor;
	
	std::map< std::vector<int> , size_t > sig2color; // takes far less memory.
	for(int idx=0; idx<signature.size(); ++idx) {
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


/* 	
	helper function for compressFactorGraph
	@params
		IN colors: colors assigned to each object.
		OUT group: output of the group as a list of objects per group.
			Each group essentially corresponds to one color.
			it must be empty or it will be made empty.
	@returns
		Nothing

*/
void constructGroupsFromColors(std::vector<int>& colors, std::vector< std::vector<size_t> >& group) {
	std::vector<int> colorToGroup(colors.size(), -1);
	size_t nColorsSeen = 0;
	group.clear();
	for(size_t vid=0; vid<colors.size(); ++vid) {
		size_t vcolor = colors[vid];
		int gid = colorToGroup[vcolor];
		if (gid==-1) {
			// this is the first time this color is being seen.
			colorToGroup[vcolor] = nColorsSeen;
			gid = nColorsSeen;
			group.push_back( std::vector<size_t>() );
			group[gid].push_back(vid);
			nColorsSeen++;
		} else {
			group[gid].push_back(vid);
		}
	}
}
/* helper function for compressFactorGraph.
	@params
		Returns true if the oldGrouping and newGrouping are different.
*/
bool detectGroupsChanged(std::vector< std::vector<size_t> >& newGrouping, std::vector< std::vector<size_t> >& oldGrouping) {
	if (newGrouping.size()!=oldGrouping.size()) { // trivial exit.
		return true;
	}
	// if we assign colors in a topological fashion, we are good.
	// if no change, then groups in newGrouping and oldGrouping are either identical or disjoint... and viceversa.
	// infact, because of the way constructGroupsFromColors() behaves, the two groups are either identical or different.

	size_t nGroups = newGrouping.size(); // ASSERT: newGrouping.size() == oldGrouping.size()
	size_t nVarsInGroup;
	for(size_t gid=0; gid<nGroups; ++gid) {
		if (newGrouping[gid].size()!=oldGrouping[gid].size()) {
			return true;
		}
		nVarsInGroup = newGrouping[gid].size();
		for(size_t vidid=0; vidid< nVarsInGroup; ++vidid) {
			if (newGrouping[gid][vidid] != oldGrouping[gid][vidid]) {
				return true;
			}
		}
	}
	return false;
}



// symmetry 1 : compressFactorGraph from counting bp ... gives us 
template<class T, class OPERATOR, class FUNCTION_TYPE_LIST, class SPACE>
std::vector< std::vector<size_t> > getVarGroupingsUsingMapSyms(
		// arguments
			const opengm::GraphicalModel<T, OPERATOR, FUNCTION_TYPE_LIST, SPACE>& gm, 
			int nMaxIterations=-1,
			double nBinsPerUnaryDimension=1
		) {	// function begins.

	// the graphical model (gm) already stores it in a nice factor-graph-eseq version.
	int iteration = 0; // used to control the number of iterations in the loop.
	size_t iterationParity = 1;
	size_t nVariables = gm.numberOfVariables();
	size_t nFactors = gm.numberOfFactors();

	


	// things related to factors
	std::vector< std::vector< std::vector<size_t> > > facGroupings(2); // one grouping for present, one for old.
		// facGroupings[0].push_back(std::vector<size_t>(nFactors,-1)); facGroupings[0].push_back(std::vector<size_t>(nFactors,-1));
		// facGroupings[1].push_back(std::vector<size_t>(nFactors,-1)); facGroupings[1].push_back(std::vector<size_t>(nFactors,-1));
	std::vector<int> facColor(nFactors, -1);
	// std::cout << "about to initialize factor colors by binning" << std::endl;
	size_t nFacColors = mapSymmetries<T, OPERATOR, FUNCTION_TYPE_LIST, SPACE>(gm, nBinsPerUnaryDimension ,facColor);
	// std::cout << "initiatelized factor colors .. nFacColors=" << nFacColors << std::endl;
	constructGroupsFromColors(facColor, facGroupings[iterationParity]);
		// ASSERT: facGroupings[iterationParity].size() == nFacColors.
	// std::cout << "got initial factor groups ... nFacGroups=" << facGroupings[iterationParity].size() << std::endl;
	// reserve space.
	std::vector< std::vector<int> > facSignature(nFactors);
	for(size_t fid=0; fid<nFactors; ++fid) {
		facSignature[fid] = std::vector<int>(1 + gm.numberOfVariables(fid),-1);
		facSignature[fid][0] = facColor[fid];
	}
	

	// things related to variables
	std::vector< std::vector< std::vector<size_t> > > varGroupings(2);
		// varGroupings[0].push_back(std::vector<size_t>(nVariables,-1));	varGroupings[0].push_back(std::vector<size_t>(nVariables,-1));
		// varGroupings[1].push_back(std::vector<size_t>(nVariables,-1));	varGroupings[1].push_back(std::vector<size_t>(nVariables,-1));
	std::vector<int> varColor(nVariables, -1);
	std::vector< std::vector<int> > varSignature(nVariables);
	// reserve space and initialize colors/signatures.
	size_t nFactorsForVid;
	for(size_t vid=0; vid<nVariables; ++vid) {
		nFactorsForVid = gm.numberOfFactors(vid);
		varSignature[vid] = std::vector<int>(1 + nFactorsForVid,-1);
		varSignature[vid][0] = 0; // all variables are unknown for now. so...
		for(size_t fidid=0; fidid<nFactorsForVid; ++fidid) {
			size_t fid = gm.factorOfVariable(vid, fidid);
			varSignature[vid][1+fidid] = facColor[fid];
		}
	}
	size_t nVarColors = assignColorBySignature(varSignature, nFacColors, varColor);
	constructGroupsFromColors(varColor, varGroupings[iterationParity]);
	// std::cout << " before we start, we end up with " << varGroupings[iterationParity].size() << " varGroups and " << nVarColors << " varColors" << std::endl;
		// ASSERT: nVarColors == varGroupings[iterationParity].size()

	// std::cout << "initialization done" << std::endl;

	iterationParity = (iterationParity+1)%2;
	bool groupsChanged = true;
	while( (groupsChanged || iteration==1) && ((nMaxIterations==-1) || (iteration<nMaxIterations)) ){ 
		std::cout << "starting iteration=" << iteration << " of" << nMaxIterations-1;
		std::cout << "\n\tnGroups=" << varGroupings[(iterationParity+1)%2].size() << std::endl;
		groupsChanged = false;
		// construct signature for each factor.
		for(size_t fid=0; fid<nFactors; ++fid) {
			facSignature[fid][0] = facColor[fid];
			int vid = -1;
			for(size_t vidid=0; vidid<gm.numberOfVariables(fid); ++vidid) {
				vid = gm.variableOfFactor(fid, vidid);
				facSignature[fid][1+vidid] = varColor[vid]; // +1 because of the 
			}
		}
		// std::cout << "got here" << std::endl;
		// recoloring factors.
		nFacColors = assignColorBySignature(facSignature, nVarColors, facColor);
		// construct  facGroups[newIdx] using colors
		// std::cout << "here too" << std::endl;
		constructGroupsFromColors(facColor, facGroupings[iterationParity]);
		// std::cout << "on iteration=" << iteration << "there were " << facGroupings[iterationParity].size() << " fac groups";
		// std::cout << "\n\t and " << nFacColors << "facColors"  << std::endl;
		// construct signature for each variable
		for(size_t vid=0; vid<nVariables; ++vid) {
			varSignature[vid][0] = varColor[vid];
			size_t fid = -1;
			for(size_t fidid=0; fidid<gm.numberOfFactors(vid); ++fidid) {
				fid = gm.factorOfVariable(vid, fidid);
				varSignature[vid][1+fidid] = facColor[fid];
			}
		}
		// recoloring variables
		// std::cout << "and then here" << std::endl;
		nVarColors = assignColorBySignature(varSignature, nFacColors, varColor);
		// construct varGroups[newIdx] using colors
		// std::cout << "here at the end" << std::endl;
		constructGroupsFromColors(varColor, varGroupings[iterationParity]);
		// std::cout << "one last thing to do" << std::endl;
		// detect if grouping has changed.
		groupsChanged = detectGroupsChanged(varGroupings[(iterationParity)%2], // old group
						varGroupings[iterationParity]); // new group.
		// update control variables
		// std::cout << "on iteration=" << iteration << "there were " << varGroupings[iterationParity].size() << " groups" << std::endl;
		iteration++;
		iterationParity = (iterationParity+1)%2;
	}

	// std::cout << "mean,min,max FactorDiff= " << totalDiff << "," << minFactorDiff << "," << maxFactorDiff << std::endl;
	return varGroupings[(iterationParity+1)%2]; // return the latest grouping.
}

/*
	getVarGroupings goes from group to list of variables. (syms)
	However, we also need something that goes from variables to group. (invsyms)
	This function produces invsyms given syms
	
	@params
		IN syms - fully populated
		OUT invsysms - initialized, and has space reserved.
	@return
		Nothing
*/
void getVar2GroupID(
	// arguments
		std::vector< std::vector<size_t> >& syms,
		std::vector<size_t>& invsysms
	) { // function begins

	size_t nObj = invsysms.size();
	size_t nGroups = syms.size();
	size_t nGroupSize = -1;
	size_t vid = -1;
	for(size_t gid=0; gid<nGroups; ++gid) {
		nGroupSize = syms[gid].size();
		for(size_t vidid=0; vidid<nGroupSize; ++vidid) {
			vid = syms[gid][vidid];
			invsysms[vid] = gid;
		}
	}
}

#endif

#endif
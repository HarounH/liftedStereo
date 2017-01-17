#ifndef _COMPRESS_FACTOR_GRAPH_
#define _COMPRESS_FACTOR_GRAPH_ 

#include <iostream>
#include <opengm/opengm.hxx>
#include <opengm/graphicalmodel/graphicalmodel.hxx>

#include <opengm/operations/minimizer.hxx>
#include <opengm/operations/adder.hxx>

#include <opengm/graphicalmodel/space/simplediscretespace.hxx>
#include <opengm/functions/potts.hxx>

#include <vector>  // for vector
#include <utility> // for pair
#include <map>


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
	TODO - test.
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


/* 	TODO: Test.
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
	TODO: test
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



/* TODO: test this.
	This function tells you if two factors in a graphical model have their label wise difference less than a threshold
	@params
		gm - the graph itself
		fid1 - factor idx 1
		fid2 - factor idx 2
		threshold - the acceptable threshold for equality... compared against l1 norm difference.
	@returns
		true if within threshold difference
		false otherwise
*/
template<class T, class OPERATOR, class FUNCTION_TYPE_LIST, class SPACE>
bool L1factorMatch(
	//arguments
		const opengm::GraphicalModel<T, OPERATOR, FUNCTION_TYPE_LIST, SPACE>& gm, 
		int fid1, int fid2,
		double threshold
	) {
	// std::cout << "start L1factorMatch" << std::endl;
	double difference = 0.0;
	
	// measure difference... l1 norm difference will do :P
	if (gm[fid1].numberOfVariables() != gm[fid2].numberOfVariables()) {
		return false;
	}
	size_t nDim = gm[fid1].numberOfVariables();
	// std::cout << "nDim=" << nDim << std::endl;
	std::vector<size_t> nLabels(nDim);
	
	size_t* iterator = new size_t[nDim]; // used to iterate over factors.
	
	for(size_t i=0; i<nDim; ++i) {
		if(gm[fid1].numberOfLabels(i)!=gm[fid2].numberOfLabels(i)) {
			return false;
		}
		nLabels[i] = gm[fid1].numberOfLabels(i);
		iterator[i] = 0;
	}

	// funky iteration method.
	bool done = false;
	int lol_iter = 0;
	while(!done) {
		// std::cout << "lol_iter=" << lol_iter++ << "abs diff=" << abs( gm[fid1].template operator()<size_t*>(iterator) - gm[fid2].template operator()<size_t*>(iterator) ) << std::endl;
		// update difference.
		// for(int i=0; i<nDim; ++i) { std::cout<<iterator[i] << " "; } std::cout << std::endl;
		double x = (double)gm[fid1].template operator()<size_t*>(iterator) - (double)gm[fid2].template operator()<size_t*>(iterator);
		// std::cout << (double)gm[fid1].template operator()<size_t*>(iterator) << " ";
		// std::cout << (double)gm[fid2].template operator()<size_t*>(iterator)  << " " ;
		// std::cout << x << std::endl;
		difference += x>0?x:-x;
		
		// update iterator.
		int carry = 1;
		
		for(int dim=0; ((dim<nDim) && (carry!=0)); ++dim) {
			// std::cout << "carry=" << carry << "iterator[dim],nLabels[dim]=" << iterator[dim] << "," << nLabels[dim] << std::endl;
			size_t lol = iterator[dim]+carry; // temporary variable.
			iterator[dim] = (lol)%nLabels[dim];
			carry = (lol)/nLabels[dim];
		}

		// check if done.
		if(carry!=0) {
			done = true;
		}
	}
	// std::cout << "end L1factorMatch" << std::endl;
	delete[] iterator;
	if (difference>maxFactorDiff) maxFactorDiff = difference;
	if(minFactorDiff < 0.0 || difference<minFactorDiff) minFactorDiff = difference;
	totalDiff = ((nDiff/(nDiff+1))*totalDiff + (difference/(nDiff+1)));
	// std::cout << "diff/threshold=" << (double)difference << "/" << (double)threshold << std::endl;
	return (difference<threshold);
}

/*	TODO test this
	This function initializes factors and gives them colors if their 
		function difference is "acceptable"
	@params
		IN gm
		OUT facColor - an assignment of colors to each factor in gm.
			the output should be initiatelized to -1 for each factor.
	@returns
		number of colors used to initialize the coloring.
	
	depends on the hash-define of function-threshold.

*/
template<class T, class OPERATOR, class FUNCTION_TYPE_LIST, class SPACE>
size_t initializeFactorColors(
	// arguments
		const opengm::GraphicalModel<T, OPERATOR, FUNCTION_TYPE_LIST, SPACE>& gm, 
		std::vector<int>& facColor,
		double funcDiffThreshold = DEFAULT_FUNCTION_L1_THRESHOLD // completely arbitrary number.
	) { // function begin
	size_t nColorsUsed = 0;
	// use a hashmap from factor's to size_t.
	std::vector<int> functionToFactor(gm.numberOfFunctions(0)+gm.numberOfFunctions(1), -1); // maps each function to the first factor that uses it.
	// std::cout << "iniFacColors START" /*".. ftf.size()=" << functionToFactor.size()*/ << std::endl;
	int funcId = -1;
	for(int fid=0; fid<facColor.size(); ++fid) {
		// if(fid == facColor.size()%10000) {
		// 	std::cout << "done with " << fid << "/" << facColor.size() << std::endl;
		// }
		// std::cout << "initializeFactors: iteration=" << fid << "/" << facColor.size() << "\r" << std::flush;
		funcId = gm[fid].functionIndex();
		// std::cout << "funcId=" << funcId  << std::endl;
		if (functionToFactor[funcId] == -1) { // first time that it is being seen.
			// see if it matches any earlier factor... 
			// std::cout << "getting in for fid=" << fid << std::endl;
			
			for(int fid2=((fid-500)>0?(fid-500):(0)); fid2<fid; ++fid2) { // FIXME make this more exhaustive.
				// std::cout << "got in here!" << std::endl;
				// does it match?
				// if yes - set the function's color, and the facColor[fid] aptly.
				if( L1factorMatch<T, OPERATOR, FUNCTION_TYPE_LIST, SPACE>(gm, fid, fid2, funcDiffThreshold) ) {
					// std::cout << "matched fid=" << fid << " with fid2=" << fid2 << std::endl;
					functionToFactor[funcId] = fid2; // doesn't exactly matter what happens here :P
					facColor[fid] = facColor[fid2];
					break;
				}
			}

			if(functionToFactor[funcId]==-1) { // nothing matched.
				// std::cout << "making new color for factor with id=" << fid << std::endl;
				functionToFactor[funcId] = fid; // first time being used.
				facColor[fid] = nColorsUsed;
				nColorsUsed++;
			}
		} else { // same function has been seen before.
			facColor[fid] = facColor[ functionToFactor[funcId] ]; // doesn't use a new color.
		}
	}
	// std::cout << "iniFacColors LAST " << std::endl;
	return nColorsUsed;

}
// symmetry 1 : compressFactorGraph from counting bp ... gives us 
template<class T, class OPERATOR, class FUNCTION_TYPE_LIST, class SPACE>
std::vector< std::vector<size_t> > getVarGroupings(
		// arguments
			const opengm::GraphicalModel<T, OPERATOR, FUNCTION_TYPE_LIST, SPACE>& gm, 
			int nMaxIterations=-1,
			double funcL1Epsilon=DEFAULT_FUNCTION_L1_THRESHOLD
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
	// std::cout << "about to initialize factor colors " << std::endl;
	size_t nFacColors = initializeFactorColors<T, OPERATOR, FUNCTION_TYPE_LIST, SPACE>(gm, facColor, funcL1Epsilon);
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

	// BRK
	// std::cout << "initialization done" << std::endl;

	iterationParity = (iterationParity+1)%2;
	bool groupsChanged = true;
	while( (groupsChanged || iteration==1) && ((nMaxIterations==-1) || (iteration<nMaxIterations)) ){ // until stop:
		std::cout << "starting iteration=" << iteration << " of " << (nMaxIterations-1);
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
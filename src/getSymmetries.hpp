#ifndef _GET_SYMMETRIES_HPP_
#define _GET_SYMMETRIES_HPP_


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

struct pair_sort_by_first_predicate { // explanatory first name.
    bool operator()(const std::pair<double,int> &left, const std::pair<double,int> &right) {
        return left.first < right.first;
    }
};

/*
	This function takes in a graphical model, and returns colors for each 
	factor in the graph.

	For practicality, we know that there are only unary and pairwise factors.
	We also know that there are three kinds of pairwise factors.
	
	Two unaries get the same color if their nRanksConsidered most likely states are the same,
	in the same order.

*/
template<class T, class OPERATOR, class FUNCTION_TYPE_LIST, class SPACE>
size_t mostLikelyFactorColorerForTSGO(
	// arguments
		const opengm::GraphicalModel<T, OPERATOR, FUNCTION_TYPE_LIST, SPACE>& gm, 
		int nRanksConsidered, // this argument should get more special as time goes on.
		std::vector<int>& colors
	) { // function begins 

	// two variables are symmetric if their unaries have similar MAP thingies.
	int nVar = gm.numberOfVariables();
	int nFac = gm.numberOfFactors();
	int nLabels= gm.numberOfLabels(0);	
	
	size_t nColorsUsed = 0;	
	std::map< std::vector<int> , int > mostLikely2Color;
	std::map< double, int > pairwiseParameter2Color;
	
	for(int fid=0; fid<nFac; ++fid) {
		// The factor is either unary or pairwise.
		if( gm[fid].dimension()==1 ) {

		} else if (gm[fid].dimension()==2) {

		} else {
			
		}
	}
	return nColorsUsed;
}

#endif
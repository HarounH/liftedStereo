#ifndef _USE_SYM_HPP_ // header guards!
#define _USE_SYM_HPP_
/**
	This file implements AlphaExpansionWithSymmetries as follows:

	It depends on the AlphaExpansion (and AlphaExpansionWithFusion) of OpenGM.
	However, given a GM and the set of symmetry groups in the GM (in the form of a vector of vector of variables (ints) )

	It really just has a function called reduceGraphUsingSymmetries(graphModel, symmetries)
*/
#include <iostream>
#include <opengm/opengm.hxx>
#include <opengm/graphicalmodel/graphicalmodel.hxx>

#include <opengm/operations/minimizer.hxx>
#include <opengm/operations/adder.hxx>

#include <opengm/graphicalmodel/space/simplediscretespace.hxx>
#include <opengm/functions/potts.hxx>

#include <vector>
/**
	This function takes a graphical model, symmetries (groups of symmetric variables) and invsyms ( which tells you which group each variable is in. )
	@params graph - a 2d grid pixel labelling graph.
	@params syms - a vector<vector<size_t>>. each vector<size_t> is a group, containing variable indices that are symmetric (hard constraint). 

	PS: This function has a lot of restrictions on the functions it can accept etc.

	WARNING: It doesn't reduce the number of edges. it is inefficient.
*/

template<class T, class OPERATOR, class FUNCTION_TYPE_LIST, class SPACE>
opengm::GraphicalModel<T, OPERATOR, FUNCTION_TYPE_LIST, SPACE> reduceGraphUsingSymmetries(
		// arguments
		const opengm::GraphicalModel<T, OPERATOR, FUNCTION_TYPE_LIST, SPACE>& graph, 
		const std::vector<std::vector<size_t>>& syms, // size=nGroups.
		const std::vector<size_t>& invsyms, // groups are labelled 0,1...nGroups 
		const double truncatedParameter1,
		const double truncatedParameter2,
		const double pottsnParameter
	) { // function starts
	size_t nGroups = syms.size();
	size_t nVariables = graph.numberOfVariables();
	size_t nLabels = graph.numberOfLabels(0); // FIXME: Number of labels of 1st variable... should be number of variables, period.
	int maxOrderPotts = graph.factorOrder();
	typedef typename opengm::GraphicalModel<T, OPERATOR, FUNCTION_TYPE_LIST, SPACE>::FunctionIdentifier funcid_type;
	
	std::vector< funcid_type > pottsIdentifiers(1+maxOrderPotts);

	opengm::GraphicalModel<T, OPERATOR, FUNCTION_TYPE_LIST, SPACE> newGraph(opengm::SimpleDiscreteSpace<size_t, size_t>(nGroups, nLabels));
	
	for(int i=1; i<=maxOrderPotts; ++i) {
		std::vector<size_t> shapeVec(i, graph.numberOfLabels(0));
		opengm::PottsNFunction<double> pn(shapeVec.begin(), shapeVec.end(), 0.0, pottsnParameter); // FIXME 50.0 is hardcoded yo
		pottsIdentifiers[i] = newGraph.addFunction(pn);
	}
	// we need to create a new graphical model
	// it has nVariables = # of groups
	// each factor is now between groups
	// 		the new factors must be constructed from factors of the graph.
	

	// step1: declare variable for new graph.

	// Step2: fix the factors 
	std::vector<bool> factorAccounted(graph.numberOfFactors(), false);
	size_t nFactorsOrig = graph.numberOfFactors();
	size_t unaryShape[] = {nLabels};
	std::vector< opengm::ExplicitFunction<double> > unaries( nGroups, opengm::ExplicitFunction<double>(unaryShape, unaryShape+1) );
	// FIXME: generalize the factors that are usable.
	// in graph, we have two kinds of factors - Unary and Potts.
	for(size_t fid=0; fid<nFactorsOrig; ++fid) {
		auto fid_factor = graph[fid];
		size_t dim = fid_factor.dimension();
		// std::cout << "dealing with fid=" << fid  << "/" << nFactorsOrig << "\r" << std::flush;
		if ( dim==1 ) { // if factor is unary, it remains unary.
			factorAccounted[fid] = true;
			size_t vid = graph.variableOfFactor(fid, 0);
			size_t gid = invsyms[vid];
			for(size_t l=0; l<nLabels; ++l) {
				size_t labels[] = {l};
				auto addValue = (fid_factor).template operator()<size_t*>(labels);
				auto presentValue = unaries[gid].template operator()<size_t*>(labels);
				unaries[gid](l) = presentValue + addValue;
			}
		} else if ( dim==2 ) { // potts function.
			size_t vid1 = graph.variableOfFactor(fid, 0); 
			size_t vid2 = graph.variableOfFactor(fid, 1);
			size_t gid1 = invsyms[vid1];
			size_t gid2 = invsyms[vid2];
			if (gid1==gid2) { // both variables are in the same group
				// modify unaries[gid] as we wish :P
				for(size_t l=0; l<nLabels; ++l) {
					size_t labels[] = {l,l};
					size_t singleLabel[] = {l};
					auto presentValue = unaries[gid1].template operator()<size_t*>(singleLabel);
					auto addValue = (fid_factor).template operator()<size_t*>(labels);
					unaries[gid1](l) = presentValue + addValue; 
				}
				factorAccounted[fid] = true;
			} else { // add it as a pairwise function.
				// simple it add it as a factor to the graph with correct variables etc.
				size_t vars[] = {gid1, gid2};
				std::sort(vars, vars+2);
				// typename opengm::GraphicalModel<T, OPERATOR, FUNCTION_TYPE_LIST, SPACE>::FunctionIdentifier 
				auto func = graph[fid].template function<(size_t)1>(); // FIXME ... need to generalize this.
				auto nfid = newGraph.addSharedFunction(func);
				newGraph.addFactor(nfid, vars, vars+2);
				factorAccounted[fid] = true;
			}
		} else { // higher order functions, unfortunately
			int nVids = graph.numberOfVariables(fid);
			std::vector<size_t> gids;
			size_t vid,gid;
			for(int vidid=0; vidid<nVids; ++vidid) {
				vid = graph.variableOfFactor(fid, vidid);
				gid = invsyms[vid];
				if ( find(gids.begin(), gids.end(), gid) == gids.end()) {
					gids.push_back(gid);
				}
			}
			if(gids.size()<2) {
				// Do Nothing, essentially ... wrong :(
			} else {
				sort(gids.begin(), gids.end());
				auto fid = pottsIdentifiers[gids.size()];
				newGraph.addFactor(fid, gids.begin(), gids.end());	
			}
		}
	}
	// Step3: add the factors to the graphical model.

	// add unaries
	for(size_t gid=0; gid< nGroups; ++gid) {
		size_t vars[] = {gid};
		auto fid = newGraph.addFunction(unaries[gid]);
		newGraph.addFactor(fid, vars, vars+1);
	}
	return newGraph; // Returns the new graph
}

template<class T, class OPERATOR, class FUNCTION_TYPE_LIST, class SPACE>
opengm::GraphicalModel<T, OPERATOR, FUNCTION_TYPE_LIST, SPACE> reduceGraphUsingSymmetriesEfficienctly(
		// arguments
		const opengm::GraphicalModel<T, OPERATOR, FUNCTION_TYPE_LIST, SPACE>& graph, 
		const std::vector<std::vector<size_t>>& syms, // size=nGroups.
		const std::vector<size_t>& invsyms, // groups are labelled 0,1...nGroups 
		const double truncatedParameter1,
		const double truncatedParameter2,
		const double pottsnParameter
	) { // function starts
	size_t nGroups = syms.size();
	size_t nVariables = graph.numberOfVariables();
	size_t nLabels = graph.numberOfLabels(0); // FIXME: Number of labels of 1st variable... should be number of variables, period.
	int maxOrderPotts = graph.factorOrder();
	typedef typename opengm::GraphicalModel<T, OPERATOR, FUNCTION_TYPE_LIST, SPACE>::FunctionIdentifier funcid_type;

	// Need to capture the unaries
	size_t unaryShape[] = {nLabels};
	std::vector< opengm::ExplicitFunction<double> > unaries( nGroups, opengm::ExplicitFunction<double>(unaryShape, unaryShape+1) );

	// Need to capture the number of pairwise per pair, since all pairwise are identitical
	std::map< std::vector<size_t> , int > pair2nTruncated;
	int nMaxPairs=-1;
	// Need to capture the number of higher order per clique since all clique potentials are identitical
	std::map< std::vector<size_t> , int > clique2nPottsn;
	int nMaxCliques=-1; // counts the maximum number of times a pottsn edge appears on a clique
	// newGraph is the reduced graph.
	opengm::GraphicalModel<T, OPERATOR, FUNCTION_TYPE_LIST, SPACE> newGraph(opengm::SimpleDiscreteSpace<size_t, size_t>(nGroups, nLabels));
	// Capture all factors appropriately.
	size_t nFactorsOrig = graph.numberOfFactors();
	for(size_t fid=0; fid<nFactorsOrig; ++fid) {
		auto fid_factor = graph[fid];
		size_t dim = fid_factor.dimension();

		if (dim==1) { // definitely a unary.
			size_t vid = graph.variableOfFactor(fid, 0);
			size_t gid = invsyms[vid];
			for(size_t l=0; l<nLabels; ++l) {
				size_t labels[] = {l};
				auto addValue = (fid_factor).template operator()<size_t*>(labels);
				auto presentValue = unaries[gid].template operator()<size_t*>(labels);
				unaries[gid](l) = presentValue + addValue;
			}
		} else if (dim==2) { // a truncatedDifference ... either becomes a unary or remains a unary.
			size_t vid1 = graph.variableOfFactor(fid, 0); 
			size_t vid2 = graph.variableOfFactor(fid, 1);
			size_t gid1 = invsyms[vid1];
			size_t gid2 = invsyms[vid2];
			if (gid1==gid2) { // both variables are in the same group
				// modify unaries[gid] as we wish :P
				for(size_t l=0; l<nLabels; ++l) {
					size_t labels[] = {l,l};
					size_t singleLabel[] = {l};
					auto presentValue = unaries[gid1].template operator()<size_t*>(singleLabel);
					auto addValue = (fid_factor).template operator()<size_t*>(labels);
					unaries[gid1](l) = presentValue + addValue; 
				}
			} else { // add it as a pairwise function.
				// simple it add it as a factor to the graph with correct variables etc.
				std::vector<size_t> vars = {gid1, gid2};
				std::sort(vars.begin(), vars.end());

				if (pair2nTruncated.find(vars)==pair2nTruncated.end()) { // not found.
					pair2nTruncated[vars] = 1;
				} else {
					pair2nTruncated[vars] += 1;
				}
				if (pair2nTruncated[vars]>nMaxPairs) nMaxPairs = pair2nTruncated[vars];
			}

		} else { // a pottsNfunction ... it will either become a unary, or remain a pottsN
			int nVids = graph.numberOfVariables(fid);
			std::vector<size_t> gids;
			size_t vid,gid;
			// construct the variables over which the new function will be (based on gids)
			for(int vidid=0; vidid<nVids; ++vidid) {
				vid = graph.variableOfFactor(fid, vidid);
				gid = invsyms[vid];
				if ( find(gids.begin(), gids.end(), gid) == gids.end()) {
					gids.push_back(gid);
				}
			}

			if(gids.size()==1) {
				for(size_t l=0; l<nLabels; ++l) {
					std::vector<size_t> labels(nVids, l);
					size_t singleLabel[] = {l};
					auto addValue = (fid_factor).template operator()< std::vector<size_t>::iterator >(labels.begin());
					auto presentValue = unaries[gid].template operator()<size_t*>(singleLabel);
					unaries[gid](l) = presentValue + addValue;
				}
			} else { // TODO
				std::sort(gids.begin(), gids.end());
				if ( clique2nPottsn.find(gids)==clique2nPottsn.end() ) {
					clique2nPottsn[gids] = 1;
				} else {
					clique2nPottsn[gids] += 1;
				}
				if (clique2nPottsn[gids]>nMaxCliques) nMaxCliques = clique2nPottsn[gids];
			}
		}
	} // handled all the factors.

	// addFactors

	// 1. unaries
	for(size_t i=0; i<unaries.size(); ++i) {
		size_t vars[]={i};
		auto fid = newGraph.addFunction(unaries[i]);
		newGraph.addFactor(fid, vars, vars+1);
	}
	// 2. all pair2nTruncated
	std::vector< funcid_type > truncated(nMaxPairs+1);
	std::vector<bool> truncatedFunctionCreated(nMaxPairs+1, false);
	typedef std::map< std::vector<size_t> , int >::iterator mapIterator;
	for(mapIterator iter=pair2nTruncated.begin(); iter!= pair2nTruncated.end(); iter++) {
		// Get the clique, the variables are the key.
		auto key = iter->first;
		int val = iter->second;
		// Figure out the function
		if (!truncatedFunctionCreated[val]) { // need to create the funcion?

			opengm::TruncatedAbsoluteDifferenceFunction<double> f(nLabels, nLabels, truncatedParameter1, val*truncatedParameter2);
			truncated[val] = newGraph.addFunction(f);
			truncatedFunctionCreated[val] = true;
		}
		// Add the function as a factor.
		newGraph.addFactor(truncated[val], key.begin(), key.end());
	}
	// 3. all clique2nPottsn
	std::vector< std::vector< funcid_type > > pottsnIdentifiers(maxOrderPotts+1, std::vector<funcid_type>(nMaxCliques+1) );
	std::vector< std::vector<bool> > pottsnCreated(maxOrderPotts+1, std::vector<bool>(nMaxCliques+1, false) );

	for(mapIterator iter=clique2nPottsn.begin(); iter!=clique2nPottsn.end(); iter++) {
		auto key = iter->first;
		int val = iter->second;
		if(pottsnCreated[key.size()][val]) {
			std::vector<size_t> shape(key.size(), nLabels);
			opengm::PottsNFunction<double> f(shape.begin(), shape.end(), 0.0, val*pottsnParameter);
			pottsnIdentifiers[key.size()][val] = newGraph.addFunction(f);
			pottsnCreated[key.size()][val] = true;
		}
		newGraph.addFactor(pottsnIdentifiers[key.size()][val], key.begin(), key.end());
	}
	std::cout << " no error yet 4\n"; // BRK 4
	return newGraph;
}
#endif
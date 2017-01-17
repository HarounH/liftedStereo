#ifndef _SYMMETRIES_2_HPP_
#define _SYMMETRIES_2_HPP_

struct sort_pred2 {
    bool operator()(const std::pair<double,int> &left, const std::pair<double,int> &right) {
        return left.first < right.first;
    }
};
void initializePairwiseColors(std::vector<int>& factorColors, png::image<png::rgb_pixel>& img) {
	int nx = img.get_width();
	int ny = img.get_height();
	double deltaf;
	int nv = nx*ny;
	for(int y=0; y<ny; ++y) {
		for(int x=0; x<nx; ++x) {
			// Rightward factor
			int vid=x + y*nx;
			if( x<(nx-1) ) {
				deltaf = absolute_color_difference(img[y][x], img[y][x+1]);
				if (deltaf<u1) {
					factorColors[nv + vid] = 0;
				} else if (deltaf < u2) {
					factorColors[nv + vid] = 1;
				} else {
					factorColors[nv + vid] = 2;
				}
			}
			// Leftward factor
			if( y<(ny-1) ) {
				deltaf = absolute_color_difference(img[y][x], img[y+1][x]);
				if (deltaf<u1) {
					factorColors[2*nv + vid] = 0;
				} else if (deltaf < u2) {
					factorColors[2*nv + vid] = 1;
				} else {
					factorColors[2*nv + vid] = 2;
				}
			}
		}
	}
}

// Returns number of colors used.
int initializeUnaryColors(std::vector<int>& factorColors, double* mcost, 
		int nx, int ny, int nLabels, int nRanksConsidered
	) {
	int nColorsUsed = 3; // 0,1, and 2 have already been used for pairwise factors.
	std::vector<double> factor(nLabels, 0.0);
	int nVariables = nx*ny;
	boost::unordered_map< std::vector<int> , int > rank2color;
	boost::unordered_map< std::vector<int> , int >::const_iterator findRanks;
	std::vector< std::pair<double, int> > potentialAndIdx(nLabels);
	std::vector<int> ranks(nRanksConsidered);
	int v;
	for(int y=0; y<ny; ++y) {
		for(int x=0; x<nx; ++x) {
			v = y*nx + x;
			for(int l=0; l<nLabels; ++l) {
				factor[l] = mcost[y*nx + x  + l*nVariables];
				potentialAndIdx[l] = std::pair<double, int>(factor[l],l);
			}
			// Get the k most likely states
			std::sort(potentialAndIdx.begin(), potentialAndIdx.end(),sort_pred2());
			for(int k=0; k<nRanksConsidered; ++k) {
				ranks[k] = potentialAndIdx[k].second;
			}
			// Now to get the color.
			findRanks = rank2color.find(ranks);
			if ( findRanks != rank2color.end() ) {
				factorColors[v] = findRanks->second;
			} else {
				rank2color[ranks] = nColorsUsed;
				factorColors[v] = nColorsUsed;
				nColorsUsed++;
			}
		}
	}
	std::cout << "Used " << (nColorsUsed-3) << " colors." << std::endl;
	return nColorsUsed;
}

void getSymmetries(
		double* mcost, png::image< png::rgb_pixel >& img, 
		int nx, int ny, int nLabels, 
		int nRanksConsidered, int nCBPIterations, 
		std::vector<std::vector<int> >& syms, std::vector<int>& invsyms
	) {
	// Step0 : Get constants
	int nVariables = nx*ny;
	int nFactors = 3*nVariables; // Some of these aren't used.
	
	// Step0.5 : Clear outputs
	invsyms.resize(nVariables, -1); // None of them is in any group.
	syms.clear();

	// Step1 : Create factor graph representation
	std::vector< std::vector<int> > variableAdjList(nVariables);
	std::vector< std::vector<int> > variableSignatures(nVariables);
	std::vector<int> variableColors(nVariables, 0);
	std::vector<std::vector<std::vector<int> > > variableGroups(2);

	// First nVariables variables are unary. Then we have rightward factors.
	// Then we have downward factors.
	std::vector< std::vector<int> > factorAdjList(nFactors);
	std::vector< std::vector<int> > factorSignatures(nFactors);
	std::vector<int> factorColors(nFactors, -1);
	std::vector<std::vector<std::vector<int> > > factorGroups(2);

	// Initializing adjList.
	for(int y=0; y<ny; ++y) {
		for(int x=0; x<nx; ++x) {
			int vid = x + y*nx;
			variableAdjList[vid].push_back(vid); // unary
			factorAdjList[vid].push_back(vid);
			if (x<(nx-1))  {
				variableAdjList[vid].push_back(nVariables + vid); // rightward
				factorAdjList[nVariables + vid].push_back(vid);
			}
			if (y<(ny-1)) {
				variableAdjList[vid].push_back(2*nVariables + vid); // downward
				factorAdjList[2*nVariables + vid].push_back(vid);
			}
			variableSignatures[vid] = std::vector<int>(1+variableAdjList[vid].size(), 0);
		}
	}
	for(int f=0; f<nFactors; ++f)  {
		if (factorAdjList[f].size()>0) {
			factorSignatures[f] = std::vector<int>(1 + factorAdjList[f].size(),-1);
		} else {
			factorSignatures[f] = std::vector<int>(1,-1);
		}
	}

	// Step2 : Initialize colors
	initializePairwiseColors(factorColors, img); // Uses three colors.
	int nFactorColors = initializeUnaryColors(factorColors, mcost, nx, ny, nLabels, nRanksConsidered);
	std::cout << "I got here 0 :D" << std::endl;
	for(int i=0; i<factorColors.size(); ++i) {
		if (factorAdjList[i].size()>0) { // factor exists
			factorSignatures[i][0] = factorColors[i];
		} else {
			factorSignatures[i][0] = -1;
		}
	}
	int newid = 0;
	std::cout << "I got here 1 :D" << std::endl;
	constructGroupsFromColors(factorColors, factorGroups[newid]);
	std::cout << "nFactorColors:" << nFactorColors << " factorGroups.size(): " << factorGroups[newid].size() << std::endl;
	
	// Need to initialize variable colors too.
	for(int v=0; v<nVariables; ++v) {
		variableSignatures[v][0] = variableSignatures[v][1] = factorColors[v];
	}
	std::cout << "I got here 2 :D" << std::endl;
	int nVariableColors = assignColorBySignature(variableSignatures, variableColors);
	std::cout << "I got here 3 :D" << std::endl;
	constructGroupsFromColors(variableColors, variableGroups[newid]);
	int cbp_iter=0;
	
	// Step3 : CBP
	while(cbp_iter<nCBPIterations) {
		newid = (newid+1)%2;
		std::cout << "Counting BP iteration #" << cbp_iter << ", " << variableGroups[ (newid+1) % 2].size() << " groups now" << std::endl;

		// Messages from variables to factors
		for(int fid=0; fid<nFactors; ++fid) {
			factorSignatures[fid][0] = factorColors[fid];
			int vid = -1;
			for(int vidid=0; vidid<factorAdjList[fid].size(); ++vidid) {
				vid = factorAdjList[fid][vidid];
				factorSignatures[fid][1+vidid] = variableColors[vid]; // +1 because of the 
			}
		}

		assignColorBySignature(factorSignatures, factorColors);// assignColors by signatures
		constructGroupsFromColors(factorColors, factorGroups[newid]); // construct groups

		// Messages from factors to variables
		for(int vid=0; vid<nVariables; ++vid) {
			variableSignatures[vid][0] = variableColors[vid];
			int fid = -1;
			for(int fidid=0; fidid<variableAdjList[vid].size(); ++fidid) {
				fid = variableAdjList[vid][fidid];
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
#endif
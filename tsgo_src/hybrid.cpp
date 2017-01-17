#ifndef _HYBRID_CPP_
#define _HYBRID_CPP_

double delta_threshold = 2.0;

// This file contains functions that support our coarse-to-fine algorithm
bool detectPlateau1(int iter, std::vector< std::vector<double> >& energies) {
	if (iter<energies[0].size()) return false; // Give the algorithm a chance man.
	int n = energies[0].size();
	// assert(n>2);
	// Calculate difference.
	std::vector<double> deltas(2,0.0);
	double temp;
	for(int lr=0; lr<2; ++lr) {
		deltas[lr] = energies[lr][(iter+1)%n] - energies[lr][(iter)%n];
		// If deltas[lr] is too huge, ditch
		if (deltas[lr]>delta_threshold) return false;
		// Else if its small then check the next one.
		for(int i=1; i<(n-1); ++i) {
			temp = energies[lr][(iter+i+1)%n] - energies[lr][(iter+i)%n];
			if ( (abs(temp-deltas[lr]) > delta_threshold) || (temp>delta_threshold) ) {
				return false;
			}
			deltas[lr] = temp;
		}
	}
	return true;
}

bool detectPlateau2(int iter, std::vector< std::vector<double> >& energies) {
	if (iter<energies[0].size()) return false; // Give the algorithm a chance man.
	int n = energies[0].size();
	// assert(n>2);
	// Calculate difference.
	std::vector<double> deltas(n-1,0.0);
	for(int lr=0; lr<2; ++lr) {
		for(int i=(deltas.size()-1); i>=0; i--) {
			// deltas[i] = energies[lr][ ] - energies[lr][ ]; // TODO
		}
		// If any of them change, avoid a switch.
	}
	return true;
}

void updateLabelsWithChangedSymmetries(
		std::vector<int>& srcLabels, 
		std::vector<std::vector<int> >& srcSyms,
		std::vector<int>& dstLabels,
		std::vector<int>& dstInvSyms
	) {
	for(int i=0; i<srcSyms.size(); ++i) { // for each group
		for(int vidx=0; vidx<srcSyms[i].size(); ++vidx) {
			// read this by writing it out. sorry.
			dstLabels[ dstInvSyms[ srcSyms[i][vidx] ] ] = srcLabels[i];
		}
	}
}

#endif
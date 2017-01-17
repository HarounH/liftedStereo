#ifndef _POST_PROCESS_HPP_
#define _POST_PROCESS_HPP_

void left_to_right_cross_checking(
	// arguments
		std::vector<size_t>& llabels,
		std::vector<size_t>& rlabels,
		std::vector<bool>& uncertain,
		int nx, int ny,
		int parity=1 // used to flip left-right to right-left.
	) {
	auto variableIndex = [&nx](const size_t x, const size_t y) {return x + nx*y;};
	int ll, rl;
	for(int y=0; y<ny; ++y) {
		for(int x=0; x<ny; ++x) {
			ll = parity*llabels[variableIndex(x,y)];
			if ((x-ll)>0 && ((x-ll) < nx)) {
				rl = rlabels[variableIndex(x-ll, y)];
				uncertain[variableIndex(x,y)] = (parity*ll != rl);
				llabels[variableIndex(x,y)] = std::min(parity*ll, rl);
			}
		}
	}
}



// Replace each pixel by the median of gaussian-bilateral weighted pixels' labels
void gaussian_median_filter(
	// arguments
		png::image< png::rgb_pixel >& img, 
		std::vector<size_t>& labels,
		std::vector<bool>& uncertain,
		size_t nLabels,
		int radius,
		double sqsigmax,
		double sqsigmaf
	) {
	std::vector<size_t> olabels(labels);
	size_t nx = img.get_width();
	size_t ny = img.get_height();
	auto variableIndex = [&nx](const size_t x, const size_t y) {return x + nx*y;};
	auto vx = [&nx](const size_t vid) { return vid%nx; };
	auto vy = [&nx](const size_t vid) { return (int)(vid/nx); };
	
	// Things we need for a gaussian filter.
	std::vector<double> gaussx(radius+1,0.0);
	for(int i=0; i<gaussx.size(); ++i) {
		gaussx[i] = exp(-0.5*(i*i)/sqsigmax);
	}
	std::vector<double> gaussf(255*3, 0.0);
	for(int i=0; i<gaussf.size(); ++i) {
		gaussf[i] = exp(-0.5*(i*i)/sqsigmaf);
	}

	int x, y, xo, yo, l, mi, mx; // mi, mx are min and max labels in the gaussian window.
	double w, sum, threshold;
	std::vector<double> histogram(nLabels, 0.0); // The weights for each label.
	
	for(int vid=0; vid<labels.size(); ++vid) {
		x = vx(vid);
		y = vy(vid);
		mi = mx = olabels[vid];
		// make a histogram.
		for(int i=0; i<histogram.size(); ++i) histogram[i] = 0.0; // clear out the histogram.

		sum=0.0;
		threshold=0.0; // used to figure out what label is the median.
		for(int i=-radius; i<=radius; ++i) {
			for(int j=-radius; j<=radius; ++j) {
				// Get label.
				xo = x+i;
				if (xo<0) xo = 0;
				if (xo>=nx) xo = nx-1;
				
				yo = y+i;
				if (yo<0) yo = 0;
				if (yo>=ny) yo = ny-1;

				l = olabels[variableIndex( xo, yo )];
				// if(uncertain[variableIndex(xo,yo)]) continue; // Don't bother with this!
				mi = (mi>l)?l:mi;
				mx = (mx<l)?l:mx;
				// Get gauss weight.
				w = gaussx[abs(i)]*gaussx[abs(j)]*gaussf[absolute_difference_color(img[y][x], img[yo][xo])];
				// Add weight to histogram
				histogram[l] += w;
				sum += w;
			}
		}
		threshold = 0.5*sum; // threshold to tell us that the median has passed us by.
		sum = 0.0;
		int median_idx = mi;
		for(int i=mi; i<=mx; ++i) {
			if (sum>= threshold) {
				median_idx = i;
				break;
			} else {
				sum += histogram[i];
			}
		}
		labels[vid] = median_idx;
	}
}

#endif
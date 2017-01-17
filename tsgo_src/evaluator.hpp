#ifndef _EVALUATOR_HPP_
#define _EVALUATOR_HPP_

double getBadPixelScore(int* dL, png::image< png::ga_pixel >& truthimg, int width, int height, int nLabels, double threshold) {
	int errors=0;
	for(int y=0; y<height; ++y) {
		for(int x=0; x<width; ++x) {
			// if error is greater threshold, it's an error.
			if (truthimg[y][x].value==0) continue;
			// std::cerr << dL[x+y*width] << " " << nLabels*double(truthimg[y][x].value)/255.0 << std::endl;
			if( abs((dL[x + y*width]) - (nLabels*double(truthimg[y][x].value)/255.0)) 
				> threshold 
				) errors++;
		}
	}
	std::cerr << "error pct. =" << double(errors)/(width*height) << std::endl;
	return double(errors)/(width*height);
}
#endif
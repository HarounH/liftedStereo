#ifndef _IMAGE_STATISTICS_HPP_
#define _IMAGE_STATISTICS_HPP_


inline double square_deltaf( png::rgb_pixel& p, png::rgb_pixel& q ) {
	return abs(p.red-q.red) + abs(p.green-q.green) + abs(p.blue-q.blue);
}

double getSquaredSigmaF(
	// arguments
		png::image< png::rgb_pixel >& img
	) {
	int nx = img.get_width();
	int ny = img.get_height();
	std::vector<double> means(3, 0.0);
	for(int y=0; y<ny; ++y) {
		for(int x=0; x<ny; ++x) {
			means[0] += img[y][x].red;
			means[1] += img[y][x].green;
			means[2] += img[y][x].blue;
		}
	}
	for(int i=0; i<means.size(); ++i) means[i] /= double(nx*ny);
	
	double squared_sigma = 0.0;
	double temp;
	for (int y = 0; y < ny; ++y) {
		for(int x = 0; x<nx; ++x) {
			temp = img[y][x].red - means[0];
			squared_sigma += (temp*temp);
			temp = img[y][x].green - means[1];
			squared_sigma += (temp*temp);
			temp = img[y][x].blue - means[2];
			squared_sigma += (temp*temp);
		}
	}
	return squared_sigma/(nx*ny);
}

double getSquaredSigmaOfLabels(
	// arguments
		png::image< png::rgb_pixel >& img,
		std::vector<size_t>& labels
	) {
	int ny = img.get_height();
	int nx = img.get_width();
	double mean = 0.0;
	for(int i=0; i<labels.size(); ++i) mean+= double(labels[i]);
	mean /= ny*nx;
	double var = 0.0;
	for(int i=0; i<labels.size(); ++i) var += ((mean-(double(labels[i])))*(mean-(double(labels[i]))));

	return (var/(ny)/(nx));
}
double getSquaredSigmaX(
	// arguments
		png::image< png::rgb_pixel >& img
	) {
	int nx = img.get_width();
	int ny = img.get_height();
	std::vector<double> means(2, 0.0);	// one for x, one for y.
	means[0] = double(nx)*0.5;
	means[1] = double(ny)*0.5;
	double squared_sigma = 0.0;
	double temp;
	for (int y = 0; y < ny; ++y) {
		for(int x = 0; x<nx; ++x) {
			temp = x - means[0];
			squared_sigma += temp*temp;
			temp = y - means[1];
			squared_sigma += temp*temp;
		}
	}
	return squared_sigma/(nx*nx*ny*ny);
}


#endif
#ifndef _UPPER_LAYER_HPP
#define _UPPER_LAYER_HPP

using namespace std;

/*
	W will be the matrix containing weights which are in bilateral filtering form.
	W is a 2d matrix, representing pixel*pixel.

	NOTE: this thing takes O(n**4) memory. Its bound to kill any machine.
*/
void constructWpq(
	// arguments
		png::image< png::rgb_pixel >& img,
		std::vector<std::vector<double> >& W
	) {
	int nx = img.get_width();
	int ny = img.get_height();
	int x1,x2, y1,y2;
	double squared_deltax, squared_deltaf;
	double sigmax2, sigmaf2;
		auto variableIndex = [&nx](const size_t x, const size_t y) {return x + nx*y;};
		auto vx = [&nx](const size_t vid) { return vid%nx; };
		auto vy = [&nx](const size_t vid) { return (int)(vid/nx); };
	for(int vid1=0; vid1<W.size(); ++vid1){
		W[vid1][vid1] = 1.0;
		for(int vid2=0; vid2<vid1; ++vid2) {		
			// change in position.
			x1 = vx(vid1);
			y1 = vy(vid1);
			x2 = vx(vid2);
			y2 = vy(vid2);
			squared_deltax = (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1);
			// change in colors
			squared_deltaf = square_deltaf(img[y1][x1], img[y2][x2]);

			W[vid1][vid2] = W[vid2][vid1] = exp( -0.5*squared_deltax*(1.0/sigmax2) -0.5*squared_deltaf*(1.0/sigmaf2) );
			
		}
	}
}

/*
	Modifies lcost and rcost accoring to the tsgo paper.
*/
void getUnariesAfterFirstLayer(
	// arguments
		png::image< png::rgb_pixel >& limg,
		png::image< png::rgb_pixel >& rimg,
		std::vector<std::vector<std::vector<double> > >& lcost,
		std::vector<std::vector<std::vector<double> > >& rcost,
		int nLabels
	) { // function begins

	int nx = limg.get_width();
	int ny = limg.get_height();
		auto variableIndex = [&nx](const size_t x, const size_t y) {return x + nx*y;};
		auto vx = [&nx](const size_t vid) { return vid%nx; };
		auto vy = [&nx](const size_t vid) { return (int)(vid/nx); };

	// Bad idea. The following 6 lines consume a lot of memory.
	// std::vector<std::vector<double> > matWlpq(nx*ny, std::vector<double>(nx*ny,0.0));
	// std::vector<std::vector<double> > matWrpq(nx*ny, std::vector<double>(nx*ny,0.0));
	// clock_t istart = clock_t();
	// constructWpq(limg, matWlpq);
	// constructWpq(rimg, matWrpq);
	// std::cout << "constructed Wpq in " << double(clock_t()-istart)/CLOCKS_PER_SEC << " s"<< endl;
	std::vector<std::vector<std::vector<double> > > newlcost(ny, std::vector<std::vector<double>>(nx, std::vector<double>(nLabels, 0.0)));
	std::vector<std::vector<std::vector<double> > > newrcost(ny, std::vector<std::vector<double>>(nx, std::vector<double>(nLabels, 0.0)));

	int vid1, vid2;
	double wl,wr;
	clock_t start = clock_t();
	double sqsigmax = getSquaredSigmaX(limg); // Should be same for rimg too.
	double sqsigmafl= getSquaredSigmaF(limg);
	double sqsigmafr= getSquaredSigmaF(rimg);
	std::cout << "Got squared sigmas in " << std::to_string(double(clock_t()-start)/CLOCKS_PER_SEC) << " s" << endl;
	double sqdeltax, sqdeltafl, sqdeltafr;
	for(int y1=0; y1<ny; ++y1) {
		for(int x1=0; x1<nx; ++x1) {
			vid1 = variableIndex(x1,y1);
			std::cout << vid1 << "/" << nx*ny << "\r";
			std::cout << std::flush;
			for(int y2=0; y2<ny; ++y2) {
				for(int x2=0; x2<nx; ++x2) {
					vid2 = variableIndex(x2,y2);
					sqdeltax = (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1);
					sqdeltafl= square_deltaf(limg[y1][x1], limg[y2][x2]);
					sqdeltafr= square_deltaf(rimg[y1][x1], rimg[y2][x2]);
					// add it in for each label.
					wl = exp(-0.5*( sqdeltax/sqsigmax + sqdeltafl/sqsigmafl));
					wr = exp(-0.5*( sqdeltax/sqsigmax + sqdeltafr/sqsigmafr));
					for(int l=0; l<nLabels; ++l) {
						newlcost[y1][x1][l] += wl*lcost[y2][x2][l];
						newrcost[y1][x1][l] += wr*rcost[y2][x2][l];
					}
				}
			}
		}
	}
	
	// Only one iteration of message passing is so weird, man.

	// TODO : transform the costs using the erf function thingy.
	lcost = newlcost;
	rcost = newrcost;
	return;
}


void doGuassianBilateralFilter(
	// arguments
		png::image< png::rgb_pixel >& img,
		std::vector<std::vector<std::vector<double> > >& cost,
		int radius, 
		double sqsigmax, double sqsigmaf
	) { // function begins

	int nx = img.get_width();
	int ny = img.get_height();
	int nLabels = cost[0][0].size();

	// make a copy.
	std::vector<std::vector<std::vector<double> > > ocost(cost);

	// initialize gaussian stuff.
	std::vector<double> gaussx(radius+1,0.0);
	for(int i=0; i<gaussx.size(); ++i) {
		gaussx[i] = exp(-0.5*(i*i)/sqsigmax);
	}
	std::vector<double> gaussf(255*3, 0.0);
	for(int i=0; i<gaussf.size(); ++i) {
		gaussf[i] = exp(-0.5*(i*i)/sqsigmaf);
	}

	int xo, yo;
	double w, tw, newcost;
	for(int y=0; y<ny; ++y) {
		for(int x=0; x<nx; ++x) {
			for(int l=0; l<nLabels; ++l) {
				newcost = 0;
				tw = 0.0;
				for(int i=-radius; i<=radius; ++i) {
					for(int j=-radius; j<=radius; ++j) {
						xo = x + i;
						if(xo<0) xo=0;
						if(xo>=nx) xo = nx-1;

						yo = y + j;
						if(yo<0) yo=0;
						if(yo>=ny) yo=ny-1;

						w = gaussx[abs(xo-x)]*gaussx[abs(yo-y)]*gaussf[abs(img[y][x].red - img[yo][xo].red)]*gaussf[abs(img[y][x].green - img[yo][xo].green)]*gaussf[abs(img[y][x].blue - img[yo][xo].blue)];
						tw += w;
						newcost += w*ocost[yo][xo][l];
					}
				}
				newcost /= tw;
				cost[y][x][l] = newcost;
			}
		}
	}
}



#endif
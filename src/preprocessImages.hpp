#ifndef _PREPROCESS_IMAGES_HPP_
#define _PREPROCESS_IMAGES_HPP_

// Systems includes... I/O stuff and PNG include. 
#include <iostream>
#include <fstream>
#include <png++/png.hpp>
#include <stdlib.h>

void getImageMeanByColor(png::image< png::rgb_pixel >& img, std::vector<double>& mean) {
	int ny = img.get_height();
	int nx = img.get_width();
	for(int i=0; i<3; ++i) {
		mean[i] = 0.0;
	}
	for(int y=0; y<ny; ++y) {
		for(int x=0; x<nx; ++x) {
			mean[0] += img[y][x].red; 
			mean[1] += img[y][x].green;
			mean[2] += img[y][x].blue;
		}
	}
	for(int i=0; i<3; ++i) {
		mean[i] /= nx;
		mean[i] /= ny;
	}
}

void getImageDeviationByColor(png::image< png::rgb_pixel >& img, std::vector<double>& mean, std::vector<double> dev) {
	int ny = img.get_height();
	int nx = img.get_width();
	for(int i=0; i<3; ++i) {
		dev[i] = 0.0;
	}
	for(int y=0; y<ny; ++y) {
		for(int x=0; x<nx; ++x) {
			dev[0] += (mean[0] - img[y][x].red)*(mean[0] - img[y][x].red);
			dev[1] += (mean[1] - img[y][x].green)*(mean[1] - img[y][x].green);
			dev[2] += (mean[2] - img[y][x].blue)*(mean[2] - img[y][x].blue);
		}
	}
	for(int i=0; i<3; ++i) {
		dev[i] /= ny;
		dev[i] /= nx;
		dev[i] = sqrt(dev[i]);
	}

}

void imageFix(png::image< png::rgb_pixel >& img,
		std::vector<double>& inmean,std::vector<double>& indev,
		std::vector<double>& outmean,std::vector<double>& outdev
	) {
	int ny = img.get_height();
	int nx = img.get_width();
	int temp;
	for(int y=0; y<ny; ++y) {
		for(int x=0; x<nx; ++x) {
			temp = (((double(img[y][x].red - inmean[0]))/indev[0])*outdev[0]) + outmean[0];
			temp = (temp<0)?0 :( temp>255? 255: temp);
			img[y][x].red = temp;

			temp = (((double(img[y][x].green - inmean[1]))/indev[1])*outdev[1]) + outmean[1];
			temp = (temp<0)?0 :( temp>255? 255: temp);
			img[y][x].green = temp;

			temp = (((double(img[y][x].blue - inmean[2]))/indev[2])*outdev[2]) + outmean[2];
			temp = (temp<0)?0 :( temp>255? 255: temp);
			img[y][x].blue = temp;
		}
	}
}

void normalize(
	// arguments
		png::image< png::rgb_pixel >& limg,
		png::image< png::rgb_pixel >& rimg
	) {
	// get their statistics.
	std::vector<double> lmean(3,0.0);
	std::vector<double> ldev(3,0.0);
	std::vector<double> rmean(3,0.0);
	std::vector<double> rdev(3,0.0);
	getImageMeanByColor(limg, lmean);
	getImageMeanByColor(rimg, rmean);

	getImageDeviationByColor(limg, lmean, ldev);
	getImageDeviationByColor(rimg, rmean, rdev);

	double threshold = 0.16;
	bool normalizationNeeded = false;
	for(int i=0; i<3; ++i) {
		if ((2*(lmean[i]-rmean[i])/(lmean[i]+rmean[i])) > threshold) {
			normalizationNeeded = true;
			break;
		} else if ((2*(ldev[i]-rdev[i])/(ldev[i]+rdev[i])) > threshold) {
			normalizationNeeded = true;
			break;
		}
	}

	// normalize them
	if(normalizationNeeded) {
		std::cout << "Need to alter mean and deviation." << std::endl;
		std::vector<double> targetmean(3, 0.0);
		std::vector<double> targetdev(3, 0.0);
		for(int i=0; i<3; ++i) {
			targetmean[i] = 0.5*(lmean[i]+rmean[i]);
			targetdev[i] = 0.5*(ldev[i]+rdev[i]);
		}
		imageFix(limg, lmean, ldev, targetmean, targetdev);
		imageFix(rimg, rmean, rdev, targetmean, targetdev);
	} else {
		std::cout << "Don't need to normalize." << std::endl;
	}
}
#endif
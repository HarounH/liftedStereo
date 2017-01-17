/*
Copyright 2011. All rights reserved.
Institute of Measurement and Control Systems
Karlsruhe Institute of Technology, Germany

This file is part of libelas.
Authors: Andreas Geiger

libelas is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either version 3 of the License, or any later version.

libelas is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
libelas; if not, write to the Free Software Foundation, Inc., 51 Franklin
Street, Fifth Floor, Boston, MA 02110-1301, USA 
*/

// basic image I/O, based on Pedro Felzenszwalb's code

#ifndef IMAGE_H
#define IMAGE_H

#include <cstdlib>
#include <climits>
#include <cstring>
#include <fstream>
#include <algorithm>
#include <png++/png.hpp>

// use imRef to access image data.
#define imRef(im, x, y) (im->access[y][x])
  
// use imPtr to get pointer to image data.
#define imPtr(im, x, y) &(im->access[y][x])

#define BUF_SIZE 256

typedef unsigned char uchar;
typedef struct { uchar r, g, b; } rgb;

inline bool operator==(const rgb &a, const rgb &b) {
  return ((a.r == b.r) && (a.g == b.g) && (a.b == b.b));
}

// image class
template <class T> class image {
public:

  // create image
  image(const int width, const int height, const bool init = false);

  // delete image
  ~image();

  // init image
  void init(const T &val);

  // deep copy
  image<T> *copy() const;
  
  // get image width/height
  int width() const { return w; }
  int height() const { return h; }
  
  // image data
  T *data;
  
  // row pointers
  T **access;
  
private:
  int w, h;
};

template <class T> image<T>::image(const int width, const int height, const bool init) {
  w = width;
  h = height;
  data = new T[w * h];  // allocate space for image data
  access = new T*[h];   // allocate space for row pointers
  
  // initialize row pointers
  for (int i = 0; i < h; i++)
    access[i] = data + (i * w);  
  
  // init to zero
  if (init)
    memset(data, 0, w * h * sizeof(T));
}

template <class T> image<T>::~image() {
  delete [] data; 
  delete [] access;
}

template <class T> void image<T>::init(const T &val) {
  T *ptr = imPtr(this, 0, 0);
  T *end = imPtr(this, w-1, h-1);
  while (ptr <= end)
    *ptr++ = val;
}


template <class T> image<T> *image<T>::copy() const {
  image<T> *im = new image<T>(w, h, false);
  memcpy(im->data, data, w * h * sizeof(T));
  return im;
}

class pnm_error {};

void pnm_read(std::ifstream &file, char *buf) {
  char doc[BUF_SIZE];
  char c;
  
  file >> c;
  while (c == '#') {
    file.getline(doc, BUF_SIZE);
    file >> c;
  }
  file.putback(c);
  
  file.width(BUF_SIZE);
  file >> buf;
  file.ignore();
}
image<rgb> *loadBuf(unsigned char *buf, int w, int h) {

  // read data
  image<rgb> *im = new image<rgb>(w, h);
  int cnt = w*h;
for(int p=0; p < cnt*3; p++) 
{
	int cl = p/cnt;
	if(cl==0){im->data[p%cnt].b = buf[p];  }
	if(cl==1){im->data[p%cnt].g = buf[p];  }
	if(cl==2){im->data[p%cnt].r = buf[p];  }

}

  return im;
}

image<rgb> *loadPGM(const char *name) {
  char buf[BUF_SIZE];
  
  // read header
  std::ifstream file(name, std::ios::in | std::ios::binary);
  pnm_read(file, buf);
  if (strncmp(buf, "P6", 2)) {
    std::cout << "ERROR: Could not read file " << name << std::endl;
    throw pnm_error();
  }

  pnm_read(file, buf);
  int width = atoi(buf);
  pnm_read(file, buf);
  int height = atoi(buf);

  pnm_read(file, buf);
  if (atoi(buf) > UCHAR_MAX) {
    std::cout << "ERROR: Could not read file " << name << std::endl;
    throw pnm_error();
  }

  // read data
  image<rgb> *im = new image<rgb>(width, height);
  file.read((char *)imPtr(im, 0, 0), width * height * sizeof(rgb));

  return im;
}


// Function by Haroun H.
image<rgb>* loadPNG(const char* fname) {
  png::image< png::rgb_pixel > img(fname);
  // Now to convert img into a image<rgb>
  int width = img.get_width();
  int height= img.get_height();
  image<rgb> *im = new image<rgb>(width, height);
  for(int y=0; y<height; ++y) {
    for(int x=0; x<width; ++x) {
      imRef(im, x, y).r = img[y][x].red;
      imRef(im, x, y).g = img[y][x].green;
      imRef(im, x, y).b = img[y][x].blue;
    }
  }
  return im;
}


void savePGM(image<rgb> *im, const char *name) {
  int width = im->width();
  int height = im->height();
  std::ofstream file(name, std::ios::out | std::ios::binary);

  file << "P6\n" << width << " " << height << "\n" << UCHAR_MAX << "\n";
  file.write((char *)imPtr(im, 0, 0), width * height * sizeof(rgb));
}

void savePGM(image<rgb> *im, std::string name) {
  int width = im->width();
  int height = im->height();
  std::ofstream file(name, std::ios::out | std::ios::binary);

  file << "P6\n" << width << " " << height << "\n" << UCHAR_MAX << "\n";
  file.write((char *)imPtr(im, 0, 0), width * height * sizeof(rgb));
}

#include <png++/png.hpp>

void ucharArray2pngImage(
	// arguments 
		png::image< png::rgb_pixel >& img,
		unsigned char * imgBuf,
		int nx, int ny
	) {
	
	for(int y=0; y<ny; ++y) {
		for(int x=0; x<nx; ++x) {
		img[y][x] = png::rgb_pixel(imgBuf[0 + x + y*nx], imgBuf[nx*ny + x + y*nx], imgBuf[nx*ny*2 + x + y*nx]);
		}
	}
}

void saveVecDisparities(
	// arguments
		std::vector<int>& disp,
		int nx, int ny, double scalefactor
	) {
	png::image< png::ga_pixel > img(nx,ny);
	for(int y=0;y<ny; ++y) {
		for(int x=0; x<nx; ++x) {
			img[y][x] = png::ga_pixel( scalefactor*disp[y*nx + x] );
		}
	}
	img.write("lol.png");
}

#endif

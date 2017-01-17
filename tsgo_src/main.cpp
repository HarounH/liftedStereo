// modified 2016-2017 -- messing around to add support for symmetries.
// modified 11/28/2013 -- added no_interpolation option via command-line
// changed 12/6/2013 to just use the ipol_gap_width parameter for this purpose

// 11/29/2013 -- changed unknown values from -10 to INFINITY when saving as .pfm

// older changes:

// modified 10/24/2013
// got rid of demo modes
// take disparity range as command line arg
// output pfms directly

// modified by DS 4/4/2013
// - added mode 'midd' to process 9 full-size pairs
// - added support for two param settings (need to recompile)
//   (1 == 'robotics', 2 == 'Middlebury', i.e. hole filling)
// - added printing of timing info
// - turned off autoscaling of disp's
#define  INFINITY (float)256*256*256*125
static const char *usage = "\n  usage: %s im0.png im1.png disp.png maxdisp nRanksConsidered nCBPIterations truth_file [no_interp=0]\n";

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

// Demo program showing how libelas can be used, try "./elas -h" for help

#include <iostream>
#include <ctime>

#include <cmath>
#include <algorithm>
#include <inttypes.h>
#include <png++/png.hpp>


int ranksConsidered=0;
int cbpIters=1;
std::string truthFile;
double threshold=1.0;

#include "evaluator.hpp"
#include "image.h"
#include "output.hpp"

#include "EDP.h"
#include "EDP.cpp"


using namespace std;


float  STD_get(uchar * I,  float inm,  int cnt ) {
	float sum = 0;
	for(int p=0; p < cnt; p++){ 
		sum +=((float)I[p]- inm)*((float)I[p]- inm);
	}
	return sqrt(sum/cnt);
}
float  MN_get(uchar * I,   int cnt )
{
	float sum = 0;
	for(int p=0; p < cnt; p++)
	{ 
		sum += I[p];
	}
	return (sum/cnt);
}
void STD_fix(uchar * I,  float inm, float ind, float outm, float outd, int cnt )
{
	for(int p=0; p < cnt; p++) {
		float vl = (float)(I[p] - inm)/ind*outd + outm;
		if(vl <0) vl =0; if(vl>255) vl = 255;
		I[p] = vl;
	}
}
float * wt_fl( int sc)
{
	int cnt = (2*sc+1)*(2*sc+1);
	float * wt = new float [cnt];
	float sum =0;
	int wcnt =0;
	float sg = (float)sc/2.125; sg *= sg;
	for(int i = -sc; i <= sc; i++)
	for(int j = -sc; j <= sc; j++)
	{ 
		float v = j*j+i*i;
		sum += wt[wcnt++] = exp( -v/2./sg);
	}
	for(int i = 0; i < cnt; i++) wt[i] /= sum; 
	return wt;
}
void Img_scl_( uchar *out, float * wt,  uchar *in, int w, int h, int xi, int yi, int sc)
{
	
	int xx = xi*sc; int yy = yi*sc;
	int sh = w*h;
	
	float sum[3] = {0,0,0};
	for(int c = 0; c <3; c++){ 
		int wcnt =0;
		for(int i = -sc/2; i <= sc/2; i++)
			for(int j = -sc/2; j <= sc/2; j++)
			{ 
				int x = xx + i; if(x>w-1) x = w-1; if(x<0) x=0;
				int y = yy + j; if(y>h-1) y =  h-1; if(y< 0) y=0;
				sum[c] += wt[wcnt++]*in[x+y*w + c*sh];
			}
			out[c] = (uchar) sum[c];
	}


}
//   void Img_scl_med_( uchar *md,  uchar *in,  char * sh_xy, int * H, int w, int h, int xi, int yi, int sc)
//{
//  
//   int xx = xi*sc; int yy = yi*sc;
//    int sh = w*h; int sh_h = sc/2*2 + 1; sh_h *= sh_h;
//  int sum[3] = {0,0,0};   
//   uchar mx[3], mn[3],  mdcr[3];
//  for(int c = 0; c <3; c++)
//  for(int i = -sc/2, cnt = 0; i <= sc/2; i++)
//  for(int j = -sc/2; j <= sc/2; j++, cnt++)
//  {
//    
//   int x = xx + i; if(x>w-1) x = w-1; if(x<0) x=0;
//   int y = yy + j; if(y>h-1) y =  h-1; if(y< 0) y=0;
//     uchar V =  in[x+y*w + c*sh];  
//   if(!cnt) {mx[c] = V; mn[c] = V;}
//   else {if(mx[c] <V) mx[c] = V;  if(mn[c] >V) mn[c] = V; }
//   H[V+ c*256]++; 
//      
//  }
//  for(int c = 0; c <3; c++){ int st = 1;
//  for(int i = mn[c]; i <=mx[c]; i++)
//  {
//    if(H[i + c*256]){sum[c] += H[i + c*256]; if(st&&sum[c]>sh_h/2){st=0; mdcr[c] = i;} H[i + c*256] = 0;}
//    
//  }}
//    float mnr;
//  for(int i = -sc/2, cnt = 0; i <= sc/2; i++)
//  for(int j = -sc/2; j <= sc/2; j++, cnt++)
//  {
//     int x = xx + i; if(x>w-1) x = w-1; if(x<0) x=0;
//   int y = yy + j; if(y>h-1) y =  h-1; if(y< 0) y=0;
//   float R = 0;  for(int c = 0; c <3; c++){float v = (mdcr[c]-in[x+y*w + c*sh]); R += v*v;}  
//   if(!cnt) {mnr = R; sh_xy[0] = i; sh_xy[1] = j; for(int c = 0; c <3; c++)md[c] = in[x+y*w + c*sh]; }
//   else   if(R<mnr){mnr = R;  sh_xy[0] = i; sh_xy[1] = j;  for(int c = 0; c <3; c++)md[c] = in[x+y*w + c*sh]; }
//      }
//
//
//
//
//}
uchar *  Img_scl(  uchar *in, int w, int h, int sc)
{
	float *wt =  wt_fl( sc/2);
	int ws = (w%sc) ? w/sc+1: w/sc;
	int hs = (h%sc) ? h/sc+1: h/sc;
	uchar rgbB[3];
	uchar * out = new uchar [ws*hs*3];
	for(int i = 0; i < ws; i++)
	for(int j = 0; j < hs; j++)
	{
		Img_scl_( rgbB,  wt,  in,  w,h, i, j,  sc);
		for(int c = 0; c < 3; c++)out[i + j*ws + c*ws*hs] = rgbB[c];
	}
	delete [] wt;
	return out;
}
uchar *  Img_scl_s(  uchar *in, int w, int h, int sc)
{
	int ws = (w%sc) ? w/sc+1: w/sc;
	int hs = (h%sc) ? h/sc+1: h/sc;
	uchar rgbB[3];
	uchar * out = new uchar [ws*hs*3];
	for(int i = 0; i < ws; i++)
		for(int j = 0; j < hs; j++)
		{
			for(int c = 0; c < 3; c++)
				out[i + j*ws + c*ws*hs] = in[i*sc + j*sc*w + c*w*h];
		}

	return out;
}
/* uchar *  Img_scl_med(  uchar *in, char *dd, int w, int h, int sc)
{
	int ws = (w%sc) ? w/sc+1: w/sc;
	int hs = (h%sc) ? h/sc+1: h/sc;
	uchar rgbB[3];
	char sh_xy[2];
	int H[256*3]; for(int i = 0; i <256*3; i++) H[i] = 0;
	uchar * out = new uchar [ws*hs*3];
	for(int i = 0; i < ws; i++)
	for(int j = 0; j < hs; j++)
	{
	 Img_scl_med_( rgbB,  in,  sh_xy, H, w, h, i, j, sc);
	  for(int c = 0; c < 3; c++)out[i + j*ws + c*ws*hs] =  rgbB[c];
	  for(int c = 0; c < 2; c++)dd[i + j*ws + c*ws*hs] = sh_xy[c];
	}
 return out;
}*/
void   Dmap_scl_(float sig, uchar* outim, uchar* inim,   float * out, float *in, int w, int h, int sc) {
	sig *=sig; 
	int ws = (w%sc) ? w/sc+1: w/sc;
	int hs = (h%sc) ? h/sc+1: h/sc;


	int sh_in = ws*hs;
	for(int i = 0; i < w; i++)
	for(int j = 0; j < h; j++)
	{
	 
		int wss = i/sc; int hss = j/sc;
		int wsp = (wss + 1 < ws)? wss +1 : wss;
		int hsp = (hss + 1 < hs)? hss +1 : hss;
		int p[4] = {wss + hss*ws, wss + hsp*ws, wsp + hss*ws, wsp  + hsp*ws};
		float wtc[4] = {0,0,0,0}, rc[4];

		for(int c = 0; c < 3; c++)
		{
			float vloutim = outim[i + j*w + c*h*w];
			for(int pi = 0; pi < 4; pi++){
				rc[pi] =(float)inim[p[pi]+ c*sh_in] - vloutim; 
				rc[pi] *= rc[pi];
				wtc[pi] += rc[pi];
			}
		}
		
		for(int pi = 0; pi < 4; pi++) wtc[pi]  =  exp(- wtc[pi]/2/sig);
		
		float wtsum = 0; 
		out[i + j*w] = 0;
		for(int pi = 0; pi < 4; pi++)
		{
			out[i + j*w] += in[p[pi]]/**wt[pi][i%sc][j%sc]*/*wtc[pi];
			wtsum += /*wt[pi][i%sc][j%sc]**/wtc[pi];
		}

		out[i+j*w] *= (wtsum)? (float)sc/wtsum : sc;
	
	}

}

void Dmap_scl(  float * out, float *in, int w, int h, int sc)
{
	
	int ws = (w%sc) ? w/sc+1: w/sc;
	int hs = (h%sc) ? h/sc+1: h/sc;
	float wt[4][4][4];
	for(int sx = 0; sx <sc; sx++)
	for(int sy = 0; sy <sc; sy++)
	{
	 wt[0][sx][sy] = (sc-sx)*(sc-sy);
	 wt[1][sx][sy] = (sc-sx)*(     sy);
	 wt[2][sx][sy] = (     sx)*(sc-sy);
	 wt[3][sx][sy] = (     sx)*(    sy);
	 float sum =wt[0][sx][sy] + wt[1][sx][sy] + wt[2][sx][sy] + wt[3][sx][sy];  
	 for(int i = 0 ; i < 4; i++)  wt[i][sx][sy] /= sum;
	}

	for(int i = 0; i < w; i++)
	for(int j = 0; j < h; j++)
	{
	 
	int wss = i/sc; int hss = j/sc;
	int wsp = (wss + 1 < ws)? wss +1 : wss;
	int hsp = (hss + 1 < hs)? hss +1 : hss;

	 out[i + j*w] = in[wss + hss*ws]*wt[0][i%sc][j%sc] + in[wss + hsp*ws]*wt[1][i%sc][j%sc] 
	 + in[wsp + hss*ws]*wt[2][i%sc][j%sc] + in[wsp  + hsp*ws]*wt[3][i%sc][j%sc] ;
	 out[i+j*w] *= sc;
	
	}

}

int Img_set(uchar * I1, uchar *I2, rgb * I1_, rgb *I2_, int cnt ) {
	float m1[3]={0,0,0}; 
	float d1[3]={0,0,0}; 
	float m2[3]={0,0,0}; 
	float d2[3]={0, 0, 0};
	
	for(int p=0; p < cnt*3; p++) 
	{
		int cl = p/cnt;
		if(cl==0){m1[0] += I1[p] = I1_[p%cnt].b; m2[0] += I2[p] = I2_[p%cnt].b; }
		if(cl==1){m1[1] += I1[p] = I1_[p%cnt].g; m2[1] += I2[p] = I2_[p%cnt].g; }
		if(cl==2){m1[2] += I1[p] = I1_[p%cnt].r; m2[2] += I2[p] = I2_[p%cnt].r; }
	}
	for(int c =0; c <3; c++){
		m1[c] /= cnt; m2[c] /= cnt;
	}
	for(int p=0; p < cnt*3; p++) 
	{
		int cl = p/cnt;
		if(cl==0){d1[0] += (I1[p]-m1[0])*(I1[p]-m1[0]); d2[0] += (I2[p]-m2[0])*(I2[p]-m2[0]); }
		if(cl==1){d1[1] += (I1[p]-m1[1])*(I1[p]-m1[1]); d2[1] += (I2[p]-m2[1])*(I2[p]-m2[1]); }
		if(cl==2){d1[2] += (I1[p]-m1[2])*(I1[p]-m1[2]); d2[2] += (I2[p]-m2[2])*(I2[p]-m2[2]); }
	}
	float thr =0.16;
	int steq=0;
	for(int c =0; c <3; c++){
		d1[c] = sqrt(d1[c]/cnt); 
		d2[c] = sqrt(d2[c]/cnt);
	}
	for(int c =0; c <3; c++){ 
		float vlm, vld;
		if((vlm = fabs(m1[c]-m2[c])*2/(m1[c]+m2[c]))> thr )steq =1;
		if((vld = fabs(d1[c]-d2[c])*2/(d1[c]+d2[c]))> thr)steq =1;
	}
	if(steq) {
		for(int c =0; c <3; c++) {
			STD_fix(&I1[cnt*c],  m1[c], d1[c], (m1[c]+m2[c])*0.5, (d1[c]+d2[c])*0.5, cnt );
			STD_fix(&I2[cnt*c],  m2[c], d2[c], (m1[c]+m2[c])*0.5, (d1[c]+d2[c])*0.5, cnt );

		}
	}
	return steq;
}
int littleendian()
{
	int intval = 1;
	uchar *uval = (uchar *)&intval;
	return uval[0] == 1;
}


// write pfm image (added by DS 10/24/2013)
// 1-band PFM image, see http://netpbm.sourceforge.net/doc/pfm.html
void WriteFilePFM(float *data, int width, int height, const char* filename, float scalefactor=1/255.0) {
	// Open the file
	FILE *stream = fopen(filename, "wb");
	if (stream == 0) {
		fprintf(stderr, "WriteFilePFM: could not open %s\n", filename);
		exit(1);
	}

	// sign of scalefact indicates endianness, see pfms specs
	if (littleendian()) {
		scalefactor = -scalefactor;
	}

	// write the header: 3 lines: Pf, dimensions, scale factor (negative val == little endian)
	fprintf(stream, "Pf\n%d %d\n%f\n", width, height, scalefactor);

	int n = width;
	// write rows -- pfm stores rows in inverse order!
	for (int y = height-1; y >= 0; y--) {
		float* ptr = data + y * width;
		// change invalid pixels (which seem to be represented as -10) to INF
		for (int x = 0; x < width; x++) {
			if (ptr[x] < 0)
			ptr[x] = INFINITY;
		}
		if ((int)fwrite(ptr, sizeof(float), n, stream) != n) {
			fprintf(stderr, "WriteFilePFM: problem writing data\n");
			exit(1);
		}
	}	
	// close file
	fclose(stream);
}

// compute disparities of pgm image input pair file_1, file_2
void process (const char* file_1, const char* file_2, const char* outfile, int maxdisp, int no_interp) 
{
   
	clock_t c0 = clock();

	// load images
	image<rgb> *I1,*I2, *Is;
	// I1 = loadPGM(file_1);
	// I2 = loadPGM(file_2);
	I1 = loadPNG(file_1);
	I2 = loadPNG(file_2);
	// string filename1(file_1);
	// string filename2(file_2);
	// savePGM(I1, (filename1.append(".pgm")));
	// savePGM(I2, (filename2.append(".pgm")));
	cout << "Loaded PNGs" << endl;
	//----------------------------
	
	// check for correct size
	if (I1->width()<=0 || I1->height() <=0 || I2->width()<=0 || I2->height() <=0 ||
		I1->width()!=I2->width() || I1->height()!=I2->height()) {
		cout << "ERROR: Images must be of same size, but" << endl;
		cout << "       I1: " << I1->width() <<  " x " << I1->height() << 
			", I2: " << I2->width() <<  " x " << I2->height() << endl;
		delete I1;
		delete I2;
		return; 
	}

	// get image width and height
	int32_t width  = I1->width();
	int32_t height = I1->height();

	// allocate memory for disparity images
	/*    const int32_t dims[3] = {width,height,width}; */// bytes per line = width
	float* D1_data = (float*)malloc(width*height*sizeof(float));
	
	// float* D2_data = (float*)malloc(width*height*sizeof(float));

	uchar *Images[2]; 
	for(int i =0; i <2; i++)  {
		Images[i] = new uchar [width*height*3];
	}
	
	int Lt = Img_set(Images[0], Images[1], I1->data, I2->data, width*height );
	Lt = 1;
	int scf =1;
	int med_f =0;

	uchar *Imsc[2];
	uchar * sc_msk;
	char * Sh_xy[2] = {NULL, NULL};
//////////////////////////// SCALE   ////////////////////////
	if(width >1500) {scf = 4;}
	if(width<=1500 && width >750) {scf = 2;}
	if(scf == 2 && maxdisp > 180) scf = 3;
	if(scf == 4 && maxdisp > 300) scf = 6;
	int ws = (width%scf) ? width/scf+1: width/scf;
	int hs = (height%scf) ? height/scf+1: height/scf;
	if(scf > 1) {
		/* if(!med_f)*/
		{
			Imsc[0] = Img_scl_s(  Images[0], width, height, scf);
			Imsc[1] = Img_scl_s(  Images[1], width, height, scf);
		}
		//else {Sh_xy[0] = new char  [ws*hs*2];
		//Sh_xy[1] = new char  [ws*hs*2];
		//Imsc[0] = Img_scl_med( Images[0], Sh_xy[0],width, height, scf);
		//Imsc[1] = Img_scl_med( Images[1], Sh_xy[1],width, height, scf);}
	}

////////////////////////////////////////////////////////////
   
	 //Is = loadBuf(Imsc[0], width/scf, height/scf);
	 //savePGM(Is, file_1);
	/////
	EDP *edp = NULL;  
	float* D1_sc = NULL;
	if(scf==1) {
		edp = new EDP(1, Images, width, height, maxdisp+1);
		edp->Size_fw = (no_interp) ? no_interp : 9;
		edp->Lt = Lt;
		if(Lt) {
			printf("1\n"); 
		}
		else {
			printf("0\n");
		}
		edp->findSlt(D1_data);
	} else {
		D1_sc = (float*)malloc(ws*hs*sizeof(float));
		edp = new EDP(1, Imsc, ws, hs, (maxdisp+1)/scf);
		edp->Size_fw = (no_interp) ? no_interp : 9;
		edp->Lt = Lt;
		if(Lt) {
			printf("1\n"); 
		} else {
			printf("0\n");
		}
		edp->Sc_out = scf;
		edp->findSlt(D1_sc);
		/*else edp->findSlt(D1_sc, Sh_xy[0], Sh_xy[1], 1./scf);*/
		float sig = 0.05* 256.;
		//Dmap_scl(  D1_data, D1_sc, width, height, scf);
		Dmap_scl_(sig, Images[0], Imsc[0], D1_data, D1_sc, width, height, scf);
	}
	
	//////
	clock_t c1 = clock();
	double secs = (double)(c1 - c0) / CLOCKS_PER_SEC;

	// save disparity image
	// for(int i=0; i< (width*height); ++i) { // For loop added by haroun.
	// 	if ((D1_data[i] - ((int)(D1_data[i]))) != 0.0) {
	// 		printf("We have a problem, houston, at i= %d, D1_data=%.2f\n", i, D1_data[i]);
	// 	}
	// }

	// FIXME haroun commented out WriteFilePFM so that the output is no longer a pfm...
	outputGrayDisparities(D1_data, width, height, outfile, 255.0/maxdisp);
	std::cout << "scf= " << scf << std::endl;
	// WriteFilePFM(D1_data, width, height, outfile, 1.0/maxdisp);

	printf("runtime: %.2fs  (%.2fs/MP)\n", secs, secs/(width*height/1000000.0));
	
	// free memory
	delete I1;
	delete I2;
	delete [] (D1_data);
	if(D1_sc) delete [] D1_sc;
	for(int i =0; i <2; i++) delete [] Images[i];
	if(edp) delete edp;
	if(Imsc[0]) delete [] Imsc[0];
	if(Imsc[1]) delete [] Imsc[1];
	if(Sh_xy[0]) delete [] Sh_xy[0];
	if(Sh_xy[1]) delete [] Sh_xy[1];
}

int main (int argc, char** argv) 
{

	if (argc < 7) {
		fprintf(stderr, usage, argv[0]);
		exit(1);
	}

	const char *file1 = argv[1];
	const char *file2 = argv[2];
	const char *outfile = argv[3];
	int maxdisp = atoi(argv[4]);
	ranksConsidered = atoi(argv[5]);
	cbpIters = atoi(argv[6]);
	truthFile = argv[7];
	
	int no_interp = 0;
	if (argc > 8) {
		no_interp = atoi(argv[8]);
	}

	process(file1, file2, outfile, maxdisp, no_interp);

	return 0;
}

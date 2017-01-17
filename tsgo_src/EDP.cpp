#include <cstdio>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cassert>
// #include <new>
#include "EDP.h"
#include <algorithm>

#include "modelCreation.hpp"
#include "infer.hpp"

#define private public
#undef private

#define MY_LFT 1
#define MY_RGT 2
#define MY_LFT_A 3
#define MY_RGT_A 4
#define MY_VIS 0
#define MY_FLT 5
#define N_PX m_nPixels
#define N_LB  m_nLabels
#define BI_WS(x,y,f)   Bi_ws[(x) + (y)*sigmai + (f)*sigmai2]
#define BI_WSx(x,y,f)   Bi_ws[(x) + (y)*sigmaix + (f)*sigmai2]
#define RNDP(x)     (int)((x) + 0.5)
#define CUTI(x,c)     ((x) < (c)) ? (x) : (c) -1
#define CUTFL(x,c)  CUTI(RNDP((x)), (c))
#define IND_HC(r,g,b,n) ((r) + (n)*((g) + (n)*b))
#define IND_HC_B(c,i,n)  ((c)==0)? (i)%(n) : (((c)==1)? ((i)/(n))%(n) : (i)/((n)*(n)))
#define FOR_I(t) for( int i =0; i < (t); i++)
#define FOR_PX_p for( int p =0; p < m_nPixels; p++)
#define FOR_DSI_p for( int p =0; p < m_nPixels*m_nLabels; p++)
#define FOR_IxP(t) for( int i =0; i < ((t)*m_nPixels); i++)

#define IND_IC(x,y,c) ((x) + (y)*m_width + (c)*m_nPixels)
#define IND_I(x,y)       ((x) + (y)*m_width)
#define IND(x,y)       ((x) + (y)*m_width)
#define TR_XX(x)      ((x)  < m_width - 1) ? (x)+1 :  m_width - 1
#define TR_XY0(x)      ((x)  > 0) ? (x)-1 : 0
#define TR_YY(y)      ((y)  < m_height - 1) ? (y)+1 :  m_height - 1
#define IF_1_L     if( lb_i1 == MY_LFT ||lb_i1 == MY_LFT_A )
#define IF_2_L     if( lb_i2 == MY_LFT ||lb_i2 == MY_LFT_A )
#define IF_1_R     if( lb_i1 == MY_RGT ||lb_i1 == MY_RGT_A )
#define IF_2_R     if( lb_i2 == MY_RGT ||lb_i2 == MY_RGT_A )
#define IF_1_V     if( lb_i1 == MY_VIS )
#define IF_2_V     if( lb_i2 == MY_VIS )
#define FOR_C1(l,x)   for( int c = 0; c<(l); c ++)im_max1[(x) +c*m_width] = rgb1[c]
#define FOR_C(l,x)   for( int c = 0; c<(l); c ++)im_max1[(x) +c*m_width] = (w1*rgb1[c] + w2*rgb2[c])/(w1+w2)
#define FOR_C2(l,x)   for( int c = 0; c<(l); c ++)im_max1[(x) +c*m_width] = rgb2[c]
#define FLR_FL(x) ((x)<0)? (((x)-(int)(x) != 0) ? (int)(x)-1 : (int)(x)) : (int)(x)
#define SLG_FL(x)  ((x)<0)?  (int)(x) :  (((x)-(int)(x) != 0) ? (int)(x)+1 : (int)(x))
#define IF_1_2_V     if( lb_i1 == MY_VIS && lb_i2 == MY_VIS)
#define Msg(y,dir)  Msg[((y)*4 +(dir))*Tau]
#define Mrg(y)  Mrg[(y)*Tau]
#define mm_D(y)  m_D[(y)*Tau]
#define MsgC(x,y,dir)  Msg[(((x)+ (y)*m_width)*4 + dir)*m_nLabels]
#define MrgC(x,y)  Mrg[((x)+ (y)*m_width)*m_nLabels]
#define Msg_add(x,y,dir)  Msg_add[(y)*m_nLabels + (x)+ m_nLabels*m_height*(dir)]
#define  LOAD_D(i, k)            for(int (l) =0; (l)<(k); (l)++) tmp_vl[(l)] = spx_D[(i) +(l)]
#define  DIV_MRG(k)            for(int (l) =0; (l)<(k); (l)++) tmp_vl[(l)] *= 0.5
#define  ADD_MS(e,n, k)      for(int (l) =0; (l)<(k); (l)++) tmp_vl[(l)] += ms[(e)][ind_img_dsi[(n)]+(l)]
#define  SUB_MS(e,n, k)      for(int (l) =0; (l)<(k); (l)++) tmp_vl[(l)] -= ms[(e)][ind_img_dsi[(n)]+(l)]
#define XLD(x)  tau_XL_d[(x)*m_nLabels]
#define XRD(x)  tau_XR_d[(x)*m_nLabels]
#define m_D(pix,l)  m_D[(pix)*m_nLabels+(l)]
#define m_V(l1,l2)  m_V[(l1)*m_nLabels+(l2)]
#define MY_INF 100000000
#define MY_SHRT_INF  10000
#define MIN(a,b)  (((a) < (b)) ? (a) : (b))
#define MAX(a,b)  (((a) > (b)) ? (a) : (b))
#define TRUNCATE_MIN(a,b) { if ((a) > (b)) (a) = (b); }
#define TRUNCATE_MAX(a,b) { if ((a) < (b)) (a) = (b); }
#define TRUNCATE TRUNCATE_MIN
#define FLOOR(a)  (((a) < 0) ? ((((int)(a)-(a))>0)? (int)(a) -1 : (int)(a)): (((-(int)(a)+(a))>0)? (int)(a) + 1 : (int)(a))
/////////////////////////////////////////////////////////////////////////////
//                  Operations on vectors (arrays of size K)               //
/////////////////////////////////////////////////////////////////////////////
inline double EDP::erf(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0) {
        sign = -1;
    }
    x = fabs(x);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
    // std::cout << "erf(" << x << ")=" << sign*y << std::endl;
    return sign*y;
}
inline void EDP::drw()
{
	for(int i = 0; i< m_width; i++ )
	{
		double x = (double)3*(i)/((double)m_width);
		/*double erv = erf(x);*/
		double ksi = i;
		double sigc = 100;
		//double erv =  fi_bf(ksi, x, sigc);

		double erv =  fi_erf(ksi,  sigc);
		if(x>2.9) {
			int a =0;
		}
		for(int j = 0; j< m_height; j++ )
		{
			double thr =( j - (double)m_height/2)/((double)m_height/2);
			for(int c =0; c<3; c++) I_ims[1][i + j*m_width + c*m_nPixels] = (erv>thr) ? 0:190;
		}
	}
}
inline double EDP::fi_bf(double ksi, double vi, double sigc)
{
	int wn =200; 
	double dlt = 1./wn; 
	double sum = 0, sumw =0;
	for(int i = -wn; i<= wn; i++ ) {
		double w = exp(-(ksi*ksi*(dlt*i-vi)*(dlt*i-vi)/2/sigc/sigc));
		sumw += w;
		sum += w*dlt*i;
	}
	return sum/sumw;
}
inline double EDP::fi_erf(double ksi, double vi, double sigc) {
	if(vi == 0 ) return 0;
	double sq2_ = sqrt(2.);  
	double mPi_ = asin(1.)*2;
	
	if(vi >  5.*sq2_*sigc/ksi+1) vi = (5.*sq2_*sigc/ksi+1);
	if(vi < -5.*sq2_*sigc/ksi-1) vi = -(5.*sq2_*sigc/ksi+1);
	
	double erf1 = erf((1.-vi)*ksi/sigc/sq2_);
	double erf2 = erf((1.+vi)*ksi/sigc/sq2_);
	double exp1 =exp(-(1.+vi)*ksi*(1.+vi)*ksi/2/sigc/sigc);
	double exp2 =exp(-(1.-vi)*ksi*(1.-vi)*ksi/2/sigc/sigc);
	double g = sqrt(2.*mPi_)*ksi*(erf1+erf2);
	double ret =(vi - ( g*vi+ 2.*sigc*(exp1-exp2))/g)/vi;

	return ret;
}

inline double  EDP::fi_erf(double ksi,  double sigc) {
	if(fabs(ksi)<0.0001) 
		ksi = 0.0001;
	double sq2_ = 1.41421356;  
	double mPi_ = 1.77245385;
	double gi = sq2_*ksi/sigc;
	double erf1 = erf(gi);
	double exp1 =exp(-gi*gi);
	double g =  2.*(1.- exp1)/(mPi_*gi*erf1);
	return g;
}

inline void EDP::fi_erf_ini() {
	sq2 = sqrt(2.);  
	mPi = asin(1.)*2;
	fi_erf_q =100; 
	double dlt = 3./fi_erf_q;
	fi_erf_tab = new double [fi_erf_q*fi_erf_q];
	for(int i = 0; i < fi_erf_q; i++ ) {
		for(int j = 0; j< fi_erf_q; j++ ) {
			fi_erf_tab[i + j*fi_erf_q]  =  fi_erf(dlt*i,  dlt*j, 1.);
		}
	}
}
inline double EDP::fi_erf_tab_f(double ksi, double vi, double sigc) {
	double ret;
	double ksin = ksi*fi_erf_q/sigc; 
	if(ksin>3*(fi_erf_q-1)) {
		ksin =3*(fi_erf_q-1);
	}
	double vin = vi*fi_erf_q;
	return ret;
}

inline void EDP::CopyVector(REAL* to, CostVal* from, int K) {
    REAL* to_finish = to + K;
    do {
		*to ++ = *from ++;
	} while (to < to_finish);
}

inline void EDP::AddVector( REAL* to,  REAL* from, int K) {
	REAL* to_finish = to + K;
    do{
		*to ++ += *from ++;
	} while (to < to_finish);
}

inline double EDP::SubtractMin(double *D, int K) {
    int k;
    double delta;
	
    delta = D[0];
    for (k=1; k<K; k++) TRUNCATE(delta, D[k]);
    for (k=0; k<K; k++) D[k] -= delta;

    return delta;
}



inline double EDP::UpdateMessageL1(REAL* M, REAL* Di_hat, int K, REAL gamma, CostVal lambda, CostVal smoothMax) {
	int k;
	EDP::REAL delta;
	delta = M[0] = gamma*Di_hat[0] - M[0];

	for (k=1; k<K; k++) {
		M[k] = gamma*Di_hat[k] - M[k];
		TRUNCATE(delta, M[k]);
		TRUNCATE(M[k], M[k-1] + lambda);
	}

	M[--k] -= delta;
	TRUNCATE(M[k], lambda*smoothMax);
	for (k--; k>=0; k--) {
		M[k] -= delta;
		TRUNCATE(M[k], M[k+1] + lambda);
		TRUNCATE(M[k], lambda*smoothMax);
	}

	return delta;
}
//inline EDP::REAL  EDP:: UpdMsgL1_(int * xr, EDP::REAL* M, EDP::REAL* Di_hat, EDP::REAL* add,int K, EDP::REAL lambda, EDP::REAL smoothMax)
//{
//    int k;
//    REAL delta = M[xr[0]]= Di_hat[0];
//    for (k=1; k<K; k++)
//	{
//		M[xr[k]] = Di_hat[k];
//	    TRUNCATE(delta, M[xr[k]]);
//	    TRUNCATE(M[xr[k]], M[xr[k-1]] + lambda);
//	}
//    k--;
//    M[xr[k]] -= delta;
//    TRUNCATE(M[xr[k]],  smoothMax);
//    for (k--; k>=0; k--)
//	{
//	    M[xr[k]] -= delta;
//	    TRUNCATE(M[xr[k]], M[xr[k+1]] + lambda);
//	    TRUNCATE(M[xr[k]], smoothMax);
//	}
//	add[0] = (M[xr[K-1]]+ lambda< smoothMax) ? M[xr[K-1]]+ lambda :smoothMax;
//    return delta;
//}
//
//
//

inline EDP::REAL EDP::UpdMsgL1_(int * xr, EDP::REAL* M, EDP::REAL* Di_hat, EDP::REAL* add,int K, EDP::REAL lambda, EDP::REAL smoothMax) {
	int k;
	REAL vlm,  vlp,   delta =  Di_hat[0]; for (k=1; k<K; k++)TRUNCATE(delta,  Di_hat[k]);
	for (k=0; k<K; k++) {
		M[xr[k]] = Di_hat[k] - delta;
		if(k){
			vlm = Di_hat[k-1] - delta + lambda; 
			TRUNCATE(M[xr[k]],  vlm);
		}
		if(k<K-1){
			vlp = Di_hat[k+1] - delta+ lambda; 
			TRUNCATE(M[xr[k]],  vlp);
		}
	}

	for (k=0; k<K; k++) {
		TRUNCATE(M[xr[k]],  smoothMax);
	}
	add[0] = (M[xr[K-1]]+ lambda< smoothMax) ? M[xr[K-1]]+ lambda :smoothMax;
	return delta;
}




inline EDP::REAL  EDP::UpdMsgL1(int * xr, EDP::REAL* M, EDP::REAL* Di_hat, int K, EDP::REAL lambda, EDP::REAL smoothMax) {
	int k;
	REAL vlm,  vlp,   delta =  Di_hat[0]; 
	for (k=1; k<K; k++)
		TRUNCATE(delta,  Di_hat[k]);
	
	for (k=0; k<K; k++) {
		M[xr[k]] = Di_hat[k] - delta;
		if(k){vlm = Di_hat[k-1] - delta + lambda; TRUNCATE(M[xr[k]],  vlm);}
		if(k<K-1){vlp = Di_hat[k+1] - delta + lambda; TRUNCATE(M[xr[k]],  vlp);}
	}

	for (k=0; k<K; k++) 
		TRUNCATE(M[xr[k]],  smoothMax);
	return delta;
}


inline EDP::REAL  EDP::UpdMsgCL1( EDP::REAL* M, EDP::REAL* Di_hat, int K, EDP::REAL lambda, EDP::REAL smoothMax) {
	int k;
	REAL vlm,  vlp,   delta =  Di_hat[0]; 
	for (k=1; k<K; k++)
		TRUNCATE(delta,  Di_hat[k]);
	
	for (k=0; k<K; k++) {
		M[k] = Di_hat[k] - delta;
		if(k){vlm = Di_hat[k-1] - delta + lambda; TRUNCATE(M[k],  vlm);}
		if(k<K-1){vlp = Di_hat[k+1] - delta + lambda; TRUNCATE(M[k],  vlp);}
	}

	for (k=0; k<K; k++) 
		TRUNCATE(M[k],  smoothMax);
    return delta;
}



//inline EDP::REAL  EDP:: UpdMsgL1(int * xr, EDP::REAL* M, EDP::REAL* Di_hat, int K, EDP::REAL lambda, EDP::REAL smoothMax)
//{
//    int k;
//    REAL delta = M[xr[0]]= Di_hat[0];
//    for (k=1; k<K; k++)
//	{
//		M[xr[k]] = Di_hat[k];
//	    TRUNCATE(delta, M[xr[k]]);
//	    TRUNCATE(M[xr[k]], M[xr[k-1]] + lambda);
//	}
//    k--;
//    M[xr[k]] -= delta;
//    TRUNCATE(M[xr[k]],  smoothMax);
//    for (k--; k>=0; k--)
//	{
//	    M[xr[k]] -= delta;
//	    TRUNCATE(M[xr[k]], M[xr[k+1]] + lambda);
//	    TRUNCATE(M[xr[k]], smoothMax);
//	}
//
//    return delta;
//}
//
//

inline void EDP::addMsg(int * xr, int * xsh, EDP::REAL* M, EDP::REAL* Mrg, int K)
{
	for (int k=0; k<K; k++) Mrg[xr[k]] += M[xsh[k]];
}
inline void EDP::addMsgC(EDP::REAL* M, EDP::REAL* Mrg, int K)
{
	for (int k=0; k<K; k++) Mrg[k] += M[k];
}

inline void  EDP::UpdMsgYL_mm(int x, int y,EDP::REAL* Di,EDP::REAL* Msg,EDP::REAL* Msg_add, EDP::REAL* Mrg)
{
	if (y >0) {
		int d_mx = (x + 1 <m_nLabels) ? x +1 : m_nLabels;
		subMsg(
			&XLD(x),
			&XLD(x),Di,&Msg(y-1,1),
			&Mrg(y),
			d_mx, 0.5);
		
		REAL  lmbd =Smthnss[(int)B_L_buf_mm_[IND_I(x,y)].y].semi*Lmb_Y; //(!G_L_buf_y[(y-1)*m_width + x]) ? this->C_gr_1*Lmbd_y : Lmbd_y;
		REAL  max =Smthnss[(int)B_L_buf_mm_[IND_I(x,y)].y].max*Lmb_Y;
		
		if(d_mx < m_nLabels)
			UpdMsgL1_(
			&XLD(x),
			&Msg(y,3),Di,
			&Msg_add(x,y,3),d_mx, lmbd, max);
		else
			UpdMsgL1(
			&XLD(x),
			&Msg(y,3),Di,d_mx, lmbd, max);
	}
}

inline void  EDP::UpdMsgYL_pp(int x, int y,EDP::REAL* Di,EDP::REAL* Msg,EDP::REAL* Msg_add, EDP::REAL* Mrg) {
	if (y < m_height-1)  {
		int d_mx = (x + 1 <m_nLabels) ? x +1 : m_nLabels;
		subMsg(
			&XLD(x),
			&XLD(x),Di,&Msg(y+1,3),
			&Mrg(y),
			d_mx, 0.5);
		REAL  lmbd =Smthnss[(int) B_L_buf_pp_[IND_I(x,y)].y].semi*Lmb_Y;// (!G_L_buf_y[y*m_width + x]) ? this->C_gr_1*Lmbd_y : Lmbd_y;
		REAL  max =Smthnss[(int) B_L_buf_pp_[IND_I(x,y)].y].max*Lmb_Y;
		if(d_mx < m_nLabels)
			UpdMsgL1_(
			&XLD(x),
			&Msg(y,1),Di,
			&Msg_add(x,y,1),d_mx, lmbd, max);
		else
			UpdMsgL1(
			&XLD(x),
			&Msg(y,1),Di,d_mx, lmbd, max);	
	}
}

inline void  EDP::UpdMsgXL_mm(int x, int y,EDP::REAL* Di,EDP::REAL* Msg,EDP::REAL* Msg_add, EDP::REAL* Mrg) {
	if (x >0) {
		int d_mx = (x + 1 <m_nLabels) ? x +1 : m_nLabels;
		if(x < m_nLabels)
			subMsg_(
			&XLD(x),
			&XLD(x-1),Di,&Msg(y,0),
			&Mrg(y),
			d_mx-1, Msg_add(x-1,y,0),
			0.5);
		else
			subMsg(
			&XLD(x),
			&XLD(x-1),Di,&Msg(y,0),
			&Mrg(y),
			d_mx, 0.5);

		REAL  lmbd =Smthnss[(int)B_L_buf_mm_[IND_I(x,y)].x].semi;//(!G_L_buf_x[y*m_width + x-1]) ? this->C_gr_1*Lmbd_y : Lmbd_y;
		REAL  max =Smthnss[(int)B_L_buf_mm_[IND_I(x,y)].x].max;
		
		if(d_mx < m_nLabels)
			UpdMsgL1_(
			&XLD(x),
			&Msg(y,2),Di,
			&Msg_add(x,y,2),d_mx, lmbd, max);
		else
			UpdMsgL1(
			&XLD(x),
			&Msg(y,2),Di,d_mx, lmbd, max);	
	}
	//----------------------------------------
}


inline void  EDP::UpdMsgXL_pp(int x, int y,EDP::REAL* Di,EDP::REAL* Msg,EDP::REAL* Msg_add, EDP::REAL* Mrg) {
	if (x < m_width-1) {
		int d_mx = (x + 1 <m_nLabels) ? x +1 : m_nLabels;
		subMsg(
			&XLD(x),
			&XLD(x+1),Di,&Msg(y,2),
			&Mrg(y),
			d_mx, 0.5);
		
		REAL  lmbd = Smthnss[(int)B_L_buf_pp_[IND_I(x,y)].x].semi;//(!G_L_buf_x[y*m_width + x]) ? this->C_gr_1*Lmbd_y : Lmbd_y;
		REAL  max =Smthnss[(int)B_L_buf_pp_[IND_I(x,y)].x].max;
		if(d_mx < m_nLabels)
			UpdMsgL1_(
				&XLD(x),
				&Msg(y,0),Di,
				&Msg_add(x,y,0),d_mx, lmbd, max);
		else
			UpdMsgL1(
				&XLD(x),
				&Msg(y,0),Di,d_mx, lmbd, max);	
	}
	//----------------------------------------
}

inline void EDP::UpdMsgC(int x, int y,EDP::REAL* Di, EDP::REAL* Msg, EDP::REAL* Mrg, int dir) {
	int sm = 0;/* double dd = (dir%2)? adapY[IND(x,y)]:adapX[IND(x,y)];*/
	if(dir<4) {
		if(dir%4 == 0)sm = B_L_buf_pp_[x + y*m_width].x;
		if(dir%4 == 2)sm = B_L_buf_mm_[x + y*m_width].x;
		if(dir%4 == 1)sm = B_L_buf_pp_[x + y*m_width].y;
		if(dir%4 == 3)sm = B_L_buf_mm_[x + y*m_width].y;
	}else{
		if(dir%4 == 0)sm = B_R_buf_pp_[x + y*m_width].x;
		if(dir%4 == 2)sm = B_R_buf_mm_[x + y*m_width].x;
		if(dir%4 == 1)sm = B_R_buf_pp_[x + y*m_width].y;
		if(dir%4 == 3)sm = B_R_buf_mm_[x + y*m_width].y;
	}
	dir = dir%4;
	REAL  max  = Smthnss[sm].max;
	REAL  lmbd   = Smthnss[sm].semi;
	//double aa  = 1. - al*dd; if(aa<0)aa=0;
	//REAL  lmbd = max*aa;
	if(dir ==0)
		if (x < m_width-1) {
			subMsgC( Di,&MsgC(x+1,y,2),  &MrgC(x,y), m_nLabels, 0.5);
			UpdMsgCL1( &MsgC(x,y,0),Di,m_nLabels, lmbd, max);	
		}
	if(dir ==1)
		if (y < m_height - 1) {
			subMsgC( Di,&MsgC(x,y+1,3),  &MrgC(x,y), m_nLabels, 0.5);
			UpdMsgCL1( &MsgC(x,y,1),Di,m_nLabels, lmbd, max);	
		}
	if(dir ==2)
		if (x > 0) {
			subMsgC( Di,&MsgC(x-1,y,0),  &MrgC(x,y), m_nLabels, 0.5);
			UpdMsgCL1( &MsgC(x,y,2),Di,m_nLabels, lmbd, max);	
		}
	if(dir ==3)
		if (y > 0)  {
			subMsgC( Di,&MsgC(x,y-1,1),  &MrgC(x,y), m_nLabels, 0.5);
			UpdMsgCL1( &MsgC(x,y,3),Di,m_nLabels, lmbd, max);	
		}
	//----------------------------------------
}

inline void  EDP::MkMrgL(int x, int y,EDP::REAL* Msg,EDP::REAL* Msg_add, EDP::REAL* Mrg)
{
	//int xinv = m_width-1 -x;
	int d_mx = (x + 1 <m_nLabels) ? x +1 : m_nLabels;
	if (x > 0) {
		if(x < m_nLabels)
			addMsg_(
				&XLD(x),
				&XLD(x-1),&Msg(y,0),
				&Mrg(y),
				d_mx-1,
				Msg_add(x-1,y,0));
		else
			addMsg(
				&XLD(x),
				&XLD(x-1),&Msg(y,0),
				&Mrg(y),
				d_mx);
	}// message (x-1,y)->(x,y)
	if (y > 0){
		addMsg(
		&XLD(x),
		&XLD(x),&Msg(y-1,1),
		&Mrg(y),
		d_mx);
	}// message (x,y-1)->(x,y)
	if (x < m_width-1) {
		addMsg(
		&XLD(x),
		&XLD(x +1),&Msg(y,2),
		&Mrg(y),
		d_mx);
	}// message (x+1,y)->(x,y)
	if (y < m_height-1) {
		addMsg(
		&XLD(x),
		&XLD(x), &Msg(y+1,3),
		&Mrg(y),
		d_mx);
	}// message (x,y+1)->(x,y)
}

inline void EDP::MkMrgC(int x, int y, EDP::REAL* Msg, EDP::REAL* Mrg) {
	if (x > 0)  {
		addMsgC( &MsgC(x-1,y,0),  &MrgC(x,y), m_nLabels); // message (x-1,y)->(x,y)
	}
	if (y > 0)  {
		addMsgC( &MsgC(x,y-1,1),  &MrgC(x,y), m_nLabels);// message (x,y-1)->(x,y)
	}
	if (x < m_width-1) {
		addMsgC( &MsgC(x+1,y,2),  &MrgC(x,y), m_nLabels);   // message (x+1,y)->(x,y)
	}
	if (y < m_height-1) {
		addMsgC( &MsgC(x,y+1,3),  &MrgC(x,y), m_nLabels);// message (x,y+1)->(x,y)
	}
}

inline void  EDP::addMsg_(int * xr, int * xsh,EDP::REAL* M, EDP::REAL* Mrg,int Km,EDP::REAL add)
{

    for (int k=0; k<Km; k++)
    	Mrg[xr[k]] += M[xsh[k]];
	Mrg[xr[Km]] += add;
}

inline void  EDP::subMsg_(int * xr,int * xsh, EDP::REAL* Di, EDP::REAL* M, EDP::REAL* Mrg,int Km,EDP::REAL add, EDP::REAL gamma)
{

    for (int k=0; k<Km; k++)Di[k]= gamma*Mrg[xr[k]] - M[xsh[k]];
	Di[Km]=gamma*Mrg[xr[Km]]-add;

}
inline void  EDP::subMsg(int * xr,int * xsh, EDP::REAL* Di, EDP::REAL* M, EDP::REAL* Mrg,int K, EDP::REAL gamma)
{

    for (int k=0; k<K; k++)Di[k]= gamma*Mrg[xr[k]] - M[xsh[k]];
	

}
inline void  EDP::subMsgC(EDP::REAL* Di, EDP::REAL* M, EDP::REAL* Mrg,int K, EDP::REAL gamma)
{
    for (int k=0; k<K; k++) {
    	Di[k]= gamma*Mrg[k] - M[k];
    }
}

/////////////////////////////////////
EDP::EDP(int nmIm, unsigned char ** i_ims, int width, int height, int nLabels ) {
	nm_Ims = nmIm;
	I_ims = i_ims;
	m_width       = width;
	m_height      = height;
	m_nLabels     = nLabels;
	m_nPixels     = width*height;
	initializeAlg();
}


EDP::~EDP()
{
	if( IG_ims[0]) delete [] IG_ims[0];
	if(IG_ims[1]) delete [] IG_ims[1];
	if(Mask_RL[0]) delete [] Mask_RL[0];
	if(Mask_RL[1]) delete [] Mask_RL[1];
	if(DspFL) delete [] DspFL;
    if(m_answer) delete[] m_answer;
	if(gtL_answer) delete[] gtL_answer;
	if (gtR_answer) delete[] gtR_answer;
	if(XL_tau) delete [] XL_tau;
	if(XR_tau) delete []  XR_tau;
	if(tau_XL_d) delete []  tau_XL_d;
	if(V_tau_bck) delete []  V_tau_bck ;
	if(R_tau_bck) delete []  R_tau_bck ;
	if(L_tau_bck) delete []  L_tau_bck ;
	if(V_tau_frw) delete []  V_tau_frw ;
	if(R_tau_frw) delete []  R_tau_frw ;
	if(L_tau_frw) delete []  L_tau_frw ;
	if(tau_XR_d) delete []  tau_XR_d;
	if(T_max_y) delete []  T_max_y;
	if (Sltn_crv) delete []  Sltn_crv ;
    if (B_L_buf_pp_) delete []  B_L_buf_pp_;
	if (B_L_buf_mm_) delete []  B_L_buf_mm_;
	if(B_R_buf_pp_) delete []  B_R_buf_pp_;
	if(B_R_buf_mm_) delete []  B_R_buf_mm_;
	if(Mrg_bf) delete [] Mrg_bf;
    if(Dst_B) delete [] Dst_B;
    if ( m_D ) delete [] m_D;
    if ( m_messages ) delete [] m_messages;
	if(fi_erf_tab) delete [] fi_erf_tab;
}






void EDP::initializeAlg()
{
	m_D = new CostVal[m_nPixels*m_nLabels];
	bufcc = new float [m_nPixels*3];///
	Itr_num = 0;
	WpH = m_width*2;
	T_max_y =          new int [m_height] ;
	this->Sltn_crv  =  new SLTN_CRV [WpH*m_height];
	//this->m_D_pnl =  new float [m_nPixels*m_nLabels];
	B_L_buf_pp_  =  new POINT_CH [m_nPixels];
	B_L_buf_mm_ = new POINT_CH [m_nPixels];
	B_R_buf_pp_  =  new POINT_CH [m_nPixels];
	B_R_buf_mm_ = new POINT_CH [m_nPixels];
	
	/////////// param add //////////////
	Trn_n = 3;
	float alp = 1;
	C_app_1 = 3;
	C_app_2  = 0.2;
	C_dsp_1  =   alp*30;
	C_dsp_2   = 30;
	C_gr_1      = 2;
	C_gr_2       =1;
	Gr_prc =0.5;
	int br_tms =1; // brich Tomasy
	int Sqrt_d =1;
	int chTr = 16;
	float mean[3] = {0, 0, 0};
	for(int p = 0; p < m_nPixels; p++)  {
		for(int c =0; c <3; c++) {
			mean[c] += I_ims[0][p + c*m_nPixels];
		}
	}
	
	for(int c =0; c <3; c++) {
		mean[c] /= m_nPixels;
	}
	
	sigma_img =0;
	for(int p = 0; p < m_nPixels; p++){
		float sum =0;
		for(int c =0; c <3; c++){
			float vl = I_ims[0][p +c*m_nPixels]- mean[c]; 
			sum += vl*vl;
		}
		sigma_img +=sum;
	}
	sigma_img = sqrt(sigma_img/m_nPixels);

	//////////////////////////////////////
	//med_fl( 1, 7, I_ims[0] );med_fl( 1, 7,  I_ims[1] );
	{ //---------------------------------------------------------------
		XL_XR_Tau(); //grid calculation
	} //---------------------------------------------------------------
	Mrg_bf = new REAL [Tau]; // working buffer fo marginal
	Dst_B  = new int [Tau];
	DspFL = new float [m_nPixels];
	Mask_RL[0]= new unsigned char [m_nPixels];
	Mask_RL[1] = new unsigned char [m_nPixels];
	IG_ims[0]  = new unsigned char [m_nPixels*6];
	IG_ims[1]  = new unsigned char [m_nPixels*6];
	GrBuf();

}


void EDP::findSlt(float * dsp) {
	int* DspM = new int [m_nPixels];
	m_answer = new int [m_nPixels];
	/////////////////////////////////////


	{ /////////////////// COST CLC
		int Sqrt_d =0; int br_tms =1;  chTr = 30;
		float * c_D = new float [m_height*Tau];
		float avmin =0;
		if(0) { // Not executed
			computeDSI_tau( I_ims[0],  I_ims[1],  br_tms,  Sqrt_d, chTr);
			avmin   = min_mean();
			memcpy(c_D, m_D, m_height*Tau*sizeof(float));
		}
		for(int i =0; i < m_height*Tau; i++) {
			m_D[i]=0;
		}
		
		Cost_mean = computeDSI_tauG( IG_ims[0],  IG_ims[1],  0,  Sqrt_d, chTr*4, 1);
		float avmingr  = min_mean();
		float gmul  = avmin/avmingr*3.5;
		
		if(0) { // not executed.
			for(int i =0; i < m_height*Tau; i++) {
				m_D[i] = c_D[i] + gmul*m_D[i] ;
			}
		}
		float thr = 0.00095;
		transform_D_( thr);
		delete [] c_D;
	}///////////////////////////////

	int itr = (N_LB<30) ? 5:8;
	png::image< png::ga_pixel > timg(truthFile);
	#if defined(HYBRID)
	hybrid_MRF__z(itr, m_answer, DspM, timg);
	#else
	modified_MRF__z( itr, m_answer, DspM, timg);
	#endif
	for(int p = 0 ; p < m_nPixels; p++) {
		dsp[p] = 0.25* m_answer[p]; // 0.25*m_answer[p]
	}
    delete [] DspM;
}

inline void EDP::SampledFlt_Add_new(float sigma_c_mul, float sigmaxy, int cl_Q, unsigned char * I_b, double *buf_to_flt, short int * K_ind,  EDP::COLOR_F * K_cls,  int nb_ch ) {
	///////////////////////////
	double sigc = sigma_c_mul*255;
	sigmaX = (m_width%sigmaix)? m_width/sigmaix+2 :m_width/sigmaix+1;
	sigmaY = (m_height%sigmaiy)? m_height/sigmaiy+2 :m_height/sigmaiy+1;
	int size1 = sigmaX*sigmaY;
	int size2 = size1*cl_Q;
	short int * map_nb = new short int [size2];
	short int * map_nbi = new short int [size2];
	double * sum_vl;
	double * sumf_vl;
	double ** sum_vl_pt;
	double ** sumf_vl_pt;
	double * sum_wt = new double [size2];
	double * sumf_wt = new double [size2];
	int * cnt_p = new int [size1];
	double mulg1 =exp(-0.5);
	double mulg2 =exp(-2.);


	//int div = 256/r;
	int n_c = 1000;
	int cc[3] ;
	double ccf[3];
	double * gss_c_wt = new double [n_c];
	double * ijK = new double [cl_Q*cl_Q];
	int *sigma_cnt = new int [cl_Q];
	double **i_sgm_pnt = new double *[cl_Q];
	int **i_q_pnt = new int * [cl_Q];
	short int * q_map = K_ind;
	double nsigm = 3;
	get_std( sigma_c_mul, I_b,  gss_c_wt,  n_c );
	//-------MAP-----------------

	//---------------------------
	int full_cnt =0;
	for(int i = 0; i < cl_Q; i++) {
		for(int j = 0; j< cl_Q; j++) {
			double dst = 0; 
			for(int c =0; c<3; c++) {
				double r = (K_cls[i].cl[c] - K_cls[j].cl[c]); 
				dst += r*r;
			} 
			dst = sqrt(dst); 
			double m = ijK[i+ j*cl_Q] = ( dst/sigc>nsigm)? 0:  gss_c_wt[round_fl(dst)]; 
			if(m!=0)
				full_cnt++;
		}
	}
	double * i_sgm = new double [full_cnt];
	int * i_q = new int [full_cnt];
	//----------------------------------------------------------------------------------------
	for(int i = 0; i < cl_Q; i++){
		if(!i){
			i_sgm_pnt[i] = i_sgm; 
			i_q_pnt[i] = i_q;
		}else{
			i_sgm_pnt[i] = &i_sgm_pnt[i-1][sigma_cnt[i-1]]; 
			i_q_pnt[i] =  &i_q_pnt[i-1][sigma_cnt[i-1]];
		}
		sigma_cnt[i] =0; 
		for(int j = 0; j< cl_Q; j++) {
			if(ijK[i+ j*cl_Q]) {
				i_sgm_pnt[i][sigma_cnt[i]] = ijK[i+ j*cl_Q];  
				i_q_pnt[i][sigma_cnt[i]] = j;
				sigma_cnt[i]++;
				/*  int a =0;*/
			}
		}
	}

	////////////////////INI filter /////////////////////////
	for(int p =0; p <size2; p++) {
		map_nb[p] = -1;
		sum_wt[p]  = 0;
		sumf_wt[p] = 0;
	}
	for(int p =0; p <size1; p++) 
		cnt_p[p] = 0;
	//---------- accumulation of f ------------------------------
	//--- wt
	for(int y=0; y<m_height; y++)
	for(int x = 0; x < m_width; x++) {
		int p = x + y*m_width;
		int xx0 = x/sigmaix, xx1 = xx0 +1; int xxm = x%sigmaix;
		int yy0 = (y/sigmaiy); int yy1 = yy0+1; int yym = y%sigmaiy;
		int pp[4] = { xx0+ sigmaX*yy0,  xx1+ sigmaX*yy0,  xx0+ sigmaX*yy1,  xx1+ sigmaX*yy1};
		int q = q_map[p];
		//for(int n =0; n <sigma_cnt[q]; n++){////////////////////////////
		/*int qn = i_q_pnt[q][n];*/
		for(int i =0; i< 4; i++)
		//int i = TBI_p0(xxm,yym);
			if(map_nb[pp[i]*cl_Q+q] == -1){
				map_nbi[pp[i]*cl_Q+cnt_p[pp[i]]] = q;  
				map_nb[pp[i]*cl_Q+q] = cnt_p[pp[i]]++;
			}//}//////////////////////////
		/* for(int i =0; i< 4; i++) { int nbq = map_nb[pp[i]*cl_Q+q]; sum_wt[pp[i]*cl_Q + nbq] += BI_WS(xxm, yym, i); }*/
	}
	//---- func
	int cnt_big = 0; 
	for(int p = 0 ; p < size1; p++) 
		cnt_big += cnt_p[p];

	sum_vl = new double [cnt_big*nb_ch];
	sumf_vl = new double [cnt_big*nb_ch];
	for(int i = 0 ; i <cnt_big*nb_ch; i++)
		sum_vl[i] =  0;
	sum_vl_pt = new double* [size1]; 
	sumf_vl_pt = new double* [size1];
	sum_vl_pt[0] = sum_vl; 
	sumf_vl_pt[0] = sumf_vl;

	for(int p = 1; p < size1; p++) {
		sum_vl_pt[p] = sum_vl_pt[p-1] + cnt_p[p-1]*nb_ch;
		sumf_vl_pt[p] = sumf_vl_pt[p-1] + cnt_p[p-1]*nb_ch;
	}

	for(int y=0; y<m_height; y++)
	for(int x = 0; x < m_width; x++) {
		int p = x + y*m_width;
		int xx0 = x/sigmaix, xx1 = xx0 +1; 
		int xxm = x%sigmaix;
		
		int yy0 = (y/sigmaiy); 
		int yy1 = yy0+1; 
		int yym = y%sigmaiy;
		int pp[4] =	{ xx0+ sigmaX*yy0,  xx1+ sigmaX*yy0,
						xx0+ sigmaX*yy1,  xx1+ sigmaX*yy1};
		int q = q_map[p];

		//for(int n =0; n <sigma_cnt[q]; n++){////////////////////////////
		//	int qn = i_q_pnt[q][n]; float wc = i_sgm_pnt[q][qn];
		for(int i =0; i< 4; i++) {
			int cntp = cnt_p[pp[i]];
			int nbq = map_nb[pp[i]*cl_Q+q]; 
			double wt = BI_WSx(xxm, yym, i);
			sum_wt[pp[i]*cl_Q + nbq] += wt;
			for(int ch = 0; ch <nb_ch; ch++){
				double vl = buf_to_flt[p + ch*m_nPixels];
				sum_vl_pt[pp[i]][nbq + cntp*ch] += wt*vl;
			}
		}
		//}///////////////////////
	}
	//-------------------------- fiter firs ---------------
	double * vlmc = new double [nb_ch];
	double * vlpc =  new double [nb_ch];
	double * vlmc2 = new double [nb_ch];
	double * vlpc2 =  new double [nb_ch];
	//-----------YY
	for(int y = 0; y < sigmaY; y ++)
	for(int x = 0; x < sigmaX; x++) {

		int p = x + sigmaX*y; 
		int cntp =  cnt_p[p]; 
		int pm = (y - 1)*sigmaX + x; 
		int pp = (y + 1)*sigmaX + x; 
		int pm2 = (y - 2)*sigmaX + x; 
		int pp2 = (y + 2)*sigmaX + x;
		for(int ci =0; ci < cntp; ci++ ) {
			int q_ci = map_nbi[p*cl_Q+ci];   
			double vlm=0, vlp =0, vlm2 =0, vlp2 =0;  
			for(int ch = 0; ch <nb_ch; ch++){
				vlmc[ch] = 0; 
				vlpc[ch] = 0; 
				vlmc2[ch] = 0; 
				vlpc2[ch] = 0;
			}
			if(y-1>=0) {
				int pm_ci = map_nb[pm*cl_Q+q_ci];
				if(pm_ci != -1) {
					vlm = sum_wt[pm*cl_Q + pm_ci];
					for(int ch = 0; ch <nb_ch; ch++)
						vlmc[ch] = sum_vl_pt[pm][pm_ci + cnt_p[pm]*ch]; 
				}
			}

			if(y+1< sigmaY) {
				int pp_ci = map_nb[pp*cl_Q+q_ci]; 
				if(pp_ci != -1){
					vlp = sum_wt[pp*cl_Q + pp_ci] ;  
					for(int ch = 0; ch <nb_ch; ch++)
						vlpc[ch] = sum_vl_pt[pp][pp_ci + cnt_p[pp]*ch]; 
				}
			}

			if(y-2>=0) {int pm_ci2 = map_nb[pm2*cl_Q+q_ci]; if(pm_ci2 != -1){vlm2 =     sum_wt[pm2*cl_Q + pm_ci2] ;  for(int ch = 0; ch <nb_ch; ch++)vlmc2[ch] = sum_vl_pt[pm2][pm_ci2 + cnt_p[pm2]*ch]; }}

			if(y+2< sigmaY) {int pp_ci2 = map_nb[pp2*cl_Q+q_ci]; if(pp_ci2 != -1){vlp2 =     sum_wt[pp2*cl_Q + pp_ci2] ;  for(int ch = 0; ch <nb_ch; ch++)vlpc2[ch] = sum_vl_pt[pp2][pp_ci2 + cnt_p[pp2]*ch]; }}

			sumf_wt[p*cl_Q + ci] = sum_wt[p*cl_Q + ci] + vlp*mulg1 + vlm*mulg1 + vlp2*mulg2 + vlm2*mulg2;
			for(int ch = 0; ch <nb_ch; ch++) 
				sumf_vl_pt[p][ci + cnt_p[p]*ch] = sum_vl_pt[p][ci + cnt_p[p]*ch] + vlpc[ch]*mulg1 + vlmc[ch]*mulg1 + vlpc2[ch]*mulg2 + vlmc2[ch]*mulg2;
		}
	}
	for(int i = 0 ; i <cnt_big*nb_ch; i++)
		sum_vl[i] =  sumf_vl[i];
	for(int p =0; p <size2; p++)  
		sum_wt[p]  =   sumf_wt[p];
	//-----------XX
	for(int y = 0; y < sigmaY; y ++)
	for(int x = 0; x < sigmaX; x++) {

		int p = x + sigmaX*y; 
		int cntp =  cnt_p[p]; 
		int pm = (y )*sigmaX + x -1; 
		int pp = (y )*sigmaX + x +1; 
		int pm2 = (y )*sigmaX + x -2; 
		int pp2 = (y )*sigmaX + x +2;
		
		for(int ci =0; ci < cntp; ci++ ) {
			int q_ci = map_nbi[p*cl_Q+ci];  
			double vlm=0, vlp =0, vlm2 =0, vlp2 =0;  
			for(int ch = 0; ch <nb_ch; ch++){
				vlmc[ch] = 0; vlpc[ch] = 0; vlmc2[ch] = 0; vlpc2[ch] = 0;
			}

			if(x-1>=0) {int pm_ci = map_nb[pm*cl_Q+q_ci];  if(pm_ci != -1){vlm =     sum_wt[pm*cl_Q + pm_ci] ;  for(int ch = 0; ch <nb_ch; ch++)vlmc[ch] = sum_vl_pt[pm][pm_ci + cnt_p[pm]*ch]; }}
			if(x+1< sigmaX) {int pp_ci = map_nb[pp*cl_Q+q_ci]; if(pp_ci != -1){vlp =     sum_wt[pp*cl_Q + pp_ci] ;  for(int ch = 0; ch <nb_ch; ch++)vlpc[ch] = sum_vl_pt[pp][pp_ci + cnt_p[pp]*ch]; }}
			if(x-2>=0) {int pm_ci2 = map_nb[pm2*cl_Q+q_ci];  if(pm_ci2 != -1){vlm2 =     sum_wt[pm2*cl_Q + pm_ci2] ;  for(int ch = 0; ch <nb_ch; ch++)vlmc2[ch] = sum_vl_pt[pm2][pm_ci2 + cnt_p[pm2]*ch]; }}
			if(x+2< sigmaX) {int pp_ci2 = map_nb[pp2*cl_Q+q_ci]; if(pp_ci2 != -1){vlp2 =     sum_wt[pp2*cl_Q + pp_ci2] ;  for(int ch = 0; ch <nb_ch; ch++)vlpc2[ch] = sum_vl_pt[pp2][pp_ci2 + cnt_p[pp2]*ch]; }}
			sumf_wt[p*cl_Q + ci] = sum_wt[p*cl_Q + ci] + vlp*mulg1 + vlm*mulg1 + vlp2*mulg2 + vlm2*mulg2;
			for(int ch = 0; ch <nb_ch; ch++) 
				sumf_vl_pt[p][ci + cnt_p[p]*ch] = sum_vl_pt[p][ci + cnt_p[p]*ch] + vlpc[ch]*mulg1 + vlmc[ch]*mulg1 + vlpc2[ch]*mulg2 + vlmc2[ch]*mulg2;
		}
	}
	for(int i = 0 ; i <cnt_big*nb_ch; i++)
		sum_vl[i] =  sumf_vl[i];
	for(int p =0; p <size2; p++)  
		sum_wt[p]  =   sumf_wt[p];
	//-----------cc
	for(int y = 0; y < sigmaY; y ++)
	for(int x = 0; x < sigmaX; x++) {

		int p = x + sigmaX*y; 
		int cntp =  cnt_p[p];
		for(int ci =0; ci < cntp; ci++ ) {
			int q_ci = map_nbi[p*cl_Q+ci];  sumf_wt[p*cl_Q + ci] = 0; for(int ch = 0; ch <nb_ch; ch++)sumf_vl_pt[p][ci + cnt_p[p]*ch]  = 0;
			for(int n =0; n <sigma_cnt[q_ci]; n++) {
			//////////////////////
				int qn = i_q_pnt[q_ci][n]; 
				int in_vl = map_nb[p*cl_Q+qn];
				if(in_vl != -1) {
					float w = i_sgm_pnt[q_ci][n];  
					sumf_wt[p*cl_Q + ci] += w*sum_wt[p*cl_Q + in_vl];
					for(int ch = 0; ch <nb_ch; ch++)
						sumf_vl_pt[p][ci + cnt_p[p]*ch] +=w*sum_vl_pt[p][in_vl + cnt_p[p]*ch];
				}
			}////////////////////////
			int a = 0;
			//for(int ch = 0; ch <nb_ch; ch++)sumf_vl_pt[p][ci + cntp*ch]  /=  (sumf_wt[p*cl_Q + ci])? sumf_wt[p*cl_Q + ci]:1;
		}
		int aa =0;
	}
	//-------------------- back to image
	for(int y=0; y<m_height; y++)
	for(int x = 0; x < m_width; x++) {
		int p = x + y*m_width;
		int xx0 = x/sigmaix, xx1 = xx0 +1; 
		int xxm = x%sigmaix;
		int yy0 = (y/sigmaiy); 
		int yy1 = yy0+1; 
		int yym = y%sigmaiy;
		int pp[4] = { xx0+ sigmaX*yy0,  xx1+ sigmaX*yy0,  xx0+ sigmaX*yy1,  xx1+ sigmaX*yy1};
		int q = q_map[p];
		for(int ch = 0; ch <nb_ch; ch++)  
			buf_to_flt[p + ch*m_nPixels]= 0;
		double wtsm = 0;
		//for(int n =0; n <sigma_cnt[q]; n++){////////////////////////////
		//int qn = i_q_pnt[q][n]; float wc = i_sgm_pnt[q][qn];
		for(int i =0; i< 4; i++) {
			int cntp = cnt_p[pp[i]];
			int nbq = map_nb[pp[i]*cl_Q+q];
			if(nbq != -1){ 
				double wt = BI_WSx(xxm, yym, i);
				wtsm += 	 sumf_wt[pp[i]*cl_Q + nbq]*wt;
				for(int ch = 0; ch <nb_ch; ch++)  
					buf_to_flt[p + ch*m_nPixels] += sumf_vl_pt[pp[i]][nbq+ cnt_p[pp[i]]*ch]*wt;
			}
		}
		for(int ch = 0; ch <nb_ch; ch++){
			buf_to_flt[p + ch*m_nPixels] /= (wtsm)? wtsm:1;		///  sumf_vl_pt[pp[i]][nbq+ cnt_p[pp[i]]*ch]*wt;
		}

	}
	delete [] vlmc;
	delete [] vlpc;
	delete [] vlmc2;
	delete [] vlpc2;

	delete [] gss_c_wt;
	delete [] ijK;
	/* delete [] q_map;*/
	delete [] sigma_cnt;
	delete [] i_sgm_pnt;
	delete [] i_q_pnt;
	delete [] i_sgm;
	delete [] i_q;
	delete [] map_nb;

	delete [] sum_vl;
	delete [] sumf_vl;
	delete [] sumf_vl_pt;
	delete [] sum_vl_pt;

	delete [] sum_wt;
	delete [] sumf_wt;
	delete [] cnt_p;
	delete [] map_nbi;
}////////////////////////////////////

	
inline void EDP::SampledFlt_Add(float sigma_c_mul, float sigmaxy, int cl_Q, unsigned char * I_b, double *buf_to_flt, short int * clr_cls,  EDP::COLOR_F * K_cls, int r, int nb_ch )
{ ///////////////////////////


	  double sigc = sigma_c_mul*255;
	  sigmaX = (m_width%sigmaix)? m_width/sigmaix+2 :m_width/sigmaix+1;
	  sigmaY = (m_height%sigmaiy)? m_height/sigmaiy+2 :m_height/sigmaiy+1;
	  int size1 = sigmaX*sigmaY;
	  int size2 = size1*cl_Q;
	  short int * map_nb = new short int [size2];
	  short int * map_nbi = new short int [size2];
	  double * sum_vl;
	   double * sumf_vl;
	   double ** sum_vl_pt;
	   double ** sumf_vl_pt;
	   double * sum_wt = new double [size2];
	    double * sumf_wt = new double [size2];
		int * cnt_p = new int [size1];
		double mulg1 =exp(-0.5);
		double mulg2 =exp(-2.);


	  int div = 256/r;
	  int n_c = 1000;
	  int cc[3] ;
	  double ccf[3];
	  double * gss_c_wt = new double [n_c];
	  double * ijK = new double [cl_Q*cl_Q];
	  int *sigma_cnt = new int [cl_Q];
	  double **i_sgm_pnt = new double *[cl_Q];
	  int **i_q_pnt = new int * [cl_Q];
	  short int * q_map = new short int [m_nPixels];
	  double nsigm = 3;
	  get_std( sigma_c_mul, I_b,  gss_c_wt,  n_c );
	  //-------MAP-----------------
						  for(int p = 0; p < m_nPixels; p++)
						  {
							  for(int c =0; c < 3; c ++){cc[c] = CUTFL((float)I_b[p + c*m_nPixels]/div, r); }
							  q_map[p] =  clr_cls[IND_HC(cc[0], cc[1], cc[2], r)];							
						  }
//------------------------------------------------------------------------------------------
						  int full_cnt =0;for(int i = 0; i < cl_Q; i++)for(int j = 0; j< cl_Q; j++)
						  {double dst = 0; for(int c =0; c<3; c++)
						  {double r = (K_cls[i].cl[c] - K_cls[j].cl[c]); dst += r*r;} dst = sqrt(dst); double m = ijK[i+ j*cl_Q] = ( dst/sigc>nsigm)? 0:  gss_c_wt[round_fl(dst)]; if(m!=0)full_cnt++;}
double * i_sgm = new double [full_cnt];
int * i_q = new int [full_cnt];
						  //----------------------------------------------------------------------------------------
						  for(int i = 0; i < cl_Q; i++)
						  {   if(!i){i_sgm_pnt[i] = i_sgm; i_q_pnt[i] = i_q;}
						       else
							   {
								   i_sgm_pnt[i] = &i_sgm_pnt[i-1][sigma_cnt[i-1]]; i_q_pnt[i] =  &i_q_pnt[i-1][sigma_cnt[i-1]];
						  }
							   sigma_cnt[i] =0; for(int j = 0; j< cl_Q; j++)
							   {
								   if(ijK[i+ j*cl_Q])
								   { i_sgm_pnt[i][sigma_cnt[i]] = ijK[i+ j*cl_Q];  i_q_pnt[i][sigma_cnt[i]] = j;
								   sigma_cnt[i]++;
								 /*  int a =0;*/
								   }
							   }
						  }

////////////////////INI filter /////////////////////////
  for(int p =0; p <size2; p++)
  {
	  	   map_nb[p] = -1;
	       sum_wt[p]  = 0;
	       sumf_wt[p] = 0;
  }
  for(int p =0; p <size1; p++) cnt_p[p] = 0;
  //---------- accumulation of f ------------------------------
  //--- wt
 for(int y=0; y<m_height; y++)
 for(int x = 0; x < m_width; x++)
 {
	 int p = x + y*m_width;
	 int xx0 = x/sigmaix, xx1 = xx0 +1; int xxm = x%sigmaix;
	 int yy0 = (y/sigmaiy); int yy1 = yy0+1; int yym = y%sigmaiy;
	 int pp[4] = { xx0+ sigmaX*yy0,  xx1+ sigmaX*yy0,  xx0+ sigmaX*yy1,  xx1+ sigmaX*yy1};
     int q = q_map[p];
	 	 //for(int n =0; n <sigma_cnt[q]; n++){////////////////////////////
			/*int qn = i_q_pnt[q][n];*/
			for(int i =0; i< 4; i++)
	        //int i = TBI_p0(xxm,yym);
			if(map_nb[pp[i]*cl_Q+q] == -1){map_nbi[pp[i]*cl_Q+cnt_p[pp[i]]] = q;  map_nb[pp[i]*cl_Q+q] = cnt_p[pp[i]]++;}//}//////////////////////////
	/* for(int i =0; i< 4; i++) { int nbq = map_nb[pp[i]*cl_Q+q]; sum_wt[pp[i]*cl_Q + nbq] += BI_WS(xxm, yym, i); }*/
 }
 //---- func
 int cnt_big = 0; for(int p = 0 ; p < size1; p++) cnt_big += cnt_p[p];

   sum_vl = new double [cnt_big*nb_ch]; sumf_vl = new double [cnt_big*nb_ch];    for(int i = 0 ; i <cnt_big*nb_ch; i++)sum_vl[i] =  0;
   sum_vl_pt = new double* [size1]; sumf_vl_pt = new double* [size1];
   sum_vl_pt[0] = sum_vl; sumf_vl_pt[0] = sumf_vl;

   for(int p = 1; p < size1; p++)
   {
	   sum_vl_pt[p] = sum_vl_pt[p-1] + cnt_p[p-1]*nb_ch;
	   sumf_vl_pt[p] = sumf_vl_pt[p-1] + cnt_p[p-1]*nb_ch;
   }

 for(int y=0; y<m_height; y++)
 for(int x = 0; x < m_width; x++)
 {
	 int p = x + y*m_width;
	 int xx0 = x/sigmaix, xx1 = xx0 +1; int xxm = x%sigmaix;
	 int yy0 = (y/sigmaiy); int yy1 = yy0+1; int yym = y%sigmaiy;
	 int pp[4] =
	 { xx0+ sigmaX*yy0,  xx1+ sigmaX*yy0,
	 xx0+ sigmaX*yy1,  xx1+ sigmaX*yy1};
     int q = q_map[p];
	
	 //for(int n =0; n <sigma_cnt[q]; n++){////////////////////////////
		//	int qn = i_q_pnt[q][n]; float wc = i_sgm_pnt[q][qn];
	 for(int i =0; i< 4; i++)
	/* int i = TBI_p0(xxm,yym); */
	 {
		 int cntp = cnt_p[pp[i]];
		 int nbq = map_nb[pp[i]*cl_Q+q]; double wt = BI_WSx(xxm, yym, i);
		  sum_wt[pp[i]*cl_Q + nbq] += wt;
		 for(int ch = 0; ch <nb_ch; ch++){
		 double vl = buf_to_flt[p + ch*m_nPixels];
		 sum_vl_pt[pp[i]][nbq + cntp*ch] += wt*vl;
		 }
	 }
	 //}///////////////////////
 }
//-------------------------- fiter firs ---------------
  double * vlmc = new double [nb_ch];
  double * vlpc =  new double [nb_ch];
    double * vlmc2 = new double [nb_ch];
  double * vlpc2 =  new double [nb_ch];
  //-----------YY
   for(int y = 0; y < sigmaY; y ++)
   for(int x = 0; x < sigmaX; x++)
   {

		int p = x + sigmaX*y; int cntp =  cnt_p[p]; int pm = (y - 1)*sigmaX + x; int pp = (y + 1)*sigmaX + x; int pm2 = (y - 2)*sigmaX + x; int pp2 = (y + 2)*sigmaX + x;
		for(int ci =0; ci < cntp; ci++ )
		{
		int q_ci = map_nbi[p*cl_Q+ci];   double vlm=0, vlp =0, vlm2 =0, vlp2 =0;  for(int ch = 0; ch <nb_ch; ch++){vlmc[ch] = 0; vlpc[ch] = 0; vlmc2[ch] = 0; vlpc2[ch] = 0;}

		if(y-1>=0)
		{int pm_ci = map_nb[pm*cl_Q+q_ci];
		if(pm_ci != -1)
		{vlm =     sum_wt[pm*cl_Q + pm_ci] ;
		for(int ch = 0; ch <nb_ch; ch++)
			vlmc[ch] = sum_vl_pt[pm][pm_ci + cnt_p[pm]*ch]; }}

		   if(y+1< sigmaY) {int pp_ci = map_nb[pp*cl_Q+q_ci]; if(pp_ci != -1){vlp =     sum_wt[pp*cl_Q + pp_ci] ;  for(int ch = 0; ch <nb_ch; ch++)vlpc[ch] = sum_vl_pt[pp][pp_ci + cnt_p[pp]*ch]; }}
		
		   if(y-2>=0) {int pm_ci2 = map_nb[pm2*cl_Q+q_ci]; if(pm_ci2 != -1){vlm2 =     sum_wt[pm2*cl_Q + pm_ci2] ;  for(int ch = 0; ch <nb_ch; ch++)vlmc2[ch] = sum_vl_pt[pm2][pm_ci2 + cnt_p[pm2]*ch]; }}

		   if(y+2< sigmaY) {int pp_ci2 = map_nb[pp2*cl_Q+q_ci]; if(pp_ci2 != -1){vlp2 =     sum_wt[pp2*cl_Q + pp_ci2] ;  for(int ch = 0; ch <nb_ch; ch++)vlpc2[ch] = sum_vl_pt[pp2][pp_ci2 + cnt_p[pp2]*ch]; }}

	   sumf_wt[p*cl_Q + ci] = sum_wt[p*cl_Q + ci] + vlp*mulg1 + vlm*mulg1 + vlp2*mulg2 + vlm2*mulg2;
	   for(int ch = 0; ch <nb_ch; ch++) sumf_vl_pt[p][ci + cnt_p[p]*ch] = sum_vl_pt[p][ci + cnt_p[p]*ch] + vlpc[ch]*mulg1 + vlmc[ch]*mulg1 + vlpc2[ch]*mulg2 + vlmc2[ch]*mulg2;
		}
   }
    for(int i = 0 ; i <cnt_big*nb_ch; i++)sum_vl[i] =  sumf_vl[i];
	for(int p =0; p <size2; p++)  sum_wt[p]  =   sumf_wt[p];
	  //-----------XX
   for(int y = 0; y < sigmaY; y ++)
   for(int x = 0; x < sigmaX; x++)
   {

		int p = x + sigmaX*y; int cntp =  cnt_p[p]; int pm = (y )*sigmaX + x -1; int pp = (y )*sigmaX + x +1; int pm2 = (y )*sigmaX + x -2; int pp2 = (y )*sigmaX + x +2;
		for(int ci =0; ci < cntp; ci++ )
		{
			int q_ci = map_nbi[p*cl_Q+ci];  double vlm=0, vlp =0, vlm2 =0, vlp2 =0;  for(int ch = 0; ch <nb_ch; ch++){vlmc[ch] = 0; vlpc[ch] = 0; vlmc2[ch] = 0; vlpc2[ch] = 0;}

		if(x-1>=0) {int pm_ci = map_nb[pm*cl_Q+q_ci];  if(pm_ci != -1){vlm =     sum_wt[pm*cl_Q + pm_ci] ;  for(int ch = 0; ch <nb_ch; ch++)vlmc[ch] = sum_vl_pt[pm][pm_ci + cnt_p[pm]*ch]; }}
		   if(x+1< sigmaX) {int pp_ci = map_nb[pp*cl_Q+q_ci]; if(pp_ci != -1){vlp =     sum_wt[pp*cl_Q + pp_ci] ;  for(int ch = 0; ch <nb_ch; ch++)vlpc[ch] = sum_vl_pt[pp][pp_ci + cnt_p[pp]*ch]; }}
		if(x-2>=0) {int pm_ci2 = map_nb[pm2*cl_Q+q_ci];  if(pm_ci2 != -1){vlm2 =     sum_wt[pm2*cl_Q + pm_ci2] ;  for(int ch = 0; ch <nb_ch; ch++)vlmc2[ch] = sum_vl_pt[pm2][pm_ci2 + cnt_p[pm2]*ch]; }}
		   if(x+2< sigmaX) {int pp_ci2 = map_nb[pp2*cl_Q+q_ci]; if(pp_ci2 != -1){vlp2 =     sum_wt[pp2*cl_Q + pp_ci2] ;  for(int ch = 0; ch <nb_ch; ch++)vlpc2[ch] = sum_vl_pt[pp2][pp_ci2 + cnt_p[pp2]*ch]; }}
	   sumf_wt[p*cl_Q + ci] = sum_wt[p*cl_Q + ci] + vlp*mulg1 + vlm*mulg1 + vlp2*mulg2 + vlm2*mulg2;
	   for(int ch = 0; ch <nb_ch; ch++) sumf_vl_pt[p][ci + cnt_p[p]*ch] = sum_vl_pt[p][ci + cnt_p[p]*ch] + vlpc[ch]*mulg1 + vlmc[ch]*mulg1 + vlpc2[ch]*mulg2 + vlmc2[ch]*mulg2;
		}
   }
    for(int i = 0 ; i <cnt_big*nb_ch; i++)sum_vl[i] =  sumf_vl[i];
	for(int p =0; p <size2; p++)  sum_wt[p]  =   sumf_wt[p];
	  //-----------cc
   for(int y = 0; y < sigmaY; y ++)
   for(int x = 0; x < sigmaX; x++)
   {

		int p = x + sigmaX*y; int cntp =  cnt_p[p];
		for(int ci =0; ci < cntp; ci++ )
		{
		int q_ci = map_nbi[p*cl_Q+ci];  sumf_wt[p*cl_Q + ci] = 0; for(int ch = 0; ch <nb_ch; ch++)sumf_vl_pt[p][ci + cnt_p[p]*ch]  = 0;
		for(int n =0; n <sigma_cnt[q_ci]; n++)
		{//////////////////////
			int qn = i_q_pnt[q_ci][n]; int in_vl = map_nb[p*cl_Q+qn];
			if(in_vl != -1)
			{
				float w = i_sgm_pnt[q_ci][n];  sumf_wt[p*cl_Q + ci] += w*sum_wt[p*cl_Q + in_vl];
				for(int ch = 0; ch <nb_ch; ch++)sumf_vl_pt[p][ci + cnt_p[p]*ch] +=w*sum_vl_pt[p][in_vl + cnt_p[p]*ch];
			}
		}////////////////////////
		int a = 0;
		//for(int ch = 0; ch <nb_ch; ch++)sumf_vl_pt[p][ci + cntp*ch]  /=  (sumf_wt[p*cl_Q + ci])? sumf_wt[p*cl_Q + ci]:1;
		}
		int aa =0;

   }
   //-------------------- back to image
 for(int y=0; y<m_height; y++)
 for(int x = 0; x < m_width; x++)
 {
	 int p = x + y*m_width;
	 int xx0 = x/sigmaix, xx1 = xx0 +1; int xxm = x%sigmaix;
	 int yy0 = (y/sigmaiy); int yy1 = yy0+1; int yym = y%sigmaiy;
	 int pp[4] = { xx0+ sigmaX*yy0,  xx1+ sigmaX*yy0,  xx0+ sigmaX*yy1,  xx1+ sigmaX*yy1};
     int q = q_map[p];
	 for(int ch = 0; ch <nb_ch; ch++)  buf_to_flt[p + ch*m_nPixels]= 0;
	 double wtsm = 0;
	  //for(int n =0; n <sigma_cnt[q]; n++){////////////////////////////
			//int qn = i_q_pnt[q][n]; float wc = i_sgm_pnt[q][qn];
	 for(int i =0; i< 4; i++)
	 {
		 int cntp = cnt_p[pp[i]];
		 int nbq = map_nb[pp[i]*cl_Q+q];
		 if(nbq != -1){ double wt =   BI_WSx(xxm, yym, i);
		 wtsm += 	 sumf_wt[pp[i]*cl_Q + nbq]*wt;
		 for(int ch = 0; ch <nb_ch; ch++)  buf_to_flt[p + ch*m_nPixels] += sumf_vl_pt[pp[i]][nbq+ cnt_p[pp[i]]*ch]*wt;}
	 }
	for(int ch = 0; ch <nb_ch; ch++)  buf_to_flt[p + ch*m_nPixels] /= (wtsm)? wtsm:1;///  sumf_vl_pt[pp[i]][nbq+ cnt_p[pp[i]]*ch]*wt;
 }
   delete [] vlmc;
   delete [] vlpc;
      delete [] vlmc2;
   delete [] vlpc2;
///////////////////////////////////////////////

	//					  for(int i =0; i< m_nPixels*3; i++ )I_ims[1][i] =0;
	//for(int i =0; i< cl_Q; i++ )
	////for(int j =0; j< cl_Q; j++ )
	//{
	//for(int n =0; n <sigma_cnt[i]; n++){  int j = i_q_pnt[i][n]; for(int c =0; c < 3; c ++)I_ims[1][i + c*m_nPixels + i_q_pnt[i][n]*m_width] =  i_sgm_pnt[i][n]*255;}
	////{  for(int c =0; c < 3; c ++)I_ims[1][i + c*m_nPixels + j*m_width] =  ijK[i+ j*cl_Q] *255;}
	// }
						
	 //------------------------------------------------------------
//	  for(int i = 0; i < cl_Q; i++)buf[i] = new float [m_nPixels*nb_ch];
//	   //---------------------------------------------------------
//	  	                  for(int ic = 0; ic < nb_ch; ic++)
//						  for(int i = 0; i < cl_Q; i++)
//	                       {
//						    for(int p = 0; p < m_nPixels; p++) buf[i][p+ic*m_nPixels] = buf_to_flt[p+ic*m_nPixels]*mask[i][p];
//						   GaussCosConv2DFst( sigma,  m_width, m_height,  &buf[i][ic*m_nPixels] );
//							}
//	  for(int i = 0; i < cl_Q; i++) GaussCosConv2DFst( sigma, m_width, m_height, mask[i]);
//      float * ksi_b = new float [m_nPixels];
//	  float * ksi_m = new float [m_nPixels]; for(int p = 0; p < m_nPixels; p++) ksi_m[p] =1;  GaussCosConv2DFst( sigma,  m_width, m_height, ksi_m);
//
///////////////////////////////////////////////////////////////////////////////
//	                          double ksi_mean = 0; double ml =  mul_tab(cl_Q);
//							  //---------------
//
//
//	  						 for(int p = 0; p < m_nPixels; p++)
//						    {
//							  int iq = q_map[p]; float sum =0;  float ds, dst = 0;
//							  for(int c =0; c < 3; c ++){ ds = buf[iq][p +c*m_nPixels]/mask[iq][p] - I_b[p + c*m_nPixels]; ds *= ds; sum += ds; }
//                             ksi_b[p] = sum;
//							 }
//							 GaussCosConv2DFst( sigma,  m_width, m_height, ksi_b);
//							 for(int p = 0; p < m_nPixels; p++){ ksi_b[p] /= ksi_m[p]; ksi_b[p] = sqrt(ksi_b[p])*ml;}
//							 for(int p = 0; p < m_nPixels; p++)
//							 {
//							 ksi_b[p] = fi_erf(ksi_b[p] , sigc);
//							 }
///////////////////////////////////////////////////////////////////////////////
//
//						 for(int p = 0; p < m_nPixels; p++)
//						  {
//
//							  for(int c =0; c < 3; c ++){ cc[c] = I_b[p + c*m_nPixels];  buf_to_flt[p+ c*m_nPixels] =0;}
//							   float w_sum =0;   							   							
//								for(int i = 0; i <cl_Q; i++)
//								{
//                                    float cv, w,  dist = 0;  int iq = q_map[p];
//
//									if(mask[i][p])
//									{for(int c =0; c < 3; c ++)
//									{
//                                   /*  float ccv = K_cls[i].cl[c];*/
//                                     cv = buf[i][p +c*m_nPixels]/mask[i][p];
//									 if(iq==i){ cv = cv*ksi_b[p] + (1- ksi_b[p])*cc[c];}
//										// cv  = cc[c] + (cv - cc[c])*mul_buf[p+c*m_nPixels];
//										//	 /*cv  = cc[c]*mul_buf[p+c*m_nPixels];*/
//										// /*cv =(mul_buf[p+c*m_nPixels])? cc[c] - (cv - cc[c])*(mul_buf[p+c*m_nPixels]): cv; if(cv>255)cv =255; if(cv<0) cv=0;*/}
//				
//									dist +=(cv - cc[c])*(cv - cc[c]);
//									ccf[c] = cv;
//									}
//									dist = sqrt(dist); w =  gss_c_wt[round_fl(dist)]*mask[i][p];
//									
//									}
//									else {w = 0;}
//									
//									w_sum += w;
//									for(int c =0; c < 3; c ++)buf_to_flt[p +c*m_nPixels] += (w)? ccf[c]*w :0;
//								}
//								for(int c =0; c < 3; c ++){ buf_to_flt[p +c*m_nPixels] /= (w_sum)? w_sum:1;
//								if(buf_to_flt[p +c*m_nPixels]<0)buf_to_flt[p +c*m_nPixels]=0; if(buf_to_flt[p +c*m_nPixels]>255)buf_to_flt[p +c*m_nPixels]=255;}
//						  }
//						  //}
   //for(int p = 0; p < m_nPixels*3; p++)buf_to_flt[p ] = mul_buf[p];
     //------------------------------------
	 delete [] gss_c_wt;
	 delete [] ijK;
	 delete [] q_map;
	 delete [] sigma_cnt;
     delete [] i_sgm_pnt;
	 delete [] i_q_pnt;
	 delete [] i_sgm;
	 delete [] i_q;
	 delete [] map_nb;

	 delete [] sum_vl;
	 delete [] sumf_vl;
	 delete [] sumf_vl_pt;
	 delete [] sum_vl_pt;
	
	 delete [] sum_wt;
	 delete [] sumf_wt;
	 delete [] cnt_p;
	 delete [] map_nbi;


	}////////////////////////////////////

	inline void EDP:: SampledFlt_Add(float sigma_c_mul, float sigma, int cl_Q, unsigned char * I_b, float *buf_to_flt, short int * clr_cls,  EDP::COLOR_F * K_cls, int r, int nb_ch )
{ ///////////////////////////


	  double sigc = sigma_c_mul*255;
	  sigmaX = (m_width%sigmai)? m_width/sigmai+2 :m_width/sigmai+1;
	  sigmaY = (m_height%sigmai)? m_height/sigmai+2 :m_height/sigmai+1;
	  int size1 = sigmaX*sigmaY;
	  int size2 = size1*cl_Q;
	  short int * map_nb = new short int [size2];
	  short int * map_nbi = new short int [size2];
	  float * sum_vl;
	   float * sumf_vl;
	   float ** sum_vl_pt;
	   float ** sumf_vl_pt;
	   float * sum_wt = new float [size2];
	    float * sumf_wt = new float [size2];
		int * cnt_p = new int [size1];
		float mulg1 =exp(-0.5);
		float mulg2 =exp(-2.);


	  int div = 256/r;
	  int n_c = 1000;
	  int cc[3] ;
	  float ccf[3];
	  float * gss_c_wt = new float [n_c];
	  float * ijK = new float [cl_Q*cl_Q];
	  int *sigma_cnt = new int [cl_Q];
	  float **i_sgm_pnt = new float *[cl_Q];
	  int **i_q_pnt = new int * [cl_Q];
	  short int * q_map = new short int [m_nPixels];
	  float nsigm = 3;
	  get_std( sigma_c_mul, I_b,  gss_c_wt,  n_c );
	  //-------MAP-----------------
						  for(int p = 0; p < m_nPixels; p++)
						  {
							  for(int c =0; c < 3; c ++){cc[c] = CUTFL((float)I_b[p + c*m_nPixels]/div, r); }
							  q_map[p] =  clr_cls[IND_HC(cc[0], cc[1], cc[2], r)];							
						  }
//------------------------------------------------------------------------------------------
						  int full_cnt =0;for(int i = 0; i < cl_Q; i++)for(int j = 0; j< cl_Q; j++)
						  {float dst = 0; for(int c =0; c<3; c++)
						  {float r = (K_cls[i].cl[c] - K_cls[j].cl[c]); dst += r*r;} dst = sqrt(dst); float m = ijK[i+ j*cl_Q] = ( dst/sigc>nsigm)? 0:  gss_c_wt[round_fl(dst)]; if(m!=0)full_cnt++;}
float * i_sgm = new float [full_cnt];
int * i_q = new int [full_cnt];
						  //----------------------------------------------------------------------------------------
						  for(int i = 0; i < cl_Q; i++)
						  {   if(!i){i_sgm_pnt[i] = i_sgm; i_q_pnt[i] = i_q;}
						       else
							   {
								   i_sgm_pnt[i] = &i_sgm_pnt[i-1][sigma_cnt[i-1]]; i_q_pnt[i] =  &i_q_pnt[i-1][sigma_cnt[i-1]];
						  }
							   sigma_cnt[i] =0; for(int j = 0; j< cl_Q; j++)
							   {
								   if(ijK[i+ j*cl_Q])
								   { i_sgm_pnt[i][sigma_cnt[i]] = ijK[i+ j*cl_Q];  i_q_pnt[i][sigma_cnt[i]] = j;
								   sigma_cnt[i]++;
								 /*  int a =0;*/
								   }
							   }
						  }

////////////////////INI filter /////////////////////////
  for(int p =0; p <size2; p++)
  {
	  	   map_nb[p] = -1;
	       sum_wt[p]  = 0;
	       sumf_wt[p] = 0;
  }
  for(int p =0; p <size1; p++) cnt_p[p] = 0;
  //---------- accumulation of f ------------------------------
  //--- wt
 for(int y=0; y<m_height; y++)
 for(int x = 0; x < m_width; x++)
 {
	 int p = x + y*m_width;
	 int xx0 = x/sigmai, xx1 = xx0 +1; int xxm = x%sigmai;
	 int yy0 = (y/sigmai); int yy1 = yy0+1; int yym = y%sigmai;
	 int pp[4] = { xx0+ sigmaX*yy0,  xx1+ sigmaX*yy0,  xx0+ sigmaX*yy1,  xx1+ sigmaX*yy1};
     int q = q_map[p];
	 	 //for(int n =0; n <sigma_cnt[q]; n++){////////////////////////////
			/*int qn = i_q_pnt[q][n];*/
			for(int i =0; i< 4; i++)
	        //int i = TBI_p0(xxm,yym);
			if(map_nb[pp[i]*cl_Q+q] == -1){map_nbi[pp[i]*cl_Q+cnt_p[pp[i]]] = q;  map_nb[pp[i]*cl_Q+q] = cnt_p[pp[i]]++;}//}//////////////////////////
	/* for(int i =0; i< 4; i++) { int nbq = map_nb[pp[i]*cl_Q+q]; sum_wt[pp[i]*cl_Q + nbq] += BI_WS(xxm, yym, i); }*/
 }
 //---- func
 int cnt_big = 0; for(int p = 0 ; p < size1; p++) cnt_big += cnt_p[p];

   sum_vl = new float [cnt_big*nb_ch]; sumf_vl = new float [cnt_big*nb_ch];    for(int i = 0 ; i <cnt_big*nb_ch; i++)sum_vl[i] =  0;
   sum_vl_pt = new float* [size1]; sumf_vl_pt = new float* [size1];
   sum_vl_pt[0] = sum_vl; sumf_vl_pt[0] = sumf_vl;

   for(int p = 1; p < size1; p++)
   {
	   sum_vl_pt[p] = sum_vl_pt[p-1] + cnt_p[p-1]*nb_ch;
	   sumf_vl_pt[p] = sumf_vl_pt[p-1] + cnt_p[p-1]*nb_ch;
   }

 for(int y=0; y<m_height; y++)
 for(int x = 0; x < m_width; x++)
 {
	 int p = x + y*m_width;
	 int xx0 = x/sigmai, xx1 = xx0 +1; int xxm = x%sigmai;
	 int yy0 = (y/sigmai); int yy1 = yy0+1; int yym = y%sigmai;
	 int pp[4] = { xx0+ sigmaX*yy0,  xx1+ sigmaX*yy0,  xx0+ sigmaX*yy1,  xx1+ sigmaX*yy1};
     int q = q_map[p];
	
	 //for(int n =0; n <sigma_cnt[q]; n++){////////////////////////////
		//	int qn = i_q_pnt[q][n]; float wc = i_sgm_pnt[q][qn];
	 for(int i =0; i< 4; i++)
	/* int i = TBI_p0(xxm,yym); */
	 {
		 int cntp = cnt_p[pp[i]];
		 int nbq = map_nb[pp[i]*cl_Q+q]; float wt = BI_WS(xxm, yym, i);
		  sum_wt[pp[i]*cl_Q + nbq] += wt;
		 for(int ch = 0; ch <nb_ch; ch++){
		 float vl = buf_to_flt[p + ch*m_nPixels];
		 sum_vl_pt[pp[i]][nbq + cntp*ch] += wt*vl;
		 }
	 }
	 //}///////////////////////
 }
//-------------------------- fiter firs ---------------
  float * vlmc = new float [nb_ch];
  float * vlpc =  new float [nb_ch];
    float * vlmc2 = new float [nb_ch];
  float * vlpc2 =  new float [nb_ch];
  //-----------YY
   for(int y = 0; y < sigmaY; y ++)
   for(int x = 0; x < sigmaX; x++)
   {

		int p = x + sigmaX*y; int cntp =  cnt_p[p]; int pm = (y - 1)*sigmaX + x; int pp = (y + 1)*sigmaX + x; int pm2 = (y - 2)*sigmaX + x; int pp2 = (y + 2)*sigmaX + x;
		for(int ci =0; ci < cntp; ci++ )
		{
		int q_ci = map_nbi[p*cl_Q+ci];   float vlm=0, vlp =0, vlm2 =0, vlp2 =0;  for(int ch = 0; ch <nb_ch; ch++){vlmc[ch] = 0; vlpc[ch] = 0; vlmc2[ch] = 0; vlpc2[ch] = 0;}

		if(y-1>=0) {int pm_ci = map_nb[pm*cl_Q+q_ci]; if(pm_ci != -1){vlm =     sum_wt[pm*cl_Q + pm_ci] ;  for(int ch = 0; ch <nb_ch; ch++)vlmc[ch] = sum_vl_pt[pm][pm_ci + cnt_p[pm]*ch]; }}

		   if(y+1< sigmaY) {int pp_ci = map_nb[pp*cl_Q+q_ci]; if(pp_ci != -1){vlp =     sum_wt[pp*cl_Q + pp_ci] ;  for(int ch = 0; ch <nb_ch; ch++)vlpc[ch] = sum_vl_pt[pp][pp_ci + cnt_p[pp]*ch]; }}
		
		   if(y-2>=0) {int pm_ci2 = map_nb[pm2*cl_Q+q_ci]; if(pm_ci2 != -1){vlm2 =     sum_wt[pm2*cl_Q + pm_ci2] ;  for(int ch = 0; ch <nb_ch; ch++)vlmc2[ch] = sum_vl_pt[pm2][pm_ci2 + cnt_p[pm2]*ch]; }}

		   if(y+2< sigmaY) {int pp_ci2 = map_nb[pp2*cl_Q+q_ci]; if(pp_ci2 != -1){vlp2 =     sum_wt[pp2*cl_Q + pp_ci2] ;  for(int ch = 0; ch <nb_ch; ch++)vlpc2[ch] = sum_vl_pt[pp2][pp_ci2 + cnt_p[pp2]*ch]; }}

	   sumf_wt[p*cl_Q + ci] = sum_wt[p*cl_Q + ci] + vlp*mulg1 + vlm*mulg1 + vlp2*mulg2 + vlm2*mulg2;
	   for(int ch = 0; ch <nb_ch; ch++) sumf_vl_pt[p][ci + cnt_p[p]*ch] = sum_vl_pt[p][ci + cnt_p[p]*ch] + vlpc[ch]*mulg1 + vlmc[ch]*mulg1 + vlpc2[ch]*mulg2 + vlmc2[ch]*mulg2;
		}
   }
    for(int i = 0 ; i <cnt_big*nb_ch; i++)sum_vl[i] =  sumf_vl[i];
	for(int p =0; p <size2; p++)  sum_wt[p]  =   sumf_wt[p];
	  //-----------XX
   for(int y = 0; y < sigmaY; y ++)
   for(int x = 0; x < sigmaX; x++)
   {

		int p = x + sigmaX*y; int cntp =  cnt_p[p]; int pm = (y )*sigmaX + x -1; int pp = (y )*sigmaX + x +1; int pm2 = (y )*sigmaX + x -2; int pp2 = (y )*sigmaX + x +2;
		for(int ci =0; ci < cntp; ci++ )
		{
			int q_ci = map_nbi[p*cl_Q+ci];  float vlm=0, vlp =0, vlm2 =0, vlp2 =0;  for(int ch = 0; ch <nb_ch; ch++){vlmc[ch] = 0; vlpc[ch] = 0; vlmc2[ch] = 0; vlpc2[ch] = 0;}

		if(x-1>=0) {int pm_ci = map_nb[pm*cl_Q+q_ci];  if(pm_ci != -1){vlm =     sum_wt[pm*cl_Q + pm_ci] ;  for(int ch = 0; ch <nb_ch; ch++)vlmc[ch] = sum_vl_pt[pm][pm_ci + cnt_p[pm]*ch]; }}
		   if(x+1< sigmaX) {int pp_ci = map_nb[pp*cl_Q+q_ci]; if(pp_ci != -1){vlp =     sum_wt[pp*cl_Q + pp_ci] ;  for(int ch = 0; ch <nb_ch; ch++)vlpc[ch] = sum_vl_pt[pp][pp_ci + cnt_p[pp]*ch]; }}
		if(x-2>=0) {int pm_ci2 = map_nb[pm2*cl_Q+q_ci];  if(pm_ci2 != -1){vlm2 =     sum_wt[pm2*cl_Q + pm_ci2] ;  for(int ch = 0; ch <nb_ch; ch++)vlmc2[ch] = sum_vl_pt[pm2][pm_ci2 + cnt_p[pm2]*ch]; }}
		   if(x+2< sigmaX) {int pp_ci2 = map_nb[pp2*cl_Q+q_ci]; if(pp_ci2 != -1){vlp2 =     sum_wt[pp2*cl_Q + pp_ci2] ;  for(int ch = 0; ch <nb_ch; ch++)vlpc2[ch] = sum_vl_pt[pp2][pp_ci2 + cnt_p[pp2]*ch]; }}
	   sumf_wt[p*cl_Q + ci] = sum_wt[p*cl_Q + ci] + vlp*mulg1 + vlm*mulg1 + vlp2*mulg2 + vlm2*mulg2;
	   for(int ch = 0; ch <nb_ch; ch++) sumf_vl_pt[p][ci + cnt_p[p]*ch] = sum_vl_pt[p][ci + cnt_p[p]*ch] + vlpc[ch]*mulg1 + vlmc[ch]*mulg1 + vlpc2[ch]*mulg2 + vlmc2[ch]*mulg2;
		}
   }
    for(int i = 0 ; i <cnt_big*nb_ch; i++)sum_vl[i] =  sumf_vl[i];
	for(int p =0; p <size2; p++)  sum_wt[p]  =   sumf_wt[p];
	  //-----------cc
   for(int y = 0; y < sigmaY; y ++)
   for(int x = 0; x < sigmaX; x++)
   {

		int p = x + sigmaX*y; int cntp =  cnt_p[p];
		for(int ci =0; ci < cntp; ci++ )
		{
		int q_ci = map_nbi[p*cl_Q+ci];  sumf_wt[p*cl_Q + ci] = 0; for(int ch = 0; ch <nb_ch; ch++)sumf_vl_pt[p][ci + cnt_p[p]*ch]  = 0;
		for(int n =0; n <sigma_cnt[q_ci]; n++)
		{//////////////////////
			int qn = i_q_pnt[q_ci][n]; int in_vl = map_nb[p*cl_Q+qn];
			if(in_vl != -1)
			{
				float w = i_sgm_pnt[q_ci][n];  sumf_wt[p*cl_Q + ci] += w*sum_wt[p*cl_Q + in_vl];
				for(int ch = 0; ch <nb_ch; ch++)sumf_vl_pt[p][ci + cnt_p[p]*ch] +=w*sum_vl_pt[p][in_vl + cnt_p[p]*ch];
			}
		}////////////////////////
		int a = 0;
		//for(int ch = 0; ch <nb_ch; ch++)sumf_vl_pt[p][ci + cntp*ch]  /=  (sumf_wt[p*cl_Q + ci])? sumf_wt[p*cl_Q + ci]:1;
		}
		int aa =0;

   }
   //-------------------- back to image
 for(int y=0; y<m_height; y++)
 for(int x = 0; x < m_width; x++)
 {
	 int p = x + y*m_width;
	 int xx0 = x/sigmai, xx1 = xx0 +1; int xxm = x%sigmai;
	 int yy0 = (y/sigmai); int yy1 = yy0+1; int yym = y%sigmai;
	 int pp[4] = { xx0+ sigmaX*yy0,  xx1+ sigmaX*yy0,  xx0+ sigmaX*yy1,  xx1+ sigmaX*yy1};
     int q = q_map[p];
	 for(int ch = 0; ch <nb_ch; ch++)  buf_to_flt[p + ch*m_nPixels]= 0;
	 float wtsm = 0;
	  //for(int n =0; n <sigma_cnt[q]; n++){////////////////////////////
			//int qn = i_q_pnt[q][n]; float wc = i_sgm_pnt[q][qn];
	 for(int i =0; i< 4; i++)
	 {
		 int cntp = cnt_p[pp[i]];
		 int nbq = map_nb[pp[i]*cl_Q+q];
		 if(nbq != -1){ float wt =   BI_WS(xxm, yym, i);
		 wtsm += 	 sumf_wt[pp[i]*cl_Q + nbq]*wt;
		 for(int ch = 0; ch <nb_ch; ch++)  buf_to_flt[p + ch*m_nPixels] += sumf_vl_pt[pp[i]][nbq+ cnt_p[pp[i]]*ch]*wt;}
	 }
	for(int ch = 0; ch <nb_ch; ch++)  buf_to_flt[p + ch*m_nPixels] /= (wtsm)? wtsm:1;///  sumf_vl_pt[pp[i]][nbq+ cnt_p[pp[i]]*ch]*wt;
 }
   delete [] vlmc;
   delete [] vlpc;
      delete [] vlmc2;
   delete [] vlpc2;
///////////////////////////////////////////////

	//					  for(int i =0; i< m_nPixels*3; i++ )I_ims[1][i] =0;
	//for(int i =0; i< cl_Q; i++ )
	////for(int j =0; j< cl_Q; j++ )
	//{
	//for(int n =0; n <sigma_cnt[i]; n++){  int j = i_q_pnt[i][n]; for(int c =0; c < 3; c ++)I_ims[1][i + c*m_nPixels + i_q_pnt[i][n]*m_width] =  i_sgm_pnt[i][n]*255;}
	////{  for(int c =0; c < 3; c ++)I_ims[1][i + c*m_nPixels + j*m_width] =  ijK[i+ j*cl_Q] *255;}
	// }
						
	 //------------------------------------------------------------
//	  for(int i = 0; i < cl_Q; i++)buf[i] = new float [m_nPixels*nb_ch];
//	   //---------------------------------------------------------
//	  	                  for(int ic = 0; ic < nb_ch; ic++)
//						  for(int i = 0; i < cl_Q; i++)
//	                       {
//						    for(int p = 0; p < m_nPixels; p++) buf[i][p+ic*m_nPixels] = buf_to_flt[p+ic*m_nPixels]*mask[i][p];
//						   GaussCosConv2DFst( sigma,  m_width, m_height,  &buf[i][ic*m_nPixels] );
//							}
//	  for(int i = 0; i < cl_Q; i++) GaussCosConv2DFst( sigma, m_width, m_height, mask[i]);
//      float * ksi_b = new float [m_nPixels];
//	  float * ksi_m = new float [m_nPixels]; for(int p = 0; p < m_nPixels; p++) ksi_m[p] =1;  GaussCosConv2DFst( sigma,  m_width, m_height, ksi_m);
//
///////////////////////////////////////////////////////////////////////////////
//	                          double ksi_mean = 0; double ml =  mul_tab(cl_Q);
//							  //---------------
//
//
//	  						 for(int p = 0; p < m_nPixels; p++)
//						    {
//							  int iq = q_map[p]; float sum =0;  float ds, dst = 0;
//							  for(int c =0; c < 3; c ++){ ds = buf[iq][p +c*m_nPixels]/mask[iq][p] - I_b[p + c*m_nPixels]; ds *= ds; sum += ds; }
//                             ksi_b[p] = sum;
//							 }
//							 GaussCosConv2DFst( sigma,  m_width, m_height, ksi_b);
//							 for(int p = 0; p < m_nPixels; p++){ ksi_b[p] /= ksi_m[p]; ksi_b[p] = sqrt(ksi_b[p])*ml;}
//							 for(int p = 0; p < m_nPixels; p++)
//							 {
//							 ksi_b[p] = fi_erf(ksi_b[p] , sigc);
//							 }
///////////////////////////////////////////////////////////////////////////////
//
//						 for(int p = 0; p < m_nPixels; p++)
//						  {
//
//							  for(int c =0; c < 3; c ++){ cc[c] = I_b[p + c*m_nPixels];  buf_to_flt[p+ c*m_nPixels] =0;}
//							   float w_sum =0;   							   							
//								for(int i = 0; i <cl_Q; i++)
//								{
//                                    float cv, w,  dist = 0;  int iq = q_map[p];
//
//									if(mask[i][p])
//									{for(int c =0; c < 3; c ++)
//									{
//                                   /*  float ccv = K_cls[i].cl[c];*/
//                                     cv = buf[i][p +c*m_nPixels]/mask[i][p];
//									 if(iq==i){ cv = cv*ksi_b[p] + (1- ksi_b[p])*cc[c];}
//										// cv  = cc[c] + (cv - cc[c])*mul_buf[p+c*m_nPixels];
//										//	 /*cv  = cc[c]*mul_buf[p+c*m_nPixels];*/
//										// /*cv =(mul_buf[p+c*m_nPixels])? cc[c] - (cv - cc[c])*(mul_buf[p+c*m_nPixels]): cv; if(cv>255)cv =255; if(cv<0) cv=0;*/}
//				
//									dist +=(cv - cc[c])*(cv - cc[c]);
//									ccf[c] = cv;
//									}
//									dist = sqrt(dist); w =  gss_c_wt[round_fl(dist)]*mask[i][p];
//									
//									}
//									else {w = 0;}
//									
//									w_sum += w;
//									for(int c =0; c < 3; c ++)buf_to_flt[p +c*m_nPixels] += (w)? ccf[c]*w :0;
//								}
//								for(int c =0; c < 3; c ++){ buf_to_flt[p +c*m_nPixels] /= (w_sum)? w_sum:1;
//								if(buf_to_flt[p +c*m_nPixels]<0)buf_to_flt[p +c*m_nPixels]=0; if(buf_to_flt[p +c*m_nPixels]>255)buf_to_flt[p +c*m_nPixels]=255;}
//						  }
//						  //}
   //for(int p = 0; p < m_nPixels*3; p++)buf_to_flt[p ] = mul_buf[p];
     //------------------------------------
	 delete [] gss_c_wt;
	 delete [] ijK;
	 delete [] q_map;
	 delete [] sigma_cnt;
     delete [] i_sgm_pnt;
	 delete [] i_q_pnt;
	 delete [] i_sgm;
	 delete [] i_q;
	 delete [] map_nb;

	 delete [] sum_vl;
	 delete [] sumf_vl;
	 delete [] sumf_vl_pt;
	 delete [] sum_vl_pt;
	
	 delete [] sum_wt;
	 delete [] sumf_wt;
	 delete [] cnt_p;
	 delete [] map_nbi;


	}////////////////////////////////////
void EDP::K_mean_Flt_Add( EDP::BYTE * I_b,  float sigma_c_mul, float sigmax,float sigmay, int cl_Q, double * buf_to_flt, int n_bufs)
{///------------------------ K_mean

	

	 sigmaix = round_fl(sigmax); if (sigmaix == 0) sigmaix = 1;
	 sigmaiy = round_fl(sigmay); if (sigmaiy == 0) sigmaiy = 1;
	 sigmai2 = sigmaix*sigmaiy;
	 Bi_ws = new double [sigmai2*4];
	 TBi_ws = new double [sigmai2*4];
      TBi_p0 = new short int [sigmai2];
	 for(int f =0; f< 4; f++)
	for(int x =0; x< sigmaix; x++)
	for(int y =0; y< sigmaiy; y++)
	{
		if(!f)Bi_ws[x + y*sigmaix + f*sigmai2] = (float)(sigmaix - x)*(sigmaiy - y)/sigmai2;
		if(f==1)Bi_ws[x + y*sigmaix + f*sigmai2] = (float)(x)*(sigmaiy - y)/sigmai2;
		if(f==2)Bi_ws[x + y*sigmaix + f*sigmai2] = (float)(y)*(sigmaix - x)/sigmai2;
		if(f==3)Bi_ws[x + y*sigmaix + f*sigmai2] = (float)(x)*(y)/sigmai2;
	}
		
	for(int x =0; x< sigmaix; x++)
	for(int y =0; y< sigmaiy; y++)
	{
		float max = Bi_ws[x + y*sigmaix ]; int min_ind = 0;
		for(int f =1; f< 4; f++)if(Bi_ws[x + y*sigmaix + f*sigmai2]>max)
		{max = Bi_ws[x + y*sigmaix + f*sigmai2]; min_ind = f;};
		TBi_p0[x + y*sigmaix] = min_ind;
	}



	 COLOR_F * K_cls = new COLOR_F [cl_Q];
     int r = 32;  int div = 256/r; int size = r*r*r;
	 short int * clr_cls  = new short int [size];
	 int * buf_to_sort  = new int [size]; memset(buf_to_sort, 0, sizeof(int)*size);
         for(int y = 0; y < m_height; y++)for(int x = 0; x < m_width; x++)
		 {
			 int vl[3] = {CUTFL((float)I_b[IND_IC(x,y,0)]/div, r ), CUTFL((float)I_b[IND_IC(x,y,1)]/div, r ),CUTFL((float)I_b[IND_IC(x,y,2)]/div, r )} ;
             buf_to_sort[IND_HC(vl[0],vl[1],vl[2], r)]++;
		 }
	 K_means(cl_Q,  K_cls, r, buf_to_sort, clr_cls);
	 //------------------------------------ test
	 //	 for(int p =0; p < m_nPixels; p++ )
	 //{   int x = p%m_width; int y = p/m_width; int cc[3];
		// int iv = 255;//(TBI_p0(x%sigmai, y%sigmai))? 128:255;
		// for(int c =0; c < 3; c ++)I_ims[1][p + c*m_nPixels ] = TBI_WS(x%sigmai, y%sigmai, c+1)*iv;
	 //}
		// return;
    SampledFlt_Add(sigma_c_mul, sigmax,  cl_Q,  I_b, buf_to_flt,  clr_cls,  K_cls, r,n_bufs );
	 delete [] buf_to_sort;
	 delete [] clr_cls;
	 delete [] K_cls;
	 delete [] Bi_ws;
	 delete [] TBi_ws;
	 delete [] TBi_p0;

	} ///------------------ end K means
void EDP::K_mean_Flt_Add_new_( int * I_b,  float sigma_c_mul, float sigmax,float sigmay, int cl_Q, double * buf_to_flt, int n_bufs)
{///------------------------ K_mean

	

	 sigmaix = round_fl(sigmax); if (sigmaix == 0) sigmaix = 1;
	 sigmaiy = round_fl(sigmay); if (sigmaiy == 0) sigmaiy = 1;
	 sigmai2 = sigmaix*sigmaiy;
	 Bi_ws = new double [sigmai2*4];
	 TBi_ws = new double [sigmai2*4];
      TBi_p0 = new short int [sigmai2];
	 for(int f =0; f< 4; f++)
	for(int x =0; x< sigmaix; x++)
	for(int y =0; y< sigmaiy; y++)
	{
		if(!f)Bi_ws[x + y*sigmaix + f*sigmai2] = (float)(sigmaix - x)*(sigmaiy - y)/sigmai2;
		if(f==1)Bi_ws[x + y*sigmaix + f*sigmai2] = (float)(x)*(sigmaiy - y)/sigmai2;
		if(f==2)Bi_ws[x + y*sigmaix + f*sigmai2] = (float)(y)*(sigmaix - x)/sigmai2;
		if(f==3)Bi_ws[x + y*sigmaix + f*sigmai2] = (float)(x)*(y)/sigmai2;
	}
		
	for(int x =0; x< sigmaix; x++)
	for(int y =0; y< sigmaiy; y++)
	{
		float max = Bi_ws[x + y*sigmaix ]; int min_ind = 0;
		for(int f =1; f< 4; f++)if(Bi_ws[x + y*sigmaix + f*sigmai2]>max)
		{max = Bi_ws[x + y*sigmaix + f*sigmai2]; min_ind = f;};
		TBi_p0[x + y*sigmaix] = min_ind;
	}



	  //{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{
	BYTE * I_b3 = new BYTE [3*m_nPixels];  for (int p =0; p< m_nPixels; p++)I_b3[p]  = I_b3[p+ m_nPixels]  = I_b3[p+ 2*m_nPixels]  = I_b[p];
short int *K_ind = new short int [m_nPixels];    //cluster`s indexes
COLOR_F *K_cls = new COLOR_F [cl_Q]; for (int q =0; q< cl_Q; q++)K_cls[q].cl[0] = K_cls[q].cl[1] = K_cls[q].cl[2] = 0; //COLOR_F *  3 float value per cluster
int * Cnt =new int [cl_Q];for (int q =0; q< cl_Q; q++)Cnt[q] = 0;
 for (int p =0; p< m_nPixels; p++)
{
	int q = K_ind[p] = I_b[p]; Cnt[q]++;
	for (int c =0; c< 3; c++) K_cls[q].cl[c]  += I_b[p];

}
 for (int q =0; q< cl_Q; q++){ for (int c =0; c< 3; c++) K_cls[q].cl[c] /= Cnt[q];}
delete [] Cnt;
//for(int p = 0; p < m_nPixels; p++)for (int c =0; c< 3; c++)
//{I_ims[1][p + c*m_nPixels] =  K_cls[K_ind[p]].cl[c];}
//return;
//}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}
    SampledFlt_Add_new(sigma_c_mul, sigmax,  cl_Q,  I_b3, buf_to_flt,  K_ind,  K_cls, n_bufs );
	 delete [] K_ind;
	 delete [] I_b3;
	 //delete [] clr_cls;
	 delete [] K_cls;
	 delete [] Bi_ws;
	 delete [] TBi_ws;
	 delete [] TBi_p0;

	} ///------------------ end K means
void EDP::ITR_BF(int iti,  EDP::BYTE * I_b,  float sigma_c_mul, float sigmax,float sigmay, int cl_Q, double * C, int n_bufs)
{
	double  *CC = new double [m_nPixels*m_nLabels];
	memcpy(CC, C, sizeof(double)*m_nPixels*m_nLabels);
	
	double  *M = new double [m_nPixels*m_nLabels];
	double  *MM = new double [m_nPixels*m_nLabels];

	for(int it =0; it < iti; it++) {
		//----------------------------------------------------------------
		UpdMsgC_Z_full( CC, M);
		memcpy(MM, M, sizeof(double)*m_nPixels*m_nLabels);
		K_mean_Flt_Add_new(  I_b,  sigma_c_mul,sigmax, sigmay,  cl_Q, MM, n_bufs) ;
		for(int p =0; p < m_nPixels*m_nLabels; p++ ) {
			CC[p] = C[p] - MM[p]+M[p];
		}
	}//----------------------------------------------------------------
	memcpy(C, CC, sizeof(double)*m_nPixels*m_nLabels);
	delete [] CC;
	delete [] M;
	delete [] MM;
}
void EDP::K_mean_Flt_Add_new( EDP::BYTE * I_b,  float sigma_c_mul, float sigmax,float sigmay, int cl_Q, double * buf_to_flt, int n_bufs) {
	///------------------------ K_mean

	sigmaix = round_fl(sigmax);
	if (sigmaix == 0) sigmaix = 1;
	
	sigmaiy = round_fl(sigmay);
	if (sigmaiy == 0) sigmaiy = 1;
	
	sigmai2 = sigmaix*sigmaiy;
	Bi_ws = new double [sigmai2*4];
	TBi_ws = new double [sigmai2*4];
	TBi_p0 = new short int [sigmai2];
	for(int f =0; f< 4; f++)
	for(int x =0; x< sigmaix; x++)
	for(int y =0; y< sigmaiy; y++) {
		if(!f) {
			Bi_ws[x + y*sigmaix + f*sigmai2] = (float)(sigmaix - x)*(sigmaiy - y)/sigmai2;
		}
		if(f==1) {
			Bi_ws[x + y*sigmaix + f*sigmai2] = (float)(x)*(sigmaiy - y)/sigmai2;
		}
		if(f==2) {
			Bi_ws[x + y*sigmaix + f*sigmai2] = (float)(y)*(sigmaix - x)/sigmai2;
		}
		if(f==3) {
			Bi_ws[x + y*sigmaix + f*sigmai2] = (float)(x)*(y)/sigmai2;
		}
	}

	for(int x =0; x< sigmaix; x++)
	for(int y =0; y< sigmaiy; y++) {
		float max = Bi_ws[x + y*sigmaix];
		int min_ind = 0;
		
		for(int f =1; f< 4; f++) {
			if(Bi_ws[x + y*sigmaix + f*sigmai2]>max) {
				max = Bi_ws[x + y*sigmaix + f*sigmai2];
				min_ind = f;
			};
		}
		TBi_p0[x + y*sigmaix] = min_ind;
	}


	//{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{
	short int *K_ind = new short int [m_nPixels];    //cluster`s indexes
	COLOR_F *K_cls; //COLOR_F *  3 float value per cluster

	BYTE * im[3] = {&I_b[0], & I_b[m_nPixels], &I_b[2*m_nPixels]}; // R G B layers of input image
	if(cl_Q == 1) {
		/////////////// CLUSTERS = 1 /////////////////
		K_cls = new COLOR_F [1]; // Average value of cluster  n => K_cls[n].cl[0] = B; K_cls[n].cl[1] = G; K_cls[n].cl[2] = R; COLOR_F *  3 float value per cluster
		for (int c =0; c< 3; c++) {
			K_cls[0].cl[c] = 0;
		}
		for(int p = 0; p < m_nPixels; p++) {
			for (int c =0; c< 3; c++) {
				K_cls[0].cl[c]  += im[c][p];
			}
			K_ind[p] = 0;
		}
		for (int c =0; c< 3; c++) {
			K_cls[0].cl[c]  /= m_nPixels;
		}
	}
	///////////////////END CL =1 ///////////////////
	{
		/////////////// CLUSTERS > 1 /////////////////
		unsigned int exp2[32]; 
		exp2[0]=1;
		for(int i = 1; i < 32; i++) {
			exp2[i] = exp2[i-1]*2;
		}
		int lg2Q = 0;
		while(exp2[lg2Q] < cl_Q  ) {
			lg2Q++;
		}
		cl_Q = exp2[lg2Q]; // first power of 2 greater than input cl_Q
		K_cls = new COLOR_F [cl_Q];// Average value of cluster  n => K_cls[n].cl[0] = B; K_cls[n].cl[1] = G; K_cls[n].cl[2] = R;
		int *H = new int [256*cl_Q];
		int *Cnt = new int [cl_Q];
		int  *Thr = new int [cl_Q];
		//----------------------------- ini K_ind --------------
		for(int p = 0; p < m_nPixels; p++) {
			K_ind[p] = 0;
		}
		for(int l = 0;  l < lg2Q ; l++) {
			//============== level of tree ===========
			BYTE * imch = im[(l+1)%3];
			memset(H, 0, sizeof(int)*256*cl_Q);
			memset(Cnt, 0, sizeof(int)*cl_Q);
			memset(Thr, 0, sizeof(int)*cl_Q);
			int lB = exp2[l];
			for(int p = 0; p < m_nPixels; p++) {
				//............... statistic.............
				int q = K_ind[p];
				H[imch[p] + q*256]++; Cnt[q]++;
			}//.........................................
			for(int ll = 0; ll < lB; ll++){
				Cnt[ll] /= 2;  int sum =0; int md =0;
				while(Cnt[ll] > sum) {
					sum += H[(md++) +  ll*256];
				}
				md--;
				Thr[ll] = md;
				Cnt[ll] = H[md + ll*256]/2;
			}
			for(int p = 0; p < m_nPixels; p++) {
				//............ dvd into 2 .....
				int q = K_ind[p];
				if(imch[p] == Thr[q] && Cnt[q] == 0) {
					K_ind[p] += lB;
				}
				if(imch[p] == Thr[q] && Cnt[q]) {
					Cnt[q]--;
				}
			    if(imch[p] >   Thr[q] ) {
			    	K_ind[p] += lB;
			    }
			}//..................................
		}//====================================
		for (int q =0; q< cl_Q; q++) {
			Cnt[q] =0;
			for (int c =0; c< 3; c++) {
				K_cls[q].cl[c] = 0;
			}
		}
		for(int p = 0; p < m_nPixels; p++) {
			int q = K_ind[p]; Cnt[q]++;
			for (int c =0; c< 3; c++) {
				K_cls[q].cl[c]  += im[c][p];
			}
		}

		for (int q =0; q< cl_Q; q++) {
			for (int c =0; c< 3; c++)
				K_cls[q].cl[c] /= Cnt[q];
		}
		delete [] H;
		delete [] Cnt;
		delete [] Thr;
	}
	//////////////////END CL > 1//////////////////
	//for(int p = 0; p < m_nPixels; p++)for (int c =0; c< 3; c++)
	//{I_ims[1][p + c*m_nPixels] =  K_cls[K_ind[p]].cl[c];}
	//return;
	//}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}
	SampledFlt_Add_new(sigma_c_mul, sigmax,  cl_Q,  I_b, buf_to_flt,  K_ind,  K_cls, n_bufs );
	delete [] K_ind;
	//delete [] clr_cls;
	delete [] K_cls;
	delete [] Bi_ws;
	delete [] TBi_ws;
	delete [] TBi_p0;

} ///------------------ end K means

void EDP::K_mean_Flt_Add_gr( EDP::BYTE * I_b,  float sigma_c_mul, float sigmax,float sigmay, int cl_Q, double * buf_to_flt, int n_bufs)
{///------------------------ K_mean

	sigmaix = round_fl(sigmax); if (sigmaix == 0) sigmaix = 1;
	sigmaiy = round_fl(sigmay); if (sigmaiy == 0) sigmaiy = 1;
	sigmai2 = sigmaix*sigmaiy;
	Bi_ws = new double [sigmai2*4];
	TBi_ws = new double [sigmai2*4];
	TBi_p0 = new short int [sigmai2];
	for(int f =0; f< 4; f++)
	for(int x =0; x< sigmaix; x++)
	for(int y =0; y< sigmaiy; y++)
	{
		if(!f)Bi_ws[x + y*sigmaix + f*sigmai2] = (float)(sigmaix - x)*(sigmaiy - y)/sigmai2;
		if(f==1)Bi_ws[x + y*sigmaix + f*sigmai2] = (float)(x)*(sigmaiy - y)/sigmai2;
		if(f==2)Bi_ws[x + y*sigmaix + f*sigmai2] = (float)(y)*(sigmaix - x)/sigmai2;
		if(f==3)Bi_ws[x + y*sigmaix + f*sigmai2] = (float)(x)*(y)/sigmai2;
	}
		
	for(int x =0; x< sigmaix; x++)
	for(int y =0; y< sigmaiy; y++) {
		float max = Bi_ws[x + y*sigmaix ]; 
		int min_ind = 0;
		for(int f =1; f< 4; f++)
			if(Bi_ws[x + y*sigmaix + f*sigmai2]>max) {
				max = Bi_ws[x + y*sigmaix + f*sigmai2]; 
				min_ind = f;
			};
		TBi_p0[x + y*sigmaix] = min_ind;
	}



	COLOR_F * K_cls = new COLOR_F [cl_Q];
	int r = 32;  int div = 256/r; int size = r*r*r;
	short int * clr_cls  = new short int [size];
	int * buf_to_sort  = new int [size]; memset(buf_to_sort, 0, sizeof(int)*size);
		for(int y = 0; y < m_height; y++)
			for(int x = 0; x < m_width; x++) {
				int vl[3] = {CUTFL((float)I_b[IND_IC(x,y,0)]/div, r ), CUTFL((float)I_b[IND_IC(x,y,1)]/div, r ),CUTFL((float)I_b[IND_IC(x,y,2)]/div, r )} ;
				buf_to_sort[IND_HC(vl[0],vl[1],vl[2], r)]++;
			}
	K_means(cl_Q,  K_cls, r, buf_to_sort, clr_cls);
	//------------------------------------ test
	//	 for(int p =0; p < m_nPixels; p++ )
	//{   int x = p%m_width; int y = p/m_width; int cc[3];
	// int iv = 255;//(TBI_p0(x%sigmai, y%sigmai))? 128:255;
	// for(int c =0; c < 3; c ++)I_ims[1][p + c*m_nPixels ] = TBI_WS(x%sigmai, y%sigmai, c+1)*iv;
	//}
	// return;
	SampledFlt_Add(sigma_c_mul, sigmax,  cl_Q,  I_b, buf_to_flt,  clr_cls,  K_cls, r,n_bufs );
	delete [] buf_to_sort;
	delete [] clr_cls;
	delete [] K_cls;
	delete [] Bi_ws;
	delete [] TBi_ws;
	delete [] TBi_p0;

} ///------------------ end K means

void EDP::K_mean_Flt_Add( EDP::BYTE * I_b,  float sigma_c_mul, float sigma, int cl_Q, float * buf_to_flt, int n_bufs) {
	///------------------------ K_mean
	sigmai = round_fl(sigma); if (sigmai == 0) sigmai = 1;
	sigmai2 = sigmai*sigmai;
	Bi_ws = new double[sigmai2*4];
	TBi_ws = new double [sigmai2*4];
	TBi_p0 = new short int [sigmai2];
	for(int f =0; f< 4; f++)
	for(int x =0; x< sigmai; x++)
	for(int y =0; y< sigmai; y++)
	{
		if(!f)Bi_ws[x + y*sigmai + f*sigmai2] = (float)(sigmai - x)*(sigmai - y)/sigmai2;
		if(f==1)Bi_ws[x + y*sigmai + f*sigmai2] = (float)(x)*(sigmai - y)/sigmai2;
		if(f==2)Bi_ws[x + y*sigmai + f*sigmai2] = (float)(y)*(sigmai - x)/sigmai2;
		if(f==3)Bi_ws[x + y*sigmai + f*sigmai2] = (float)(x)*(y)/sigmai2;
	}
		
	for(int x =0; x< sigmai; x++)
	for(int y =0; y< sigmai; y++)
	{
		float max = Bi_ws[x + y*sigmai ]; int min_ind = 0; for(int f =1; f< 4; f++)if(Bi_ws[x + y*sigmai + f*sigmai2]>max){max = Bi_ws[x + y*sigmai + f*sigmai2]; min_ind = f;};
		TBi_p0[x + y*sigmai] = min_ind;
	}

	//  float rr[4];	
	//for(int x =0; x< sigmai; x++)
	//for(int y =0; y< sigmai; y++)
	//{
	//	rr[0] = exp(-((float)x*x + y*y)/(2*sigmai2));
	//	rr[1] = exp(-((float)(sigmai-x)*(sigmai-x) + y*y)/(2*sigmai2));
	//	rr[2] = exp(-((float)(sigmai-y)*(sigmai-y) + x*x)/(2*sigmai2));
	//	rr[3] = exp(-((float)(sigmai-x)*(sigmai-x) + (sigmai-y)*(sigmai-y))/(2*sigmai2));
 //       float sum = 0; for(int f =0; f< 4; f++) sum += rr[f];
	//	for(int f =0; f< 4; f++) Bi_ws[x + y*sigmai + f*sigmai2] = rr[f]/sum;
	//}
	//  	 for(int f =0; f< 4; f++)
	//for(int x =0; x< sigmai; x++)
	//for(int y =0; y< sigmai; y++)
	//{
	//	if(!f)Bi_ws[x + y*sigmai + f*sigmai2] = (float)(sigmai - x)*(sigmai - y)/sigmai2;
	//	if(f==1)Bi_ws[x + y*sigmai + f*sigmai2] = (float)(x)*(sigmai - y)/sigmai2;
	//	if(f==2)Bi_ws[x + y*sigmai + f*sigmai2] = (float)(y)*(sigmai - x)/sigmai2;
	//	if(f==3)Bi_ws[x + y*sigmai + f*sigmai2] = (float)(x)*(y)/sigmai2;
	//}
	// float rr[4]; float sqt = (float)sigmai2/2; float diag = sqrt((float)sigmai2*2);
	//for(int x =0; x< sigmai; x++)
	//for(int y =0; y< sigmai; y++)
	//{
	//	rr[0] = sqrt((float)x*x + y*y);
	//	rr[1] = sqrt((float)(sigmai-x)*(sigmai-x) + y*y);
	//	rr[2] = sqrt((float)(sigmai-y)*(sigmai-y) + x*x);
	//	rr[3] = sqrt((float)(sigmai-x)*(sigmai-x) + (sigmai-y)*(sigmai-y));
	//	int p0 = TBi_p0[x + y*sigmai] = (rr[0]<rr[3]) ? 0:3;
	//	for(int f =0; f< 4; f++){
	//		if(!f){ float pr = (diag + rr[1] + rr[2])/2;
	//			TBi_ws[x + y*sigmai + f*sigmai2] = sqrt(pr*(pr-rr[1])*(pr - rr[2])*(pr - diag))/sqt;}
	//	if(f==1){  float pr = ((float)sigmai + rr[p0] + rr[2])/2;
	//			TBi_ws[x + y*sigmai + f*sigmai2] = sqrt(pr*(pr-rr[p0])*(pr - rr[2])*(pr - sigmai))/sqt;}
	//	if(f==2){  float pr = ((float)sigmai + rr[p0] + rr[1])/2;
	//			TBi_ws[x + y*sigmai + f*sigmai2] = sqrt(pr*(pr-rr[p0])*(pr - rr[1])*(pr - sigmai))/sqt;}
	//	if(f==3){ float pr = (diag + rr[1] + rr[2])/2;
	//			TBi_ws[x + y*sigmai + f*sigmai2] = sqrt(pr*(pr-rr[1])*(pr - rr[2])*(pr - diag))/sqt;}
	//	}
	//}
  //float test = TBI_WS(73%sigmai, 52%sigmai, 0);
	 //float test1 = TBI_WS(73%sigmai, 52%sigmai, 1);
	 //float test2 = TBI_WS(73%sigmai, 52%sigmai, 2);
	 //float test3 = TBI_WS(73%sigmai, 52%sigmai, 3);
	 //float a = test1+test2+test3;
	COLOR_F * K_cls = new COLOR_F [cl_Q];
	int r = 32;  int div = 256/r; int size = r*r*r;
	short int * clr_cls  = new short int [size];
	int * buf_to_sort  = new int [size]; memset(buf_to_sort, 0, sizeof(int)*size);
	for(int y = 0; y < m_height; y++)
		for(int x = 0; x < m_width; x++) {
			int vl[3] = {CUTFL((float)I_b[IND_IC(x,y,0)]/div, r ), CUTFL((float)I_b[IND_IC(x,y,1)]/div, r ),CUTFL((float)I_b[IND_IC(x,y,2)]/div, r )} ;
			buf_to_sort[IND_HC(vl[0],vl[1],vl[2], r)]++;
		}
	K_means(cl_Q,  K_cls, r, buf_to_sort, clr_cls);
	//------------------------------------ test
	//	 for(int p =0; p < m_nPixels; p++ )
	//{   int x = p%m_width; int y = p/m_width; int cc[3];
	// int iv = 255;//(TBI_p0(x%sigmai, y%sigmai))? 128:255;
	// for(int c =0; c < 3; c ++)I_ims[1][p + c*m_nPixels ] = TBI_WS(x%sigmai, y%sigmai, c+1)*iv;
	//}
	// return;
	SampledFlt_Add(sigma_c_mul, sigma,  cl_Q,  I_b, buf_to_flt,  clr_cls,  K_cls, r,n_bufs );
	delete [] buf_to_sort;
	delete [] clr_cls;
	delete [] K_cls;
	delete [] Bi_ws;
	delete [] TBi_ws;
	delete [] TBi_p0;

} ///------------------ end K means

void EDP::MRF__(int itr, int * dL, int * dR)
{
	 //float sigma_c_mul =  Cnst_2; float sigma =  Cnst_1; int cl_Q =300;/*100*(1+iii);*/ int n_bufs =m_nLabels;
    float * c_D = new float [m_height*Tau]; memcpy(c_D, m_D, m_height*Tau*sizeof(float));
	REAL * Mrgnl = new REAL [m_height*Tau];
	//float * m_cost = new float [m_nPixels*this->m_nLabels];
	//copy_m_D( m_cost, 0);
	//prf_cost_inv( m_cost);	
	//K_mean_Flt_Add(I_ims[0],  sigma_c_mul, sigma, cl_Q, m_cost, n_bufs);
	//this->Cost_mean = copy_m_D_b( m_cost, 0);
	//delete [] m_cost;
	Make_Gr_inf(0, 3,12); 
	Smthnss[0].max = Lmbd_trn*Cost_mean; 
	Smthnss[0].semi  =  Lmbd_trn*Cost_mean/4;
    TRWS_L(itr,Mrgnl, dL);
	//TRWS_R(itr,Mrgnl, dR);
	//
	//for(int i =0; i<3; i++){
	//	  MkDspCorr_pp(dL, dR, 0);
	//      TRWS_L(itr,Mrgnl, dL);
	//	  memcpy(m_D, c_D, m_height*Tau*sizeof(float));
	//      MkDspCorr_pp(dL, dR, 1);
	//      TRWS_R(itr,Mrgnl, dR);
	//      memcpy(m_D, c_D, m_height*Tau*sizeof(float));}
	//MkDspCorr_RL(  dL, dR);
	//for(int i=0; i<2; i++)
	//{
	//	int r_w = 5;  float cl_prc = 0.5;
	//Gss_wnd_( r_w,  cl_prc, dL, 0, m_nLabels );
	//Gss_wnd_( r_w,  cl_prc, dR, 1, m_nLabels );
	//}
 //   MkDspCorr_RL(  dL, dR);

	delete [] Mrgnl;
	delete [] c_D;
		
}
void EDP::mk_BF_cost( float sigma_c_mul, float sigmax,float sigmay, int cl_Q, double **intr_cst, int n_bufs) {
	for(int lr = 0; lr<2;lr++) {
		BYTE * I_b = I_ims[lr];
		copy_m_D( intr_cst[lr], lr);
		//OneD_filterX(I_b, m_cost, 3, 0.01);
		/*OneD_filter(I_b, m_cost, sigma, sigma_c_mul);*/
		K_mean_Flt_Add(I_b,  sigma_c_mul, sigmax,sigmay, cl_Q, intr_cst[lr], n_bufs);
	}
}

void EDP::MRF__z(int itr, int * dL, int * dR) {
	Smthnss = new DELTA_D [5];
	float c_lbl = 10, occ =0;
	Size_fw = 18;
	float sigma_c_mul = 0.75/255;
	float sigmax = Size_fw;
	float sigmay = Size_fw ;
	int cl_Q = (N_LB <30) ? 2000 : 4000;/*100*(1+iii);*/
	int n_bufs =m_nLabels;
	Lmbd_trn = 0.5;

	REAL * Mrgnl = new REAL [m_nPixels*this->m_nLabels];
	double *intr_cst[2];
	double * m_cost = new double [m_nPixels*this->m_nLabels];
	//////////////////////////////////////////////////////////////////////
	//float mns[3] = {0,0,0};  float sgs[3] ={0,0,0}; FOR_PX_p { for(int c=0; c <3; c++)mns[c] += I_ims[0][p+ c*N_PX];}
	// for(int c=0; c <3; c++)mns[c] /= N_PX;  FOR_PX_p { for(int c=0; c <3; c++)sgs[c] += (I_ims[0][p+ c*N_PX]-mns[c])*(I_ims[0][p+ c*N_PX]-mns[c]);}
	//  for(int c=0; c <3; c++)sgs[c] = sqrt(sgs[c]/N_PX);
	//////////////////////////////////////////////////////////////////
	
	if(N_LB>30) {
		Make_Gr_inf_pp(0.011,15,  7, 15);
	} else {
		Make_Gr_inf(0,  7, 15);
	}
	
	float Lm_tr[3] = {2.5*0.125, 1.25*0.125,1.*0.125};
	float ml[3] = {1./6, 1./6, 1./6};

	for(int q =0; q<3; q++){
		Smthnss[q].max = Lm_tr[q];
		Smthnss[q].semi = Lm_tr[q]*ml[q];
	}
	//////////////////////////////////////////////////////////
	for(int lr = 0; lr<2;lr++) {
		int * rez =  (lr)? dR : dL;
		BYTE * I_b = I_ims[lr];
		copy_m_D(m_cost, lr);
		int thr = 5;
		float alp =0.;
		
		Make_Gr_fl_buf(  thr, lr,m_cost, N_LB, alp );
		/*if(Sc_out == 1)*/
		K_mean_Flt_Add_new(I_b, sigma_c_mul, sigmax,sigmay, cl_Q, m_cost, n_bufs);
		#if defined ORIGINAL
			TRWS_CST(itr,Mrgnl, m_cost, rez,lr);
		#elif defined AE
			infer(std::max(itr,100), m_cost, rez, lr, I_b, m_height, m_width, N_LB);
		#elif defined SYM_AE
			int nRanksConsidered = ranksConsidered;
			int nCBPIterations = cbpIters;
			symmetries_infer(std::max(itr,100), m_cost, rez, lr, I_b, m_height, m_width, N_LB, nRanksConsidered, nCBPIterations);
		#endif
		if(!lr){
			occ  = sigma_img_get(rez);
			if (occ <= c_lbl) {
				lr = 2;
			}
		}
	} // end for(lr=0...)

	{	
		sigma_img_get(dL);
	
		if(occ > c_lbl) {
			MkDspCorr_RL(  dL, dR);
		}

		{
			if(occ > c_lbl)
				for(int i=0; i<1; i++) {
					int r_w = 4;
					float cl_prc = sigma_img;//0.23;
					Gss_wnd_( r_w,  cl_prc, dL, 0, m_nLabels );
					Gss_wnd_( r_w,  cl_prc, dR, 1, m_nLabels );
				}
		}

		{
			if(occ > c_lbl) {
				MkDspCorr_RL(  dL, dR);
			}
			int r_w1 = 4;
			float cl_prc1 = sigma_img;//0.23;
			Gss_wnd_( r_w1,  cl_prc1, dL, 0, m_nLabels );
		}

		Outlr_reg_filter(0.001, 3, dL, 0, 10, 0.01);

		if(0)  {
			Intermed(dL, dR);
		}
		int r = 4;//(N_LB >16) ? 4:2;
		if(occ > c_lbl) {
			CorrectSub_new(dL, dR, /*intr_cst[0],*/ r, 1);
		} else {
			for(int i =0; i < m_nPixels; i++ ){
				dL[i] *= 1;
			}
		}
	}
	//----------------------------------------------
	delete [] m_cost;
	delete [] Mrgnl;
}


// This function attempts to get error per iteration along with other details.
void EDP::modified_MRF__z(int itr, int * dL, int * dR, png::image< png::ga_pixel >& truthimg) {
	using namespace std;
	clock_t start,end; // timers.
	typedef opengm::AlphaExpansionFusion<Model, opengm::Minimizer> AEFInferType;
	std::cout << "Entered modified_MRF__z" << std::endl;
	AEFInferType::Parameter params(1); // we're gonna step through this.
	// std::vector< AEFInferType::TimingVisitorType > visitors(2, AEFInferType::TimingVisitorType());	
	std::vector<AEFInferType> aes; // Empty for now.

	// Initialize variables
	std::vector<Model> gms(2); // gm[0] is LR, gm[1] is RL
	std::vector< std::vector< std::vector<int> > > syms(2);
	std::vector< std::vector<int> > invsyms(2);
	std::vector< std::vector<int> > veclabels(2);
	int nRanksConsidered = ranksConsidered;
	int nCBPIterations = cbpIters;
	int nx=m_width;
	int ny=m_height;
	int nLabels=this->m_nLabels;
	// ASSERT: assumes that global variables truth_file and scale are correctly set.
	// STEP : Initialize m_cost
	std::cerr << "nLabels=" << nLabels << std::endl;
	Smthnss = new DELTA_D [5];
	float c_lbl = 10, occ =0;
	Size_fw = 18;
	float sigma_c_mul = 0.75/255;
	float sigmax = Size_fw;
	float sigmay = Size_fw ;
	int cl_Q = (N_LB <30) ? 2000 : 4000;/*100*(1+iii);*/
	int n_bufs =m_nLabels;
	Lmbd_trn = 0.5;
	REAL * Mrgnl = new REAL [m_nPixels*this->m_nLabels];
	double *intr_cst[2];
	double * m_cost = new double [m_nPixels*this->m_nLabels];
	if(N_LB>30) {
		Make_Gr_inf_pp(0.011,15,  7, 15);
	} else {
		Make_Gr_inf(0,  7, 15);
	}
	float Lm_tr[3] = {2.5*0.125, 1.25*0.125,1.*0.125};
	float ml[3] = {1./6, 1./6, 1./6};
	for(int q =0; q<3; q++){
		Smthnss[q].max = Lm_tr[q];
		Smthnss[q].semi = Lm_tr[q]*ml[q];
	} 

	
	// STEP : Inference
	// Upper layer - KMeans
	png::image< png::rgb_pixel > img(nx, ny);
	for(int lr = 0; lr<2;lr++) {
		int * rez =  (lr)? dR : dL;
		BYTE * I_b = I_ims[lr];
		copy_m_D(m_cost, lr);
		int thr = 5;
		float alp =0.;
		
		Make_Gr_fl_buf(  thr, lr,m_cost, N_LB, alp );
		/*if(Sc_out == 1)*/
		K_mean_Flt_Add_new(I_b, sigma_c_mul, sigmax,sigmay, cl_Q, m_cost, n_bufs);
		
		// Need RVO here to speed it up. FIXME: use Return-by-value optimization.
		ucharArray2pngImage(img, I_b, nx, ny);
		#if defined AE
		
		start = clock();
		gms[lr] = constructGraphicalModel(img, nLabels, m_cost); // TODO
		end = clock();
		std::cout << "AE model construction took " << double(end-start)/CLOCKS_PER_SEC << std::endl;
		
		#elif defined SYM_AE
		
		std::cout << "SYM Sym execution started " << nRanksConsidered << " " << nCBPIterations << std::endl;
		start = clock();
		getSymmetries(m_cost, img, nx, ny, nLabels, nRanksConsidered, nCBPIterations, syms[lr], invsyms[lr]);
		end = clock();
		std::cout << "SYM symmetry finding took " << double(end-start)/CLOCKS_PER_SEC << std::endl;
		
		std::cout << "lr=" << lr << "has " << syms.size() << " groups" << std::endl;
		start = clock();
		gms[lr] = constructGraphicalModel(img, nLabels, m_cost, syms[lr], invsyms[lr]); // Using symmetries
		end = clock();
		std::cout << "SYM model construction took " << double(end-start)/CLOCKS_PER_SEC << std::endl;
		
		#endif
		
		aes.push_back( AEFInferType(gms[lr], params) );
		veclabels[lr] = std::vector<int>(gms[lr].numberOfVariables(),0);

	} // end for(lr=0...)

	// Lower layer - MRF
	// Save labelling from last time - LR and RL
	AEFInferType::TimingVisitorType visitorl;
	AEFInferType::TimingVisitorType visitorr;
	for(int iter=0; iter<100; ++iter) {
		// We don't really understand the magic that happens in the original MRF__z
		// hence, we'll just use the entire TSGO algorithm here...
		std::cout << "Starting iteration=" << iter << std::endl;
		for(int lr = 0; lr<2;lr++) {
			int * rez =  (lr)? dR : dL;

			// Do a single step of inference on gms[lr]
			aes[lr].setStartingPoint( veclabels[lr].begin() );
			std::cout << "lr=" << lr << std::endl;
			if(lr==0) aes[lr].infer(visitorl);
			else aes[lr].infer(visitorr);
			
			aes[lr].arg(veclabels[lr]);

			// Copy the result into veclabels[lr], rez
			#if defined AE
			copyVecLabels2Rez(veclabels[lr], rez);
			#elif defined SYM_AE
			copyVecLabels2Rez(veclabels[lr], invsyms[lr], rez);
			#endif
			// Important: Removed this block because it messes with iterations.
			// if(!lr){
			// 	occ  = sigma_img_get(rez);
			// 	if (occ <= c_lbl) {
			// 		lr = 2;
			// 	}
			// }

		} // end for(lr=0...)
		// Postprocessing
		{
			sigma_img_get(dL);
		
			if(occ > c_lbl) {
				MkDspCorr_RL(  dL, dR);
			}

			{
				if(occ > c_lbl)
					for(int i=0; i<1; i++) {
						int r_w = 4;
						float cl_prc = sigma_img;//0.23;
						Gss_wnd_( r_w,  cl_prc, dL, 0, m_nLabels );
						Gss_wnd_( r_w,  cl_prc, dR, 1, m_nLabels );
					}
			}

			{
				if(occ > c_lbl) {
					MkDspCorr_RL(  dL, dR);
				}
				int r_w1 = 4;
				float cl_prc1 = sigma_img;//0.23;
				Gss_wnd_( r_w1,  cl_prc1, dL, 0, m_nLabels );
			}

			Outlr_reg_filter(0.001, 3, dL, 0, 10, 0.01);

			if(0)  {
				Intermed(dL, dR);
			}
			int r = 4;//(N_LB >16) ? 4:2;
			if(occ > c_lbl) {
				CorrectSub_new(dL, dR, /*intr_cst[0],*/ r, 1);
			} else {
				for(int i =0; i < m_nPixels; i++ ){
					dL[i] *= 1;
				}
			}
		}
		double badpixelscore;
		// if(iter==99)
		badpixelscore = getBadPixelScore(dL, truthimg, nx, ny, nLabels, threshold);
		std::cout << "On step " << iter << " error was " << badpixelscore << std::endl;
	}
	// STEP : cleanup
	delete [] m_cost;
	delete [] Mrgnl;
}

#include "hybrid.cpp"
void EDP::hybrid_MRF__z(int itr, int * dL, int * dR, png::image< png::ga_pixel >& truthimg) {
	// Step0 : Initializations
	using namespace std;
	clock_t start,end; // timers.
	typedef opengm::AlphaExpansionFusion<Model, opengm::Minimizer> AEFInferType;
	std::cout << "Entered hybrid_MRF__z" << std::endl;
	AEFInferType::Parameter params(1); // we're gonna step through this.
	std::vector<AEFInferType> aes; // Empty for now.
	int nx=m_width;
	int ny=m_height;
	int nLabels=this->m_nLabels;
	std::cerr << "nLabels=" << nLabels << std::endl;
	Smthnss = new DELTA_D [5];
	float c_lbl = 10, occ =0;
	Size_fw = 18;
	float sigma_c_mul = 0.75/255;
	float sigmax = Size_fw;
	float sigmay = Size_fw ;
	int cl_Q = (N_LB <30) ? 2000 : 4000;/*100*(1+iii);*/
	int n_bufs =m_nLabels;
	Lmbd_trn = 0.5;
	REAL * Mrgnl = new REAL [m_nPixels*this->m_nLabels];
	double *intr_cst[2];
	double * m_cost = new double [m_nPixels*this->m_nLabels];
	if(N_LB>30) {
		Make_Gr_inf_pp(0.011,15,  7, 15);
	} else {
		Make_Gr_inf(0,  7, 15);
	}
	float Lm_tr[3] = {2.5*0.125, 1.25*0.125,1.*0.125};
	float ml[3] = {1./6, 1./6, 1./6};
	for(int q =0; q<3; q++){
		Smthnss[q].max = Lm_tr[q];
		Smthnss[q].semi = Lm_tr[q]*ml[q];
	}

	// Model 
	typedef std::pair<int, int> hyperparameter_t;
	typedef std::vector< std::vector<int> > sym_t;
	typedef std::vector<int> invsym_t;
	typedef std::vector<int> labels_t;
	std::vector< hyperparameter_t > hyperparameters = {std::make_pair(1,0), /*std::make_pair(1,1),*/ std::make_pair(2,0), std::make_pair(3,0)};
	std::vector< std::map<hyperparameter_t , Model> > vecgms(2, std::map<hyperparameter_t , Model>());
	std::vector< std::map<hyperparameter_t , sym_t> > vecsyms(2, std::map<hyperparameter_t , sym_t>());
	std::vector< std::map<hyperparameter_t , invsym_t> > vecinvsyms(2, std::map<hyperparameter_t , invsym_t>());
	std::vector< std::map<hyperparameter_t , labels_t> > veclabels(2, std::map<hyperparameter_t , labels_t>());
	std::vector< std::map<hyperparameter_t, AEFInferType*> > vecaesptr(2, std::map<hyperparameter_t , AEFInferType*>());
	
	// Step1 : upper mrf - kmeans?
	png::image< png::rgb_pixel > img(nx, ny);
	for(int lr = 0; lr<2;lr++) {
		int * rez =  (lr)? dR : dL;
		BYTE * I_b = I_ims[lr];
		copy_m_D(m_cost, lr);
		int thr = 5;
		float alp =0.;
		
		Make_Gr_fl_buf(  thr, lr,m_cost, N_LB, alp );
		/*if(Sc_out == 1)*/
		K_mean_Flt_Add_new(I_b, sigma_c_mul, sigmax,sigmay, cl_Q, m_cost, n_bufs);
		
		// Need RVO here to speed it up. FIXME: use Return-by-value optimization.
		ucharArray2pngImage(img, I_b, nx, ny);
		
		// TODO - create different graphical models with symmetries
		
		for(int i=0; i<hyperparameters.size(); ++i) {
			// Get syms, invsyms
			vecsyms[lr][hyperparameters[i]] = std::vector< std::vector<int> >();
			vecinvsyms[lr][hyperparameters[i]] = std::vector<int>();
			getSymmetries(m_cost, img, nx, ny, nLabels, hyperparameters[i].first, hyperparameters[i].second, vecsyms[lr][hyperparameters[i]], vecinvsyms[lr][hyperparameters[i]]);
			vecgms[lr][hyperparameters[i]] = constructGraphicalModel(img, nLabels, m_cost, vecsyms[lr][hyperparameters[i]], vecinvsyms[lr][hyperparameters[i]]);
			// std::cerr << "invsyms size:" << vecinvsyms[lr][hyperparameters[i]].size() << "\n\n\n\n\n" << endl;
			// create inference thing and labels.
			veclabels[lr][hyperparameters[i]] = std::vector<int>(vecgms[lr][hyperparameters[i]].numberOfVariables());
			vecaesptr[lr][hyperparameters[i]] = new AEFInferType(vecgms[lr][hyperparameters[i]], params); // TODO - need a copy function.
		
			std::cerr << "lr=" << lr << std::endl;
			std::cerr << "hyperparam=" << hyperparameters[i].first << "," << hyperparameters[i].second << std::endl;
			std::cerr << "   gm size=" << vecgms[lr][hyperparameters[i]].numberOfVariables() << "v," << vecgms[lr][hyperparameters[i]].numberOfFactors() << "f" << endl;
		}
	} // end for(lr=0...)
	std::cout << "Models created. Now starting inference" << std::endl;
	// Step3 : inference
	AEFInferType::TimingVisitorType visitorl;
	AEFInferType::TimingVisitorType visitorr;
	hyperparameter_t active_hyperparameter = std::make_pair(1,0); // start coarse
	int active_hyperparameter_idx = 0;
	std::vector< std::vector<double> > energies(2, std::vector<double>(5,0.0));

	for(int iter=0; iter<100; ++iter) {
		// We don't really understand the magic that happens in the original MRF__z
		// hence, we'll just use the entire TSGO algorithm here...
		
		// IMPORTANT - create a function that does this.
		if(active_hyperparameter_idx<(hyperparameters.size()-1) && (detectPlateau1(iter, energies))) {
			// TODO - check if you need to switch.
			active_hyperparameter_idx++;
			std::cout << "Switching hyperparameters from ";
			std::cout << active_hyperparameter.first << " " << active_hyperparameter.second << " to ";
			active_hyperparameter = hyperparameters[active_hyperparameter_idx];
			std::cout << active_hyperparameter.first << " " << active_hyperparameter.second << std::endl;
			// copy labels over correctly from veclabels[lr][active_hyperparameter] to veclabels[lr][new_hyperparameter]
			for(int lr=0; lr<2; ++lr){
				updateLabelsWithChangedSymmetries(veclabels[lr][hyperparameters[active_hyperparameter_idx-1]], vecsyms[lr][hyperparameters[active_hyperparameter_idx-1]], veclabels[lr][active_hyperparameter], vecinvsyms[lr][active_hyperparameter]);
			}
		}
		std::cout << "Starting iteration=" << iter << " using parameters=(" << active_hyperparameter.first << "," << active_hyperparameter.second << ")" << std::endl;
		for(int lr = 0; lr<2;lr++) {
			int * rez =  (lr)? dR : dL;
			// Do a single step of inference on gms[lr]
			vecaesptr[lr][active_hyperparameter]->setStartingPoint( veclabels[lr][active_hyperparameter].begin() );
			std::cout << "lr: " << lr << "  active_hyperparameter: " << active_hyperparameter.first << "," << active_hyperparameter.second << std::endl;
			if(lr==0) vecaesptr[lr][active_hyperparameter]->infer(visitorl);
			else vecaesptr[lr][active_hyperparameter]->infer(visitorr);
			vecaesptr[lr][active_hyperparameter]->arg(veclabels[lr][active_hyperparameter]);
			
			// Copy the result into veclabels[lr], rez
			copyVecLabels2Rez(veclabels[lr][active_hyperparameter], vecinvsyms[lr][active_hyperparameter], rez);
			// save energy
			energies[lr][iter%(energies[lr].size())] = vecgms[lr][active_hyperparameter].evaluate(veclabels[lr][active_hyperparameter]);
		} // end for(lr=0...)
		
		// Postprocessing
		{
			sigma_img_get(dL);
		
			if(occ > c_lbl) {
				MkDspCorr_RL(  dL, dR);
			}

			{
				if(occ > c_lbl)
					for(int i=0; i<1; i++) {
						int r_w = 4;
						float cl_prc = sigma_img;//0.23;
						Gss_wnd_( r_w,  cl_prc, dL, 0, m_nLabels );
						Gss_wnd_( r_w,  cl_prc, dR, 1, m_nLabels );
					}
			}

			{
				if(occ > c_lbl) {
					MkDspCorr_RL(  dL, dR);
				}
				int r_w1 = 4;
				float cl_prc1 = sigma_img;//0.23;
				Gss_wnd_( r_w1,  cl_prc1, dL, 0, m_nLabels );
			}

			Outlr_reg_filter(0.001, 3, dL, 0, 10, 0.01);

			if(0)  {
				Intermed(dL, dR);
			}
			int r = 4;//(N_LB >16) ? 4:2;
			if(occ > c_lbl) {
				CorrectSub_new(dL, dR, /*intr_cst[0],*/ r, 1);
			} else {
				for(int i =0; i < m_nPixels; i++ ){
					dL[i] *= 1;
				}
			}
		}
		double badpixelscore;
		// if(iter==99)
		badpixelscore = getBadPixelScore(dL, truthimg, nx, ny, nLabels, threshold);
		std::cout << "On step " << iter << " error was " << badpixelscore << std::endl;
	}
	// Step4: clean up
	delete [] m_cost;
	delete [] Mrgnl;
}

void EDP::FilterGrDsp( unsigned char * sup, int *outc)
{

  m_G0 = new POINT_G [m_nPixels];
  float *out = new float  [m_nPixels];
	//--------------------------------------------------------------
   double  mean;

		for(int y = 0; y< m_height ; y++)
		for(int x = 0; x< m_width ; x++)
		{
		m_G0[IND(x,y)].m =0;
         if(x == m_width -1)m_G0[IND(x,y)].x =0;
		 else                           m_G0[IND(x,y)].x  =outc[IND(x+1,y)] - outc[IND(x,y)] ;
		 if(y == m_height -1)m_G0[IND(x,y)].y =0;
		 else                           m_G0[IND(x,y)].y  = outc[IND(x,y+1) ] - outc[IND(x,y)] ;
	
		}
		for(int y = 0; y< m_height ; y++)
		for(int x =  m_width/2; x< m_width/2 +20 ; x++)
		{
		m_G0[IND(x,y)].y =0; m_G0[IND(x,y)].x = 0;
		}

		mean =0; for(int p = 0; p < m_nPixels; p++) mean += (outc[p]); mean /= m_nPixels;
		MkRecFromG( out, mean, 0);
		for(int p = 0; p < m_nPixels; p++){ if(out[p]>m_nLabels-1)out[p]=m_nLabels-1;  if(out[p]<0)out[p]=0; outc[p] = (abs(round_fl(out[p])-out[p])>2)? round_fl(out[p]) : outc[p] ;}


//-----------------///
		 //for(int p = 0; p < m_nPixels; p++)  outc[p]= B_L_buf_pp_[p].x*m_nLabels/2;
delete [] m_G0;
delete [] out;
}
void   EDP:: FilterGrDsp( unsigned char * sup, int *outc, int thr)
{

  m_G0 = new POINT_G [m_nPixels];
  float *out = new float  [m_nPixels];
	//--------------------------------------------------------------
   double  mean;

		for(int y = 0; y< m_height ; y++)
		for(int x = 0; x< m_width ; x++)
		{
		m_G0[IND(x,y)].m =0;
         if(x == m_width -1)m_G0[IND(x,y)].x =0;
		 else                           m_G0[IND(x,y)].x  =outc[IND(x+1,y)] - outc[IND(x,y)] ;
		 if(y == m_height -1)m_G0[IND(x,y)].y =0;
		 else                           m_G0[IND(x,y)].y  = outc[IND(x,y+1) ] - outc[IND(x,y)] ;
	
		}
		for(int y = 0; y< m_height ; y++)
		for(int x = 0; x< m_width ; x++)
		{
	
			if(x != m_width -1)
			{
				int sum =0; for(int c=0; c <3; c++)sum += abs(sup[IND_IC(x+1,y,c)] - sup[IND_IC(x,y,c)]) ;
				if(sum <= thr)m_G0[IND(x,y)].x = 0;
			}
            if(y != m_width -1)
			{
				int sum =0; for(int c=0; c <3; c++)sum += abs(sup[IND_IC(x,y+1,c)] - sup[IND_IC(x,y,c)]) ;
				if(sum <= thr)m_G0[IND(x,y)].y = 0;
			}
		}

		mean =0; for(int p = 0; p < m_nPixels; p++) mean += (outc[p]); mean /= m_nPixels;
		MkRecFromG( out, mean, 0);
		for(int p = 0; p < m_nPixels; p++){ if(out[p]>255)out[p]=255;  if(out[p]<0)out[p]=0; outc[p] =round_fl(out[p]);}


//-----------------///
		 //for(int p = 0; p < m_nPixels; p++)  outc[p]= B_L_buf_pp_[p].x*m_nLabels/2;
delete [] m_G0;
delete [] out;
}
void   EDP:: FilterGrC( unsigned char * sup, unsigned char *outc, int thr)
{

  m_G0 = new POINT_G [m_nPixels];
  float *out = new float  [m_nPixels];
	//--------------------------------------------------------------
   double  mean;

		for(int y = 0; y< m_height ; y++)
		for(int x = 0; x< m_width ; x++)
		{
		m_G0[IND(x,y)].m =0;
         if(x == m_width -1)m_G0[IND(x,y)].x =0;
		 else                           m_G0[IND(x,y)].x  =outc[IND(x+1,y)] - outc[IND(x,y)] ;
		 if(y == m_height -1)m_G0[IND(x,y)].y =0;
		 else                           m_G0[IND(x,y)].y  = outc[IND(x,y+1) ] - outc[IND(x,y)] ;
	
		}
		for(int y = 0; y< m_height ; y++)
		for(int x = 0; x< m_width ; x++)
		{
	
			if(x != m_width -1)
			{
				int sum =0; for(int c=0; c <3; c++)sum += abs(sup[IND_IC(x+1,y,c)] - sup[IND_IC(x,y,c)]) ;
				if(sum <= thr)m_G0[IND(x,y)].x = 0;
			}
            if(y != m_width -1)
			{
				int sum =0; for(int c=0; c <3; c++)sum += abs(sup[IND_IC(x,y+1,c)] - sup[IND_IC(x,y,c)]) ;
				if(sum <= thr)m_G0[IND(x,y)].y = 0;
			}
		}

		mean =0; for(int p = 0; p < m_nPixels; p++) mean += (outc[p]); mean /= m_nPixels;
		MkRecFromG( out, mean, 0);
		for(int p = 0; p < m_nPixels; p++){ if(out[p]>255)out[p]=255;  if(out[p]<0)out[p]=0; outc[p] =round_fl(out[p]);}


//-----------------///
		 //for(int p = 0; p < m_nPixels; p++)  outc[p]= B_L_buf_pp_[p].x*m_nLabels/2;
delete [] m_G0;
delete [] out;
}
void   EDP:: FilterBF_C( unsigned char *outc,float sig, float sigc, int cl_Q)
{


  double *out = new double  [m_nPixels*3]; FOR_IxP(3)out[i]=outc[i];
	//--------------------------------------------------------------

        K_mean_Flt_Add_new( outc, sigc, sig,sig, cl_Q, out, 3);

		for(int p = 0; p < m_nPixels*3; p++){ if(out[p]>255)out[p]=255;  if(out[p]<0)out[p]=0; outc[p] =round_fl(out[p]);}


//-----------------///

delete [] out;
}
void   EDP:: FilterGrC( unsigned char * sup, double *outc, int thr)
{

  m_G0 = new POINT_G [m_nPixels];
  float *out = new float  [m_nPixels];
	//--------------------------------------------------------------
   double  mean;

		for(int y = 0; y< m_height ; y++)
		for(int x = 0; x< m_width ; x++)
		{
		m_G0[IND(x,y)].m =0;
         if(x == m_width -1)m_G0[IND(x,y)].x =0;
		 else                           m_G0[IND(x,y)].x  =outc[IND(x+1,y)] - outc[IND(x,y)] ;
		 if(y == m_height -1)m_G0[IND(x,y)].y =0;
		 else                           m_G0[IND(x,y)].y  = outc[IND(x,y+1) ] - outc[IND(x,y)] ;
	
		}
		for(int y = 0; y< m_height ; y++)
		for(int x = 0; x< m_width ; x++)
		{
	
			if(x != m_width -1)
			{
				int sum =0; for(int c=0; c <3; c++)sum += abs(sup[IND_IC(x+1,y,c)] - sup[IND_IC(x,y,c)]) ;
				if(sum <= thr)m_G0[IND(x,y)].x = 0;
			}
            if(y != m_width -1)
			{
				int sum =0; for(int c=0; c <3; c++)sum += abs(sup[IND_IC(x,y+1,c)] - sup[IND_IC(x,y,c)]) ;
				if(sum <= thr)m_G0[IND(x,y)].y = 0;
			}
		}

		mean =0; for(int p = 0; p < m_nPixels; p++) mean += (outc[p]); mean /= m_nPixels;
		MkRecFromG( out, mean, 0);
		for(int p = 0; p < m_nPixels; p++){ outc[p] =out[p];}


//-----------------///
		 //for(int p = 0; p < m_nPixels; p++)  outc[p]= B_L_buf_pp_[p].x*m_nLabels/2;
delete [] m_G0;
delete [] out;
}
void EDP::PrjectSlt(int st_intrp, int *outL, int *outR)
{

	 float dirI;
	 int cmp_im;
	 if(this->nm_Ims-3){int n_intr = ((this->nm_Ims-3)+1)/2*2-1;
	 /*float i_sp  =  1./(n_intr-3); */float i_sp;
	 if(n_intr ==1){i_sp =  0.5; dirI = 0.5;cmp_im = 2;}
	 if(n_intr ==3){i_sp =  0.5; dirI = 0; cmp_im =3;}
	 if(n_intr >3 ){i_sp =  1./(n_intr-3);dirI = - i_sp; cmp_im = (n_intr-3)/2+3;}
	
	if(1)for(int im =0; im<n_intr; dirI += i_sp, im++ )
	 for(int i = 0; i< m_height; i++) ProjectImCrvNew(i, this->I_ims[2+im], dirI);
     //int r =2;
     //for(int im =0; im<n_intr;  im++ )color_med_( r,this->I_ims[2+im] );

	 }
	if(0){for(int i = 0; i< m_height; i++)ProjectDCrv(i, 1, outR);
	                       for(int i = 0; i< m_height; i++)ProjectDCrv(i, 0, outL);}
	


}

void EDP::PrjectSltG( int *outL)
{

	 float dirI;
	 int cmp_im;
	 int n_intr;
	 if(this->nm_Ims-3){n_intr = ((this->nm_Ims-3)+1)/2*2-1;
	 /*float i_sp  =  1./(n_intr-3); */float i_sp;
	 if(n_intr ==1){i_sp =  0.5; dirI = 0.5;cmp_im = 2;}
	 if(n_intr ==3){i_sp =  0.5; dirI = 0; cmp_im =3;}
	 if(n_intr >3 ){i_sp =  1./(n_intr-3);dirI = - i_sp; cmp_im = (n_intr-3)/2+3;}
	
	 for(int im =0; im<n_intr; dirI += i_sp, im++ )
	 for(int i = 0; i< m_height; i++)
	 ProjectImCrv(i, this->I_ims[2+im], dirI);

	 }

	
	  for(int i = 0; i< m_height; i++)ProjectDCrv(i, 0, outL);
     memcpy(I_ims[n_intr], I_ims[cmp_im], m_nPixels*3);
}

void EDP::CorrectSub(int *outL, int *outR, int dh, int sc)
{

	
	int itr =4;  float strt_dsp =0.5;  subpix_cst_tune = 0.25; float f2 = 1;  int dpth=dh; int ilr =0;
	for(int i =0; i<m_nPixels; i++)DspFL[i] = outL[i];
	SubpixelAcc( itr,strt_dsp, f2, dpth, ilr);
	this->L_occ(outL);
	for(int i =0; i<m_nPixels; i++){ outL[i] =  (!Mask_RL[0][i])? round_fl( DspFL[i]*sc) : outL[i]*sc;
	if(outL[i]<0)outL[i]=0; if(outL[i]>(m_nLabels-1)*sc )outL[i]=(m_nLabels-1)*sc ;}


	
}

void EDP::CorrectSub_new(int *outL, int *outR, /*double * m_cost,*/ int r,  int sc) {

	subpix_cnst = 1024;
	for(int i =0; i<4; i++ )
		fine_BC[i] = new float [subpix_cnst];

	for(int i=0; i < subpix_cnst; i++) {
		float x =(float)i/subpix_cnst;
		fine_BC[0][i] =  0.5*(x-1)*(x-2)*(x+1); //f(0)
		fine_BC[1][i] = -0.5*(x)*(x-2)*(x+1); // f(1)		
		fine_BC[3][i] = -1./6.*(x)*(x-2)*(x-1); // f(-1)
		fine_BC[2][i] =  1./6.*(x)*(x+1)*(x-1); //f(2)
	}
	///////////////////////////////////////////////////////////
	
	int * rez = new int [N_PX];
	int dm = r +1;
	FOR_PX_p {
		outL[p] *= r;
	}
	for(int i =0; i < 1; i++){ // really, what?
		srand( 27644437 );
		FOR_PX_p{
			int rdl  = (abs(rand()))%dm - r/2;
			rez[p] =  outL[p] + rdl;
			
			if (rez[p]<0) rez[p] = 0;
			
			if(rez[p]>=N_LB*r) rez[p] = N_LB*r -1;
		}
		int r_w = 11;
		float cl_prc = (N_LB> 16)? ((N_LB>20)? 3.5 : 2.5):0.5;//0.23;
		Gss_wnd_ML( r_w,  cl_prc ,outL, rez, m_nLabels*r , r);
		FOR_PX_p {
			outL[p] = rez[p];
		}
	}

	FOR_PX_p {
		outL[p] = rez[p];
	}

	delete [] rez;
	for(int i =0; i<4; i++ )
		delete [] fine_BC[i] ;
}

void EDP::Intermed(int *outL, int *outR)
{

	this->Cost_mean = 1000;
	m_D_pR = new float [m_height*Tau];
	m_D_pL = new float [m_height*Tau];
   	this->MkPnlCst_RL(outL, outR);
	Mk_DRL(outL, outR,1);
	//this->drwCost(231);
	int r_w; float cl_prc;
	//-------------------------------------
	
	int rr = 8;     float clp = 0.1;
	//filterDSI_tau(I_ims[0],    rr,    clp);
	for(int i = 0; i< m_height; i++) { findSb(i); PaintCrv(i);}
	PrjectSlt(1,  outL, outR);
	
	delete []  m_D_pR ;
    delete []  m_D_pL;
}

float EDP::L_occ(int *outL )
{

	FOR_PX_p {Mask_RL[0][p] = 0; }
		for(int y =0; y<m_height; y++)
		{ int dp;
	    for(int x =0; x< m_width; x++)
	    { int d = outL[x+ y*m_width];
		  if(d > x ) {Mask_RL[0][x+ y*m_width] = 1;}
		  else
		  {
			  int lim = (m_width - x < N_LB  - 1 - d)  ?  m_width - x : N_LB  - 1 - d ;
			  for(int dd =1; dd < lim; dd++)if(dp=outL[x+ y*m_width + dd] == d+dd){Mask_RL[0][x+ y*m_width] = 1; dd = lim+1;}
		  }

		}}
	
       float ret = 0;
		FOR_PX_p{ret += I_ims[2][p + N_PX] = I_ims[2][p + 2*N_PX] = I_ims[2][p] = Mask_RL[0][p]*255;}
		return ret/N_PX/255;
}

void EDP::GT_gnr(int *outL , int *outR)
{
	Cost_mean = 1000;
	int sc = 1;
	for(int i =0; i<m_nPixels; i++){outL[i] /=sc; outR[i] /=sc;}
	m_D_pR = new float [m_height*Tau];
	m_D_pL = new float [m_height*Tau];
   	this->MkPnlCst_RL(outL, outR);
	Mk_DRL(outL, outR, 1);
	
	//-------------------------------------
	
	for(int i = 0; i< m_height; i++) { findSb(i); PaintCrv(i);}
	PrjectSlt(1,  outL, outR);
	MskRL_Corr(outL, 0);//RL seg
	MskRL_Corr(outR, 1);//RL seg
  /*  this->drwCost(231);*/
	for(int i =0; i<m_nPixels; i++){outL[i] *=sc; outR[i] *=sc;}

	delete []  m_D_pR ;
    delete []  m_D_pL;
}
EDP::REAL EDP::getErrGT()
{
	int sum =0; int cnt =0;
	int sumn =0; int cntn =0;
	for(int i =0; i < m_nPixels; i++ )
	{
		int diff = abs(gtL_answer[i]-m_answer[i]/**out_sc*/);
		
						if(gtL_answer[i])
						{	
							if(!Mask_RL[0][i]){cntn++; if(diff>1){sumn++; /*m_answer[i] = 0;*/ } }
							cnt++;
							if(diff>1){sum++; /*m_answer[i] = 0;*/}
							
						}
	}
 //err_nonnoc = (float)sumn*100/cntn;
 return (((REAL)sum*100.)/cnt);
}

EDP::REAL EDP::getErrGTI()
{
	 int cmp_im;  float dirI;
	 double sum =0; int cnt =0;
	 if(this->nm_Ims-3){int n_intr = ((this->nm_Ims-3)+1)/2*2-1;
	 /*float i_sp  =  1./(n_intr-3); */float i_sp;
	 if(n_intr ==1){i_sp =  0.5; dirI = 0.5;cmp_im = 2;}
	 if(n_intr ==3){i_sp =  0.5; dirI = 0; cmp_im =3;}
	 if(n_intr >3 ){i_sp =  1./(n_intr-3);dirI = - i_sp; cmp_im = (n_intr-3)/2+3;}
	 unsigned char * im1 = I_ims[cmp_im]; unsigned char * im2 =  I_ims[nm_Ims-1];
	
	for(int i =0; i < m_nPixels*3; i++ )
	{
		int x = (i%m_nPixels)%m_width;
		if(x>m_nLabels/2&& x<m_width - m_nLabels/2){ sum += (im1[i]-im2[i])*(im1[i]-im2[i]); cnt++;}

		//im1[i]=  abs(im1[i]-im2[i]);

	}}
	 else cnt=1;
	double mse = sqrt(sum/(cnt));
	mse = log10((double)255/mse)*20;
 return mse;
	 }

	  EDP::REAL EDP::getErrGTI(float *im1, float * im2)
{
	 int cmp_im;  float dirI;
	 double sum =0; int cnt =0;

	
	for(int i =0; i < m_nPixels*3; i++ )
	{
		int x = (i%m_nPixels)%m_width;
		{ sum += (im1[i]-im2[i])*(im1[i]-im2[i]); cnt++;}


	}
	double mse = sqrt(sum/(cnt));
	mse = log10((double)255/mse)*20;
 return mse;
	 }

void EDP::XL_XR_Tau()
{
///// XX_TAU ///////////////////////////////
	int sub_cnst = m_nLabels*(m_nLabels -1)/2;
	Tau = m_width*m_nLabels - sub_cnst;
	XL_tau  = new short int [Tau];
	XR_tau  = new short int [Tau];
	tau_XL_d = new int [m_width*m_nLabels];
	tau_XR_d = new int [m_width*m_nLabels];
	//---------------------------
	for(int i =0; i <m_width*m_nLabels; i++) {
		tau_XL_d[i] = tau_XR_d[i] =-1;
	}
	//----------------------------
	V_tau_bck  =  new int [Tau];
	R_tau_bck  =  new int [Tau];
	L_tau_bck  =  new int [Tau];
	V_tau_frw  =  new int [Tau];
	R_tau_frw  =  new int [Tau];
	L_tau_frw  =  new int [Tau];
//________________________________________
	XL_tau[0] = 0;    XR_tau[0] = 0;   tau_XL_d[0]= 0 ;						   tau_XR_d[0]= 0 ;
	XL_tau[1] = 1;    XR_tau[1] = 0;   tau_XL_d[m_nLabels +1]= 1;    tau_XR_d[1]=  1;
	for(int t=2; t<Tau; t++) {
		int dl, dr;
		int d;
		int xl = XL_tau[t-1], xr = XR_tau[t-1];
		//----------- find dl,dr tau then xl(tua+1)
		if( xl-xr > 0) {
			dl= 0;    dr=1;
		}else {
			dl=1;
			dr = (xr<m_nLabels-1) ? -xr : 2-m_nLabels;
		}
		xl += dl;
		xr += dr;
		d = xl-xr;

		XL_tau[t] = xl;
		XR_tau[t] = xr;
		tau_XL_d[ xl*m_nLabels + d]= t;
		tau_XR_d[ xr*m_nLabels + d]= t;
     }

//////////////////////
	for(int t=0; t<Tau; t++) {
		int d;
		int xl = XL_tau[t], xr = XR_tau[t];
		d = xl-xr;
		V_tau_bck[t]  = (xl-1<0)? -1:((xr-1<0)? -1:tau_XL_d[d+ (xl-1)*m_nLabels]);
		V_tau_frw[t] = (xl+1>=m_width)? -1:((xr+1>=m_width)? -1:tau_XL_d[d+ (xl+1)*m_nLabels]);
		L_tau_bck[t]  = (xl-1<0)? -1:((d<1)? -1:tau_XL_d[d-1+ (xl-1)*m_nLabels]);
		L_tau_frw[t] = (xl+1>=m_width)? -1:((d>=m_nLabels-1)? -1:tau_XL_d[d+1+ (xl+1)*m_nLabels]);
		R_tau_bck[t]  = (xr-1<0)? -1:((d>=m_nLabels-1)? -1:tau_XR_d[d+1+ (xr-1)*m_nLabels]);
		R_tau_frw[t] = (xr+1>=m_width)? -1:((d<1)? -1:tau_XR_d[d-1 + (xr+1)*m_nLabels]);
	}


}
void EDP::drwCost(int yy)
{

	for(int i =0; i<m_nPixels*3; i++)
		I_ims[2][i]=100 + ((i/m_nPixels)%3)*30;

	if(1) {
		for(int i =0; i<Tau; i++) {
			int shf = m_height*m_width;
			int x = XL_tau[i];
			int y = x - XR_tau[i];
			x -= y/2;
			//int vl = m_D[i+yy*Tau];
			y = m_height - y-1; /*if(y<0) y= 0;*/
			int vl = (int)m_D[yy*Tau+i];//pow(m_D[yy*Tau+i], (float)0.7);
			int vl1 =(m_D[yy*Tau+i] < 255) ? m_D[yy*Tau+i]: 255; int vl2 =(vl/3 < 255) ? vl/3: 255; int vl3 =(vl/6 < 255) ? vl/6: 255;//(pow((float)vl,(float)0.45)*21); if(vl1>255) vl1=255;
			for(int c =0; c<3; c++ )this->I_ims[2][x + y*m_width + shf*c] =(c<2) ? ((c<1)? vl3:vl2) : vl1;
		}
	}
	//for(int i =0; i<m_nPixels*3; i++)I_ims[2][i]=255;
	if(1) {
		for(int i =0; i<Tau; i++) {
			int shf = m_height*m_width;
			int x = XL_tau[i];
			int y = x - XR_tau[i];
			x -= y/2;
			//int vl = m_D[i+yy*Tau];
			y = m_height - y-1; /*if(y<0) y= 0;*/
			int vl1,vl2, vl3, vlL= (int)m_D_pL[yy*Tau+i]; int vlR= (int)m_D_pR[yy*Tau+i];
			//if(vlL<Cost_mean) {vl1 = 0; vl2=(int)((float)vlL/Cost_mean*70); vl3=((float)vlL/Cost_mean*255);} else {vl1 =(int)((float)vlL/Cost_mean*127); vl2=(int)((float)vlL/Cost_mean*127); vl3=(int)((float)vlL/Cost_mean*127);}
			//f(vlR<Cost_mean) {vl1 = 0; vl2=50; vl3=220;} else {vl1 =(int)(vlL/Cost_mean*127); vl2=(int)(vlL/Cost_mean*127); vl3=(int)(vlL/Cost_mean*127);}
			if(vlL<Cost_mean/2||vlR<Cost_mean/2){ 
				if(vlL<Cost_mean/2){vl1 = 0; vl2=200; vl3=200;} 
				if(vlR<Cost_mean/2){vl1 = 0; vl2=0; vl3=255;}
				for(int c =0; c<3; c++ )
					this->I_ims[2][x + y*m_width + shf*c] =(c<2) ? ((c<1)? vl3:vl2) : vl1;
			}
		}
	}

	if(0) {
		for(int i =0; i<T_max_y[yy]; i++){
			int shf = m_height*m_width;
			int x = XL_tau[Sltn_crv[WpH*yy+i].tau];
			int d =Sltn_crv[WpH*yy+i].d;
			int y =x - XR_tau[Sltn_crv[WpH*yy+i].tau];
			x -= y/2;
			y = m_height - y-1; /*if(y<0) y= 0;*/
			int v1 = 0, v2 = 255, v3 = 0;
			for(int c =0; c<3; c++ )
				this->I_ims[2][x + y*m_width + shf*c] =(c<2) ? ((c<1)? v1:v2) : v3;
		}
	}


}
void EDP::PaintCrv(int yy )
{
	Col_outl_thr = 120;
	int crv_lth = this->T_max_y[yy];
	SLTN_CRV * Slcv = &this->Sltn_crv[yy*WpH];
	unsigned char * im_l = &this->I_ims[0][yy*m_width];
	unsigned char * im_r = &this->I_ims[1][yy*m_width];
	int * seg_i = new int [crv_lth];
	unsigned char * seg_lb = new unsigned char [2*m_width];
	unsigned char * seg_lb2 = new unsigned char [2*m_width];
	short int * seg_mx_d  = new short  [2*m_width];
	short int * seg_mx_d2  = new short  [2*m_width];
	int seg_ii=0;
	for(int i =0; i<crv_lth; i++)
	{
		unsigned char lb_i = (Slcv[i].lb)%8;
		int tau_i = Slcv[i].tau;
		int xl = XL_tau[tau_i];
		int xr = XR_tau[tau_i];
		int d = xl-xr;
        // Color
		int sumv = 0;  Slcv[i].rgb[3] =0;
		for( int c = 0; c<3; c ++)
		{
		if(lb_i == MY_VIS ) {sumv += abs(im_l[xl + c*m_nPixels] - im_r[xr+ c*m_nPixels]); Slcv[i].rgb[c] = (im_l[xl + c*m_nPixels] + im_r[xr+ c*m_nPixels])/2;}
		if(lb_i==MY_LFT||lb_i== MY_LFT_A)     Slcv[i].rgb[c] = im_l[xl + c*m_nPixels] ;
		if(lb_i==MY_RGT||lb_i== MY_RGT_A)   Slcv[i].rgb[c] = im_r[xr + c*m_nPixels] ;
		}
		if(lb_i == MY_VIS )Slcv[i].rgb[3] =(sumv<255)? sumv: 255;
         // Color end
		//segment fill
		if(i == 0){
		seg_i[i]=seg_ii;
		seg_lb[seg_ii] = lb_i;
		seg_mx_d[seg_ii]=d;
		}
		else {
			if((Slcv[i-1].lb)%8 != lb_i)
			{
		    seg_i[i]=++seg_ii;
		    seg_lb[seg_ii] = lb_i;
		    seg_mx_d[seg_ii]=d;
			}
			else
			{
		    seg_i[i]=seg_ii;
			if(lb_i==MY_RGT_A||lb_i== MY_LFT_A)seg_mx_d[seg_ii]=(seg_mx_d[seg_ii]>d)? seg_mx_d[seg_ii]:d ;
			if(lb_i==MY_RGT||lb_i== MY_LFT)seg_mx_d[seg_ii]=(seg_mx_d[seg_ii]<d)? seg_mx_d[seg_ii]:d ;
			}
		}

	}
//---------------------------------------------------
for(int i = 0; i<=seg_ii;i++)
{
	int im = (i>0) ? i-1:i;
	int ip =(i<seg_ii)? i+1:i;
	int lb_i =  seg_lb[i];
	if(lb_i==MY_VIS)
	{ seg_mx_d2[i] = seg_mx_d[i]; seg_lb2[i] = lb_i%8;}
	if(lb_i==MY_RGT_A||lb_i== MY_LFT_A)
	{
		seg_mx_d2[i]=(seg_mx_d[im]> seg_mx_d[i])? seg_mx_d[im] :seg_mx_d[i];
		seg_mx_d2[i] =(seg_mx_d2[i]> seg_mx_d[ip])? seg_mx_d2[i] :seg_mx_d[ip];
		seg_lb2[i] = lb_i%8;
	}
	if(lb_i==MY_RGT||lb_i== MY_LFT)
	{
		seg_mx_d2[i]=(seg_mx_d[im]<seg_mx_d[i])? seg_mx_d[im] :seg_mx_d[i];
		seg_mx_d2[i] =(seg_mx_d2[i]< seg_mx_d[ip])? seg_mx_d2[i] :seg_mx_d[ip];
		seg_lb2[i] = lb_i%8;
     if(lb_i==MY_RGT&&((seg_lb[ip]%8== MY_LFT)||(seg_lb[im]%8== MY_LFT))) seg_lb2[i] = lb_i%8+8;
	 if(lb_i==MY_LFT&&((seg_lb[ip]%8== MY_RGT)||(seg_lb[im]%8== MY_RGT))) seg_lb2[i] = lb_i%8+8;
	}
}
for(int i =0; i<crv_lth; i++)
{
      Slcv[i].lb = seg_lb2[seg_i[i]];
	  Slcv[i].d= seg_mx_d2[seg_i[i]];
}
		/////////////////////////////
		delete [] seg_mx_d2;
	    delete [] seg_i;
	    delete [] seg_lb;
	    delete [] seg_mx_d;
		delete [] seg_lb2;

}
void EDP::ProjectImCrv(int yy, unsigned char *img, float dir)
{
	int crv_lth = this->T_max_y[yy];
	SLTN_CRV * Slcv = &this->Sltn_crv[yy*WpH];
	unsigned char * img_st = &img[yy*m_width];
	unsigned char * im_max1 =  new unsigned char [m_width*4]; memset(im_max1,0,m_width*4);
	short int  * i_d_max1 =  new short int  [m_width];  for(int i =0; i<m_width; i++){i_d_max1[i] = -1; }


//----------------------------------------------------
	if(dir==0||dir==1){
	for(int i =0; i<crv_lth; i++)
	{
		unsigned char lb_i = (Slcv[i].lb)%8;
		int tau_i = Slcv[i].tau;
		int xl = XL_tau[tau_i];
		int xr = XR_tau[tau_i];
		int d = xl-xr;
		//////////// status LFT ///////////////
		if(dir==0)
		{
			int x_im =xl;
			if (i_d_max1[x_im] <=d)
			{
				i_d_max1[x_im]=d;
				for( int c = 0; c<4; c ++)im_max1[x_im +c*m_width] = Slcv[i].rgb[c];
			}
		 }
      //////////// status RGT /////////////
		if(dir==1)
		{
			int x_im =xr;
			if (i_d_max1[x_im] <= d)
			{
				i_d_max1[x_im]=d;
				for( int c = 0; c<4; c ++)im_max1[x_im +c*m_width] = Slcv[i].rgb[c];
			}
		 }

        /////////////////////////////////////////////////////
	}}
	
	//---------------
	else {
	for(int i =0; i<crv_lth-1; i++)
	{
		//unsigned char lb_i1 = (Slcv[i].lb)%8;
		//unsigned char lb_i2 = (Slcv[i+1].lb)%8;
		unsigned char *rgb1, *rgbt1 = &Slcv[i].rgb[0];
		unsigned char *rgb2 ,*rgbt2 = &Slcv[i+1].rgb[0];
		int tau_i1 = Slcv[i].tau;
		int tau_i2 = Slcv[i+1].tau;
		int xl1 = XL_tau[tau_i1];
		int xr1 = XR_tau[tau_i1];
		int xl2 = XL_tau[tau_i2];
		int xr2 = XR_tau[tau_i2];
		int d1 =  xl1-xr1;
		int d2 =  xl1-xr2;
		 int dm = (xl1-xr1> xl2-xr2) ? xl1-xr1: xl2-xr2;
//-------------------------------------------------------------
		if(i==0)for( int c = 0; c<4; c ++){im_max1[c*m_width] = rgbt1[c];}
		if(i==crv_lth-2)for( int c = 0; c<4; c ++){im_max1[c*m_width +m_width-1] = rgbt2[c];}
			//
			float  x_imf1, xf1 = (float)xl1- dir*d1;
			float  x_imf2, xf2 = (float)xl2- dir*d2;
			float df1, df2;
			if(xf1<=xf2){x_imf1 = xf1; x_imf2 = xf2; rgb1 =rgbt1; rgb2 = rgbt2; df1=d1; df2 = d2;}
			else            {x_imf1 = xf2; x_imf2 = xf1; rgb1 =rgbt2; rgb2 = rgbt1;df1=d2; df2 = d1;}
			int xim1 =(x_imf1!=x_imf2)? (int)x_imf1+1 : (int)x_imf1;
			//int xim2 = (int)x_imf2;
			//int n = xim2 - xim1;
			/*int add =(x_imf1!= x_imf2)? 0:1;*/	
			float wcc;
			if(dir<0||dir>0)wcc = 1;
			if(dir>0&&dir<1)wcc=0.5;
			if(dir==0.5)wcc=0;
			for(int ii =xim1; ii<=x_imf2+wcc; ii++)
			{
				int xx = ii;
				if (xx>=0&&xx<m_width){
					
				float w1 = fabs((float)x_imf2 - xx); float  w2 = fabs((float)x_imf1 - xx);
				if(w1+w2==0)
				{float vl = (df1 + df2)/2; if (i_d_max1[xx] <=round_fl(vl)){ i_d_max1[xx]=round_fl(vl);
					for( int c = 0; c<4; c ++)
					{float vl = ((float)(rgb1[c] + rgb2[c])/2); im_max1[xx +c*m_width] = (vl>255) ? 255: round_fl(vl);}}
				}
				else {float vl = (w1*df1 + w2*df2)/(w1+w2); if (i_d_max1[xx] <=round_fl(vl)){ i_d_max1[xx]=round_fl(vl);
					for( int c = 0; c<4; c ++)
					{float vl =(w1*rgb1[c] + w2*rgb2[c])/(w1+w2) ; im_max1[xx +c*m_width] = (vl>255) ? 255: round_fl(vl); } }
				}
					
					
				
				}
			}	
	}}

	for(int i =0; i<m_width; i++) for( int c = 0; c<4; c ++) img_st[i + c*m_nPixels] = im_max1[i +c*m_width];

		/////////////////////////////
		delete [] im_max1;
        delete [] i_d_max1;


}
void EDP::ProjectImCrvNew(int yy, unsigned char *img, float dir) {
	int crv_lth = this->T_max_y[yy];
	SLTN_CRV * Slcv = &this->Sltn_crv[yy*WpH];
	unsigned char * img_st = &img[yy*m_width];
	unsigned char * im_max1 =  new unsigned char [m_width*4]; memset(im_max1,0,m_width*4);
	short int  * i_d_max1 =  new short int  [m_width];  for(int i =0; i<m_width; i++){i_d_max1[i] = -1; }


//----------------------------------------------------
	if(dir==0||dir==1){
	for(int i =0; i<crv_lth; i++)
	{
		unsigned char lb_i = (Slcv[i].lb)%8;
		int tau_i = Slcv[i].tau;
		int xl = XL_tau[tau_i];
		int xr = XR_tau[tau_i];
		int d = xl-xr;
		//////////// status LFT ///////////////
		if(dir==0)
		{
			int x_im =xl;
			if (i_d_max1[x_im] <=d)
			{
				i_d_max1[x_im]=d;
				for( int c = 0; c<4; c ++)im_max1[x_im +c*m_width] = Slcv[i].rgb[c];
			}
		 }
      //////////// status RGT /////////////
		if(dir==1)
		{
			int x_im =xr;
			if (i_d_max1[x_im] <= d)
			{
				i_d_max1[x_im]=d;
				for( int c = 0; c<4; c ++)im_max1[x_im +c*m_width] = Slcv[i].rgb[c];
			}
		 }

        /////////////////////////////////////////////////////
	}}
	
	//---------------
	else {
	for(int i =0; i<crv_lth-1; i++)
	{
		unsigned char lb_i1 = (Slcv[i].lb)%8;
		unsigned char lb_i2 = (Slcv[i+1].lb)%8;
		unsigned char *rgb1   = &Slcv[i].rgb[0];
		unsigned char *rgb2  = &Slcv[i+1].rgb[0];
		int tau_i1 = Slcv[i].tau;
		int tau_i2 = Slcv[i+1].tau;
		int xl1 = XL_tau[tau_i1];
		int xr1 = XR_tau[tau_i1];
		int xl2 = XL_tau[tau_i2];
		int xr2 = XR_tau[tau_i2];
		int d_c1 = Slcv[i].d;
		int d_c2 = Slcv[i+1].d;
		int d1 =   xl1-xr1;
		int d2 =   xl2-xr2;
		 //int dm = (xl1-xr1> xl2-xr2) ? xl1-xr1: xl2-xr2;
//-------------------------------------------------------------
			float  x_imf1= (float)xl1- dir*d1;
			float  x_imf2= (float)xl2- dir*d2;
			int st_p, p, fn_p;
			float w1, w2;

		  if(dir>=0&&dir<=1)
		  {
			  IF_1_L  {d1 = d_c1;  x_imf1= (float)xl1- dir*d1; }
			  IF_2_L  {d2 = d_c2; x_imf2= (float)xl2- dir*d2; }
			  IF_1_R  {x_imf1= (float)xl1 + d_c1 - d1- dir*d_c1; d1 = d_c1;}
			  IF_2_R { x_imf2= (float)xl2 + d_c2 - d2 - dir*d_c2; d2 = d_c2; }
		  }
           if(d1==d2)//////////////////// d1 == d2
		   {
			               	if(x_imf1 -(int) x_imf1 != 0)
                           {	
							p = (int)x_imf2;
							if(p>=0&&p<m_width){ w1 = fabs(x_imf2 - p); w2 = fabs(x_imf1 - p);
							IF_1_2_V { if(i_d_max1[p] <m_nLabels + d1){i_d_max1[p] = m_nLabels + d1; FOR_C(4, p);} }
							else {if(i_d_max1[p] < d1){i_d_max1[p] = d1; FOR_C(4, p);}}}
							}
							else
						   {	
							st_p = (int)x_imf1; fn_p = (int)x_imf2;
							for(p=st_p; p<=fn_p; p++)
							if(p>=0&&p<m_width){ w1 = fabs(x_imf2 - p); w2 = fabs(x_imf1 - p);
							IF_1_2_V { if(i_d_max1[p] <m_nLabels + d1){i_d_max1[p] = m_nLabels + d1; FOR_C(4, p);} }
							else {if(i_d_max1[p] < d1){i_d_max1[p] = d1; FOR_C(4, p);}}}
							}
		   }////////////////////////// end d1 ==d2
		   else //////////////////// d1 == d2
		   {			
               st_p = FLR_FL(x_imf1);  fn_p = FLR_FL(x_imf2);
			   int inc = (fn_p - st_p >=0) ? 1:-1; int dmin = (d1<d2) ? d1:d2;

							for(p=st_p; p<=fn_p; p+=inc)
							if(p>=0&&p<m_width){
								                                       w1 = fabs(x_imf2 - p); w2 = fabs(x_imf1 - p);
							if(p == st_p){ IF_1_V  { if(i_d_max1[p]<m_nLabels + d1){i_d_max1[p] = m_nLabels + d1;  if(w1+w2)FOR_C(4, p); else FOR_C1(4, p);} }
							                         else       { if(i_d_max1[p]< d1){i_d_max1[p] = d1;  if(w1+w2)FOR_C(4, p); else  FOR_C1(4, p);} }  }
							if(p == fn_p){ IF_2_V  { if(i_d_max1[p]<m_nLabels + d2){i_d_max1[p] = m_nLabels + d2; if(w1+w2)FOR_C(4, p); else FOR_C2(4, p);} }
							                         else        { if(i_d_max1[p]< d2){i_d_max1[p] = d2;  if(w1+w2)FOR_C(4, p); else  FOR_C2(4, p);} }  }
							if(p!=fn_p&&p!=st_p){  if(i_d_max1[p]<dmin){i_d_max1[p] = dmin;  FOR_C(4, p);  }}
							                                        }

		   }////////////////////////// end d1 ==d2

	}}

	for(int i =0; i<m_width; i++) for( int c = 0; c<4; c ++) img_st[i + c*m_nPixels] = im_max1[i +c*m_width];

		/////////////////////////////
		delete [] im_max1;
        delete [] i_d_max1;


}
void EDP::ProjectImCrvNew_(int yy, unsigned char *img, float dir) {
	int crv_lth = this->T_max_y[yy];
	SLTN_CRV * Slcv = &this->Sltn_crv[yy*WpH];
	unsigned char * img_st = &img[yy*m_width];
	unsigned char * im_max1 =  new unsigned char [m_width*4]; memset(im_max1,0,m_width*4);
	short int  * i_d_max1 =  new short int  [m_width];  for(int i =0; i<m_width; i++){i_d_max1[i] = -1; }


//----------------------------------------------------
	if(dir==0||dir==1){
	for(int i =0; i<crv_lth; i++)
	{
		unsigned char lb_i = (Slcv[i].lb)%8;
		int tau_i = Slcv[i].tau;
		int xl = XL_tau[tau_i];
		int xr = XR_tau[tau_i];
		int d = xl-xr;
		//////////// status LFT ///////////////
		if(dir==0)
		{
			int x_im =xl;
			if (i_d_max1[x_im] <=d)
			{
				i_d_max1[x_im]=d;
				for( int c = 0; c<4; c ++)im_max1[x_im +c*m_width] = Slcv[i].rgb[c];
			}
		 }
      //////////// status RGT /////////////
		if(dir==1)
		{
			int x_im =xr;
			if (i_d_max1[x_im] <= d)
			{
				i_d_max1[x_im]=d;
				for( int c = 0; c<4; c ++)im_max1[x_im +c*m_width] = Slcv[i].rgb[c];
			}
		 }

        /////////////////////////////////////////////////////
	}}
	
	//---------------
	else {
	for(int i =0; i<crv_lth-1; i++)
	{
		unsigned char lb_i1 = (Slcv[i].lb)%8;
		unsigned char lb_i2 = (Slcv[i+1].lb)%8;
		unsigned char *rgb1   = &Slcv[i].rgb[0];
		unsigned char *rgb2  = &Slcv[i+1].rgb[0];
		int tau_i1 = Slcv[i].tau;
		int tau_i2 = Slcv[i+1].tau;
		int xl1 = XL_tau[tau_i1];
		int xr1 = XR_tau[tau_i1];
		int xl2 = XL_tau[tau_i2];
		int xr2 = XR_tau[tau_i2];
		int d_c1 = Slcv[i].d;
		int d_c2 = Slcv[i+1].d;
		int d1 =   xl1-xr1;
		int d2 =   xl1-xr2;
		 //int dm = (xl1-xr1> xl2-xr2) ? xl1-xr1: xl2-xr2;
//-------------------------------------------------------------
			float  x_imf1= (float)xl1- dir*d1;
			float  x_imf2= (float)xl2- dir*d2;

		  if(dir>=0&&dir<=1)
		  {
			  if( lb_i1 == MY_LFT ||lb_i1 == MY_LFT_A ){d1 = d_c1;  x_imf1= (float)xl1- dir*d1; }
			  if( lb_i2 == MY_LFT ||lb_i2 == MY_LFT_A ) {d2 = d_c2; x_imf2= (float)xl2- dir*d2; }
			  if( lb_i1 == MY_RGT ||lb_i1 == MY_RGT_A ) {x_imf1= (float)xl1 + d_c1 - d1- dir*d_c1; d1 = d_c1;}
			  if( lb_i2 == MY_RGT ||lb_i2 == MY_RGT_A ) { x_imf2= (float)xl2 + d_c2 - d2 - dir*d_c2; d2 = d_c2; }
		  }

			int st_p = (x_imf1<0)? (int)x_imf1-1 : (int)x_imf1;
			int fn_p = (x_imf2<0)? (int)x_imf2-1 : (int)x_imf2;
			int inc = (fn_p - st_p >=0) ? 1:-1;
            int dsp;
			for(int ii =st_p; ii<=fn_p; ii+=inc)
			{
				int xx = ii;
				if (xx>=0&&xx<m_width){
					
			     if(ii == st_p)
				 {
                  if( lb_i1 == MY_VIS){ i_d_max1[xx]=m_nLabels; for( int c = 0; c<4; c ++)im_max1[xx +c*m_width] = rgb1[c]; }
				  else	{if( i_d_max1[xx] < (dsp = (d1 < d2) ? d1 : d2)){  i_d_max1[xx]=dsp; for( int c = 0; c<4; c ++)im_max1[xx +c*m_width] = rgb1[c];} }
				 }
				  if(ii == fn_p)
				 {
					 if( lb_i2 == MY_VIS){if( i_d_max1[xx] < m_nLabels){ i_d_max1[xx]=m_nLabels; for( int c = 0; c<4; c ++)im_max1[xx +c*m_width] = rgb2[c]; }}
				  else	{if( i_d_max1[xx] < (dsp = (d1 < d2) ? d1 : d2)){  i_d_max1[xx]=dsp; for( int c = 0; c<4; c ++)im_max1[xx +c*m_width] = rgb2[c];} }
				 }
				 if( ii != fn_p && ii != st_p )
				 {
                  if( lb_i1 == MY_VIS && lb_i2 == MY_VIS){ i_d_max1[xx]=m_nLabels; for( int c = 0; c<4; c ++)im_max1[xx +c*m_width] = (rgb2[c]+  rgb2[c])/2 ; }
				  else	{if( i_d_max1[xx] < (dsp = (d1 < d2) ? d1 : d2)){  i_d_max1[xx]=dsp; for( int c = 0; c<4; c ++)im_max1[xx +c*m_width] = (rgb2[c]+  rgb2[c])/2 ;} }
				 }
				}///////////////// xx>=0&&xx<m_width
			} //////intr	
	}}

	for(int i =0; i<m_width; i++) for( int c = 0; c<4; c ++) img_st[i + c*m_nPixels] = im_max1[i +c*m_width];

		/////////////////////////////
		delete [] im_max1;
        delete [] i_d_max1;


}
void EDP::ProjectImCrvB(int yy, unsigned char *img, float dir) {
	int crv_lth = this->T_max_y[yy];
	SLTN_CRV * Slcv = &this->Sltn_crv[yy*WpH];
	unsigned char * img_st = &img[yy*m_width];
	unsigned char * im_max1 =  new unsigned char [m_width*4]; memset(im_max1,0,m_width*4);
	short int  * i_d_max1 =  new short int  [m_width];  for(int i =0; i<m_width; i++){i_d_max1[i] = -1; }

//----------------------------------------------------
	if(dir==0||dir==1){
	for(int i =0; i<crv_lth; i++)
	{
		unsigned char lb_i = (Slcv[i].lb)%8;
		int tau_i = Slcv[i].tau;
		int xl = XL_tau[tau_i];
		int xr = XR_tau[tau_i];
		int d = xl-xr;
		//////////// status LFT ///////////////
		if(dir==0)
		{
			int x_im =xl;
			if (i_d_max1[x_im] <=d)
			{
				i_d_max1[x_im]=d;
				for( int c = 0; c<4; c ++)im_max1[x_im +c*m_width] = Slcv[i].rgb[c];
			}
		 }
      //////////// status RGT /////////////
		if(dir==1)
		{
			int x_im =xr;
			if (i_d_max1[x_im] <= d)
			{
				i_d_max1[x_im]=d;
				for( int c = 0; c<4; c ++)im_max1[x_im +c*m_width] = Slcv[i].rgb[c];
			}
		 }

        /////////////////////////////////////////////////////
	}}
	
	//---------------
	else {
	for(int i =0; i<crv_lth-1; i++)
	{
		//unsigned char lb_i1 = (Slcv[i].lb)%8;
		//unsigned char lb_i2 = (Slcv[i+1].lb)%8;
		unsigned char *rgb1, *rgbt1 = &Slcv[i].rgb[0];
		unsigned char *rgb2 ,*rgbt2 = &Slcv[i+1].rgb[0];
		int tau_i1 = Slcv[i].tau;
		int tau_i2 = Slcv[i+1].tau;
		int xl1 = XL_tau[tau_i1];
		int xr1 = XR_tau[tau_i1];
		int xl2 = XL_tau[tau_i2];
		int xr2 = XR_tau[tau_i2];
		int d1 =  xl1-xr1;
		int d2 =  xl1-xr2;
		 int dm = (xl1-xr1> xl2-xr2) ? xl1-xr1: xl2-xr2;
//-------------------------------------------------------------
		if(i==0)for( int c = 0; c<4; c ++){im_max1[c*m_width + m_width-1] = rgbt1[c];}
		if(i==crv_lth-2)for( int c = 0; c<4; c ++){im_max1[c*m_width ] = rgbt2[c];}
			//
			float  x_imf1, xf1 = (float)xl1- dir*d1;
			float  x_imf2, xf2 = (float)xl2- dir*d2;
			if(xf1<xf2){x_imf1 = xf1; x_imf2 = xf2; rgb1 =rgbt1; rgb2 = rgbt2;}
			else            {x_imf1 = xf2; x_imf2 = xf1; rgb1 =rgbt2; rgb2 = rgbt1;}
			int xim1 = (int)x_imf1;
			int xim2 = (int)x_imf2;
			int n = xim2 - xim1;
				
			for(int ii =0; ii<n; ii++)
			{
				int xx = xim1 +1+ ii;
				if (xx>=0&&xx<m_width){
				if (i_d_max1[xx] <=dm){ i_d_max1[xx]=dm;
				float w1 = fabs(x_imf2 - xx); float  w2 = fabs(x_imf1 - xx);
				for( int c = 0; c<4; c ++){float vl = (w1+w2 >0)? (w1*rgb1[c] + w2*rgb2[c])/(w1+w2) : 0; im_max1[xx +c*m_width] = (vl>255) ? 255: round_fl(vl); }}
				}
			}	
	}}

	for(int i =0; i<m_width; i++) for( int c = 0; c<4; c ++) img_st[i + c*m_nPixels] = im_max1[i +c*m_width];

		/////////////////////////////
		delete [] im_max1;
        delete [] i_d_max1;


}
void EDP::ProjectDCrv(int yy, int dir, int * out)
{
	int crv_lth = this->T_max_y[yy];
	SLTN_CRV * Slcv = &this->Sltn_crv[yy*WpH];
	Label * img_st = &out[yy*m_width];
	unsigned char * maska = &Mask_RL[dir][yy*m_width];
	short int  * im_max1 =  new short int  [m_width];
	short int  * i_d_max1 =  new short int  [m_width];   for(int i =0; i<m_width; i++){i_d_max1[i] = -1;    }
//----------------------------------------------------
	if(dir==0||dir==1){
	for(int i =0; i<crv_lth; i++)
	{
		unsigned char lb_i = (Slcv[i].lb);
		int tau_i = Slcv[i].tau;
		int xl = XL_tau[tau_i];
		int xr = XR_tau[tau_i];
		int d = xl-xr;
		//////////// status LFT ///////////////
		if(dir==0)
		{
			int x_im =xl;
			if (i_d_max1[x_im] <=d)
			{
				i_d_max1[x_im]=d;
				im_max1[x_im] = Slcv[i].d;
				maska[x_im] = lb_i;
			}
		 }
      //////////// status RGT /////////////
		if(dir==1)
		{
			int x_im =xr;
			if (i_d_max1[x_im] <= d)
			{
				i_d_max1[x_im]=d;
				im_max1[x_im] = Slcv[i].d;
				maska[x_im] = lb_i;
			}
		 }

        /////////////////////////////////////////////////////
	}}
	
	//---------------
	else {
	for(int i =0; i<crv_lth-1; i++)
	{
		unsigned char lb_i1 = (Slcv[i].lb)/8;
		unsigned char lb_i2 = (Slcv[i+1].lb)/8;
		int rgb1, rgbt1 = Slcv[i].d;
		int rgb2 ,rgbt2 = Slcv[i+1].d;
		int tau_i1 = Slcv[i].tau;
		int tau_i2 = Slcv[i+1].tau;
		int xl1 = XL_tau[tau_i1];
		int xr1 = XR_tau[tau_i1];
		int xl2 = XL_tau[tau_i2];
		int xr2 = XR_tau[tau_i2];
		int d1 =  xl1-xr1;
		int d2 =  xl1-xr2;
		float df1, df2;
		 int dm = (xl1-xr1> xl2-xr2) ? xl1-xr1: xl2-xr2;
//-------------------------------------------------------------
		if(i==0)for( int c = 0; c<3; c ++){im_max1[0] = rgbt1;}
		if(i==crv_lth-2)for( int c = 0; c<3; c ++){im_max1[m_width-1] = rgbt2;}
			//
			float  x_imf1, xf1 = (float)xl1- dir*d1;
			float  x_imf2, xf2 = (float)xl2- dir*d2;
			if(xf1<xf2){x_imf1 = xf1; x_imf2 = xf2; rgb1 =rgbt1; rgb2 = rgbt2; df1=d1; df2 = d2;}
			else            {x_imf1 = xf2; x_imf2 = xf1; rgb1 =rgbt2; rgb2 = rgbt1;df1=d2; df2 = d1;}
			int xim1 = (int)x_imf1;
			int xim2 = (int)x_imf2;
			int n = xim2 - xim1;
				
			for(int ii =0; ii<n; ii++)
			{
				int xx = xim1 +1+ ii;
				if (xx>=0&&xx<m_width){
				if (i_d_max1[xx] <=dm){ i_d_max1[xx]=dm;
				float w1 = fabs(x_imf2 - xx); float  w2 = fabs(x_imf1 - xx);
				{float vl = (w1+w2 >0)? (w1*rgb1 + w2*rgb2)/(w1+w2) : 0; im_max1[xx] =  round_fl(vl); maska[xx] = lb_i1*255;}}
				}
			}	
	}}

	for(int i =0; i<m_width; i++)  img_st[i] = im_max1[i];

		/////////////////////////////
		delete [] im_max1;
	    delete [] i_d_max1;




}
EDP::REAL EDP::computeDSI_tau(unsigned char *im1,   unsigned char *im2,
		int birchfield,		// use Birchfield/Tomasi costs
		int squaredDiffs,	// use squared differences
		int truncDiffs		// truncated differences
	) { // function starts
	int nColors =3;
    int worst_match = nColors * (squaredDiffs ? 255 * 255 : 255);
    // truncation threshold - NOTE: if squared, don't multiply by nColors (Eucl. dist.)
    int maxsumdiff = squaredDiffs ? truncDiffs * truncDiffs : nColors * abs(truncDiffs);
    // value for out-of-bounds matches
    int badcost = std::min(worst_match, maxsumdiff);
    REAL sum_ret =0;
    int dsiIndex=0;
    for (int y = 0; y < m_height; y++) {
		for (int t = 0; t <Tau; t++)  {//tau
			int yw0 = y*m_width;
			int x =    XL_tau[t];
			int x2 =  XR_tau[t];
			int xp =(x<m_width-1)? x+1:x;
	        int x2p =(x2<m_width-1)? x2+1:x2;
			int d = x-x2;
			int dsiValue;
			int sumdiff = 0;
			int sumdfb[4] = {0,0,0,0};
			for (int b = 0; b < nColors; b++) { // color
				int yw= yw0+b*this->m_nPixels;
				int diff;  int dfb[4];
				if (birchfield)  {
					dfb[0]= abs(im1[yw + x] +im1[yw + xp] - im2[yw+x2]- im2[yw+x2p])/2;
					dfb[1]= abs(im1[yw + x]  - im2[yw+x2]);
					dfb[2]= abs((im1[yw + x] +im1[yw + xp])/2 - im2[yw+x2]);
					dfb[3]= abs(im1[yw + x] - (im2[yw+x2]- im2[yw+x2p])/2);
					for(int q =0;q<4;q++) {
						sumdfb[q] += (squaredDiffs )? dfb[q]*dfb[q]:dfb[q];
					}
				} else {
					diff = abs(im1[yw + x] - im2[yw+x2]);
					sumdiff += (squaredDiffs ? diff * diff : diff);
				}
			}//color

			// truncate diffs
			if (birchfield) {
				sumdiff =std::min(std::min(sumdfb[2], sumdfb[3]) , std::min(sumdfb[0], sumdfb[1]));
			}
			dsiValue = std::min(sumdiff, maxsumdiff);
			/*if(t==0||t==Tau-1) dsiValue=0;*/
			////
			sum_ret += m_D[dsiIndex++] = dsiValue;
		
		} ///tau
    }
	return sum_ret/(m_height*Tau);
}
EDP::REAL EDP::computeDSI_tauG(unsigned char *im1,   unsigned char *im2,
		int birchfield,       // use Birchfield/Tomasi costs
		int squaredDiffs,     // use squared differences
		int truncDiffs,
		float alph)       // truncated differences
{
    int nColors =6;
    int worst_match = nColors/2 * (squaredDiffs ? 255 * 255 : 255);
    // truncation threshold - NOTE: if squared, don't multiply by nColors (Eucl. dist.)
    int maxsumdiff = squaredDiffs ? truncDiffs * truncDiffs : nColors/2* abs(truncDiffs);
    // value for out-of-bounds matches
    int badcost = std::min(worst_match, maxsumdiff);
    REAL sum_ret =0;
    int dsiIndex = 0;
    for (int y = 0; y < m_height; y++) {
		for (int t = 0; t <Tau; t++) {//tau
			int yw0 = y*m_width;
			int x = XL_tau[t];
			int x2 = XR_tau[t];
			int xp =(x<m_width-1)? x+1:x;
	        int x2p =(x2<m_width-1)? x2+1:x2;
			int d = x-x2;
			int dsiValue;
			int sumdiff = 0;
			int sumdfb[4] = {0,0,0,0};
			for (int b = 0; b < nColors; b++) {
				int yw= yw0+b*this->m_nPixels;
				int diff;  int dfb[4];
				if (birchfield) {
					dfb[0]= abs(im1[yw + x] +im1[yw + xp] - im2[yw+x2]- im2[yw+x2p])/2;
					dfb[1]= abs(im1[yw + x]  - im2[yw+x2]);
					dfb[2]= abs((im1[yw + x] +im1[yw + xp])/2 - im2[yw+x2]);
					dfb[3]= abs(im1[yw + x] - (im2[yw+x2]- im2[yw+x2p])/2);
					for(int q =0;q<4;q++)
						sumdfb[q] += (squaredDiffs )? dfb[q]*dfb[q]:dfb[q];
				} else  {
					diff = abs(im1[yw + x] - im2[yw+x2]);
					sumdiff += (squaredDiffs ? diff * diff : diff);
				}
			}//color

			// truncate diffs
			if (birchfield) {
				sumdiff =std::min(std::min(sumdfb[2], sumdfb[3]) , std::min(sumdfb[0], sumdfb[1]));
			}
			sumdiff  *= alph;
			dsiValue = std::min(sumdiff, maxsumdiff);
			/*if(t==0||t==Tau-1) dsiValue=0;*/
			////
			sum_ret += m_D[dsiIndex] += ( m_D[dsiIndex]+ dsiValue>maxsumdiff) ? maxsumdiff-m_D[dsiIndex]: dsiValue;
			dsiIndex++;
		
		} ///tau
    }
	return sum_ret/(m_height*Tau);
}

void EDP::filterDSI_tau(unsigned char *im,    int r_w,     float cl_prc) {
	float * c_D = new float [Tau*m_height];
/*	float * c_D_pL = new float [Tau*m_height];
	float * c_D_pR = new float [Tau*m_height];*/ 
	for (int t = 0; t <Tau*m_height; t++) {/*c_D_pR[t] = m_D_pR[t]; c_D_pL[t] = m_D_pL[t]; */c_D[t] = m_D[t]; }
	//----------------------------
	unsigned char rgb[3];
	float rgb_m[3] ={0,0,0};
	for(int i= 0; i< m_nPixels;i ++)
	for(int c= 0; c< 3;c ++)
	{	rgb_m[c] +=im[i+m_nPixels*c];}
	for(int c= 0; c< 3;c ++)
	{	rgb_m[c] /= m_nPixels;}
	float vll, std =0;
	for(int i= 0; i< m_nPixels;i ++)
	{ vll =0;
	for(int c= 0; c< 3;c ++)
	{ float vl = (im[i+m_nPixels*c]-rgb_m[c]); vl *=vl; vll +=vl;}
	std += vll;
	}
	int cl_fc = round_fl(sqrt(std/m_nPixels)*cl_prc);
	///////////////////////////////////////////////
	float *gss_r_wt = new float [m_height];
	float *gss_c_wt = new float [1000];
	float sgm_r = r_w*0.5;
	int r_x =0;
	int dm = (2*r_w+1)*(r_x*2+1);
	float * cst = new float [dm];
	//float * cstL = new float [dm];
	//float * cstR = new float [dm];
	float * wgt = new float [dm];


	for(int i=0; i<m_height; i++ )gss_r_wt[i]= exp(-(i*i)/(2.*sgm_r*sgm_r));
	for(int i=0; i<1000; i++ )gss_c_wt[i]= exp(-(i*i)/(2.*cl_fc*cl_fc));

    ///////////////////////////////////////////////////////


    for (int y = 0; y < m_height; y++) {
	for (int t = 0; t <Tau; t++)
	{//tau
		
		int x =    XL_tau[t];
		int d = x - XR_tau[t];
		int ii=0; int st =0;
		
		unsigned char vl0[3] = {im[y*m_width + x], im[y*m_width + x + m_nPixels],  im[y*m_width + x + m_nPixels*2]};
		for(int i =-r_w; i<=r_w; i++)
		for(int j =-r_x; j<=r_x; j++)
		{
			int y_i = y+i; int x_j = x+j;
			
			if(y_i>=0&&y_i<m_height &&x_j>=0&&x_j<m_width)
			{
			int tt = this->tau_XL_d[x_j*m_nLabels+d];
			if(tt!=-1){
			cst[ii] = m_D[y_i*Tau+ tt]; /*cstL[ii] = m_D_pL[y_i*Tau+ tt]; cstR[ii] = m_D_pR[y_i*Tau+ tt]; *//*if(cst[ii] < Cost_mean/2) st = 1;  */
			ii++;
			}
			}
		}
		if(1){
		ii=0;
		for(int i =-r_w; i<=r_w; i++)
		for(int j =-r_x; j<=r_x; j++)
		{
			int y_i = y+i; int x_j = x+j; int ind_t = x_j+y_i*m_width;
			int rd = round_fl(sqrt((float)i*i+j*j));
			if(y_i>=0&&y_i<m_height &&x_j>=0&&x_j<m_width)
			{
			int tt = this->tau_XL_d[x_j*m_nLabels+d];
			if(tt!=-1){
			int c_t = (abs(vl0[0]-im[ind_t]) + abs(vl0[1]-im[ind_t+ m_nPixels]) + abs(vl0[2]-im[ind_t+ m_nPixels*2]));
			if(c_t>999)c_t = 999;
			wgt[ii] = gss_c_wt[c_t]*gss_r_wt[rd];
			ii++;}
			}
		}
		float cst_p = 0;
		//float cst_pL = 0;
		//float cst_pR = 0;
		float wt_p =0;
		for(int i =0; i<ii; i++){cst_p +=cst[i]*wgt[i];/* cst_pL +=cstL[i]*wgt[i]; cst_pR +=cstR[i]*wgt[i];*/  wt_p += wgt[i];}
		c_D[y*Tau + t] = cst_p/wt_p ;/* c_D_pL[y*Tau + t] = cst_pL/wt_p ; c_D_pR[y*Tau + t] = cst_pR/wt_p ; *///(c_D[y*Tau + t] !=Cost_mean*2)? cst_p/wt_p : Cost_mean*2;
		}
		else c_D[y*Tau + t] = Cost_mean;
	//////////////////////////////////////

	        } ///tau

    }
	for (int t = 0; t <Tau*m_height; t++) {m_D[t] = c_D[t]; /*m_D_pL[t] = c_D_pL[t]; m_D_pR[t] = c_D_pR[t];*/}
delete[] c_D;
//delete[] c_D_pL;
//delete[] c_D_pR;
delete [] gss_r_wt;
delete [] gss_c_wt;
delete [] cst;
//delete [] cstL;
//delete [] cstR;
delete [] wgt;
}

double EDP::Cost_mean_vl() {
	
    double * sums = new double [m_nLabels]; memset(sums, 0, sizeof(double)*m_nLabels);
	int * cnt = new int [m_nLabels]; memset(cnt, 0, sizeof(int)*m_nLabels);
    for (int y = 0; y < m_height; y++) {
	for (int t = 0; t <Tau; t++) {//tau
			int xl = XL_tau[t];
			int xr = XR_tau[t];
			int d = xl-xr;
			sums[d] += m_D[y*Tau + t];
			cnt[d] ++;
		} ///tau
	}
	double min = sums[0]/cnt[0];
	double max = sums[0]/cnt[0];
	for(int i =1;  i<m_nLabels; i++) {
		min = (min< sums[i]/cnt[i]) ? min: sums[i]/cnt[i]; 
		max = (max> sums[i]/cnt[i]) ? max: sums[i]/cnt[i];
	}
	delete [] sums;
	delete [] cnt;
	return min;
}

EDP::REAL EDP::TRWS_L(int itr, EDP:: REAL * Mrg, int *DspM) {
	REAL max =0, min=0;
	
	REAL * Msg = new REAL [m_height*4*Tau];
	memset(Msg, 0, m_height*Tau*4*sizeof(REAL));
	
	REAL * Msg_add = new REAL [m_nLabels*m_height*4];
	memset(Msg_add, 0, 4*m_height*m_nLabels*sizeof(REAL));
	
	REAL * Di = new REAL [m_nLabels];
	for(int it =0;it<itr;it++) {
		/////////////////// ALGORITHM //////////////////////////
		for(int p= 0; p<m_height*Tau; p++) {
			Mrg[p]= m_D[p];
		}

		for(int y =0; y<m_height; y++)
		for(int x =0; x<m_width; x ++) {
			MkMrgL(x, y, Msg, Msg_add, Mrg);
			UpdMsgXL_pp(x, y, Di,Msg, Msg_add, Mrg);
			UpdMsgYL_pp(x, y, Di,Msg, Msg_add, Mrg);
		}// end dir 0 ////////////////

		for(int p= 0; p<m_height*Tau; p++) {
			Mrg[p]= m_D[p];
		}
		for(int y =m_height-1; y >=0; y--)
		for(int x =m_width-1; x >=0;  x--) {
			MkMrgL(x, y, Msg, Msg_add, Mrg);
			UpdMsgXL_mm(x, y, Di,Msg, Msg_add, Mrg);
			UpdMsgYL_mm(x, y, Di,Msg, Msg_add, Mrg);
		}// end dir ////////////////
		for(int p= 0; p<m_height*Tau; p++) {
			Mrg[p]= m_D[p];
		}
		for(int y =0; y<m_height; y++)
		for(int x =m_width-1; x >=0;  x--) {
			MkMrgL(x, y, Msg, Msg_add, Mrg);
			UpdMsgXL_mm(x, y, Di,Msg, Msg_add, Mrg);
			UpdMsgYL_pp(x, y, Di,Msg, Msg_add, Mrg);
		}// end dir ////////////////
		for(int p= 0; p<m_height*Tau; p++) {
			Mrg[p]= m_D[p];
		}
		for(int y =m_height-1; y >=0; y--)
		for(int x =0; x<m_width; x ++) {
			MkMrgL(x, y, Msg, Msg_add, Mrg);
			UpdMsgXL_pp(x, y, Di,Msg, Msg_add, Mrg);
			UpdMsgYL_mm(x, y, Di,Msg, Msg_add, Mrg);
		}// end dir ////////////////
	}/////////end itr ///////////
	//////////////////////////////////
	//////////////////////////////////
	for(int p= 0; p<m_height*Tau; p++){
		Mrg[p]= m_D[p];
	}

	for(int y =0; y<m_height; y++)
	for(int x =0; x<m_width; x ++) {
		int d_mrg;
		int d_mx = (x + 1 <m_nLabels) ? x +1 : m_nLabels;
		MkMrgL(x, y, Msg, Msg_add, Mrg);
		REAL *mrg = &Mrg(y);
		int *xld = &XLD(x);
		DspM[y*m_width+x]=0;
		REAL mrg_min = mrg[xld[0]];
		for(int d = 1; d< d_mx; d++) {
			if(mrg_min > mrg[xld[d]]) {
				d_mrg=DspM[y*m_width+x]=d;
				mrg_min = mrg[xld[d]];
			}
		}
		for(int d = 0; d< d_mx; d++) {
			mrg[xld[d]] -=mrg_min;
		}

	}// end dir ////////////////
	//for(int p= 0; p<m_height*Tau; p++)Mrg[p] -= min;
	//////////////////////
	delete [] Msg;
	delete [] Di;
	delete [] Msg_add;
	return max-min;
}
void EDP::transform_D(float thr,  float pw)
{
	float avmin =0;
	for(int y =0; y<m_height; y++)
	for(int x =0; x<m_width; x ++) {
		int d_mrg;
		int d_mx = (x + 1 <m_nLabels) ? x +1 : m_nLabels;
		float *mrg = &m_D[y*Tau];
		int *xld = &XLD(x);
		float mrg_min = mrg[xld[0]];
		for(int d = 1; d< d_mx; d++)
			if(mrg_min > mrg[xld[d]]) {
				mrg_min = mrg[xld[d]];
			}
		avmin += mrg_min;

	}// end dir ////////////////
	avmin /= m_nPixels;  thr = thr*avmin;

	for(int y =0; y<m_height; y++)
	for(int x =0; x<m_width; x ++) {
		int d_mrg;
		int d_mx = (x + 1 <m_nLabels) ? x +1 : m_nLabels;
		float *mrg = &mm_D(y);
		int *xld = &XLD(x);
		for(int d = 0; d< d_mx; d++){
			float ex = mrg[xld[d]]/thr;
			ex = pow(ex, pw);
			mrg[xld[d]] = 1. - exp(-ex);
		}
	}// end dir ////////////////
}

void EDP::view_D(int st, double * cst) {
	double min, sum = 0;
	float *vw = new float [N_PX];
	for(int x = 0; x < m_width; x++)
	{
		int i0 = x + st*m_width;
		int d =0; sum +=  vw[x] = min = cst[i0];
		for(int l =1; l < N_LB; l++){
			if(cst[i0 + l*N_PX] <0 || cst[i0 + l*N_PX] >1) {
				int a =0;
			}
			sum += vw[x + l*m_width] = cst[i0 + l*N_PX]; if(min > cst[i0 + l*N_PX]){d =l; min = cst[i0 + l*N_PX];}
		}
		vw[x + d*m_width] = -1;
	}
	sum /= m_width*N_LB;
	for(int x = 0; x < m_width; x++) {
		for(int l =0; l < N_LB; l++) {
			float vl =  vw[x + l*m_width]*200/sum;
		if(vl < 0) 
			for(int c = 0; c < 3; c++)
				I_ims[2][x + l*m_width + c*N_PX] = (c==2) ? 255: 64;
		else 
			for(int c = 0; c < 3; c++)
				I_ims[2][x + l*m_width + c*N_PX] = (vl >255) ? 255: (int)vl;
	}


	}
	delete [] vw;
}

void EDP::view_D(int st, float * cst) {
	double min, sum = 0;
	float *vw = new float [N_PX];
	for(int x = 0; x < m_width; x++) {
		int i0 = x + st*m_width;
		int d =0; sum +=  vw[x] = min = cst[i0];
		for(int l =1; l < N_LB; l++){
			sum += vw[x + l*m_width] = cst[i0 + l*N_PX];
			if(min > cst[i0 + l*N_PX]){
				d =l;
				min = cst[i0 + l*N_PX];
			}
		}
		vw[x + d*m_width] = -1;
	}
	sum /= m_width*N_LB;

	for(int x = 0; x < m_width; x++) {
		for(int l =0; l < N_LB; l++){
			float vl =  vw[x + l*m_width]*200/sum;
			if(vl < 0) {
				for(int c = 0; c < 3; c++) {
					I_ims[2][x + l*m_width + c*N_PX] = (c==2) ? 255: 64;
				}
			} else {
				for(int c = 0; c < 3; c++) {
					I_ims[2][x + l*m_width + c*N_PX] = (vl >255) ? 255: (int)vl;
				}
			}
		}
	}
	delete [] vw;
}



void EDP::transform_D_(float thr) {
	float mm  = min_mean();
	float avmin =0;
	int cnt =0;
	for(int y =0; y<m_height; y++)
	for(int x =0; x<m_width; x ++) {
		int d_mrg;
		int d_mx = (x + 1 <m_nLabels) ? x +1 : m_nLabels;
		float *mrg = &m_D[y*Tau];
		int *xld = &XLD(x);
		
		float mrg_min = mrg[xld[0]];
		for(int d = 1; d< d_mx; d++) {
			if(mrg_min > mrg[xld[d]]) {
				mrg_min = mrg[xld[d]];
			}
		}

		for(int d = 0; d< d_mx; d++){
			float dl =(mrg[xld[d]]  - mrg_min);
			avmin += dl;
			cnt ++;
		}
	}// end dir ////////////////

	avmin /= cnt;  //thr = thr*avmin;
	
	thr = thr*(avmin - mm)*(avmin - mm);
	
	if(thr> 2.5) // what even
		thr *=0.91;
	
	for(int y =0; y<m_height; y++)
	for(int x =0; x<m_width; x ++){
		int d_mrg;
		int d_mx = (x + 1 <m_nLabels) ? x +1 : m_nLabels;
		float *mrg = &mm_D(y);
		int *xld = &XLD(x);
		for(int d = 0; d< d_mx; d++){
			float arg =  ((mrg[xld[d]] - mm)/mm);
			mrg[xld[d]] = (1. + erf(arg*thr))/2;
			//float ex = mrg[xld[d]]/thr; ex = pow(ex, pw);  mrg[xld[d]] = 1. - exp(-ex);
		}
	}// end dir ////////////////
}

float EDP::min_mean() {
	float avmin =0;
	for(int y =0; y<m_height; y++)
	for(int x =0; x<m_width; x ++) {
		int d_mrg;
		int d_mx = (x + 1 <m_nLabels) ? x +1 : m_nLabels;
		float *mrg = &m_D[y*Tau];
		int *xld = &XLD(x);
		
		float mrg_min = mrg[xld[0]];
		
		for(int d = 1; d< d_mx; d++) {
			if(mrg_min > mrg[xld[d]]) {
				mrg_min = mrg[xld[d]];
			}
		}
		avmin += mrg_min;

	}// end dir ////////////////
	avmin /= m_nPixels;
	return avmin;
}

void EDP::TRWS_CST(int itr, EDP:: REAL * Mrg, double *m_cst, int *DspM, int lr) {
	REAL * Msg = new REAL [m_nPixels*m_nLabels*4];
	memset(Msg, 0, m_nPixels*m_nLabels*4*sizeof(REAL));
	
	REAL * Di = new REAL [m_nLabels];

	for(int it =0;it<itr;it++) {
		/////////////////// ALGORITHM //////////////////////////
		for(int p= 0; p<m_nPixels; p++) {
			for(int l= 0; l<m_nLabels; l++) {
				Mrg[p*m_nLabels+l]= m_cst[p+l*m_nPixels];
			}
		}

		for(int y =0; y<m_height; y++) {
			for(int x =0; x<m_width; x ++) {
				MkMrgC(x, y, Msg, Mrg);
				UpdMsgC(x, y, Di,Msg, Mrg, 0+lr*4);
				UpdMsgC(x, y, Di,Msg, Mrg, 1+lr*4);
			}// end dir 0 ////////////////
		}
		
		for(int p= 0; p<m_nPixels; p++) {
			for(int l= 0; l<m_nLabels; l++) {
				Mrg[p*m_nLabels+l]= m_cst[p+l*m_nPixels];
			}
		}

		for(int y =m_height-1; y >=0; y--) {
			for(int x =m_width-1; x >=0;  x--) {
				MkMrgC(x, y, Msg, Mrg);
				UpdMsgC(x, y, Di,Msg, Mrg, 2+lr*4);
				UpdMsgC(x, y, Di,Msg, Mrg, 3+lr*4);
			}// end dir ////////////////
		}

		for(int p= 0; p<m_nPixels; p++) {
			for(int l= 0; l<m_nLabels; l++) {
				Mrg[p*m_nLabels+l]= m_cst[p+l*m_nPixels];
			}
		}

		for(int y =0; y<m_height; y++)
			for(int x =m_width-1; x >=0;  x--) {
				MkMrgC(x, y, Msg, Mrg);
				UpdMsgC(x, y, Di,Msg, Mrg, 2+lr*4);
				UpdMsgC(x, y, Di,Msg, Mrg, 1+lr*4);
			}// end dir ////////////////

		for(int p= 0; p<m_nPixels; p++)
		for(int l= 0; l<m_nLabels; l++)
			Mrg[p*m_nLabels+l]= m_cst[p+l*m_nPixels];

		for(int y =m_height-1; y >=0; y--)
			for(int x =0; x<m_width; x ++) {
				MkMrgC(x, y, Msg, Mrg);
				UpdMsgC(x, y, Di,Msg, Mrg, 0+lr*4);
				UpdMsgC(x, y, Di,Msg, Mrg, 3+lr*4);
			}// end dir ////////////////
	}
	/////////end itr ///////////
	//////////////////////////////////
	//////////////////////////////////
	for(int p= 0; p<m_nPixels; p++) {
		for(int l= 0; l<m_nLabels; l++) {
			Mrg[p*m_nLabels+l]= m_cst[p+l*m_nPixels];
		}
	}

	for(int y =0; y<m_height; y++) {
		for(int x =0; x<m_width; x ++) {
		 MkMrgC(x, y, Msg, Mrg);
		 DspM[y*m_width+x]=m_nLabels-1; 
		 REAL mrg_min = Mrg[(y*m_width+x)*m_nLabels+ m_nLabels-1];
		 for(int d = m_nLabels-2; d>=0; d--) {
		 	if(mrg_min > Mrg[(y*m_width+x)*m_nLabels+ d]) {
		 		DspM[y*m_width+x]=d;
		 		mrg_min = Mrg[(y*m_width+x)*m_nLabels+ d];
		 	}
		 }
	}
	//for(int d = 0; d< m_nLabels; d++)Mrg[(y*m_width+x)*m_nLabels+ d] -=mrg_min;

	}// end dir ////////////////
	//////////////////////
	delete [] Msg;
	delete [] Di;

}

void EDP::adap_d( double *m_cst, int *DspM, int lr)
{
	//float * e = new float [m_nPixels];
	/*for(int x= 0; x <m_nPixels; x++)e[x]=1;*/
	this->min_cost(m_cst, DspM);
	int nint = 50;
	double hgr[3][3] = {1,1,1,1,1,1,1,1,1};
	double * hcs = new double [nint]; for(int x= 0; x <nint; x++)hcs[x] =1;
    double * cst_lt = new double [nint];
    Gss_wnd_( 4,  0.1,DspM, lr, m_nLabels );
	double grcnt[3] = {0,0,0};
	for(int x= 0; x <m_width; x++)
    for(int y= 0; y <m_height; y++)
	{
	int xp = (x < m_width-1) ?  x+1: x;
	int yp = (y < m_height-1) ? y+1: y;
	int ixp = IND(xp,y); int iyp = IND(x,yp); int i = IND(x,y);
	double mnc = m_cst[i ]; for(int l = 1; l < m_nLabels; l++) if(mnc > m_cst[i + l*m_nPixels]) mnc = m_cst[i + l*m_nPixels];
	 for(int l = 0; l < m_nLabels; l++)m_cst[i + l*m_nPixels] -= mnc;
	int d = DspM[i];
	int csti = (int)(m_cst[i + d*m_nPixels]*(double)nint); if (csti >=nint) csti = nint-1;
	hcs[csti]++;
    int ddx = abs(DspM[i]- DspM[ixp]);
	int ddy = abs(DspM[i]- DspM[iyp]);
	int grnx = (lr) ? B_R_buf_pp_[i].x :B_L_buf_pp_[i].x; grcnt[grnx] ++;
	int grny = (lr) ? B_R_buf_pp_[i].y :B_L_buf_pp_[i].y; grcnt[grny] ++;
	if(grnx == 0)
	{
		if(ddx == 0)hgr[0][0]++;
		if(ddx == 1)hgr[0][1]++;
		if(ddx >   1)hgr[0][2]++;
	}
	if(grny == 0)
	{
		if(ddy == 0)hgr[0][0]++;
		if(ddy == 1)hgr[0][1]++;
		if(ddy >   1)hgr[0][2]++;
	}
	if(grnx == 1)
	{
		if(ddx == 0)hgr[1][0]++;
		if(ddx == 1)hgr[1][1]++;
		if(ddx >   1)hgr[1][2]++;
	}
	if(grny == 1)
	{
		if(ddy == 0)hgr[1][0]++;
		if(ddy == 1)hgr[1][1]++;
		if(ddy >   1)hgr[1][2]++;
	}
	if(grnx == 2)
	{
		if(ddx == 0)hgr[2][0]++;
		if(ddx == 1)hgr[2][1]++;
		if(ddx >   1)hgr[2][2]++;
	}
	if(grny == 2)
	{
		if(ddy == 0)hgr[2][0]++;
		if(ddy == 1)hgr[2][1]++;
		if(ddy >   1)hgr[2][2]++;
	}
	}
	int sum =0; for(int pp = 0; pp< nint; pp++) sum += hcs[pp];
	double lP = log((double)m_nPixels); double min = hcs[0] = lP - log((double)hcs[0]+1) ;
	for(int x= 1; x <nint; x++){hcs[x] = lP - log((double)hcs[x]+1) ; if(min > hcs[x]) min = hcs[x];}
	for(int x= 0; x <nint; x++){ hcs[x] -= min; }
	for(int i =0; i <3; i++)hgr[i][1] += hgr[i][2];
	double mins[3] = { hgr[0][0] = -log(hgr[0][0] / grcnt[0]), hgr[1][0] = -log(hgr[1][0] / grcnt[1]),hgr[2][0] = -log(hgr[2][0] / grcnt[2])     };
	for(int i =0; i <3; i++)for(int j =1; j <3; j++) { hgr[i][j] = -log(hgr[i][j] / grcnt[i]); if(mins[i] > hgr[i][j]) mins[i] = hgr[i][j]; }
	for(int i =0; i <3; i++)for(int j =0; j <3; j++) { hgr[i][j] -= mins[i]; }
	//for(int p = 0; p < m_nPixels; p++)
	//{
	//	for(int d = 0; d < m_nLabels; d++)
	//	{int im = (int)(m_cst[p + d*m_nPixels]*(double)nint); if (im >= nint) im = nint - 1; int ip = (im+1 < nint)? im+1: im;
	//	double w2 = m_cst[p + d*m_nPixels]*(double)nint - im;  w2 = w2 - (int)w2;
 //       m_cst[p + d*m_nPixels] = hcs[im]*(1. - w2) + hcs[ip]*w2;
	//	}
	//}
	for(int q =0; q<3; q++){Smthnss[q].max = hgr[q][2]/lP; Smthnss[q].semi  =  hgr[q][2]/6/lP;}
	//for(int x= 0; x <m_width; x++)
 //   for(int y= 0; y <m_height; y++)
	//{
	//int xp = (x < m_width-1) ? x+1: x;
	//int yp = (y < m_height-1) ? y+1: y;
	//int ixp = IND(xp,y); int iyp = IND(x,yp); int i = IND(x,y);
 //   adapX[i] = abs(DspM[i]- DspM[ixp]);
	//adapY[i] = abs(DspM[i]- DspM[iyp]);
	//}
	//this->GaussCosConv2DFst(5, m_width, m_height, adapX);
	//this->GaussCosConv2DFst(5, m_width, m_height, adapY);
	//this->GaussCosConv2DFst(5, m_width, m_height, e);
	//for(int x= 0; x <m_nPixels; x++){adapX[x] /=e[x]; adapY[x] /=e[x];}
	//adapx =adapy = 0;

 //   for(int x= 0; x <m_nPixels; x++)
	//{
	//	adapy +=I_ims[2][x] =I_ims[2][x+m_nPixels]= adapY[x]*8;
	//	adapx += I_ims[2][x+m_nPixels*2]= adapX[x]*8;
	//}
	//adapx /=m_nPixels;
	//adapy /=m_nPixels;
	//delete [] e;
	delete [] hcs;
	delete [] cst_lt;
}
    void EDP:: comb_cst( double *out_cst, double ** in_cst, int lr)
{
	double * bs = (lr) ? in_cst[1] : in_cst[0];
	double * dp = (lr) ? in_cst[0] : in_cst[1];
	for(int x= 0; x <m_width; x++)
    for(int y= 0; y <m_height; y++)
	{
     int p = IND(x,y);
	 for( int l = 0; l < m_nLabels; l++)
	 {
	   int xp;
	   if(lr==0){xp = (x-l >= 0) ? x -l: -1;}
	   if(lr==1){xp = (x + l < m_width ) ? x +l: -1;}
	   out_cst[IND_IC(x,y,l)] = (xp+1) ? bs[IND_IC(x,y,l)] * sqrt(dp[IND_IC(xp,y,l)]) :  bs[IND_IC(x,y,l)] ;
	 }

	}
}

 void EDP:: TRWS_CST(int itr, EDP:: REAL * Mrg, float *m_cst, int *DspM, int lr)
{

REAL * Msg = new REAL [m_nPixels*m_nLabels*4]; memset(Msg, 0, m_nPixels*m_nLabels*4*sizeof(REAL));
REAL * Di = new REAL [m_nLabels];
for(int it =0;it<itr;it++)
{/////////////////// ALGORITHM //////////////////////////
for(int p= 0; p<m_nPixels; p++)
for(int l= 0; l<m_nLabels; l++)
Mrg[p*m_nLabels+l]= m_cst[p+l*m_nPixels];

for(int y =0; y<m_height; y++)
for(int x =0; x<m_width; x ++)
{
 MkMrgC(x, y, Msg, Mrg);
 UpdMsgC(x, y, Di,Msg, Mrg, 0+lr*4);
 UpdMsgC(x, y, Di,Msg, Mrg, 1+lr*4);
}// end dir 0 ////////////////
for(int p= 0; p<m_nPixels; p++)
for(int l= 0; l<m_nLabels; l++)
Mrg[p*m_nLabels+l]= m_cst[p+l*m_nPixels];

for(int y =m_height-1; y >=0; y--)
for(int x =m_width-1; x >=0;  x--)
{
 MkMrgC(x, y, Msg, Mrg);
 UpdMsgC(x, y, Di,Msg, Mrg, 2+lr*4);
 UpdMsgC(x, y, Di,Msg, Mrg, 3+lr*4);
}// end dir ////////////////
for(int p= 0; p<m_nPixels; p++)
for(int l= 0; l<m_nLabels; l++)
Mrg[p*m_nLabels+l]= m_cst[p+l*m_nPixels];

for(int y =0; y<m_height; y++)
for(int x =m_width-1; x >=0;  x--)
{
MkMrgC(x, y, Msg, Mrg);
 UpdMsgC(x, y, Di,Msg, Mrg, 2+lr*4);
 UpdMsgC(x, y, Di,Msg, Mrg, 1+lr*4);
}// end dir ////////////////
for(int p= 0; p<m_nPixels; p++)
for(int l= 0; l<m_nLabels; l++)
Mrg[p*m_nLabels+l]= m_cst[p+l*m_nPixels];

for(int y =m_height-1; y >=0; y--)
for(int x =0; x<m_width; x ++)
{
MkMrgC(x, y, Msg, Mrg);
 UpdMsgC(x, y, Di,Msg, Mrg, 0+lr*4);
 UpdMsgC(x, y, Di,Msg, Mrg, 3+lr*4);
}// end dir ////////////////
}/////////end itr ///////////
//////////////////////////////////
//////////////////////////////////
for(int p= 0; p<m_nPixels; p++)
for(int l= 0; l<m_nLabels; l++)
Mrg[p*m_nLabels+l]= m_cst[p+l*m_nPixels];

for(int y =0; y<m_height; y++)
for(int x =0; x<m_width; x ++)
{
 MkMrgC(x, y, Msg, Mrg);
 DspM[y*m_width+x]=m_nLabels-1; REAL mrg_min = Mrg[(y*m_width+x)*m_nLabels+ m_nLabels-1];
 for(int d = m_nLabels-2; d>=0; d--)
 if(mrg_min > Mrg[(y*m_width+x)*m_nLabels+ d])
 {DspM[y*m_width+x]=d; mrg_min = Mrg[(y*m_width+x)*m_nLabels+ d];}
//for(int d = 0; d< m_nLabels; d++)Mrg[(y*m_width+x)*m_nLabels+ d] -=mrg_min;

}// end dir ////////////////
//////////////////////
delete [] Msg;
delete [] Di;

}
inline void  EDP::copy_m_D( float * m_cost, int dr){
//Cnst_3= 0;

	if(!dr)
		for(int y =0; y<m_height; y++)
		for(int x =m_width-1; x>=0; x --) {
		int d_mx = (x + 1 <m_nLabels) ? x +1 : m_nLabels;

		int *xld = &XLD(x);
		for(int d = 0; d< m_nLabels; d++)
		if(d< d_mx)/*Cnst_3 += */m_cost[y*m_width+x + d*m_nPixels] = m_D[y*Tau+xld[d]];
		else /*Cnst_3 += */m_cost[y*m_width+x + d*m_nPixels] = m_cost[y*m_width + x + 1 + d - d_mx + d*m_nPixels];
		}
	else 
		for(int y =0; y<m_height; y++)
		for(int x =0 ; x< m_width; x ++)
		{
		int xinv = m_width-1 -x;
		int d_mx = (xinv + 1 <m_nLabels) ? xinv +1 : m_nLabels;
		int *xld = &XRD(x);
		for(int d = 0; d< m_nLabels; d++)
		if(d< d_mx)m_cost[y*m_width+x + d*m_nPixels] = m_D[y*Tau+xld[d]];
		else m_cost[y*m_width+x + d*m_nPixels] = m_cost[y*m_width + x - 1 - d + d_mx + d*m_nPixels];
		}
}
inline void  EDP::copy_m_D( double * m_cost, int dr) {
//Cnst_3= 0;

	if(!dr)
		for(int y =0; y<m_height; y++)
		for(int x =m_width-1; x>=0; x --) {
			int d_mx = (x + 1 <m_nLabels) ? x +1 : m_nLabels;

			int *xld = &XLD(x);
			for(int d = 0; d< m_nLabels; d++)
				if(d< d_mx)/*Cnst_3 += */
					m_cost[y*m_width+x + d*m_nPixels] = m_D[y*Tau+xld[d]];
				else /*Cnst_3 += */
					m_cost[y*m_width+x + d*m_nPixels] = m_cost[y*m_width + x + 1 + d - d_mx + d*m_nPixels];
		}
	else 
		for(int y =0; y<m_height; y++)
		for(int x =0 ; x< m_width; x ++) {
			int xinv = m_width-1 -x;
			int d_mx = (xinv + 1 <m_nLabels) ? xinv +1 : m_nLabels;
			int *xld = &XRD(x);
			for(int d = 0; d< m_nLabels; d++)
				if(d< d_mx)
					m_cost[y*m_width+x + d*m_nPixels] = m_D[y*Tau+xld[d]];
				else 
					m_cost[y*m_width+x + d*m_nPixels] = m_cost[y*m_width + x - 1 - d + d_mx + d*m_nPixels];
		}
}
inline double  EDP::copy_m_D_b( float * m_cost, int dr)
{
	double ret =0;
	if(!dr)
		for(int y =0; y<m_height; y++)
		for(int x =m_width-1; x>=0; x --) {
			int d_mx = (x + 1 <m_nLabels) ? x +1 : m_nLabels;
			int *xld = &XLD(x);
			for(int d = 0; d< m_nLabels; d++)
			if(d< d_mx)
				ret += m_D[y*Tau+xld[d]] = m_cost[y*m_width+x + d*m_nPixels];
		}
	else 
		for(int y =0; y<m_height; y++)
		for(int x =0 ; x< m_width; x ++) {
			int xinv = m_width-1 -x;
			int d_mx = (xinv + 1 <m_nLabels) ? xinv +1 : m_nLabels;
			int *xld = &XRD(x);
			for(int d = 0; d< m_nLabels; d++)
				if(d< d_mx)
					ret += m_D[y*Tau+xld[d]]  = m_cost[y*m_width+x + d*m_nPixels];
		}

	return ret/(m_height*Tau);

}
inline double  EDP::mean_bf( float * m_cost, int pixs)
{
double ret =0; for(int i = 0; i< pixs; i++) ret += m_cost[i]; return (ret/pixs);

}
inline void  EDP::cross_chkL( float * m_cost, int * dL, int * dR)
{

	for(int i =0; i < m_nPixels; i++ )
	{
	int x = i%m_width, y = i/m_width;
	int d = dL[i]; int dr = (x-d>=0)? dR[i-d]:m_nLabels;

	if(d==dr) dL[i] =d;  else dL[i] =m_nLabels;
	}
cross_chkL_cst( m_cost, dL);

}
		  	  	  inline void  EDP:: cross_chkL_cst( float * m_cost, int * dL)
{
	for(int i =0; i < m_nPixels*m_nLabels; i++ ) m_cost[i] = 1;
	for(int i =0; i < m_nPixels; i++ )
	{
	int x = i%m_width, y = i/m_width;
	if(dL[i] != m_nLabels) m_cost[i+dL[i]*m_nPixels] = 0;
	/*else for(int d =0; d < m_nLabels; d++ ) m_cost[i+d*m_nPixels] =0;*/
	}


}
				      inline void  EDP:: min_cost( double * m_cost, int * dsp)
{

for(int y =0; y<m_height; y++)
for(int x =m_width-1; x>=0; x --)
{

 dsp[y*m_width+x]=0; float mrg_min = m_cost[y*m_width+x ];
 for(int d = 1; d< m_nLabels; d++)
 if(mrg_min >m_cost[y*m_width+x + d*m_nPixels])
 {dsp[y*m_width+x]=d; mrg_min = m_cost[y*m_width+x + d*m_nPixels];}

}

}
		
    inline void  EDP:: min_cost( float * m_cost, int * dsp)
{

for(int y =0; y<m_height; y++)
for(int x =m_width-1; x>=0; x --)
{

 dsp[y*m_width+x]=0; float mrg_min = m_cost[y*m_width+x ];
 for(int d = 1; d< m_nLabels; d++)
 if(mrg_min >m_cost[y*m_width+x + d*m_nPixels])
 {dsp[y*m_width+x]=d; mrg_min = m_cost[y*m_width+x + d*m_nPixels];}

}

}
	    inline void  EDP:: min_cost( double* m_cost, int * dsp, int lbs, int sc)
{

for(int y =0; y<m_height; y++)
for(int x =m_width-1; x>=0; x --)
{

 dsp[y*m_width+x]=0; double mrg_min = m_cost[y*m_width+x ];
 for(int d = 1; d< lbs; d++)
 if(mrg_min >m_cost[y*m_width+x + d*m_nPixels])
 {dsp[y*m_width+x]=d; mrg_min = m_cost[y*m_width+x + d*m_nPixels];}

}
FOR_PX_p {dsp[p] *= sc;}

}
		    inline void EDP:: prf_cost( float * m_cost)
{
	float div = 3.5;
int *mrg_min  = new int [m_nPixels*3];
for(int y =0; y<m_height; y++)
for(int x =m_width-1; x>=0; x --)
{

 mrg_min[y*m_width+x ] = m_cost[y*m_width+x ];
 for(int d = 1; d< m_nLabels; d++)
 if(mrg_min[y*m_width+x ] >m_cost[y*m_width+x + d*m_nPixels])
 {mrg_min[y*m_width+x ] = m_cost[y*m_width+x + d*m_nPixels];}

}
for(int i =0; i<m_nPixels; i++){mrg_min[i + m_nPixels] =mrg_min[i]; mrg_min[i + m_nPixels*2] =1;}
SetSort(m_nPixels,&mrg_min[m_nPixels], &mrg_min[m_nPixels*2] );
float thr = mrg_min[mrg_min[(int)(m_nPixels/div) + m_nPixels*2]];
for(int i =0; i<m_nPixels; i++){ float sum=0;
for(int d = 0; d< m_nLabels; d++)sum += m_cost[i + d*m_nPixels] = (m_cost[i + d*m_nPixels]<thr)? 1:0;
//if(sum)for(int d = 0; d< m_nLabels; d++)m_cost[i+ d*m_nPixels] /=sum;
}
delete [] mrg_min;
}
	    inline void EDP:: prf_cost_inv( float * m_cost)
{
//	float div = 3.3;//Cnst_1;
int *mrg_min  = new int [m_nPixels*3];
int *dd = &mrg_min[m_nPixels];
for(int y =0; y<m_height; y++)
for(int x =m_width-1; x>=0; x --)
{

 mrg_min[y*m_width+x ] = m_cost[y*m_width+x ]; dd[y*m_width+x] =0;
 for(int d = 1; d< m_nLabels; d++)
 if(mrg_min[y*m_width+x ] >m_cost[y*m_width+x + d*m_nPixels])
 {mrg_min[y*m_width+x ] = m_cost[y*m_width+x + d*m_nPixels];  dd[y*m_width+x] =d; }
}
float ss = 0; for(int i =0; i<m_nPixels; i++)ss += mrg_min[i]; ss /=m_nPixels;
//Cnst_2 = ss;
//for(int i =0; i<m_nPixels; i++){mrg_min[i + m_nPixels] =mrg_min[i]; mrg_min[i + m_nPixels*2] =1;}
//SetSort(m_nPixels,&mrg_min[m_nPixels], &mrg_min[m_nPixels*2] );
float pw = (m_nLabels > 20)? 3:1.5;

float thr = ss*1.32; //352;//352;//mrg_min[mrg_min[(int)((float)m_nPixels/div) + m_nPixels*2]];
for(int i =0; i<m_nPixels; i++){ float sum=0;
//for(int d = 0; d< m_nLabels; d++)sum += m_cost[i + d*m_nPixels] = (d==dd[i])? 0:1;
//for(int d = 0; d< m_nLabels; d++)sum += m_cost[i + d*m_nPixels] = (m_cost[i + d*m_nPixels]<thr)? 0:1;
for(int d = 0; d< m_nLabels; d++){ float ddw =1;/*(d <= m_nLabels/10)? 2-d*10/m_nLabels:1;*/ float ex =ddw*m_cost[i + d*m_nPixels]/thr; ex = pow(ex, pw); m_cost[i + d*m_nPixels] = 1. - exp(-ex);}
//if(sum)for(int d = 0; d< m_nLabels; d++)m_cost[i+ d*m_nPixels] /=sum;
}
delete [] mrg_min;

}
			    inline float EDP:: mean_min_sm( float * m_cost)
{

int *mrg_min  = new int [m_nPixels*2];
int *dd = &mrg_min[m_nPixels];
for(int y =0; y<m_height; y++)
for(int x =m_width-1; x>=0; x --)
{

 mrg_min[y*m_width+x ] = m_cost[y*m_width+x ]; dd[y*m_width+x] =0;
 for(int d = 1; d< m_nLabels; d++)
 if(mrg_min[y*m_width+x ] >m_cost[y*m_width+x + d*m_nPixels])
 {mrg_min[y*m_width+x ] = m_cost[y*m_width+x + d*m_nPixels];  dd[y*m_width+x] =d; }
}
float ss = 0; for(int i =0; i<m_nPixels; i++)ss += mrg_min[i]; ss /=m_nPixels;

delete [] mrg_min;
return ss;
}
			    inline void EDP:: prf_cost_( float * m_cost)
{
float *cs = new float [m_nPixels*m_nLabels];
for(int y =0; y<m_height; y++)
for(int x =m_width-1; x>=0; x --)
 for(int d = 0; d< m_nLabels; d++)
{
int ic = x + y*m_width + d*m_nPixels;
int xp =     (x+1 < m_width)? x+1 : x; int xm =     (x- 1 >= 0)? x-1 : x;
int dp =     (d+1 < m_nLabels)? d+1 : d; int dm =     (d - 1 >= 0)? d -1 : d;
int icp1  = xp + y*m_width + dp*m_nPixels;
int icm1 = xm + y*m_width + dm*m_nPixels;
int icp2  = xp + y*m_width + dm*m_nPixels;
int icm2 = xm + y*m_width + dp*m_nPixels;
float cv = m_cost[ic];
float cv1 = (m_cost[icp1]+ m_cost[icm1])/2;
float cv2 = (m_cost[icp2]+ m_cost[icm2])/2;
if(cv1<cv) cv = cv1;
if(cv2<cv) cv = cv2;
cs[ic] = cv;
}
memcpy(m_cost, cs, m_nPixels*m_nLabels*sizeof(float));
delete [] cs;
}
			    inline void EDP:: prf_cost_inv_( float * m_cost)
{
	int sig =1;
	int * rz = new int [m_nPixels];
	for(int l = 0; l <m_nLabels; l++)GaussCosConv2DFst(sig, m_width,m_height, &m_cost[m_nPixels*l]);

this->min_cost(m_cost, rz);
for(int i =0; i<m_nPixels; i++){ int rzv = rz[i];
for(int d = 0; d< m_nLabels; d++) m_cost[i + d*m_nPixels] = (d==rzv)? 0:1;

}
delete [] rz;

}
	    inline void  EDP:: max_cost( float * m_cost, int * dsp)
{

for(int y =0; y<m_height; y++)
for(int x =m_width-1; x>=0; x --)
{

 dsp[y*m_width+x]=0; float mrg_max = m_cost[y*m_width+x ];
 for(int d = 1; d< m_nLabels; d++)
 if(mrg_max <m_cost[y*m_width+x + d*m_nPixels])
 {dsp[y*m_width+x]=d; mrg_max = m_cost[y*m_width+x + d*m_nPixels];}

}

}
			    inline void  EDP:: max_cost_do( float * m_cost)
{

for(int y =0; y<m_height; y++)
for(int x =m_width-1; x>=0; x --)
{

 int dsp=0; float mrg_max = m_cost[y*m_width+x ];
 for(int d = 1; d< m_nLabels; d++)
 if(mrg_max <m_cost[y*m_width+x + d*m_nPixels])
 {dsp=d; mrg_max = m_cost[y*m_width+x + d*m_nPixels];}
 for(int d = 0; d< m_nLabels; d++) m_cost[y*m_width+x + d*m_nPixels] = (d==dsp) ? 1:0;
}

}
	    inline void  EDP:: med_cost( float * m_cost, int * dsp)
{

for(int y =0; y<m_height; y++)
for(int x =m_width-1; x>=0; x --)
{
	 dsp[y*m_width+x]=0;
float sum =0;  for(int d = 0; d< m_nLabels; d++)sum += m_cost[y*m_width+x + d*m_nPixels];
float thr = sum/2;   sum =0;  for(int d = 0; d< m_nLabels; d++){sum += m_cost[y*m_width+x + d*m_nPixels]; if(sum> thr){dsp[y*m_width+x]=d;  d = m_nLabels;}}


}

}
 EDP::REAL EDP:: findS(int y)
{
REAL Sum_opt =0;
Msg_X_pp = new REAL [Tau];
SLTN_CRV * Slt = &this->Sltn_crv[y*WpH];
char *dS = new char [Tau];

CostVal * pnt_yD = &m_D[y*Tau];
float * pnt_yD_pL = &this->m_D_pL[y*Tau];
float * pnt_yD_pR = &this->m_D_pR[y*Tau];
//REAL * My_Z = new REAL [Tau]; memset( My_Z, 0,  Tau*sizeof(REAL));
//REAL * My_pp = (y<m_height-1) ? &this->Msg_Y_mm[(y+1)*Tau] : My_Z;
//REAL * My_mm = (y>0) ? &this->Msg_Y_pp[(y-1)*Tau] : My_Z;
Msg_X_pp[0] =0;//My_pp[0] +My_mm[0];
dS[0] =0;// L
/////////////////// direct///////////////////////
int tau_bk;
for(int t =1; t<Tau; t++)
{
	//REAL Msg_y = My_pp[t] +My_mm[t];
REAL vl[3] = {MY_INF, MY_INF, MY_INF};
//-------- Visible dS[t] =0
if((tau_bk = V_tau_bck[t])>=0) vl[0]  = /*Msg_y */+ Msg_X_pp[tau_bk] + ( (pnt_yD[tau_bk]<pnt_yD[t]) ? pnt_yD[tau_bk]: pnt_yD[t]); //Msg_X_pp[tau_bk] + (pnt_yD[tau_bk]+pnt_yD[t])/2;
//-------- L  dS[t] =1
if((tau_bk = L_tau_bck[t])>=0)vl[1]  = /*Msg_y */+ Msg_X_pp[tau_bk] +  ( (pnt_yD_pR[tau_bk]<pnt_yD_pR[t]) ? pnt_yD_pR[tau_bk]: pnt_yD_pR[t]);
//-------- R  dS[t] =2
if((tau_bk = R_tau_bck[t])>=0)vl[2]  = /*Msg_y*/ + Msg_X_pp[tau_bk] +   ( (pnt_yD_pL[tau_bk]<pnt_yD_pL[t]) ? pnt_yD_pL[tau_bk]: pnt_yD_pL[t]);// (pnt_yD_pnl[tau_bk]+pnt_yD_pnl[t])/2;
/////////////////////// Choice/////////////////
int ds_ch =1; REAL min_vl = vl[1];
if(vl[0]<=min_vl){ min_vl = vl[0]; ds_ch=0;};
if(vl[2]<=min_vl)ds_ch=2;
dS[t] = ds_ch;
Msg_X_pp[t] = vl[ds_ch];
/////////////////////////////////////////////
}


///ini solution zero
if(dS[Tau-1]==0){Sum_opt = Msg_X_pp[Tau-1]  ; Slt[0].lb=0;}
else { Sum_opt = Msg_X_pp[Tau-1]; Slt[0].lb = 4; }
Slt[0].tau = Tau-1;
T_max_y[y] = 1;
//________________
int tau_i;
if(dS[Tau-1] == 0 )tau_i = V_tau_bck[Tau-1];
if(dS[Tau-1] == 1 )tau_i = L_tau_bck[Tau-1];
if(dS[Tau-1] == 2 )tau_i = R_tau_bck[Tau-1];
//________________
int st =1;  while(st)
{
int tt, lbl;

Slt[T_max_y[y]].tau = tau_i;
tau_bk = Slt[T_max_y[y]-1].tau;

if(dS[tau_bk] !=1 && dS[tau_i]!=2)lbl = 0;
if(dS[tau_bk] ==1 && dS[tau_i]==2)lbl = 5;
if(dS[tau_bk]==1 && dS[tau_i]!=2)
{
if(XR_tau[tau_i]!=0 )lbl = 1;
else lbl = 3;
}
if(dS[tau_bk] !=1 && dS[tau_i]==2)
{
if(XL_tau[tau_i]!=m_width-1 )lbl = 2;
else lbl = 4;
}
Slt[T_max_y[y]].lb =    lbl;

if(tau_i)
{
if(dS[tau_i] == 0 )tt = V_tau_bck[tau_i];
if(dS[tau_i] == 1 )tt = L_tau_bck[tau_i];
if(dS[tau_i] == 2 )tt = R_tau_bck[tau_i];
tau_i = tt;
}
else st = 0;
T_max_y[y]++;
}
int ii=0;
for(int i = 0; i <T_max_y[y]; i++)
{
if(Slt[i].lb!=5) {  Slt[ii].lb = Slt[i].lb;  Slt[ii].tau = Slt[i].tau; ii++; }
else
{
int a=0;
}
}



T_max_y[y]= ii;
//----------------------------------
delete [] dS;
delete [] Msg_X_pp;
//delete [] My_Z;
return Sum_opt;

}
 EDP::REAL EDP:: findSb(int y)
{

REAL Sum_opt =0;
Msg_X_mm = new REAL [Tau];
SLTN_CRV * Slt = &this->Sltn_crv[y*WpH];
char *dSb = new char [Tau];
CostVal * pnt_yD = &m_D[y*Tau];
float * pnt_yD_pL = &this->m_D_pL[y*Tau];
float * pnt_yD_pR = &this->m_D_pR[y*Tau];
//REAL * My_Z = new REAL [Tau]; memset( My_Z, 0, Tau*sizeof(REAL));
//REAL * My_pp = (y<m_height-1) ? &this->Msg_Y_mm[(y+1)*Tau] : My_Z;
//REAL * My_mm = (y>0) ? &this->Msg_Y_pp[(y-1)*Tau] : My_Z;
Msg_X_mm[Tau-1] = 0;//My_pp[Tau-1] +My_mm[Tau-1];
dSb[Tau-1] =0;// L
int tau_fr;
/////////////////// direct///////////////////////
for(int t =Tau-2; t>=0; t--)
{
	//REAL Msg_y = My_pp[t] +My_mm[t];
REAL vl[3] = {MY_INF, MY_INF, MY_INF};
//-------- Visible dS[t] =0
if((tau_fr = V_tau_frw[t])>=0) vl[0] = Msg_X_mm[tau_fr] + ((pnt_yD[tau_fr] < pnt_yD[t]) ?  pnt_yD[tau_fr] : pnt_yD[t]);                              //    (pnt_yD[tau_fr] + pnt_yD[t])/2;
//-------- L  dS[t] =1
if((tau_fr = L_tau_frw[t])>=0)vl[1]  = Msg_X_mm[tau_fr] + ((pnt_yD_pR[tau_fr] < pnt_yD_pR[t]) ?  pnt_yD_pR[tau_fr] : pnt_yD_pR[t]);   //(pnt_yD_pR[tau_fr] + pnt_yD_pR[t])/2;
//-------- R  dS[t] =2
if((tau_fr = R_tau_frw[t])>=0)vl[2]  = Msg_X_mm[tau_fr] +  ((pnt_yD_pL[tau_fr] < pnt_yD_pL[t]) ?  pnt_yD_pL[tau_fr] : pnt_yD_pL[t]);   //( pnt_yD_pL[tau_fr] + pnt_yD_pL[t])/2;
/////////////////////// Choice/////////////////
int ds_ch =2; REAL min_vl = vl[2];
if(vl[0]<=min_vl){ min_vl = vl[0]; ds_ch=0;};
if(vl[1]<=min_vl)ds_ch=1;
dSb[t] = ds_ch;
Msg_X_mm[t] = vl[ds_ch];

}


///ini solution zero
if(dSb[0]==0){Sum_opt = Msg_X_mm[0] ; Slt[0].lb=0;  }
else { Sum_opt = Msg_X_mm[0];  Slt[0].lb= 3; }
Slt[0].tau = 0;
T_max_y[y] = 1;
//________________
int tau_i;
if(dSb[0] == 0 )tau_i = V_tau_frw[0];
if(dSb[0] == 1 )tau_i = L_tau_frw[0];
if(dSb[0] == 2 )tau_i = R_tau_frw[0];
//---------------------------------------------------
int st =1;  while(st)
{
int tt, lbl;

Slt[T_max_y[y]].tau = tau_i;
tau_fr = Slt[T_max_y[y]-1].tau;

if(dSb[tau_i] !=1 && dSb[tau_fr]!=2)lbl = 0;
if(dSb[tau_i] ==1 && dSb[tau_fr]==2)lbl = 5;
if(dSb[tau_i] ==1 && dSb[tau_fr]!=2)
{
if(XR_tau[tau_i]!=0 )lbl = 1;
else lbl = 3;
}
if(dSb[tau_i] !=1 && dSb[tau_fr]==2)
{
if(XL_tau[tau_i]!=m_width-1 )lbl = 2;
else lbl = 4;
}

Slt[T_max_y[y]].lb =    lbl;

if(tau_i<Tau-1)
{
if(dSb[tau_i] == 0 )tt = V_tau_frw[tau_i];
if(dSb[tau_i] == 1 )tt = L_tau_frw[tau_i];
if(dSb[tau_i] == 2 )tt = R_tau_frw[tau_i];
tau_i = tt;
}
else st = 0;
T_max_y[y]++;
}
int ii=0;
for(int i = 0; i <T_max_y[y]; i++)
{
if(Slt[i].lb!=5) {  Slt[ii].lb = Slt[i].lb;  Slt[ii].tau = Slt[i].tau; ii++; }
else
{
int a=0;
}
}
T_max_y[y]= ii;
//----------------------------------
delete [] dSb;
delete [] Msg_X_mm;
//delete [] My_Z;
return Sum_opt;
}

void EDP:: Make_G_lw_hi( unsigned char * G_buf,  unsigned char * I_buf, int dir, float thr)
{
	    float * dir_pp = new float [m_nPixels*3];
		float * dir_mm = new float [m_nPixels*3];
		int r =5; float c_r = 0.5; int dr = 1;
	Make_G_fromBiF( dir_pp,  I_ims[0], r, c_r, dr);
		 r =5;  c_r = 0.5;  dr = 3;
	Make_G_fromBiF( dir_mm,  I_ims[0], r, c_r, dr);
	for(int i =0; i< m_width; i++)
	for(int j =0; j< m_height; j++)
	{		 float sum =0;  if(j<m_height - 1)for(int c =0; c< 3; c++)sum  += fabs(dir_pp[IND_IC(i,j,c)]- dir_mm[IND_IC(i,j+1,c)]);
	sum = (sum<21)? 255 :0;
             for(int c =0; c< 3; c++)  I_ims[2][IND_IC(i,j,c)] =sum;
	}
int h_size = 256*3;
int Histo[256*3];  memset(Histo, 0, sizeof(int)*h_size);
int * g_c  =  new int [m_nPixels];
///////////////////////////////////////////////////////////
if(dir ==0)for(int i=0; i<m_nPixels; i++)
{int x = i%m_width; int ip = (x<m_width-1)? i+1:i;
g_c[i]=0;  for(int c=0; c<3; c++)g_c[i] += abs(I_buf[ip+c*m_nPixels]-I_buf[i+c*m_nPixels]);}
//
//if(dir ==2)for(int i=0; i<m_nPixels; i++)
//{int x = i%m_width; int im = (x>0)? i-1:i;
//g_c[i]=0;  for(int c=0; c<3; c++)g_c[i] += abs(I_buf[im+c*m_nPixels]-I_buf[i+c*m_nPixels]);}
////
if(dir ==1)for(int i=0; i<m_nPixels; i++)
{int y = i/m_width; int ip = (y<m_height-1)? i+m_width:i;
g_c[i]=0;  for(int c=0; c<3; c++)g_c[i] += abs(I_buf[ip+c*m_nPixels]-I_buf[i+c*m_nPixels]);}
//
//if(dir ==3)for(int i=0; i<m_nPixels; i++)
//{int y = i/m_width; int ip = (y<m_height-1)? i+m_width:i;
//g_c[i]=0;  for(int c=0; c<3; c++)g_c[i] += abs(I_buf[ip+c*m_nPixels]-I_buf[i+c*m_nPixels]);}
/////////////////////////////////////////////////////////////////
for(int i=0; i<m_nPixels; i++) Histo[g_c[i]]++;
int his_nmb =0; for(int i=0;i<h_size;i++)his_nmb += Histo[i];     int thr_st = ((float)his_nmb*thr);
int i_h=0; int st=1; int sum=0; while (st){sum+= Histo[i_h++];if(sum>thr_st)st=0; } int i_thr= i_h;
////////////////////////////////////////////////////////
for(int i=0;i<m_nPixels;i++)  G_buf[i]=(g_c[i]>=i_thr)? 1:0;
delete [] g_c;
delete [] dir_pp;
delete [] dir_mm;
}
void EDP:: Make_G_lw_hi_( unsigned char * G_buf,  unsigned char * I_buf, int dir, float thr)
{
int h_size = 256*3;
int Histo[256*3];  memset(Histo, 0, sizeof(int)*h_size);
int * g_c  =  new int [m_nPixels];
///////////////////////////////////////////////////////////
if(dir ==0)for(int i=0; i<m_nPixels; i++)
{int x = i%m_width; int ip = (x<m_width-1)? i+1:i;
g_c[i]=0;  for(int c=0; c<3; c++)g_c[i] += (abs(I_buf[ip+c*m_nPixels]-I_buf[i+c*m_nPixels])); g_c[i] =sqrt ((float)g_c[i]);}

if(dir ==1)for(int i=0; i<m_nPixels; i++)
{int y = i/m_width; int ip = (y<m_height-1)? i+m_width:i;
g_c[i]=0;  for(int c=0; c<3; c++)g_c[i] += abs(I_buf[ip+c*m_nPixels]-I_buf[i+c*m_nPixels]);  g_c[i] =sqrt ((float)g_c[i]);}
/////////////////////////////////////////////////////////////////
float sum =0; for(int i=0; i<m_nPixels; i++) sum += g_c[i]*g_c[i] ;
sum = sqrt(sum/m_nPixels/2);
////////////////////////////////////////////////////////
for(int i=0;i<m_nPixels;i++)  G_buf[i]= ((float)C_gr_1*g_c[i]/sum>C_gr_1)? 0 :  (BYTE)(C_gr_1 - (float)C_gr_1*g_c[i]/sum);
delete [] g_c;
}
void EDP:: Make_G_fromBiF( float * fltb,  unsigned char * I_buf, int r_w, float cl_prc, int dir)
{

   unsigned char rgb[3];
	float rgb_m[3] ={0,0,0};
	for(int i= 0; i< m_nPixels;i ++)
	for(int c= 0; c< 3;c ++)
	{	rgb_m[c] += I_buf[i+m_nPixels*c];}
	for(int c= 0; c< 3;c ++)
	{	rgb_m[c] /= m_nPixels;}
	float vll, std =0;
	for(int i= 0; i< m_nPixels;i ++)
	{ vll =0;
	for(int c= 0; c< 3;c ++)
	{ float vl = (I_buf[i+m_nPixels*c]-rgb_m[c]); vl *=vl; vll +=vl;}
	std += vll;
	}
	int cl_fc = round_fl(sqrt(std/m_nPixels)*cl_prc);
	///////////////////////////////////////////////
	int dm = 2*r_w+1;
	int size_w = dm*dm;


	/////////////////////////////////////////
	float *gss_r_wt = new float [m_width+m_height];
	float *gss_c_wt = new float [1000];
	float sgm_r = r_w*0.5;

	for(int i=0; i<m_width+m_height; i++ )gss_r_wt[i]= exp(-(i*i)/(2.*sgm_r*sgm_r));
	for(int i=0; i<1000; i++ )gss_c_wt[i]= exp(-(i*i)/(2.*cl_fc*cl_fc));

    ///////////////////////////////////////////////////////
	for(int i=0; i<m_nPixels; i++)BiF_wnd( &fltb[i], i, dir,  I_buf, r_w,gss_c_wt, gss_r_wt);
	
///////////////////////////////////////////////////////
delete [] gss_r_wt;
delete [] gss_c_wt;

}
void EDP:: Make_Inflc_Ptn( POINT_F * fltb,  unsigned char * I_buf, int r_w, float cl_prc, int dir )
{

   unsigned char rgb[3];
   int dir_x, dir_y; if(dir ==0) { dir_x = 0 ; dir_y=1; }
                                if(dir ==2) { dir_x = 2;  dir_y=3; }
	float rgb_m[3] ={0,0,0};
	for(int i= 0; i< m_nPixels;i ++)
	for(int c= 0; c< 3;c ++)
	{	rgb_m[c] += I_buf[i+m_nPixels*c];}
	for(int c= 0; c< 3;c ++)
	{	rgb_m[c] /= m_nPixels;}
	float vll, std =0;
	for(int i= 0; i< m_nPixels;i ++)
	{ vll =0;
	for(int c= 0; c< 3;c ++)
	{ float vl = (I_buf[i+m_nPixels*c]-rgb_m[c]); vl *=vl; vll +=vl;}
	std += vll;
	}
	double STD_IMG = sqrt(std/m_nPixels);
	int cl_fc = round_fl(sqrt(std/m_nPixels)*cl_prc);
	///////////////////////////////////////////////
	int dm = 2*r_w+1;
	int size_w = dm*dm;


	/////////////////////////////////////////
	float *gss_r_wt = new float [m_width+m_height];
	float *gss_c_wt = new float [1000];
	float sgm_r = r_w*0.5;

	for(int i=0; i<m_width+m_height; i++ )gss_r_wt[i]= exp(-(i*i)/(2.*sgm_r*sgm_r));
	for(int i=0; i<1000; i++ )gss_c_wt[i]= exp(-(i*i)/(2.*cl_fc*cl_fc));

    ///////////////////////////////////////////////////////
	for(int i=0; i<m_nPixels; i++)
	{ int x = i%m_width; int y = i/m_width;
     if(x > r_w && x < m_width - r_w && y > r_w && y < m_height - r_w)
	{fltb[i].x = BiF_wnd_( i, dir_x,  I_buf, r_w,gss_c_wt, gss_r_wt); fltb[i].y = BiF_wnd_( i, dir_y,  I_buf, r_w,gss_c_wt, gss_r_wt);}
	 else
	 {fltb[i].x = 0.99; fltb[i].y =  0.99;}
	}
	//for(int i=0; i<m_nPixels; i++){fltb[i].x = (fltb[i].x>INFL_TR) ? fltb[i].x:INFL_TR; fltb[i].y = (fltb[i].y>INFL_TR) ? fltb[i].y:INFL_TR;}
	
///////////////////////////////////////////////////////
delete [] gss_r_wt;
delete [] gss_c_wt;

}
void EDP:: BFfilter(unsigned char * inout, float sigma, float sig_cl, int m_q)
{
float * bf = new float [3*m_nPixels];
       for(int i =0; i < m_nPixels*3; i++ )bf[i] = inout[i];
	   K_mean_Flt_Add(inout, sig_cl, sigma, m_q, bf, 3);
	   for(int i =0; i < m_nPixels*3; i++ ){ int r = round_fl(bf[i]); inout[i] = (r<0)? 0:(r>255)? 255: r;}

          for(int i =0; i < m_nPixels*3; i++ )I_ims[2][i] = bf[i];
delete [] bf;
}
void EDP:: BFfilter(int * inout, float sigma, float sig_cl, int m_q, int lr)
{
double * bf = new double [m_nPixels];
unsigned char * sup = (lr) ? I_ims[1] : I_ims[0];
unsigned char * supp = new BYTE [m_nPixels*3];
for(int i =0; i < m_nPixels*3; i++ )supp[i] = inout[i%N_PX];//(i >=m_nPixels) ?  sup[i] : inout[i];
       for(int i =0; i < m_nPixels; i++ )bf[i] = inout[i];
	   K_mean_Flt_Add_new(supp, sig_cl, sigma, sigma, m_q, bf, 1);
	   for(int i =0; i < m_nPixels; i++ ){ int r = round_fl(bf[i]); inout[i] = (r<0)? 0:(r>255)? 255: r;}

          /*for(int i =0; i < m_nPixels*3; i++ )I_ims[2][i] = bf[i];*/
delete [] bf;
delete [] supp;
}

void EDP::Make_Gr_inf () {
	int H[1000];
	for(int i =0; i <1000; i++)
		H[i]=0;
	for(int i =0; i < m_nPixels; i++ ) {
		int x = i%m_width; int y = i/m_width;
		int xpp = (x+1<m_width)? x+1:x; int ypp = (y+1<m_height)? y+1:y;
		int iypp = ypp*m_width + x;
		int ixpp = y*m_width + xpp;
		int dx = 0;
		for(int c=0;c <3; c++) {
			dx += abs(I_ims[0][i + c*m_nPixels] - I_ims[0][ixpp + c*m_nPixels]);
		}
		int dy = 0;
		for(int c=0;c <3; c++) {
			dy += abs(I_ims[0][i + c*m_nPixels] - I_ims[0][iypp + c*m_nPixels]);
		}
		H[dx]++; H[dy]++;
	}
	int first = m_nPixels*2/3; int scnd = m_nPixels*4/3;
	int th1 = -1, th2 = -1; int sum = 0;
	for(int i =0; i <1000; i++){
		sum += H[i];
		if(th1 < 0) {
			if(sum >= first)
				th1 =i;
		}
		if(th2 < 0) {
			if(sum >= scnd){
				th2 =i; i=1001;
			}
		}
	}
	//Cnst_1 = th1; Cnst_2 = th2;
	Make_Gr_inf(0, th1, th2);
}

void EDP::Make_Gr_inf_pp(float thr_ar, int thr,  int thr1, int thr2) {
	// make copies of images.
	BYTE * li = new BYTE [m_nPixels*3];
	for(int p =0; p <3*m_nPixels; p++){
		li[p] = I_ims[0][p];
	}
	BYTE * ri = new BYTE [m_nPixels*3];
	for(int p=0; p <3*m_nPixels; p++){
		ri[p] = I_ims[1][p];
	}
	//===
    int * WMaskR = new int [m_nPixels];
	int * WMaskL = new int [m_nPixels];
	WMask = new int [m_nPixels];
	ResInd   = new POINT [m_nPixels];

	//--------------------------------
	GetClrMask(thr, li); 
	FOR_PX_p {
		if(p>m_width&& p< N_PX -m_width){
			int a1 = WMask[p] - WMask[p+1];
			int a2 = WMask[p] - WMask[p-1];
			int a3 = WMask[p] - WMask[p+m_width];
			int a4 = WMask[p] - WMask[p-m_width];
			WMaskL[p] = ((!a1)&&(!a2)&&(!a3)&&(!a4))? ResInd[WMask[p]+1].x :1;
		}
	}
	GetClrMask(thr2, ri);
	FOR_PX_p {
		if(p>m_width&& p< N_PX -m_width){
			int a1 = WMask[p] - WMask[p+1];
			int a2 = WMask[p] - WMask[p-1];
			int a3 = WMask[p] - WMask[p+m_width];
			int a4 = WMask[p] - WMask[p-m_width];
			WMaskR[p] = ((!a1)&&(!a2)&&(!a3)&&(!a4))? ResInd[WMask[p]+1].x :1;
		}
	}

	for(int i =0; i < m_nPixels; i++ ){
		B_L_buf_pp_[i].x = 2; B_L_buf_mm_[i].x = 2;
		B_L_buf_pp_[i].y = 2; B_L_buf_mm_[i].y = 2;
		B_R_buf_pp_[i].x = 2; B_R_buf_mm_[i].x = 2;
		B_R_buf_pp_[i].y = 2; B_R_buf_mm_[i].y = 2;
	}
	for(int i =0; i < m_nPixels; i++ ) {
		int x = i%m_width; int y = i/m_width;
		int xpp = (x+1<m_width)? x+1:x; int ypp = (y+1<m_height)? y+1:y;
		int iypp = ypp*m_width + x;
		int ixpp = y*m_width + xpp;
		//int xmm = (x >0)? x-1:x; int ymm = (y>0)? y-1:y;
		//int iymm = ymm*m_width + x;
		//   int ixmm = y*m_width + xmm;

		int dlx =0; int dly =0; int drx =0; int dry =0;
		for(int c = 0; c< 3; c++) {
			int dlxc = abs(li[i + c*m_nPixels]- li[ixpp + c*m_nPixels]);
			int dlyc = abs(li[i+ c*m_nPixels]- li[iypp + c*m_nPixels]);
			int drxc = abs(ri[i + c*m_nPixels]- ri[ixpp + c*m_nPixels]);
			int dryc = abs(ri[i + c*m_nPixels]- ri[iypp + c*m_nPixels]);
			dlx +=dlxc; dly += dlyc; drx += drxc; dry += dryc;
		}
      ////////////// L
		B_L_buf_pp_[i].x = 2;
		B_L_buf_mm_[ixpp].x = 2;
		if(dlx<thr2) {
			B_L_buf_pp_[i].x = 1;
			B_L_buf_mm_[ixpp].x = 1;
		}
		if(dlx<thr1) {
			B_L_buf_pp_[i].x = 0;
			B_L_buf_mm_[ixpp].x = 0;
		}
		if(WMaskL[i] > thr_ar*N_PX)	{
			B_L_buf_pp_[i].x = 0;
			B_L_buf_mm_[ixpp].x = 0;
		}
		B_L_buf_pp_[i].y = 2;
		B_L_buf_mm_[iypp].y = 2;
		if(dly<thr2) {
			B_L_buf_pp_[i].y = 1;
			B_L_buf_mm_[iypp].y = 1;
		}
		if(dly<thr1) {
			B_L_buf_pp_[i].y = 0;
			B_L_buf_mm_[iypp].y = 0;
		}
		if(WMaskL[i]  > thr_ar*N_PX) {
			B_L_buf_pp_[i].y = 0;
			B_L_buf_mm_[iypp].y = 0;
		}
	 //     ////////////// R
		B_R_buf_pp_[i].x = 2;
		B_R_buf_mm_[ixpp].x = 2;
		if(drx<thr2) {
	 		B_R_buf_pp_[i].x = 1;
			B_R_buf_mm_[ixpp].x = 1;
		}
		if(drx<thr1) {
			B_R_buf_pp_[i].x = 0;
			B_R_buf_mm_[ixpp].x = 0;
		}
		if(WMaskR[i]   > thr_ar*N_PX) {
			B_R_buf_pp_[i].x = 0;
			B_R_buf_mm_[ixpp].x = 0;
		}
	 //
		B_R_buf_pp_[i].y = 2;
		B_R_buf_mm_[iypp].y = 2;
		if(dry<thr2) {
			B_R_buf_pp_[i].y = 1;
			B_R_buf_mm_[iypp].y = 1;
		}
		if(dry<thr1) {
			B_R_buf_pp_[i].y = 0;
			B_R_buf_mm_[iypp].y = 0;
		}
		if(WMaskR[i]   > thr_ar*N_PX) {
			B_R_buf_pp_[i].y = 0;
			B_R_buf_mm_[iypp].y = 0;
		}
	 ////////////////////////

	}
	delete [] ri;
	delete [] li;
	delete [] WMaskR;
	delete [] WMaskL;
}

void EDP::Outlr_reg_filter(float thr_ar, int thr, int * rez, int lr, int r_w, float sig_c) {
	BYTE * rezb = new BYTE [m_nPixels*3]; // one BYTE for each color of each pixel.
	for(int p =0; p <3*m_nPixels; p++){
		rezb[p] = rez[p%N_PX]; // wut
	}
	BYTE * ol_msk = new BYTE [N_PX]; // one for each pixel.
	FOR_PX_p {
		ol_msk[p] =0;
	}
	
	//===
	WMask = new int [m_nPixels];
	ResInd = new POINT [m_nPixels];

	//--------------------------------
	GetClrMask(thr, rezb);
	float cof = thr_ar* N_PX;
	int cnt = 0;
	
	FOR_PX_p {
		int x = p%m_width;
		float ar = (float)ResInd[WMask[p]+1].x/cof;
		cnt +=	ol_msk[p] = (x<N_LB)?( (ar < 4) ? 0:1 ):( (ar < 1) ? 0:1 );//I_ims[2][p] = I_ims[2][p + N_PX] = I_ims[2][p + 2*N_PX] =(ar < 255 )? ar:255;
	}
	cnt = N_PX - cnt;
	/////
	float *gss_r_wt = new float [m_width+m_height];
	float *gss_c_wt = new float [1000];
	float sgm_r = r_w*0.5;
	float cl_fc = sig_c *255;
	for(int i=0; i<m_width+m_height; i++ )
		gss_r_wt[i]= exp(-(i*i)/(2.*sgm_r*sgm_r));
	for(int i=0; i<1000; i++ )
		gss_c_wt[i]= exp(-(i*i)/(2.*cl_fc*cl_fc));
	//----------------------------------
	cnt =3;
	int st =1;
	while (st&&cnt){
		st = 0;
		cnt--;
		FOR_PX_p{
			if(!ol_msk[p]) {
				int out = mask_BF_med(p, ol_msk, rez, r_w, gss_c_wt,  gss_r_wt, lr);
				if(out <0) {
					st = 1;
				} else {
					rez[p] = out;
					ol_msk[p] =1;
				}
			}
		};
	}//------------------------------------
	delete [] gss_r_wt;
	delete [] gss_c_wt;
	/////
	delete [] rezb;
	delete [] ol_msk;
}

void EDP::Make_Gr_inf(int Ql, int thr1, int thr2)
{
	BYTE * li = new BYTE [m_nPixels*3]; for(int p =0; p <3*m_nPixels; p++){li[p] = I_ims[0][p];}
	BYTE * ri =new BYTE [m_nPixels*3]; for(int p=0; p <3*m_nPixels; p++){ri[p] = I_ims[1][p];}
//===
	/*int * WMaskR = new int [m_nPixels];
	int * WMaskL = new int [m_nPixels];
	WMask = new int [m_nPixels];
	ResInd   = new POINT [m_nPixels];
	float thr_ar = 0.01;*/
	//--------------------------------
	/*GetClrMask(thr2, li); FOR_PX_p
	{if(p>m_width&& p< N_PX -m_width){
		int a1 = WMask[p] - WMask[p+1];
		int a2 = WMask[p] - WMask[p-1];
		int a3 = WMask[p] - WMask[p+m_width];
		int a4 = WMask[p] - WMask[p-m_width];
		WMaskL[p] = ((!a1)&&(!a2)&&(!a3)&&(!a4))? ResInd[WMask[p]+1].x :1;
	}
	}
	GetClrMask(thr2, ri); FOR_PX_p
	{if(p>m_width&& p< N_PX -m_width){
		int a1 = WMask[p] - WMask[p+1];
		int a2 = WMask[p] - WMask[p-1];
		int a3 = WMask[p] - WMask[p+m_width];
		int a4 = WMask[p] - WMask[p-m_width];
		WMaskR[p] = ((!a1)&&(!a2)&&(!a3)&&(!a4))? ResInd[WMask[p]+1].x :1;
	}
	}*/
//===
//===
	//int Hs[800];for(int i =0; i <800; i++)Hs[i] =0;//----------
	//for(int i =0; i < m_nPixels; i++ )
	//   {
	//   int x = i%m_width; int y = i/m_width;
	//   int xpp = (x+1<m_width)? x+1:x; int ypp = (y+1<m_height)? y+1:y;
	//   int iypp = ypp*m_width + x;
 //      int ixpp = y*m_width + xpp;
	//   int dlx =0; int dly =0; int drx =0; int dry =0;
	//   for(int c = 0; c< 3; c++)
	//   {
 //      int dlxc = abs(li[i + c*m_nPixels]- li[ixpp + c*m_nPixels]);
	//   int dlyc = abs(li[i + c*m_nPixels]- li[iypp + c*m_nPixels]);
	//   //int drxc = abs(ri[i + c*m_nPixels]- ri[ixpp + c*m_nPixels]);
	//   //int dryc = abs(ri[i + c*m_nPixels]- ri[iypp + c*m_nPixels]);
	//   dlx +=dlxc; dly += dlyc; //drx += drxc; dry += dryc;
	//   }
 //     Hs[dlx]++; Hs[dly]++;

	//}
	//int ii =0; int cnt =0; int st =1; int thr = 2*m_nPixels*0.5; while(st){cnt += Hs[ii]; if(cnt > thr)st=0; else ii++;}

//---------
	//Cnst_1 =ii;
	//for( int c =0; c < 3; c++){FilterGrC(I_ims[0],&li[c*m_nPixels], Cnst_1); FilterGrC(I_ims[1], &ri[c*m_nPixels], Cnst_1);}
	//FilterBF_C( ri, 3, 0.05, 1050); FilterBF_C( li, 3, 0.05, 1050);
       for(int i =0; i < m_nPixels; i++ )
	   {B_L_buf_pp_[i].x = 2; B_L_buf_mm_[i].x = 2;
	   B_L_buf_pp_[i].y = 2; B_L_buf_mm_[i].y = 2;
	   B_R_buf_pp_[i].x = 2; B_R_buf_mm_[i].x = 2;
	   B_R_buf_pp_[i].y = 2; B_R_buf_mm_[i].y = 2; }
       for(int i =0; i < m_nPixels; i++ )
	   {
	   int x = i%m_width; int y = i/m_width;
	   int xpp = (x+1<m_width)? x+1:x; int ypp = (y+1<m_height)? y+1:y;
	   int iypp = ypp*m_width + x;
       int ixpp = y*m_width + xpp;
	   //	   int xmm = (x >0)? x-1:x; int ymm = (y>0)? y-1:y;
	   //int iymm = ymm*m_width + x;
    //   int ixmm = y*m_width + xmm;

	   int dlx =0; int dly =0; int drx =0; int dry =0;
	   for(int c = 0; c< 3; c++)
	   {
       int dlxc = abs(li[i + c*m_nPixels]- li[ixpp + c*m_nPixels]);
	   int dlyc = abs(li[i+ c*m_nPixels]- li[iypp + c*m_nPixels]);
	   int drxc = abs(ri[i + c*m_nPixels]- ri[ixpp + c*m_nPixels]);
	   int dryc = abs(ri[i + c*m_nPixels]- ri[iypp + c*m_nPixels]);
	   dlx +=dlxc; dly += dlyc; drx += drxc; dry += dryc;
	   }
      ////////////// L
	  B_L_buf_pp_[i].x = 2;
	  B_L_buf_mm_[ixpp].x = 2;
	  if(dlx<thr2)
	 {B_L_buf_pp_[i].x = 1;
	  B_L_buf_mm_[ixpp].x = 1;}
	 if(dlx<thr1)
	 {B_L_buf_pp_[i].x = 0;
	  B_L_buf_mm_[ixpp].x = 0;}
	  //if(WMaskL[i] > thr_ar*N_PX)
	  //{B_L_buf_pp_[i].x = 0;
	  //B_L_buf_mm_[ixpp].x = 0;}
	 //
	  B_L_buf_pp_[i].y = 2;
	  B_L_buf_mm_[iypp].y = 2;
	  if(dly<thr2)
	 {B_L_buf_pp_[i].y = 1;
	  B_L_buf_mm_[iypp].y = 1;}
	 if(dly<thr1)
	 {B_L_buf_pp_[i].y = 0;
	  B_L_buf_mm_[iypp].y = 0;}
	 //if(WMaskL[i]  > thr_ar*N_PX)
	 //{B_L_buf_pp_[i].y = 0;
	 // B_L_buf_mm_[iypp].y = 0;}
	 //     ////////////// R
	  B_R_buf_pp_[i].x = 2;
	  B_R_buf_mm_[ixpp].x = 2;
	  if(drx<thr2)
	 {B_R_buf_pp_[i].x = 1;
	  B_R_buf_mm_[ixpp].x = 1;}
	 if(drx<thr1)
	 {B_R_buf_pp_[i].x = 0;
	  B_R_buf_mm_[ixpp].x = 0;}
	 //if(WMaskR[i]   > thr_ar*N_PX)
	 //{B_R_buf_pp_[i].x = 0;
	 // B_R_buf_mm_[ixpp].x = 0;}
	 //
	  B_R_buf_pp_[i].y = 2;
	  B_R_buf_mm_[iypp].y = 2;
	  if(dry<thr2)
	 {B_R_buf_pp_[i].y = 1;
	  B_R_buf_mm_[iypp].y = 1;}
	 if(dry<thr1)
	 {B_R_buf_pp_[i].y = 0;
	  B_R_buf_mm_[iypp].y = 0;}
	 //if(WMaskR[i]   > thr_ar*N_PX) //if(dry<thr1)
	 //{B_R_buf_pp_[i].y = 0;
	 // B_R_buf_mm_[iypp].y = 0;}
	 //(ResInd[WMask[p]+1].x  > 0.01*N_PX)
	 ////////////////////////

	}
delete [] ri;
delete [] li;
//delete [] WMaskR;
//delete [] WMaskL;
}
inline int EDP::round_fl(float vl) {
	return (vl >= 0 ? (int)(vl+0.5) : (int)(vl -0.5));
}
inline int EDP::round_fl(double vl) {
	return (vl >= 0 ? (int)(vl+0.5) : (int)(vl -0.5));
}
/////////////////////////////////////7
void EDP::MkPnlCst_RL(int *LMp, int *RMp) {
	float cst_mul = Cost_mean*PCst_mlt/m_nLabels/m_nLabels;
	int div = 8;
//////////////RMp////////////////

	for(int y =0; y<m_height; y++)
	for(int x =0; x<m_width; x ++) {
		int xinv = m_width-1 -x;
		int d_mx = (xinv + 1 <m_nLabels) ? xinv +1 : m_nLabels;
		float *pcst = &m_D_pR[y*Tau];
		int *xld = &XRD(x);
		//-----------------------------------
		for(int d = 0; d< d_mx; d++)
			pcst[xld[d]]= Cost_mean;//*(m_nLabels*2-d)/m_nLabels;
		if(x>0) {
			int d_x = RMp[y*m_width+x];
			int d_xg = RMp[y*m_width+x-1];
			if(d_xg < d_x) {
				float cst = 0;// PCst_mlt/(d_x- d_xg);
				for(int d = d_x; d>=d_xg; d--)
					pcst[xld[d]] *= cst;
			}
		} else {
			int d_x = RMp[y*m_width+x];
			if(0 < d_x) {
			float cst = 0;//PCst_mlt/(d_x);
			for(int d = d_x; d>= 0; d--)
				pcst[xld[d]] *=  cst;
			}
		}
	}
// end dir ////////////////
	for(int y =0; y<m_height; y++)
	for(int x =0; x<m_width; x ++)
	{
		int d_mx = (x + 1 <m_nLabels) ? x +1 : m_nLabels;
		float *pcst = &m_D_pL[y*Tau];
		int *xld = &XLD(x);
		//-----------------------------------
		for(int d = 0; d< d_mx; d++)
			pcst[xld[d]]=  Cost_mean;//*(m_nLabels*2-d)/m_nLabels;
		if(x<m_width-1) {
			int d_x = LMp[y*m_width+x];
			int d_xg = LMp[y*m_width+x+1];
			if(d_xg < d_x) {
				float cst = 0;//PCst_mlt/(d_x- d_xg);
				for(int d = d_x; d>=d_xg; d--) {
					pcst[xld[d]] *=  cst;
				}
			}
		}else {
			int d_x = LMp[y*m_width+x];
			if(0 < d_x) {
				float cst = 0;//PCst_mlt/(d_x);
				for(int d = d_x; d>=0; d--) {pcst[xld[d]] *=  cst;}
			}
		}
	}
// end dir ////////////////
}
void EDP::MkPnlCst_G(int *LMp)
{
	Cost_mean = 1000;
	//////////////RMp////////////////

	for(int y =0; y<m_height; y++)
	for(int x =0; x<m_width; x ++)
	{
		int xinv = m_width-1 -x;
		int d_mx = (xinv + 1 <m_nLabels) ? xinv +1 : m_nLabels;
		float *pcst = &m_D_pR[y*Tau];
		int *xld = &XRD(x);
		//-----------------------------------
		for(int d = 0; d< d_mx; d++)pcst[xld[d]]= Cost_mean;//*(m_nLabels*2-d)/m_nLabels;
		if(x>0)
		{

			int d_x0 = LMp[y*m_width+x];
			int d_x   = (x-d_x0 >=0)? LMp[y*m_width+x-d_x0]: LMp[y*m_width];
			int d_xg0 = LMp[y*m_width+x-1];
			int d_xg   = (x-d_xg0 >=0)? LMp[y*m_width+x-d_xg0]: LMp[y*m_width];
			if(d_xg < d_x)
			{
				float cst = 0;// PCst_mlt/(d_x- d_xg);
				for(int d = d_x; d>=d_xg; d--)
					pcst[xld[d]] *= cst;
			}
		} else {
			int d_x0 = LMp[y*m_width+x];
			int d_x   = (x-d_x0 >=0)? LMp[y*m_width+x-d_x0]: LMp[y*m_width];
			if(0 < d_x)
			{
				float cst = 0;//PCst_mlt/(d_x);
				for(int d = d_x; d>= 0; d--)
					pcst[xld[d]] *=  cst;
			}
		}

	}
	// end dir ////////////////
	for(int y =0; y<m_height; y++)
	for(int x =0; x<m_width; x ++)
	{
		int d_mx = (x + 1 <m_nLabels) ? x +1 : m_nLabels;
		float *pcst = &m_D_pL[y*Tau];
		int *xld = &XLD(x);
		//-----------------------------------
		for(int d = 0; d< d_mx; d++)pcst[xld[d]]=  Cost_mean;//*(m_nLabels*2-d)/m_nLabels;
		if(x<m_width-1)
		{
			int d_x = LMp[y*m_width+x];
			int d_xg = LMp[y*m_width+x+1];
			if(d_xg < d_x)
			{
				float cst = 0;//PCst_mlt/(d_x- d_xg);
				for(int d = d_x; d>=d_xg; d--) {
					pcst[xld[d]] *=  cst;
				}
			}
		} else {
			int d_x = LMp[y*m_width+x];
			if(0 < d_x)
			{
				float cst = 0;//PCst_mlt/(d_x);
				for(int d = d_x; d>=0; d--) {
					pcst[xld[d]] *=  cst;
				}
			}
		}

	}
// end dir ////////////////
}
/*
	Takes in the disparity for left image and disparity for right image.
	Fixes the two disparity maps.

	@params IN/OUT inL left disparity
	@params IN/OUT inR right disparity

	@returns None - instead, modifies inL, inR to be the "corrected" disparities.
*/
void EDP::MkDspCorr_RL(int *inL, int *inR) {

	int* outL = new int [m_nPixels];
	for(int y =0; y<m_nPixels; y++) {
		outL[y] = inL[y];
	}
	int* outR = new int [m_nPixels];
	for(int y =0; y<m_nPixels; y++){
		outR[y] = inR[y];
	}

	//////////////RMp////////////////

	for(int y =0; y<m_height; y++)
	for(int x =0; x<m_width; x ++) {
		int ind = x + y*m_width; // index of pixel.
		int dr = inR[ind]; // disparity in inR

		int indL =(x+dr<m_width) ? x + dr + y*m_width: m_width-1 + y*m_width; // index of shifted pixel.

		if(x+dr>m_nLabels && x+dr<m_width) { // fixing occlusion, essentially.
			int dl = inL[indL]; // disparity of shifted pixel in other image.
			outR[ind] = (dl<dr) ? dl:dr;
		} else {
			outR[ind] = dr;
		}
		
	}
	// end dir ////////////////

	// same as for loop above, but in other direction.
	for(int y =0; y<m_height; y++)
	for(int x =0; x<m_width; x ++) {
		int ind = x + y*m_width;
		int dl = inL[ind];
		int indR = x - dl + y*m_width;
		if(m_width-1 - x +dl >m_nLabels && x - dl >=0 ) {
			int dr = inR[indR];
			outL[ind] = (dr <dl) ? dr:dl;
		} else {
			outL[ind] = dl;
		}
	}

	memcpy( inL, outL, sizeof(int)*m_nPixels);
	memcpy( inR, outR, sizeof(int)*m_nPixels);
	delete [] outL;
	delete [] outR;
	// end dir ////////////////
}

void EDP::MkDspCorr_pp(int *inL, int *inR, int iLR)
{
	if(iLR ==1) {
		for(int y =0; y<m_height; y++)
		for(int x =0; x<m_width; x ++)
		{
			int xinv = m_width-1 -x;
			int d_mx = (xinv + 1 <m_nLabels) ? xinv +1 : m_nLabels;
			CostVal * mrg = &m_D[y*Tau];
			int *xrd = &XRD(x);
			int d_mrg = inR[y*m_width+x];
			int dl = (x+d_mrg<m_width) ? inL[y*m_width+x+d_mrg] : -1;
			float mul = (dl!=-1)?  ( (dl==d_mrg) ? 1 : (abs(dl-d_mrg)>1) ? 0: 0.7) : 0.5;
			for(int d = 0; d< d_mx; d++)
				mrg[xrd[d]] *= mul;
		}
	} else {
		for(int y =0; y<m_height; y++)
		for(int x =0; x<m_width; x ++)
		{
			int d_mx = (x + 1 <m_nLabels) ? x +1 : m_nLabels;
			CostVal * mrg = &m_D[y*Tau];
			int *xld = &XLD(x);
			int d_mrg = inL[y*m_width+x];
			int dr = (x-d_mrg>=0) ? inR[y*m_width+x-d_mrg] : -1;
			float mul = (dr!=-1)?  ( (dr==d_mrg) ? 1 : (abs(dr-d_mrg)>1) ? 0: 0.7) : 0.5;
			for(int d = 0; d< d_mx; d++)
				mrg[xld[d]] *= mul;

		}
	}
	//memcpy( inL, outL, sizeof(int)*m_nPixels);
	//memcpy( inR, outR, sizeof(int)*m_nPixels);
	//delete [] outL;
	//delete [] outR;
	// end dir ////////////////
}
void EDP::MskRL_Corr( int *in, int iRL)
{

	for(int y =0; y<m_height; y++)
	{
		for(int x =0; x<m_width; x ++)
		{
			int ind = x + y*m_width;
			int indyb = x+  (y-1)*m_width;
			int indxb = x -1 + (y)*m_width;
			if(this->Mask_RL[iRL][ind]/8)
			{
				if(y>0) in[ind] = (in[ind]<in[indyb])? in[ind] : in[indyb];
				if(x>0) in[ind] = (in[ind]<in[indxb])? in[ind] : in[indxb];
			}
		}
		for(int x =m_width-1; x>=0; x --)
		{
			int ind = x + y*m_width;
			int indyb = x+  (y-1)*m_width;
			int indxb = x +1 + (y)*m_width;
			if(this->Mask_RL[iRL][ind]/8)
			{
				if(y>0) in[ind] = (in[ind]<in[indyb])? in[ind] : in[indyb];
				if(x< m_width-1) in[ind] = (in[ind]<in[indxb])? in[ind] : in[indxb];
			}
		}

	}/////////////////////////////////////////
	for(int y =m_height-1; y>=0; y--)
	{
		for(int x =0; x<m_width; x ++)
		{
			int ind = x + y*m_width;
			int indyb = x+  (y+1)*m_width;
			int indxb = x -1 + (y)*m_width;
			if(this->Mask_RL[iRL][ind]/8)
			{
				if(y<m_height-1) in[ind] = (in[ind]<in[indyb])? in[ind] : in[indyb];
				if(x>0) in[ind] = (in[ind]<in[indxb])? in[ind] : in[indxb];
			}
		}
		for(int x =m_width-1; x>=0; x --)
		{
			int ind = x + y*m_width;
			int indyb = x+  (y+1)*m_width;
			int indxb = x +1 + (y)*m_width;
			if(this->Mask_RL[iRL][ind]/8)
			{
				if(y<m_height-1) in[ind] = (in[ind]<in[indyb])? in[ind] : in[indyb];
				if(x< m_width-1) in[ind] = (in[ind]<in[indxb])? in[ind] : in[indxb];
			}
		}

	}/////////////////////////////////////////
}

EDP::REAL EDP::TRWS_R(int itr, EDP:: REAL * Mrg, int *DspM)
{
	REAL max =0, min=0;
	REAL * Msg = new REAL [m_height*4*Tau]; memset(Msg, 0, m_height*Tau*4*sizeof(REAL));
	REAL * Msg_add = new REAL [m_nLabels*m_height*4]; memset(Msg_add, 0, 4*m_height*m_nLabels*sizeof(REAL));
	REAL * Di = new REAL [m_nLabels];
	for(int it =0;it<itr;it++)
	{/////////////////// ALGORITHM //////////////////////////
		for(int p= 0; p<m_height*Tau; p++)
			Mrg[p]= m_D[p];

		for(int y =0; y<m_height; y++)
			for(int x =0; x<m_width; x ++)
			{
				MkMrgR(x, y, Msg, Msg_add, Mrg);
				UpdMsgXR_pp(x, y, Di,Msg, Msg_add, Mrg);
				UpdMsgYR_pp(x, y, Di,Msg, Msg_add, Mrg);
			}// end dir 0 ////////////////

		for(int p= 0; p<m_height*Tau; p++)
			Mrg[p]= m_D[p];
		
		for(int y =m_height-1; y >=0; y--)
			for(int x =m_width-1; x >=0;  x--)
			{
				MkMrgR(x, y, Msg, Msg_add, Mrg);
				UpdMsgXR_mm(x, y, Di,Msg, Msg_add, Mrg);
				UpdMsgYR_mm(x, y, Di,Msg, Msg_add, Mrg);
			}// end dir ////////////////
		for(int p= 0; p<m_height*Tau; p++)
			Mrg[p]= m_D[p];
		for(int y =0; y<m_height; y++)
			for(int x =m_width-1; x >=0;  x--)
			{
				MkMrgR(x, y, Msg, Msg_add, Mrg);
				UpdMsgXR_mm(x, y, Di,Msg, Msg_add, Mrg);
				UpdMsgYR_pp(x, y, Di,Msg, Msg_add, Mrg);
			}// end dir ////////////////
		for(int p= 0; p<m_height*Tau; p++)
			Mrg[p]= m_D[p];
		for(int y =m_height-1; y >=0; y--)
			for(int x =0; x<m_width; x ++)
			{
				MkMrgR(x, y, Msg, Msg_add, Mrg);
				UpdMsgXR_pp(x, y, Di,Msg, Msg_add, Mrg);
				UpdMsgYR_mm(x, y, Di,Msg, Msg_add, Mrg);
			}// end dir ////////////////
	}/////////end itr ///////////
	//////////////////////////////////
	//////////////////////////////////
	for(int p= 0; p<m_height*Tau; p++)
		Mrg[p]= m_D[p];

	for(int y =0; y<m_height; y++)
		for(int x =0; x<m_width; x ++)
		{
			int xinv = m_width-1 -x;
			int d_mrg;
			int d_mx = (xinv + 1 <m_nLabels) ? xinv +1 : m_nLabels;
			MkMrgR(x, y, Msg, Msg_add, Mrg);
			REAL *mrg = &Mrg(y);
			int *xld = &XRD(x);
			DspM[y*m_width+x]=0; REAL mrg_min = mrg[xld[0]];
			for(int d = 1; d< d_mx; d++) {
				if(mrg_min > mrg[xld[d]]){
					d_mrg=DspM[y*m_width+x]=d; mrg_min = mrg[xld[d]];
				}
			}
			//if(!y&&!x)min = mrg_min;
			//else TRUNCATE(min, mrg_min);
			//for(int d = 0; d< d_mx; d++)
			// mrg[xld[d]] =(d-d_mrg)? Cost_mean: 0;
		}// end dir ////////////////
	//for(int p= 0; p<m_height*Tau; p++)Mrg[p] -= min;
	//////////////////////
	delete [] Msg;
	delete [] Di;
	delete [] Msg_add;
	return max-min;
}
void EDP::Mk_DG(int * DspL)
{

Cost_mean = 1000;
for(int y =0; y<m_height; y++)
for(int x =0; x<m_width; x ++)
{
 int d_mx = (x + 1 <m_nLabels) ? x +1 : m_nLabels;
 CostVal * mrg = &m_D[y*Tau];
 int *xld = &XLD(x);
  int d_mrg = DspL[y*m_width+x];

for(int d = 0; d< d_mx; d++)
{
 int ds = abs(d-d_mrg);
 float cst = (ds<1)? 0:Cost_mean;
mrg[xld[d]] = cst;
}
}
//--------------------------------------------
for(int y =0; y<m_height; y++)
{//YYYYYYYYYYYYYYYYYYYYYYY
CostVal * mrg = &m_D[y*Tau];
for(int i =0; i<Tau; i++)Dst_B[i] = -1;

for(int d = m_nLabels-1; d>0; d--)
for(int x =0; x<m_width; x ++)
{
 int d_mx = (x + 1 <m_nLabels) ? x +1 : m_nLabels;
 if(d<d_mx){

 int *xld = &XLD(x);
 int *xrdp = &XLD(x+1);
 int *xrdm = &XLD(x-1);
 //------------------
 if(d!=m_nLabels-1)
 {
  int bs_l = (d!=d_mx-1) ? Dst_B[xld[d+1]] : -1;
  int bs_r = (x<m_width-1) ? Dst_B[xrdp[d+1]]: bs_l;
  if(bs_l != -1||bs_r != -1)
  {
   Dst_B[xld[d]]=1;
   mrg[xld[d]] =Cost_mean*2;
  }
  else Dst_B[xld[d]]=(mrg[xld[d]]==0)? 0 : Dst_B[xld[d]];

 }
//----------------------
 else  Dst_B[xld[d]]=(mrg[xld[d]]==0)? 0 : Dst_B[xld[d]];

 }///d_mx
}////XXXXXXXXXXXXXXXXXXXX

}////////////YYYYYYYYYYYYYYY

for(int y =0; y<m_height; y++)
{//YYYYYYYYYYYYYYYYYYYYYYY
CostVal * mrg = &m_D[y*Tau];
CostVal * mrg_pL = &this->m_D_pL[y*Tau];
CostVal * mrg_pR = &this->m_D_pR[y*Tau];
for(int d = m_nLabels-1; d>0; d--)
for(int x =0; x<m_width; x ++)
{
 int d_mx = (x + 1 <m_nLabels) ? x +1 : m_nLabels;
 if(d<d_mx){

 int *xld = &XLD(x);
 int *xldp = &XLD(x+1);
 int *xldm = &XLD(x-1);
 //------------------
 if(d!=m_nLabels-1)
 {
  CostVal vl0 = mrg[xld[d]];
  CostVal vlp = (x<m_width-1) ? mrg[xldp[d]]:0;
  CostVal vlm = (x>0&&d<d_mx-1) ? mrg[xldm[d]]:0;
  if(vl0>Cost_mean+10&&vlp>Cost_mean+10&&vlm>Cost_mean+10)
  {
   mrg_pL[xld[d]] =Cost_mean;
   mrg_pR[xld[d]] =Cost_mean;
   /*mrg[xld[d]] =Cost_mean;*/
  }


 }
//----------------------


 }///d_mx
}////XXXXXXXXXXXXXXXXXXXX
}////////////YYYYYYYYYYYYYYY

for(int y =0; y<m_height; y++)
{//YYYYYYYYYYYYYYYYYYYYYYY
CostVal * mrg = &m_D[y*Tau];

for(int d = m_nLabels-1; d>0; d--)
for(int x =0; x<m_width; x ++)
{
 int d_mx = (x + 1 <m_nLabels) ? x +1 : m_nLabels;
 if(d<d_mx){

 int *xld = &XLD(x);

 //------------------
 if(d!=m_nLabels-1)
 {
  CostVal vl0 = mrg[xld[d]];

  if(vl0>Cost_mean+10)
  {

   mrg[xld[d]] =Cost_mean;
  }


 }
//----------------------


 }///d_mx
}////XXXXXXXXXXXXXXXXXXXX
}////////////YYYYYYYYYYYYYYY
}
void EDP::Mk_DRL(int * DspL,  int *DspR, int thr)
{
if(thr<1) thr=1;
for(int y =0; y<m_height; y++)
for(int x =0; x<m_width; x ++)
{
 int xinv = m_width-1 -x;
 int d_mx = (xinv + 1 <m_nLabels) ? xinv +1 : m_nLabels;
 CostVal * mrg = &m_D[y*Tau];
 int *xld = &XRD(x);
  int d_mrg = DspR[y*m_width+x];

for(int d = 0; d< d_mx; d++)
{
 int ds = abs(d-d_mrg);
 float cst = (ds<thr)? ( (ds==0)? 0: Cost_mean*ds/thr):Cost_mean;
 mrg[xld[d]] =cst;
}
}

for(int y =0; y<m_height; y++)
for(int x =0; x<m_width; x ++)
{
 int d_mx = (x + 1 <m_nLabels) ? x +1 : m_nLabels;
 CostVal * mrg = &m_D[y*Tau];
 int *xld = &XLD(x);
  int d_mrg = DspL[y*m_width+x];

for(int d = 0; d< d_mx; d++)
{
 int ds = abs(d-d_mrg);
 float cst = (ds<thr)? ( (ds==0)? 0: Cost_mean*ds/thr):Cost_mean;
mrg[xld[d]] = (mrg[xld[d]]>cst)? mrg[xld[d]]:cst;
}

if(thr>1){ double min = mrg[xld[0]]; int min_ind =0;
for(int d = 1; d< d_mx; d++)if(mrg[xld[d]]< min) {mrg[xld[d]]= min; min_ind =d;}
for(int d = 0; d< d_mx; d++)mrg[xld[d]] = (min_ind ==d)? 0:Cost_mean;
}
}
//--------------------------------------------
for(int y =0; y<m_height; y++)
{//YYYYYYYYYYYYYYYYYYYYYYY
CostVal * mrg = &m_D[y*Tau];
for(int i =0; i<Tau; i++)Dst_B[i] = -1;

for(int d = m_nLabels-1; d>0; d--)
for(int x =0; x<m_width; x ++)
{
 int d_mx = (x + 1 <m_nLabels) ? x +1 : m_nLabels;
 if(d<d_mx){

 int *xld = &XLD(x);
 int *xrdp = &XLD(x+1);
 int *xrdm = &XLD(x-1);
 //------------------
 if(d!=m_nLabels-1)
 {
  int bs_l = (d!=d_mx-1) ? Dst_B[xld[d+1]] : -1;
  int bs_r = (x<m_width-1) ? Dst_B[xrdp[d+1]]: bs_l;
  if(bs_l != -1||bs_r != -1)
  {
   Dst_B[xld[d]]=1;
   mrg[xld[d]] =Cost_mean*2;
  }
  else Dst_B[xld[d]]=(mrg[xld[d]]==0)? 0 : Dst_B[xld[d]];

 }
//----------------------
 else  Dst_B[xld[d]]=(mrg[xld[d]]==0)? 0 : Dst_B[xld[d]];

 }///d_mx
}////XXXXXXXXXXXXXXXXXXXX

}////////////YYYYYYYYYYYYYYY

//////////////////////
///////////////////////

for(int y =0; y<m_height; y++)
{//YYYYYYYYYYYYYYYYYYYYYYY
CostVal * mrg = &m_D[y*Tau];
CostVal * mrg_pL = &this->m_D_pL[y*Tau];
CostVal * mrg_pR = &this->m_D_pR[y*Tau];
for(int d = m_nLabels-1; d>0; d--)
for(int x =0; x<m_width; x ++)
{
 int d_mx = (x + 1 <m_nLabels) ? x +1 : m_nLabels;
 if(d<d_mx){

 int *xld = &XLD(x);
 int *xldp = &XLD(x+1);
 int *xldm = &XLD(x-1);
 //------------------
 if(d!=m_nLabels-1)
 {
  CostVal vl0 = mrg[xld[d]];
  CostVal vlp = (x<m_width-1) ? mrg[xldp[d]]:0;
  CostVal vlm = (x>0&&d<d_mx-1) ? mrg[xldm[d]]:0;
  if(vl0>Cost_mean+10&&vlp>Cost_mean+10&&vlm>Cost_mean+10)
  {
   mrg_pL[xld[d]] =Cost_mean;
   mrg_pR[xld[d]] =Cost_mean;
   /*mrg[xld[d]] =Cost_mean;*/
  }


 }
//----------------------


 }///d_mx
}////XXXXXXXXXXXXXXXXXXXX
}////////////YYYYYYYYYYYYYYY

for(int y =0; y<m_height; y++)
{//YYYYYYYYYYYYYYYYYYYYYYY
CostVal * mrg = &m_D[y*Tau];

for(int d = m_nLabels-1; d>0; d--)
for(int x =0; x<m_width; x ++)
{
 int d_mx = (x + 1 <m_nLabels) ? x +1 : m_nLabels;
 if(d<d_mx){

 int *xld = &XLD(x);

 //------------------
 if(d!=m_nLabels-1)
 {
  CostVal vl0 = mrg[xld[d]];

  if(vl0>Cost_mean+10)
  {

   mrg[xld[d]] =Cost_mean;
  }


 }
//----------------------


 }///d_mx
}////XXXXXXXXXXXXXXXXXXXX
}////////////YYYYYYYYYYYYYYY
}
double  EDP:: Mean_tr(double *max)
{
for(int y =0; y<m_height; y++)
for(int x =0; x<m_width; x ++)
{
 int xinv = m_width-1 -x;
 int d_mx = (xinv + 1 <m_nLabels) ? xinv +1 : m_nLabels;
 CostVal * mrg = &m_D[y*Tau];
 int *xld = &XRD(x);
 double min =  mrg[xld[0]];
 for(int d = 1; d< d_mx; d++) min =(min < mrg[xld[d]]) ? min : mrg[xld[d]];
 for(int d = 0; d< d_mx; d++) mrg[xld[d]] -=min;
}

for(int y =0; y<m_height; y++)
for(int x =0; x<m_width; x ++)
{
 int d_mx = (x + 1 <m_nLabels) ? x +1 : m_nLabels;
 CostVal * mrg = &m_D[y*Tau];
 int *xld = &XLD(x);
 double min =  mrg[xld[0]];
 for(int d = 1; d< d_mx; d++) min =(min < mrg[xld[d]]) ? min : mrg[xld[d]];
 for(int d = 0; d< d_mx; d++) mrg[xld[d]] -=min;
}
double sum =0; max[0] = 0;
for(int i =0; i<m_height*Tau; i++){ sum   += m_D[i]; max[0] = (max[0] > m_D[i]) ? max[0] : m_D[i];}
return sum/(m_height*Tau);
}
inline void  EDP::UpdMsgYR_mm(int x, int y,EDP::REAL* Di,EDP::REAL* Msg,EDP::REAL* Msg_add, EDP::REAL* Mrg)
{



if (y >0)
{
 int xinv = m_width-1 -x;
 int d_mx = (xinv + 1 <m_nLabels) ? xinv +1 : m_nLabels;
 subMsg(
 &XRD(x),
 &XRD(x),Di,&Msg(y-1,1),
 &Mrg(y),
 d_mx, 0.5);
 REAL  lmbd =  Smthnss[(int)B_R_buf_mm_[IND_I(x,y)].y].semi*Lmb_Y;  //(!G_R_buf_y[(y-1)*m_width + x]) ? this->C_gr_1*Lmbd_y : Lmbd_y;
  REAL max =    Smthnss[(int)B_R_buf_mm_[IND_I(x,y)].y].max*Lmb_Y;
 if(d_mx < m_nLabels)
 UpdMsgL1_(
 &XRD(x),
 &Msg(y,3),Di,
 &Msg_add(xinv,y,3),d_mx, lmbd, max);
 else
 UpdMsgL1(
 &XRD(x),
 &Msg(y,3),Di,d_mx, lmbd,max);	
}
}

inline void  EDP::UpdMsgYR_pp(int x, int y,EDP::REAL* Di,EDP::REAL* Msg,EDP::REAL* Msg_add, EDP::REAL* Mrg)
{



if (y < m_height-1)
{
 int xinv = m_width-1 -x;
 int d_mx = (xinv + 1 <m_nLabels) ? xinv +1 : m_nLabels;
 subMsg(
 &XRD(x),
 &XRD(x),Di,&Msg(y+1,3),
 &Mrg(y),
 d_mx, 0.5);
 REAL  lmbd = Smthnss[(int)B_R_buf_pp_[IND_I(x,y)].y].semi*Lmb_Y;// (!G_R_buf_y[y*m_width + x]) ? this->C_gr_1*Lmbd_y : Lmbd_y;
 REAL   max  = Smthnss[(int)B_R_buf_pp_[IND_I(x,y)].y].max*Lmb_Y;
 if(d_mx < m_nLabels)
 UpdMsgL1_(
 &XRD(x),
 &Msg(y,1),Di,
 &Msg_add(xinv,y,1),d_mx, lmbd, max);
 else
 UpdMsgL1(
 &XRD(x),
 &Msg(y,1),Di,d_mx, lmbd, max);	
}
}
inline void  EDP::UpdMsgXR_mm(int x, int y,EDP::REAL* Di,EDP::REAL* Msg,EDP::REAL* Msg_add, EDP::REAL* Mrg)
{

 if (x >0)
{
 int xinv = m_width-1 -x;
 int d_mx = (xinv + 1 <m_nLabels) ? xinv +1 : m_nLabels;
 subMsg(
 &XRD(x),
 &XRD(x-1),Di,&Msg(y,0),
 &Mrg(y),
 d_mx, 0.5);
 REAL  lmbd = Smthnss[(int)B_R_buf_mm_[IND_I(x,y)].x].semi;//(!G_R_buf_x[y*m_width + x-1]) ? this->C_gr_1*Lmbd_y : Lmbd_y;
 REAL  max = Smthnss[(int)B_R_buf_mm_[IND_I(x,y)].x].max;
 if(d_mx < m_nLabels)
 UpdMsgL1_(
 &XRD(x),
 &Msg(y,2),Di,
 &Msg_add(xinv,y,2),d_mx, lmbd, max);
 else
 UpdMsgL1(
 &XRD(x),
 &Msg(y,2),Di,d_mx, lmbd, max);	
}
//----------------------------------------

}


inline void  EDP::UpdMsgXR_pp(int x, int y,EDP::REAL* Di,EDP::REAL* Msg,EDP::REAL* Msg_add, EDP::REAL* Mrg)
{

 if (x < m_width-1)
{
  int xinv = m_width-1 -x;
 int d_mx = (xinv + 1 <m_nLabels) ? xinv +1 : m_nLabels;
 if(xinv < m_nLabels)
 subMsg_(
 &XRD(x),
 &XRD(x+1),Di,&Msg(y,2),
 &Mrg(y),
 d_mx-1, Msg_add(xinv-1,y,2),
 0.5);
 else
 subMsg(
 &XRD(x),
 &XRD(x+1),Di,&Msg(y,2),
 &Mrg(y),
 d_mx, 0.5);
 REAL  lmbd = Smthnss[(int)B_R_buf_pp_[IND_I(x,y)].x].semi;//(!G_R_buf_x[y*m_width + x]) ? this->C_gr_1*Lmbd_y : Lmbd_y;
 REAL  max = Smthnss[(int) B_R_buf_pp_[IND_I(x,y)].x].max;
 if(d_mx < m_nLabels)
 UpdMsgL1_(
 &XRD(x),
 &Msg(y,0),Di,
 &Msg_add(xinv,y,0),d_mx, lmbd,max);
 else
 UpdMsgL1(
 &XRD(x),
 &Msg(y,0),Di,d_mx, lmbd, max);	
}
//----------------------------------------

}

inline void  EDP::MkMrgR(int x, int y,EDP::REAL* Msg,EDP::REAL* Msg_add, EDP::REAL* Mrg)
{
 int xinv = m_width-1 -x;
 int d_mx = (xinv + 1 <m_nLabels) ? xinv +1 : m_nLabels;
 if (x > 0)
 {

 addMsg(
 &XRD(x),
 &XRD(x-1),&Msg(y,0),
 &Mrg(y),
 d_mx);
 }// message (x-1,y)->(x,y)
if (y > 0)
{
 addMsg(
 &XRD(x),
 &XRD(x),&Msg(y-1,1),
 &Mrg(y),
 d_mx);
}// message (x,y-1)->(x,y)
if (x < m_width-1)
{
 if(xinv < m_nLabels)
 addMsg_(
 &XRD(x),
 &XRD(x+1),&Msg(y,2),
 &Mrg(y),
 d_mx-1,
 Msg_add(xinv-1,y,2));
 else
 addMsg(
 &XRD(x),
 &XRD(x +1),&Msg(y,2),
 &Mrg(y),
 d_mx);
}// message (x+1,y)->(x,y)
if (y < m_height-1)
{
 addMsg(
 &XRD(x),
 &XRD(x), &Msg(y+1,3),
 &Mrg(y),
 d_mx);
}// message (x,y+1)->(x,y)


}

void   EDP::Gss_wnd_pp( int r_w, float cl_prc, int * DspM, int iLR, int o_s )
{
	unsigned char rgb[3];
	float rgb_m[3] ={0,0,0};
	for(int i= 0; i< m_nPixels;i ++)
	for(int c= 0; c< 3;c ++)
	{	rgb_m[c] += I_ims[iLR][i+m_nPixels*c];}
	for(int c= 0; c< 3;c ++)
	{	rgb_m[c] /= m_nPixels;}
	float vll, std =0;
	for(int i= 0; i< m_nPixels;i ++)
	{ vll =0;
	for(int c= 0; c< 3;c ++)
	{ float vl = (I_ims[iLR][i+m_nPixels*c]-rgb_m[c]); vl *=vl; vll +=vl;}
	std += vll;
	}
	int cl_fc = round_fl(sqrt(std/m_nPixels)*cl_prc);
	///////////////////////////////////////////////
	int dm = 2*r_w+1;
	int size_w = dm*dm;
	int * rez = new int [m_nPixels];
	float * Histo = new float [o_s];
	/////////////////////////////////////////
	float *gss_r_wt = new float [m_width+m_height];
	float *gss_c_wt = new float [1000];
	float sgm_r = r_w*0.5;

	for(int i=0; i<m_width+m_height; i++ )gss_r_wt[i]= exp(-(i*i)/(2.*sgm_r*sgm_r));
	for(int i=0; i<1000; i++ )gss_c_wt[i]= exp(-(i*i)/(2.*cl_fc*cl_fc));
    ///////////////////////////////////////////////
	int * map = new int [m_nPixels];
	 int max = 0;
     for(int y=0; y<m_height; y++)
	 {
      int st_rl =0; int cnt = 0; for(int x=m_width-1; x>=0; x--)
	  {
		  if(!st_rl){if(DspM[IND(x,y)] >= x){st_rl =1; cnt =1;} map[IND(x,y)] =0;}
		  else { max = (max> cnt) ? max : cnt;/* DspM[IND(x,y)] =*/ map[IND(x,y)] =cnt;  cnt++;}

	  }
	 }
    ///////////////////////////////////////////////////////
	for(int it =1; it <=max; it++){
	 for(int y=0; y<m_height; y++)
	for(int x=m_nLabels; x>=0; x--)
		if(map[IND(x,y)] == it){ int min,max;
	    Histo_min_max_pp(IND(x,y), &min, &max, DspM,  r_w,  Histo, map);
		if(min != -1){ if(min!=max)rez[IND(x,y)]=gss_wnd_pp( IND(x,y), min, max,  DspM, r_w,gss_c_wt, gss_r_wt, iLR, Histo, map);
		else rez[IND(x,y)] = min;}  else rez[IND(x,y)]  = DspM[IND(x,y)];
	  map[IND(x,y)]=0;
	 }  else rez[IND(x,y)]  = DspM[IND(x,y)];
     for(int y=0; y<m_height; y++)
	for(int x=m_nLabels; x>=0; x--)  DspM[IND(x,y)] = rez[IND(x,y)];
	}
	
///////////////////////////////////////////////////////
//for(int i= 0; i< m_nPixels;i ++)  DspM[i] = rez[i];
delete [] map;
delete [] gss_r_wt;
delete [] gss_c_wt;
delete [] rez;
delete [] Histo;


}
void   EDP:: OneD_filter(unsigned char * RI, float * cst, int r, float sig_c)
{
int rr = r*2; int w_sz = 2*rr + 1;
float * w = new float [w_sz];
float * c_cst = new float [m_nPixels*m_nLabels];
memcpy(c_cst, cst,m_nPixels*m_nLabels*sizeof(float));
float sg_c = 255.*sig_c; sg_c = 1./(sg_c*sg_c*2);
float sg_sp = 1./(r*r*2);
for(int x =0; x < m_width;  x++)
for(int y =0; y < m_height; y++)
{
    float sum_w = 0;
	BYTE cc[3] = {RI[IND(x,y)],RI[IND_IC(x,y,1)],RI[IND_IC(x,y,2)]};
	for(int ri = - rr, wi=0; ri <= rr; ri++, wi++)
	{
		int yy = y + ri; if(yy<0)yy=0; if(yy>=m_height)yy=m_height-1;
		float ds =0;
		float r2 = ri*ri;
		if(ri)for(int c =0; c <3; c++)
		{float dif = RI[IND_IC(x,yy,c)]-cc[c]; ds +=dif*dif; }
		w[wi] = exp(-ds*sg_c - r2*sg_sp); sum_w += w[wi];
	}
	for(int l =0; l <m_nLabels; l++){float sum =0;
	for(int ri = - rr, wi=0; ri <= rr; ri++, wi++)
	{
		int yy = y + ri; if(yy<0)yy=0; if(yy>=m_height)yy=m_height-1;
		
		sum += c_cst[IND_IC(x,yy,l)]*w[wi];
	}
    cst[IND_IC(x,y,l)]= sum/sum_w;
	}

}

delete [] w;
delete [] c_cst;
}
void   EDP:: OneD_filter(unsigned char * RI, double * cst, int r, float sig_c)
{
int rr = r*2; int w_sz = 2*rr + 1;
double * w = new double [w_sz];
double * c_cst = new double [m_nPixels*m_nLabels];
memcpy(c_cst, cst,m_nPixels*m_nLabels*sizeof(double));
double sg_c = 255.*sig_c; sg_c = 1./(sg_c*sg_c*2);
double sg_sp = 1./(r*r*2);
for(int x =0; x < m_width;  x++)
for(int y =0; y < m_height; y++)
{
    float sum_w = 0;
	BYTE cc[3] = {RI[IND(x,y)],RI[IND_IC(x,y,1)],RI[IND_IC(x,y,2)]};
	for(int ri = - rr, wi=0; ri <= rr; ri++, wi++)
	{
		int yy = y + ri; if(yy<0)yy=0; if(yy>=m_height)yy=m_height-1;
		double ds =0;
		double r2 = ri*ri;
		if(ri)for(int c =0; c <3; c++)
		{double dif = RI[IND_IC(x,yy,c)]-cc[c]; ds +=dif*dif; }
		w[wi] = exp(-ds*sg_c - r2*sg_sp); sum_w += w[wi];
	}
	for(int l =0; l <m_nLabels; l++){double sum =0;
	for(int ri = - rr, wi=0; ri <= rr; ri++, wi++)
	{
		int yy = y + ri; if(yy<0)yy=0; if(yy>=m_height)yy=m_height-1;
		
		sum += c_cst[IND_IC(x,yy,l)]*w[wi];
	}
    cst[IND_IC(x,y,l)]= sum/sum_w;
	}

}

delete [] w;
delete [] c_cst;
}
void   EDP:: OneD_filterX(unsigned char * RI, float * cst, int r, float sig_c)
{
int rr = r*2; int w_sz = 2*rr + 1;
float * w = new float [w_sz];
float * c_cst = new float [m_nPixels*m_nLabels];
memcpy(c_cst, cst,m_nPixels*m_nLabels*sizeof(float));
float sg_c = 255.*sig_c; sg_c = 1./(sg_c*sg_c*2);
float sg_sp = 1./(r*r*2);
for(int x =0; x < m_width;  x++)
for(int y =0; y < m_height; y++)
{
    float sum_w = 0;
	BYTE cc[3] = {RI[IND(x,y)],RI[IND_IC(x,y,1)],RI[IND_IC(x,y,2)]};
	for(int ri = - rr, wi=0; ri <= rr; ri++, wi++)
	{
		int xx = x + ri; if(xx<0)xx=0; if(xx>=m_width)xx=m_width-1;
		float ds =0;
		float r2 = ri*ri;
		if(ri)for(int c =0; c <3; c++)
		{float dif = RI[IND_IC(xx,y,c)]-cc[c]; ds +=dif*dif; }
		w[wi] = exp(-ds*sg_c - r2*sg_sp); sum_w += w[wi];
	}
	for(int l =0; l <m_nLabels; l++){float sum =0;
	for(int ri = - rr, wi=0; ri <= rr; ri++, wi++)
	{
		int xx = x + ri; if(xx<0)xx=0; if(xx>=m_width)xx=m_width-1;
		
		sum += c_cst[IND_IC(xx,y,l)]*w[wi];
	}
    cst[IND_IC(x,y,l)]= sum/sum_w;
	}

}

delete [] w;
delete [] c_cst;
}
void   EDP::get_std( double cl_r, unsigned char * buf,  double * gss_c_wt, int n_c )
{
	unsigned char rgb[3];
	double rgb_m[3] ={0,0,0};

	double cl_fc = cl_r*255;//round_fl(sqrt(std/m_nPixels)*cl_r);
	///////////////////////////////////////////////
	for(int i=0; i<n_c; i++ )gss_c_wt[i]= exp(-(i*i)/(2.*cl_fc*cl_fc));
    ///////////////////////////////////////////////
}

void   EDP::get_std( float cl_r, unsigned char * buf,  float * gss_c_wt, int n_c )
{
	unsigned char rgb[3];
	float rgb_m[3] ={0,0,0};
	//for(int i= 0; i< m_nPixels;i ++)
	//for(int c= 0; c< 3;c ++)
	//{	rgb_m[c] += buf[i+m_nPixels*c];}
	//for(int c= 0; c< 3;c ++)
	//{	rgb_m[c] /= m_nPixels;}
	//float vll, std =0;
	//for(int i= 0; i< m_nPixels;i ++)
	//{ vll =0;
	//for(int c= 0; c< 3;c ++)
	//{ float vl = (buf[i+m_nPixels*c]-rgb_m[c]); vl *=vl; vll +=vl;}
	//std += vll;
	//}
	float cl_fc = cl_r*255;//round_fl(sqrt(std/m_nPixels)*cl_r);
	///////////////////////////////////////////////
	for(int i=0; i<n_c; i++ )gss_c_wt[i]= exp(-(i*i)/(2.*cl_fc*cl_fc));
    ///////////////////////////////////////////////
}
void EDP::Gss_wnd_ML( int r_w, float cl_fc, int * sup,  int * DspM, int o_s, int rg ) {
	
	float * cst = new float [N_PX]; float sigcst = 16;
	FOR_PX_p {
		/*I_ims[2][p + 2*N_PX] = I_ims[2][p + N_PX] = I_ims[2][p] = */
		cst[p] = (exp(-cst_sub_pix(  (float)(DspM[p])/rg, p, I_ims[0], I_ims[1])/sigcst));
	}
	///////////////////////////////////////////////
	int dm = 2*r_w+1;
	int size_w = dm*dm;
	int * rez = new int [m_nPixels];
	float * Histo = new float [o_s];
	/////////////////////////////////////////
	float *gss_r_wt = new float [m_width+m_height];
	float *gss_c_wt = new float [1000];
	float sgm_r = r_w*0.5;

	for(int i=0; i<m_width+m_height; i++ )
		gss_r_wt[i]= exp(-(i*i)/(2.*sgm_r*sgm_r));
	
	for(int i=0; i<1000; i++ )
		gss_c_wt[i] =/* (i> 2*cl_fc)? 0 : */exp(-(i*i)/(2.*cl_fc*cl_fc));

	///////////////////////////////////////////////////////
	for(int i=0; i<m_nPixels; i++){ int min,max;
		Histo_min_max(i, &min, &max, DspM,  r_w,  Histo);
		if(min!=max)rez[i]=gss_wnd_ML( i, min, max,  sup, DspM,  r_w, gss_c_wt, gss_r_wt, cst, Histo);
		else rez[i] = min;
	}
	
	///////////////////////////////////////////////////////
	for(int i= 0; i< m_nPixels;i ++)
		DspM[i] = rez[i];
	delete [] gss_r_wt;
	delete [] gss_c_wt;
	delete [] rez;
	delete [] Histo;
	delete [] cst;

}
/*
	Gaussian Window?
	@params IN r_w
*/
void EDP::Gss_wnd_( int r_w, float cl_prc, int * DspM, int iLR, int o_s ) {
	unsigned char rgb[3];

	int cl_fc = cl_prc*255;//round_fl(sqrt(std/m_nPixels)*cl_prc);
	///////////////////////////////////////////////
	int dm = 2*r_w+1;
	int size_w = dm*dm;
	int * rez = new int [m_nPixels];
	float * Histo = new float [o_s];
	/////////////////////////////////////////
	float *gss_r_wt = new float [m_width+m_height];
	float *gss_c_wt = new float [1000];
	float sgm_r = r_w*0.5;

	for(int i=0; i<m_width+m_height; i++ )
		gss_r_wt[i]= exp(-(i*i)/(2.*sgm_r*sgm_r));
	
	for(int i=0; i<1000; i++ )
		gss_c_wt[i]= exp(-(i*i)/(2.*cl_fc*cl_fc));

    ///////////////////////////////////////////////////////
	for(int i=0; i<m_nPixels; i++){
		int min,max;
		Histo_min_max(i, &min, &max, DspM,  r_w,  Histo);
		if(min!=max) {
			rez[i]=gss_wnd_( i, min, max,  DspM, r_w,gss_c_wt, gss_r_wt, iLR, Histo);
		} else {
			rez[i] = min;
		}
	}

	///////////////////////////////////////////////////////
	for(int i= 0; i< m_nPixels;i ++) { // modify output.
		DspM[i] = rez[i];
	}
	delete [] gss_r_wt;
	delete [] gss_c_wt;
	delete [] rez;
	delete [] Histo;

}

void EDP::BF_1D( int r_w, float cl_prc,int iLR, float * fb )
{
	int cl_fc = cl_prc*255;//round_fl(sqrt(std/m_nPixels)*cl_prc);
	///////////////////////////////////////////////
	int dm = 2*r_w+1;
	int size_w = dm*dm;
	float  * rez = new float [m_nPixels];
	/////////////////////////////////////////
	float *gss_r_wt = new float [m_width+m_height];
	float *gss_c_wt = new float [1000];
	float sgm_r = r_w*0.5;

	for(int i=0; i<m_width+m_height; i++ )gss_r_wt[i]= exp(-(i*i)/(2.*sgm_r*sgm_r));
	for(int i=0; i<1000; i++ )gss_c_wt[i]= exp(-(i*i)/(2.*cl_fc*cl_fc));

    ///////////////////////////////////////////////////////
	for(int i=0; i<m_nPixels; i++) rez[i]= BF_1D(i, fb,  r_w, gss_c_wt, gss_r_wt,  iLR);
	
	///////////////////////////////////////////////////////
	for(int i= 0; i< m_nPixels;i ++)  fb[i] = rez[i];
	delete [] gss_r_wt;
	delete [] gss_c_wt;
	delete [] rez;
}

void EDP::med_fl( int r_w, int r_c,  unsigned char* inb)
{
	unsigned char * out = new unsigned char [m_nPixels*3];
	int np = (2*r_w+1) * (2*r_w+1);
	r_c *= 2*r_c;
	unsigned char * hs = new unsigned char [np*3];
	float * dsts = new float [np];
	float * w = new float [np];


	for(int p =0 ; p < m_nPixels; p ++)
	{
		int x = p%m_width;  int y  = p/m_width;
		int inp =0;
		for(int i=-r_w; i<=r_w; i++)
		for(int j=-r_w; j<=r_w; j++)
		{
			int xx=  x +i;
			if(xx<0)  xx = 0;
			if(xx>=m_width)  xx = m_width-1;
			
			int yy=  y +j;
			if(yy<0)  yy = 0;
			if(yy>=m_height) yy = m_height-1;
			
			int ind_t = xx + yy*m_width;
			for(int c = 0; c <3; c++)
				hs[inp*3 +c] = inb[ind_t + c*m_nPixels];
			inp++;
		}
		for(int j= 0; j < np; j++)
		{
			w[j] = 0;
			for(int c = 0; c <3; c++)
				w[j] +=(hs[np/2*3+c] - hs[j*3+c])*(hs[np/2*3+c] - hs[j*3+c]);
			w[j]  = exp(-w[j]/r_c);
		}
		for(int i =0; i < np; i++) {
			dsts[i] = 0;
			for(int j= 0; j < np; j++)
				if(i!=j) {
					float dst = 0;
					for(int c = 0; c <3; c++)
						dst +=(hs[i*3+c] - hs[j*3+c])*(hs[i*3+c] - hs[j*3+c])*w[j];
					dsts[i] += sqrt(dst);
				}
		}
		float min = dsts[0];
		int i_min =0;
		for(int i =1; i < np; i++)
			if(min > dsts[i]){
				min = dsts[i];
				i_min = i;
			};
		for(int c = 0; c <3; c++)
			out[x+ y*m_width + c*m_nPixels] = hs[i_min*3+c];
	}
	
	///////////////////////////////////////////////////////
	for(int i= 0; i< m_nPixels*3;i ++)
		inb[i] = out[i];
	delete [] w;
	delete [] out;
	delete [] hs;
	delete [] dsts;

}
 int EDP::Gss_wnd_Z( int r_w, float cl_prc, int * DspM, int iLR, int o_s, float div )
{
	int ret = 0;
	unsigned char rgb[3];
	float rgb_m[3] ={0,0,0};
	for(int i= 0; i< m_nPixels;i ++)
	for(int c= 0; c< 3;c ++)
	{	rgb_m[c] += I_ims[iLR][i+m_nPixels*c];}
	for(int c= 0; c< 3;c ++)
	{	rgb_m[c] /= m_nPixels;}
	float vll, std =0;
	for(int i= 0; i< m_nPixels;i ++)
	{ vll =0;
	for(int c= 0; c< 3;c ++)
	{ float vl = (I_ims[iLR][i+m_nPixels*c]-rgb_m[c]); vl *=vl; vll +=vl;}
	std += vll;
	}
	int cl_fc = cl_prc*255;//round_fl(sqrt(std/m_nPixels)*cl_prc);
	///////////////////////////////////////////////
	int dm = 2*r_w+1;
	int size_w = dm*dm;
	int * rez = new int [m_nPixels];
	float * Histo = new float [o_s];
	/////////////////////////////////////////
	float *gss_r_wt = new float [m_width+m_height];
	float *gss_c_wt = new float [1000];
	float sgm_r = r_w*0.5;

	for(int i=0; i<m_width+m_height; i++ )gss_r_wt[i]= exp(-(i*i)/(2.*sgm_r*sgm_r));
	for(int i=0; i<1000; i++ )gss_c_wt[i]= exp(-(i*i)/(2.*cl_fc*cl_fc));
     for(int i= 0; i< m_nPixels;i ++)   rez[i] = DspM[i];
	 ///////////////////////////////////////////////////////
	 float thrw =0;
for(int i=-r_w; i<=r_w; i++)
for(int j=-r_w; j<=r_w; j++)
{
	
	int rd = round_fl(sqrt((float)i*i+j*j));
	thrw += gss_c_wt[cl_fc]*gss_r_wt[rd];
}

    ///////////////////////////////////////////////////////
	for(int i=0; i<m_nPixels; i++)if(DspM[i]==o_s){ int min,max;
	    Histo_min_max_Z(i, &min, &max, DspM,  r_w,  Histo, o_s);
        if(min==max&& min ==o_s)rez[i]= o_s;
		else rez[i] = gss_wnd_Z( i, min, max,  DspM, r_w,gss_c_wt, gss_r_wt, iLR, Histo, thrw/div);
	}
	
///////////////////////////////////////////////////////
for(int i= 0; i< m_nPixels;i ++)  if((DspM[i] = rez[i])==o_s) ret++;
delete [] gss_r_wt;
delete [] gss_c_wt;
delete [] rez;
delete [] Histo;

return ret;
}

void inline EDP::Histo_min_max(int ind, int *min, int *max,  int * buf0, int r_w, float* Histo)
{
	int dm = 2*r_w+1;
	int size_w = dm*dm;
	//--------------------------------------
	int ig= ind%m_width; int jg= ind/m_width; // that's the center of the window.
	///////////////////////////////////////////////////////////////////////////
	min[0] = 100000;
	max[0] = 0;
	for(int i=-r_w; i<=r_w; i++)
	for(int j=-r_w; j<=r_w; j++)
	{	
		int x=  ig+i;
		if(x<0) x=0;
		if(x>=m_width) x = m_width-1;
		int y= jg+j;
		if(y<0) y=0;
		if(y>=m_height) y = m_height-1;
		int ind_t = x+y*m_width;
		{
			int l = buf0[ind_t];
			max[0] = (max[0]> l ) ? max[0]:l;
			min[0] = (min[0]< l) ? min[0]:l;
		}
	}
	for(int i=min[0]; i <= max[0]; i++){
		Histo[i]=0;
	}
}

void inline EDP::Histo_min_max_Z(int ind, int *min, int *max,  int * buf0, int r_w, float* Histo, int Z)
{
	int dm = 2*r_w+1;
	int size_w = dm*dm;
	//--------------------------------------
	int ig= ind%m_width; int jg= ind/m_width;
	///////////////////////////////////////////////////////////////////////////
	min[0] = 100000;
	max[0] = 0;
	int st =0;
	for(int i=-r_w; i<=r_w; i++)
	for(int j=-r_w; j<=r_w; j++)
	{ 	
		int x=  ig+i; if(x<0) x=0; if(x>=m_width) x = m_width-1;
		int y= jg+j;     if(y<0) y=0; if(y>=m_height) y = m_height-1;  int ind_t = x+y*m_width;
		int l = buf0[ind_t];
		if(l!=Z){
			max[0] = (max[0]> l ) ? max[0]:l;  min[0] = (min[0]< l) ? min[0]:l;
			st =1;
		}
		
	}
	if(!st){
		max[0] = Z; min[0] = Z;
		return;
	}
	for(int i=min[0]; i <= max[0]; i++){
		Histo[i]=0;
	}
}

void inline EDP::Histo_min_max_pp(int ind, int *min, int *max,  int * buf0, int r_w, float* Histo, int *map)
{
	int dm = 2*r_w+1;
	int size_w = dm*dm;
	//--------------------------------------
	int ig= ind%m_width; int jg= ind/m_width;


///////////////////////////////////////////////////////////////////////////
	min[0] = 100000;
	max[0] = 0;
	for(int i=0; i<=r_w; i++)
	for(int j=-r_w; j<=r_w; j++)
	{ 	
		int x=  ig+i;
		if(x<0) x=0;
		if(x>=m_width) x = m_width-1;
		int y= jg+j;
		if(y<0) y=0;
		if(y>=m_height) y = m_height-1;
		int ind_t = x+y*m_width;
		{
			int l = buf0[ind_t]; /*int st = map[ind_t];*/ /*if(!st)*/
			{
				max[0] = (max[0]> l ) ? max[0]:l;  min[0] = (min[0]< l) ? min[0]:l;
			}
		}
	}
	if(min[0] >max[0])
		min[0] = max[0] = -1;
	if(min[0]!=-1)
		for(int i=min[0]; i <= max[0]; i++){
			Histo[i]=0;
		}
}

int inline EDP::mask_BF_med(int ind,  unsigned char * msk, int * buf0, int r_w, float * gss_c_wt,  float * gss_r_wt, int iRL)
{
		float *H = new float [N_LB]; memset(H, 0, N_LB*sizeof(float));
		int dm = 2*r_w+1;
		int size_w = dm*dm;
		//--------------------------------------

		unsigned char * im = this->I_ims[iRL];
		int ig= ind%m_width; int jg= ind/m_width;


	///////////////////////////////////////////////////////////////////////////

	for(int i=-r_w; i<=r_w; i++)
	for(int j=-r_w; j<=r_w; j++)
	{
		
		int x=  ig+i; if(x<0) x=0; if(x>=m_width) x = m_width-1;
		int y= jg+j;     if(y<0) y=0; if(y>=m_height) y = m_height-1;  int ind_t = x+y*m_width;
		int rd = round_fl(sqrt((float)i*i+j*j));
						if(msk[ind_t]){
						int l = buf0[ind_t];
						int c_t = (abs(im[ind]-im[ind_t]) + abs(im[ind+ m_nPixels]-im[ind_t+ m_nPixels]) + abs(im[ind+ m_nPixels*2]-im[ind_t+ m_nPixels*2]));
						if(c_t>=1000)c_t=999;
						H[l] += gss_c_wt[c_t]*gss_r_wt[rd];
						}


	}

	/////////////////////////////////////////////////
	float sum=0;

	for(int i=0; i <N_LB; i++){ sum +=H[i];}
	if(sum == 0) { delete [] H ; return -1;}
	float thr =sum*0.5; sum=0; int med_ind = 0;
	for(int i=0; i <N_LB; i++){ sum +=H[i]; if(sum >= thr){med_ind = i; i  = N_LB;} }
	delete [] H;
	return  med_ind;
}
int inline EDP::gss_wnd_ML(int ind, int min, int max,  int * im, int * buf0, int r_w, float * gss_c_wt,float * gss_r_wt, float *cst, float* Histo)
{
	int dm = 2*r_w+1;
	int size_w = dm*dm;
	//--------------------------------------

	int ig= ind%m_width; int jg= ind/m_width;

	///////////////////////////////////////////////////////////////////////////

	for(int i=-r_w; i<=r_w; i++)
	for(int j=-r_w; j<=r_w; j++)
	{
		
		int x=  ig+i; if(x<0) x = -1; if(x>=m_width) x =  -1;
		int y= jg+j;     if(y<0) y= -1; if(y>=m_height) y = -1;  int ind_t = (x+1 && y+1)? x+y*m_width : -1;
		int rd = round_fl(sqrt((float)i*i+j*j));
		if(ind_t+1)
		{
	
			int l = buf0[ind_t];
			int c_t  =  abs(im[ind]-buf0[ind_t]);
			if(c_t>=1000)c_t=999;
			Histo[l] += gss_c_wt[c_t]*gss_r_wt[rd]*cst[ind_t];
	
		}
	}

	/////////////////////////////////////////////////
	float sum=0;
	for(int i = min; i <= max; i++){ sum +=Histo[i];}
	float thr =sum*0.5; sum=0; int med_ind =min;
	for(int i=min; i <= max; i++){ sum +=Histo[i]; if(sum>=thr){med_ind = i; i = max+1;}}
	return  med_ind;
	/////////////////////////////////////////////////
}

int inline EDP::gss_wnd_(int ind, int min, int max,  int * buf0, int r_w, float * gss_c_wt,float * gss_r_wt, int iRL, float* Histo) {
	int dm = 2*r_w+1;
	int size_w = dm*dm;
	//--------------------------------------
	unsigned char * im = this->I_ims[iRL];
	int ig= ind%m_width; int jg= ind/m_width;

	///////////////////////////////////////////////////////////////////////////
	int ii=0;
	for(int i=-r_w; i<=r_w; i++)
		for(int j=-r_w; j<=r_w; j++)
		{

			int x=  ig+i;
			if(x<0) x=0;
			if(x>=m_width) x = m_width-1;
			
			int y= jg+j;
			if(y<0) y=0;
			if(y>=m_height) y = m_height-1;

			int ind_t = x+y*m_width;
			int rd = round_fl(sqrt((float)i*i+j*j));
			int l = buf0[ind_t]; // disparity at pixel ind_t
			int c_t = (	abs(im[ind]-im[ind_t]) +
						abs(im[ind+ m_nPixels]-im[ind_t+ m_nPixels]) +
						abs(im[ind+ m_nPixels*2]-im[ind_t+ m_nPixels*2])
					); // absolute difference in colors.
			
			if(c_t>=1000)
				c_t=999;
			
			Histo[l] += gss_c_wt[c_t]*gss_r_wt[rd];
			ii++;
		}

	/////////////////////////////////////////////////
	float sum=0;
	for(int i=min; i <= max; i++){
		sum +=Histo[i];
	}
	float thr =sum*0.5;
	sum=0;
	int med_ind =min, med_indm=min;
	for(int i=min; i <= max; i++){
		sum +=Histo[i];
		if(sum>=thr) {
			med_ind = i;
			i = max+1;
		} else {
			if(Histo[i])
				med_indm=i;
		}
	}
	if(med_ind==min)
		return min;
	
	return med_ind;
	//if(sum-thr < thr + Histo[med_ind] - sum) return  med_ind;
	//else return  med_indm;
	/////////////////////////////////////////////////
}


float inline EDP::BF_1D(int ind,   float * buf0, int r_w, float * gss_c_wt,float * gss_r_wt, int iRL)
{
	int dm = 2*r_w+1;
	//--------------------------------------

	unsigned char * im = this->I_ims[iRL];
	int ig= ind%m_width; int y= ind/m_width;


///////////////////////////////////////////////////////////////////////////
float sum =0, sumw=0;
for(int i=-r_w; i<=r_w; i++)
{
	
	int x=  ig+i; if(x<0) x=0; if(x>=m_width) x = m_width-1;
    int ind_t = x+y*m_width;
	int rd = abs(i);
	float l = buf0[ind_t];
	int c_t = round_fl(sqrt((float) (im[ind]-im[ind_t])*(im[ind]-im[ind_t])
		+ (float) (im[ind + m_nPixels]-im[ind_t+ m_nPixels])*(im[ind+ m_nPixels]-im[ind_t+ m_nPixels])
		              + (float) (im[ind + m_nPixels*2]-im[ind_t+ m_nPixels*2])*(im[ind+ m_nPixels*2]-im[ind_t+ m_nPixels*2]) ));
	if(c_t>=1000)c_t=999;
	 float wt = gss_c_wt[c_t]*gss_r_wt[rd]; sum += wt*l; sumw +=wt;

}
return sum/sumw;
/////////////////////////////////////////////////
}
int inline EDP::gss_wnd_Z(int ind, int min, int max,  int * buf0, int r_w, float * gss_c_wt,float * gss_r_wt, int iRL, float* Histo, float thrw)
{
	int dm = 2*r_w+1;
	int size_w = dm*dm;
	//--------------------------------------

	unsigned char * im = this->I_ims[iRL];
	int ig= ind%m_width; int jg= ind/m_width;


///////////////////////////////////////////////////////////////////////////
int ii=0;
for(int i=-r_w; i<=r_w; i++)
for(int j=-r_w; j<=r_w; j++)
{
	
	int x=  ig+i; if(x<0) x=0; if(x>=m_width) x = m_width-1;
	int y= jg+j;     if(y<0) y=0; if(y>=m_height) y = m_height-1;  int ind_t = x+y*m_width;
	int rd = round_fl(sqrt((float)i*i+j*j));
	int l = buf0[ind_t];
	int c_t = (abs(im[ind]-im[ind_t]) + abs(im[ind+ m_nPixels]-im[ind_t+ m_nPixels]) + abs(im[ind+ m_nPixels*2]-im[ind_t+ m_nPixels*2]));
	if(c_t>=1000)c_t=999;
	if(l!=m_nLabels)Histo[l] += gss_c_wt[c_t]*gss_r_wt[rd];
    ii++;

}

/////////////////////////////////////////////////
float sum=0;
for(int i=min; i <= max; i++){ sum +=Histo[i];}
if(sum < thrw) return m_nLabels;
float thr =sum*0.5; sum=0; int med_ind =min, med_indm=min;
for(int i=min; i <= max; i++){ sum +=Histo[i]; if(sum>=thr){med_ind = i; i = max+1;} else {if(Histo[i]) med_indm=i;}}
if(med_ind==min) return min;
return  med_ind;
//if(sum-thr < thr + Histo[med_ind] - sum) return  med_ind;
//else return  med_indm;
/////////////////////////////////////////////////
}
int inline EDP::gss_wnd_pp(int ind, int min, int max,  int * buf0, int r_w, float * gss_c_wt,float * gss_r_wt, int iRL, float* Histo,  int* map)
{
	int dm = 2*r_w+1;
	int size_w = dm*dm;
	//--------------------------------------

	unsigned char * im = this->I_ims[iRL];
	int ig= ind%m_width; int jg= ind/m_width;


///////////////////////////////////////////////////////////////////////////
int ii=0;
for(int i=0; i<=r_w; i++)
for(int j=-r_w; j<=r_w; j++)
{
	
	int x=  ig+i; if(x<0) x=0; if(x>=m_width) x = m_width-1;
	int y= jg+j;     if(y<0) y=0; if(y>=m_height) y = m_height-1;  int ind_t = x+y*m_width;
	int rd = round_fl(sqrt((float)i*i+j*j));
	int l = buf0[ind_t];/* int st = map [ind_t];*/
	/*if(!st)*/{
	int c_t = (abs(im[ind]-im[ind_t]) + abs(im[ind+ m_nPixels]-im[ind_t+ m_nPixels]) + abs(im[ind+ m_nPixels*2]-im[ind_t+ m_nPixels*2]));
	if(c_t>=1000)c_t=999;
	Histo[l] += gss_c_wt[c_t]*gss_r_wt[rd];
    ii++;
	}
}

/////////////////////////////////////////////////
float sum=0;
for(int i=min; i <= max; i++){ sum +=Histo[i];}
float thr =sum*0.5; sum=0; int med_ind =min, med_indm=min;
for(int i=min; i <= max; i++){ sum +=Histo[i]; if(sum>=thr){med_ind = i; i = max+1;} else {if(Histo[i]) med_indm=i;}}
if(med_ind==min) return min;
return  med_ind;
//if(sum-thr < thr + Histo[med_ind] - sum) return  med_ind;
//else return  med_indm;
/////////////////////////////////////////////////
}
float inline EDP::gss_wnd_c(int ind,    unsigned char * im, int r_w, float * gss_c_wt,   float * gss_r_wt, float *buf_f)
{
	int dm = 2*r_w+1;
	int size_w = dm*dm;
	//--------------------------------------

	int ig= ind%m_width; int jg= ind/m_width;


///////////////////////////////////////////////////////////////////////////
double sum = 0; double sum_w =0;
for(int i=-r_w; i<=r_w; i++)
for(int j=-r_w; j<=r_w; j++)
{
	int st =1;
	int x=  ig+i; if(x<0){ x=0; st=0;} if(x>=m_width) {x = m_width-1; st =0;}
	int y= jg+j;     if(y<0){ y=0; st=0; }if(y>=m_height) {y = m_height-1; st =0;}  int ind_t = x+y*m_width;
	int rd = round_fl(sqrt((float)i*i+j*j));
	
	if(st && rd <= r_w){
	float sum_c =0; for(int c = 0; c <3; c++)sum_c += (im[ind+c*m_nPixels]-im[ind_t+c*m_nPixels])*(im[ind+c*m_nPixels]-im[ind_t+c*m_nPixels]);
	int c_t = round_fl(sqrt((float)sum_c));
	if(c_t>=1000)c_t=999;
	double w =  gss_c_wt[c_t]*gss_r_wt[rd];
    sum += w*buf_f[ind_t]; sum_w += w;
	}
}


return  sum/sum_w;
}
void inline EDP::BiF_wnd( float *ret, int ind, int dir,   unsigned char * im, int r_w, float * gss_c_wt, float * gss_r_wt)
{
	int size_w = (r_w+1)*(2*r_w+1);
	double sum[3] = {0,0,0};
	double sum_wt =0;
	//--------------------------------------

	int ig= ind%m_width; int jg= ind/m_width;
	int st_x, fn_x, st_y, fn_y;
	if(dir == 0){st_x = -r_w; fn_x = 0; st_y = -r_w; fn_y = r_w;}
	if(dir == 1){st_x = -r_w; fn_x = r_w; st_y = -r_w; fn_y = 0;}
	if(dir == 2){st_x =  0; fn_x = r_w; st_y = -r_w; fn_y = r_w;}
	if(dir == 3){st_x = -r_w; fn_x = r_w; st_y = 0 ; fn_y = r_w;}

///////////////////////////////////////////////////////////////////////////
for(int i=st_x; i<=fn_x; i++)
for(int j=st_y; j<=fn_y; j++)
{
	
	int x= ig+i;     if(x<0) x=0; if(x>=m_width) x = m_width-1;
	int y= jg+j;     if(y<0) y=0; if(y>=m_height) y = m_height-1;
	int ind_t = x+y*m_width;
	int rd = round_fl(sqrt((float)i*i+j*j));
	int c_t =0;
	for(int c =0; c<3; c++) c_t += abs(im[ind + c*m_nPixels]-im[ind_t + c*m_nPixels]);
	for(int c =0; c<3; c++)  sum[c]   += (double) gss_c_wt[c_t]*gss_r_wt[rd]*im[ind_t + c*m_nPixels];
	sum_wt += (double)gss_c_wt[c_t]*gss_r_wt[rd];
}
for(int c =0; c<3; c++) ret[c*m_nPixels]  = sum[c]/sum_wt;
}
float inline EDP::BiF_wnd_( int ind, int dir,   unsigned char * im, int r_w, float * gss_c_wt, float * gss_r_wt)
{
	 float ret; int x0, y0, ind0;
	int size_w = (r_w+1)*(2*r_w+1);
	double sum[3] = {0,0,0};
	double sum_wt =0;
	double sum_g =0;
	//--------------------------------------

	int ig= ind%m_width; int jg= ind/m_width;
	int st_x, fn_x, st_y, fn_y;
	if(dir == 0){st_x = -r_w; fn_x = 0; st_y = -r_w; fn_y = r_w;   x0 = TR_XX(ig);   y0 = jg;}
	if(dir == 1){st_x = -r_w; fn_x = r_w; st_y = -r_w; fn_y = 0;   x0 = ig; y0 =TR_YY(jg);}
	if(dir == 2){st_x =  0; fn_x = r_w; st_y = -r_w; fn_y = r_w;   x0 = TR_XY0(ig);   y0 = jg;}
	if(dir == 3){st_x = -r_w; fn_x = r_w; st_y = 0; fn_y = r_w;    x0 = ig;   y0 = TR_XY0(jg);}
    ind0 = IND_I(x0, y0);
///////////////////////////////////////////////////////////////////////////

for(int i=st_x; i<=fn_x; i++)
for(int j=st_y; j<=fn_y; j++)
{
	
	int x= ig+i;     if(x<0) x=0; if(x>=m_width) x = m_width-1;
	int y= jg+j;     if(y<0) y=0; if(y>=m_height) y = m_height-1;
	int ind_t = x+y*m_width;
	int rd = round_fl(sqrt((float)i*i+j*j));
	int c_t =0;
	for(int c =0; c<3; c++) c_t += abs(im[ind0 + c*m_nPixels]-im[ind_t + c*m_nPixels]);
	sum_wt += (double)gss_c_wt[c_t]*gss_r_wt[rd];
	sum_g    += gss_r_wt[rd];
}
return sum_wt/sum_g;
}
void  inline EDP::color_med(int ind, int r_w, unsigned char * out_vl, unsigned char * im )
{
	int dm = 2*r_w+1;
	int size_w = dm*dm;
	//--------------------------------------

	float * tmp =new float [3]; for(int c =0; c <3; c++)tmp[c]=0;
	float * sumw =new float [3]; for(int c =0; c <3; c++)sumw[c]=0;
	int ig= ind%m_width; int jg= ind/m_width;


///////////////////////////////////////////////////////////////////////////
int ii=0;
for(int i=-r_w; i<=r_w; i++)
for(int j=-r_w; j<=r_w; j++)
{
	float rdw = exp(-(float)(i*i+j*j)/2);

	
	int x= ig+i; int y= jg+j; int ind_t = ((x>=0)&&(x<m_width)&&(y>=0)&&(y<m_height))?  x+y*m_width: -1;	
if((ind_t+1))	
{
    float w =rdw/(1+ im[ind_t+3*m_nPixels]);
	for(int c =0; c <3; c++){tmp[c] += w*im[ind_t+c*m_nPixels]; sumw[c] += w;}

}

}

/////////////////////////////////////////////////
//int min, min_i;
//for(int i = 0; i<ii; i++)
//{ int sum =0;
//for(int j = 0; j<ii; j++)for(int c= 0; c<3; c++)
//{sum += fabs(tmp[i+c*size_w] -tmp[j+c*size_w])/(1+tmp[j+3*size_w]);}
//if(!i){min = sum; min_i =0;}
//else {if(sum<min){min = sum; min_i =i;}}
//}
for(int c= 0; c<3; c++) out_vl[c*m_nPixels] = (tmp[c]/sumw[c]);

delete [] tmp;
delete [] sumw;
}
void  inline EDP::color_med_( int r_w,unsigned char * im )
{

	unsigned char * tmp =new unsigned char [m_nPixels*4]; memcpy(  tmp, im, m_nPixels*4);
	int cnt =1; /*while (cnt){  cnt =0;*/
	for(int i = 0; i<m_nPixels; i++)if(im[i + 3*m_nPixels]>120){cnt++; color_med(i,  r_w, &tmp[i], im );
	/*memcpy( im, tmp, m_nPixels*4);*/
	/*}*/
   }

    memcpy( im, tmp, m_nPixels*4);
delete [] tmp;
}
///////////////////////////////////////////////
void EDP::SubpixelAcc(int itr, float strt_dsp, float pw2mul, int depth, int iRL)
	{ 	

       subpix_cnst = 1024;
        for(int i =0; i<4; i++ ) fine_BC[i] = new float [subpix_cnst];
		
		for(int i=0; i < subpix_cnst; i++)
		{
        float x =(float)i/subpix_cnst;
         fine_BC[0][i] =  0.5*(x-1)*(x-2)*(x+1); //f(0)
		 fine_BC[1][i] = -0.5*(x)*(x-2)*(x+1); // f(1)		
		 fine_BC[3][i] = -1./6.*(x)*(x-2)*(x-1); // f(-1)
		 fine_BC[2][i] =  1./6.*(x)*(x+1)*(x-1); //f(2)
		}
	
    ///////////////////////////////////////////////////////////////////
	 for(int i=0; i<depth;i++) get_dsp_subpix( 1,  strt_dsp/pow(2,pw2mul*i), itr, iRL);
	///////////////////////////////////////////////////////////
 for(int i =0; i<4; i++ ) delete [] fine_BC[i];

}
void EDP::get_dsp_subpix(int r, float st, int itr, int iLR)
	{
		
		Lbs_S = (2*r+1);
		spx_D =    new float [m_nPixels* Lbs_S];
		//--------------------------------
		cost_subpix_(r, st, iLR);
	    get_dsp_sbpx_mrf( itr*4);

	//////////////////////////////
	delete [] spx_D;
}

void EDP:: cost_subpix(int r, float stp, int iRL)
{

int dm = r*2+1;
Lbs_S = dm;
Lbs_st = stp;
Lbs_S_map = new float [Lbs_S];
for(int i=0; i<Lbs_S; i++) Lbs_S_map[i] =(float)(i-r)*stp;
float* spx_D_c = new float [m_nPixels*Lbs_S];
//-------------------------------------------------------------------------------------------------
for(int i=0; i< m_nPixels; i++)
//if(!Mask_RL[iRL][i])
{if(!iRL){cst_subpix_pix(  &spx_D[i*Lbs_S], i,  I_ims[0], I_ims[1]);/* cst_subpix_pixG(  &spx_D_c[i*Lbs_S], i,  IG_ims[0], IG_ims[1]);*/}
else { cst_subpix_pix(  &spx_D[i*Lbs_S], i,  I_ims[1], I_ims[0]); /*cst_subpix_pixG(  &spx_D_c[i*Lbs_S], i,  IG_ims[1], IG_ims[0]); */}  }
//else {for(int l =0;l<Lbs_S; l++)spx_D[i*Lbs_S+l]=0;}
//for(int i=0; i< m_nPixels*dm; i++) spx_D[i] = (spx_D[i]+spx_D_c[i]/2< 3*chTr) ? spx_D[i]+spx_D_c[i]/2: 3*chTr;
Cost_mean = 0; for(int i=0; i< m_nPixels*dm; i++) Cost_mean += spx_D[i]; Cost_mean /=  m_nPixels*dm;
float dv=(Lbs_st<1)? 1. +log(1./Lbs_st):  1. +log(Lbs_st);
Cost_mean *=  subpix_cst_tune/sqrt(dv);
delete [] spx_D_c;
}
void EDP:: cost_subpix_(int r, float stp, int iRL)
{

int dm = r*2+1;
Lbs_S = dm;
Lbs_st = stp;
Lbs_S_map = new float [Lbs_S];
for(int i=0; i<Lbs_S; i++) Lbs_S_map[i] =(float)(i-r)*stp;
//-------------------------------------------------------------------------------------------------
FOR_PX_p
{
int pl = Lbs_S*p;
//for(int l = 0; l < Lbs_S; l++)spx_D[pl + l] = 1 . - exp(-cst_sub_pix(  DspFL[p] + Lbs_S_map[l], p, I_ims[0], I_ims[1])/16);

}

}
void  inline EDP:: cst_subpix_pix(  float *out_buf, int ind0, EDP::BYTE * buf0, EDP::BYTE * buf1)
{
int x = ind0%m_width; int y = ind0/m_width;
 int indp[4];

float  dsp0 = DspFL[ind0];
//-----------------------------------------------------
for(int ii=0; ii<Lbs_S; ii++)
{
float x1= dsp0 + Lbs_S_map[ii] + x ;
if((x1 >= 0)&&(x1 < m_width))
{
int x1_0     =   (int)x1;
float m_x1= x1-x1_0;
int x1_1    = ( m_x1==0) ? x1_0 : x1_0 +1;  x1_1    =  (  x1_1 <m_width) ? x1_1 : m_width-1;
//------------------------------------------------------------------------------------------------
if((x1_1==x1_0))                     { int ind1 = y*m_width+x1_0;
    out_buf [ii] = 0;
	int ind0_ = ind0;
    for(int ic=0; ic<3; ic++, ind0_ += m_nPixels, ind1 += m_nPixels) out_buf [ii] +=abs(buf0[ind0_]-buf1[ind1]);
	out_buf [ii] = ( out_buf [ii] < 3*chTr ) ? out_buf [ii] : 3*chTr;}
else                                            {
indp[0] = x1_0 + y*m_width;
indp[1] = (x1_0+1<m_width) ? x1_0+1 + y*m_width: m_width-1  + y*m_width;
indp[2] = (x1_0+2<m_width) ? x1_0+2 + y*m_width: m_width-1  + y*m_width;
indp[3] = (x1_0>0) ? x1_0-1 + y*m_width:  y*m_width;
int i_m_x1 = (int)(m_x1*subpix_cnst); if (i_m_x1>=subpix_cnst) i_m_x1 = subpix_cnst-1;
	out_buf [ii] = 0;
	if(ST_Cub){
	for(int ic=0, is=0; ic<3; ic++, is+=m_nPixels)
	{ float dst =0; for(int iw=0; iw<4; iw++)
	dst +=  fine_BC[iw][i_m_x1]*buf1[indp[iw]+is];
	out_buf [ii]  += fabs(dst -buf0[ind0+is]);}
	out_buf [ii] =   ( out_buf [ii] < 3*chTr ) ? out_buf [ii] : 3*chTr;
	}
else
{	for(int ic=0, is=0; ic<3; ic++, is+=m_nPixels)
	{ float dst =0;
	dst =  ((1.-m_x1)*buf1[indp[0]+is] + (m_x1)*buf1[indp[1]+is]);
	out_buf [ii]  += fabs(dst -buf0[ind0+is]);}
	out_buf [ii] = ( out_buf [ii] < 3*chTr ) ? out_buf [ii] : 3*chTr;
	}
}
}
else 3*chTr;
}
}

float  inline EDP::cst_sub_pix(float dsp0, int ind0, EDP::BYTE * buf0, EDP::BYTE * buf1) {
int x = ind0%m_width; int y = ind0/m_width; float ret = 0; int indp[4];
//-----------------------------------------------------
float x1= (float)x - dsp0;
if((x1 >= 0) && (x1 < m_width))
{
			int x1_0     =   (int)x1;
			float m_x1= x1-x1_0;
			int x1_1    = ( m_x1==0) ? x1_0 : x1_0 +1;  x1_1    =  (  x1_1 <m_width) ? x1_1 : m_width-1;
			indp[0] = x1_0 + y*m_width;
			indp[1] = (x1_0+1<m_width) ? x1_0+1 + y*m_width: m_width-1  + y*m_width;
			indp[2] = (x1_0+2<m_width) ? x1_0+2 + y*m_width: m_width-1  + y*m_width;
			indp[3] = (x1_0>0) ? x1_0-1 + y*m_width:  y*m_width;
			int i_m_x1 = (int)(m_x1*subpix_cnst); if (i_m_x1>=subpix_cnst) i_m_x1 = subpix_cnst-1;
			for(int ic=0, is=0; ic<3; ic++, is+=m_nPixels)
				{
				float dst =0; float b0 = buf0[ind0+is];  for(int iw=0; iw<4; iw++) dst +=  fine_BC[iw][i_m_x1]*buf1[indp[iw]+is];
				ret  += fabs(dst -b0);
               }
}

return ret;
}

void  inline EDP::cst_subpix_pixG(  float *out_buf, int ind0, EDP::BYTE * buf0, EDP::BYTE * buf1) {
	int x = ind0%m_width; int y = ind0/m_width;
	int indp[4];

	float  dsp0 = DspFL[ind0];
	//-----------------------------------------------------
	for(int ii=0; ii<Lbs_S; ii++) {
		float x1= dsp0 + Lbs_S_map[ii] + x ;
		if((x1 >= 0)&&(x1 < m_width)) {
			int x1_0     =   (int)x1;
			float m_x1= x1-x1_0;
			int x1_1    = ( m_x1==0) ? x1_0 : x1_0 +1;  x1_1    =  (  x1_1 <m_width) ? x1_1 : m_width-1;
			//------------------------------------------------------------------------------------------------
			if((x1_1==x1_0)) { 
				int ind1 = y*m_width+x1_0;
				out_buf [ii] = 0;
				int ind0_ = ind0;
				for(int ic=0; ic<6; ic++, ind0_ += m_nPixels, ind1 += m_nPixels) 
					out_buf [ii] +=abs(buf0[ind0_]-buf1[ind1]);
				out_buf [ii] = ( out_buf [ii] < 3*chTr ) ? out_buf [ii] : 3*chTr;
			}else{
				indp[0] = x1_0 + y*m_width;
				indp[1] = (x1_0+1<m_width) ? x1_0+1 + y*m_width: m_width-1  + y*m_width;
				indp[2] = (x1_0+2<m_width) ? x1_0+2 + y*m_width: m_width-1  + y*m_width;
				indp[3] = (x1_0>0) ? x1_0-1 + y*m_width:  y*m_width;
				int i_m_x1 = (int)(m_x1*subpix_cnst); 
				if (i_m_x1>=subpix_cnst) 
					i_m_x1 = subpix_cnst-1;
				out_buf [ii] = 0;
				if(ST_Cub){
					for(int ic=0, is=0; ic<6; ic++, is+=m_nPixels) { 
						float dst =0; 
						for(int iw=0; iw<4; iw++)
							dst +=  fine_BC[iw][i_m_x1]*buf1[indp[iw]+is];
						out_buf [ii]  += fabs(dst -buf0[ind0+is]);
					}
					out_buf [ii] = ( out_buf [ii] < 3*chTr ) ? out_buf [ii] : 3*chTr;
				} else { 
					for(int ic=0, is=0; ic<6; ic++, is+=m_nPixels) { 
						float dst =0;
						dst =  ((1.-m_x1)*buf1[indp[0]+is] + (m_x1)*buf1[indp[1]+is]);
						out_buf [ii]  += fabs(dst -buf0[ind0+is]);
					}
					out_buf [ii] = ( out_buf [ii] < 3*chTr ) ? out_buf [ii] : 3*chTr;
					}
			}
		} else {
			out_buf[ii]=3*chTr;
		}
	}
}

void EDP::get_dsp_sbpx_mrf(  int iti)
{
	REAL * ms[4];
	for(int i=0;i<4;i++) {
		ms[i]=new REAL [m_nPixels*Lbs_S];  
		memset( ms[i], 0, m_nPixels*Lbs_S*sizeof(REAL));
	}
	//-------------------------------------------------------------------
	dsp_sp_mrf( iti, ms);
	//-------------------------------------------------------------------
	for(int i=0;i<4;i++) 
		delete [] ms[i];
								
}
void EDP::dsp_sp_mrf(int itr,  EDP::REAL **ms)
{
	
int    ind_img_dsi[4];										int    ind_img_dsp[4];
int  i_dir_0[2]={0, m_width-1};					int  i_dir_inc[2]={1,-1};
int  j_dir_0[2]={0, m_height-1};					int  j_dir_inc[2]={1,-1};

REAL* tmp_vl= new REAL [Lbs_S];   REAL *tmp_vl_c= new REAL[Lbs_S];
/////////////////////////////////////////////////////////////
int pp=m_width*Lbs_S;
//-----------------------------
	POINT dir_for_4[4];
	dir_for_4[0].x=0;  dir_for_4[0].y=0;
	dir_for_4[1].x=0;  dir_for_4[1].y=1;
	dir_for_4[2].x=1;  dir_for_4[2].y=0;
	dir_for_4[3].x=1;  dir_for_4[3].y=1;
//-------------------------------
	int i_dir_ind=0;   int j_dir_ind=0;
/////////////////////////////////////////////////////////////

for(int iti=0;iti<itr;iti++)
{//==========ITERATION===============

	int iti_d4=iti/4; int iti_m4=iti%4;  int iti_m4_i=(4-iti%4)%4;
	{ i_dir_ind=dir_for_4[iti_m4].x;  j_dir_ind=dir_for_4[iti_m4].y;}
//---------------------------
	int ind_img=0;
    int ind_dsi0=0;
	int add_cnst_x = i_dir_inc[i_dir_ind];
	int start_cnst_x = i_dir_0[i_dir_ind];
	int add_cnst_y =  j_dir_inc[j_dir_ind];
	int start_cnst_y = j_dir_0[j_dir_ind];

//////////IMAGE DO /////////////////////////////
for(int j0=0, j=start_cnst_y; j0<m_height; j0++, j+=add_cnst_y)
for(int i0=0, i=start_cnst_x; i0<m_width; i0++, i+=add_cnst_x )
{  //FOR X Y 0000000000000000000000000000000000000
     ind_img=(j*m_width+i);     ind_dsi0=ind_img*Lbs_S;
	////-----------------------------------------------
	ind_img_dsi[2]=(i>0)? ind_dsi0-Lbs_S:-1;//X----
    ind_img_dsi[0]= (i<m_width-1)? ind_dsi0+Lbs_S:-1;//X+++++
    ind_img_dsi[3]= (j>0)? ind_dsi0-pp:-1;
    ind_img_dsi[1]= (j<m_height-1)? ind_dsi0+pp:-1;
								//----------------ITR SUM--------------------
								for(int ixy=0; ixy<2; ixy++)
								{  int ind_edg; int st_do=1;
								if(ixy){
								if(!i_dir_ind){ind_edg=0; }
								else{ ind_edg=2;} }
								else {if(!j_dir_ind){ind_edg=1; }
								else{ ind_edg=3; }}		  	
								int ind_edg_inv=(ind_edg+2)%4;
								if((ind_img_dsi[ind_edg]+1))
								{//////start

								 LOAD_D(ind_dsi0, Lbs_S); 							
								//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
								if(st_do--)
								{for(int ind_nbr=0; ind_nbr<4; ind_nbr++)
								{ int ind_edg_nbr=(ind_nbr+2)%4;
								if(ind_img_dsi[ind_nbr]+1){ ADD_MS(ind_edg_nbr,ind_nbr, Lbs_S);}}
								DIV_MRG(Lbs_S); memcpy(tmp_vl_c, tmp_vl, sizeof(REAL)*Lbs_S);}
								else {memcpy(tmp_vl, tmp_vl_c, sizeof(REAL)*Lbs_S);}
								//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
							SUB_MS(ind_edg_inv,ind_edg, Lbs_S);						
							float dsp_fs=DspFL[ind_img_dsi[ind_edg]/Lbs_S];
	                        float dsp_st=DspFL[ind_img];
						    MSG_UP(dsp_st, dsp_fs, tmp_vl, &ms[ind_edg][ind_dsi0],Cost_mean);
						}///end
	
						}//---------------ITR SUM--------------
} //FOR X Y END
}//==========ITERATION end===============

  Add_DspFL(tmp_vl, ms);

//----------------------------
delete tmp_vl, tmp_vl_c;
}
inline void EDP::MSG_UP_(float dsp_st, float dsp_fs, EDP::REAL * outb, EDP::REAL * outb_c, EDP::REAL mul)
{
REAL vl;
mul /= Lbs_st;
float st_pp, fs_pp;
REAL min_out = 1000000000;
for( int d_fs=0; d_fs<Lbs_S;d_fs++)
{
REAL min = 1000000000;
fs_pp = dsp_fs + Lbs_S_map[d_fs];
 for( int d_st=0; d_st<Lbs_S;d_st++)
 {
	 st_pp = dsp_st + Lbs_S_map[d_st];
	 float dif=  fabs(st_pp-fs_pp);
	 float add_c = (dif)*mul;  if((vl=outb[d_st]+ add_c)<min) min=vl;
 }	 	
 outb_c[d_fs]=min; min_out = (min_out<min)? min_out : min;
}
for(int i =0; i<Lbs_S; i++ )outb_c[i] -= min_out;
}


inline void EDP::MSG_UP(float dsp_st, float dsp_fs, EDP::REAL * outb, EDP::REAL * outb_c, EDP::REAL mul)
{
REAL vl;
mul /= Lbs_st;
float st_pp, fs_pp;
REAL min_out = 1000000000;
for( int d_fs=0; d_fs<Lbs_S;d_fs++)
{
REAL min = 1000000000;
fs_pp = dsp_fs + Lbs_S_map[d_fs];
 for( int d_st=0; d_st<Lbs_S;d_st++)
 {
	 st_pp = dsp_st + Lbs_S_map[d_st];
	 float dif=  fabs(st_pp-fs_pp);
	 float add_c = (dif)*mul;  if((vl=outb[d_st]+ add_c)<min) min=vl;
 }	 	
 outb_c[d_fs]=min; min_out = (min_out<min)? min_out : min;
}
for(int i =0; i<Lbs_S; i++ )outb_c[i] -= min_out;
}

void  EDP::Add_DspFL(EDP::REAL *tmp_vl, EDP::REAL ** ms)
{
int ind_img_dsi[4];
int pp=m_width*Lbs_S;
for(int i=0; i<m_width; i++)
for(int j=0; j<m_height; j++)
{
    int  ind_img=(j*m_width+i);     int  ind_dsi0=ind_img*Lbs_S;
	////-----------------------------------------------
	ind_img_dsi[2]=(i>0)? ind_dsi0-Lbs_S:-1;//X----
    ind_img_dsi[0]= (i<m_width-1)? ind_dsi0+Lbs_S:-1;//X+++++
    ind_img_dsi[3]= (j>0)? ind_dsi0-pp:-1;
    ind_img_dsi[1]= (j<m_height-1)? ind_dsi0+pp:-1;
	//----------------ITR SUM--------------------
    	 LOAD_D(ind_dsi0, Lbs_S); 		
for(int ind_nbr=0; ind_nbr<4; ind_nbr++)
{
 int ind_edg_nbr=(ind_nbr+2)%4;
 if(ind_img_dsi[ind_nbr]+1)ADD_MS(ind_edg_nbr,ind_nbr, Lbs_S);
}
	int g = /*m_answer[ind_img] = */  Get_min(tmp_vl);
	DspFL[ind_img] += this->Lbs_S_map[g];  /*  if(DspFL[ind_img] >0)DspFL[ind_img] = 0;*/

}
}
inline int EDP::Get_min(EDP::REAL *buf)
{
	REAL min=buf[0];  int dsp_min=0;
	for(int d=1; d<Lbs_S; d++)
	{if(buf[d]<min){min=buf[d];dsp_min=d;}}
     return dsp_min;
}

void EDP::GrBuf() {
	float * f_im_g_f = new float [m_nPixels*6];
	float * s_im_g_f = new float [m_nPixels*6];
	for(int i=0;i<3;i++){
		grBuf(&f_im_g_f[i*2*m_nPixels],&I_ims[0][i*m_nPixels]);
		grBuf(&s_im_g_f[i*2*m_nPixels],&I_ims[1][i*m_nPixels]);
	}

	for(int c=0;c<6;c++){
		for(int i=0;i<m_nPixels;i++) {
			f_im_g_f[i + c*m_nPixels] = (f_im_g_f[i + c*m_nPixels] <0)? -sqrt(-f_im_g_f[i + c*m_nPixels]): sqrt(f_im_g_f[i + c*m_nPixels]);
			s_im_g_f[i + c*m_nPixels] = (s_im_g_f[i + c*m_nPixels] <0)? -sqrt(-s_im_g_f[i + c*m_nPixels]): sqrt(s_im_g_f[i + c*m_nPixels]);
		}
	}
	float mul = Get_mul_gr( f_im_g_f ,I_ims[0]);
	
	for(int c=0;c<6;c++) {
		for(int i=0;i<m_nPixels;i++) {
			f_im_g_f[i + c*m_nPixels] *= mul;
			s_im_g_f[i + c*m_nPixels] *=  mul;
			IG_ims[0][i + c*m_nPixels] = (f_im_g_f[i + c*m_nPixels] + 128<0) ? 0: ( (f_im_g_f[i + c*m_nPixels] + 128>255) ? 255:f_im_g_f[i + c*m_nPixels] + 128) ;
			IG_ims[1][i + c*m_nPixels] = (s_im_g_f[i + c*m_nPixels] + 128<0) ? 0: ( (s_im_g_f[i + c*m_nPixels] + 128>255) ? 255:s_im_g_f[i + c*m_nPixels] + 128) ;
		}
	}
	/////////////////////////////////
	delete [] f_im_g_f;
	delete [] s_im_g_f;
}


void EDP::grBuf( float* buf_out, EDP::BYTE* buf_in) {
	int ImgX = m_width, ImgY = m_height,
	ImgSh = m_nPixels;
	float mul_sq2 = 1./sqrt((float)2);
	for(int i=0;i<ImgSh;i++){
		int x= i%ImgX; int y=i/ImgX;
		int xp =(x+1<ImgX) ? x+1:ImgX-1;
		int yp =(y+1<ImgY) ? y+1:ImgY-1;
		int xm =(x-1>=0) ? x-1:0;
		int ym =(y-1>=0) ? y-1:0;
		/////////////////////////////////////////////////
		float dx = (buf_in[xp+ImgX*y]-buf_in[xm+ImgX*y]);
		float dy = (buf_in[x+ImgX*yp]-buf_in[x+ImgX*ym]);
		float dx1 = mul_sq2*(buf_in[xp+ImgX*ym]-buf_in[xm+ImgX*yp]);
		float dy1 = mul_sq2*(buf_in[xp+ImgX*yp]-buf_in[xm+ImgX*ym]);

		buf_out[i]=(dx+(dx1+dy1)*mul_sq2)*4;
		buf_out[i+ImgSh]=(dy+(-dx1+dy1)*mul_sq2)*4;
	}
}

void EDP::grBuf_( float* buf_out, EDP::BYTE* buf_in) {
	int ImgX = m_width, ImgY = m_height, ImgSh = m_nPixels;
	float mul_sq2 = 1./sqrt((float)2);
	for(int i=0;i<ImgSh;i++){
		int x= i%ImgX; int y=i/ImgX;
		int xp =(x+1<ImgX) ? x+1:ImgX-1;
		int yp =(y+1<ImgY) ? y+1:ImgY-1;
		int xm =(x-1>=0) ? x-1:0;
		int ym =(y-1>=0) ? y-1:0;
		/////////////////////////////////////////////////
		float dx = (buf_in[xp+ImgX*y]-buf_in[xm+ImgX*y]);
		float dy = (buf_in[x+ImgX*yp]-buf_in[x+ImgX*ym]);
		float dx1 = mul_sq2*(buf_in[xp+ImgX*ym]-buf_in[xm+ImgX*yp]);
		float dy1 = mul_sq2*(buf_in[xp+ImgX*yp]-buf_in[xm+ImgX*ym]);
		buf_out[i]=dx*8;//(dx+(dx1+dy1)*mul_sq2)*4;
		buf_out[i+ImgSh]= dx*8;
	}//(dy+(-dx1+dy1)*mul_sq2)*4;}
}
int EDP::get_lb(int p, float * buf_lb, int N) {
	int ret = 0;
	float min = buf_lb[p];
	for(int i =1; i < N; i ++) {
		if(min > buf_lb[p + i*m_nPixels]){
			min = buf_lb[p + i*m_nPixels];
			ret =i;
		}
	}
	return ret;
}
void EDP::ROTH_lb(int n_lb, float * buf_lb, unsigned char * I_b, int * answer) {
	int r = 32;
	int div = 256/r;
	int size = r*r*r;
	double logmin =  1./((double)m_nPixels*m_nPixels);
	double  * buf2srt = new double[size];
	for(int p = 0; p < size; p++) {
		buf2srt[p] = 0;
	}
	for(int y = 0; y < m_height; y++)for(int x = 0; x < m_width; x++){
		double w = (answer[IND(x,y)]==n_lb)? 1:0;
		int vl[3] = {CUTFL((float)I_b[IND_IC(x,y,0)]/div, r ), CUTFL((float)I_b[IND_IC(x,y,1)]/div, r ),CUTFL((float)I_b[IND_IC(x,y,2)]/div, r )} ;
		buf2srt[IND_HC(vl[0],vl[1],vl[2], r)]    +=  w ;
	}
	float sum =0;
	for(int p = 0; p < size; p++) {
		sum += buf2srt[p];
	}
	buf2srt[0] = -log(buf2srt[0]/sum + logmin);
	float min = buf2srt[0];

	for(int p = 1; p < size; p++) {
		buf2srt[p] = -log(buf2srt[p]/sum + logmin);
		min = (min < buf2srt[p] ) ? min: buf2srt[p];
	}
	for(int y = 0; y < m_height; y++)for(int x = 0; x < m_width; x++) {
		int vl[3] = {CUTFL((float)I_b[IND_IC(x,y,0)]/div, r ), CUTFL((float)I_b[IND_IC(x,y,1)]/div, r ),CUTFL((float)I_b[IND_IC(x,y,2)]/div, r )} ;
		buf_lb[IND(x,y)] = buf2srt[IND_HC(vl[0],vl[1],vl[2], r)] ;
	}

	delete [] buf2srt;
}


void EDP::getCst_lft(int n_lb, float * buf_lb,  float thr)
{
	unsigned char * I_b1 = I_ims[0]; unsigned char * I_b2 = I_ims[1];
	unsigned char * G_b1 = this->IG_ims[0]; unsigned char * G_b2 = IG_ims[1];
  for(int y = 0; y < m_height; y++)for(int x = m_width-1; x >=0; x--)
		 {
			 int x2 = x - n_lb; float cst, cstg;
			 if(x2>=0) {
				 cst = 0; cstg =0;
				 for(int c = 0 ; c < 3; c++){float vl = I_b1[IND_IC(x,y,c)] - I_b2[IND_IC(x2,y,c)]; cst += vl*vl; }
				 for(int c = 0 ; c < 6; c++){float vl = G_b1[IND_IC(x,y,c)] - G_b2[IND_IC(x2,y,c)]; cstg += vl*vl; }
				 buf_lb[IND(x,y)] = (cst+cstg< thr) ? cst +cstg : thr;  }
			 else { buf_lb[IND(x,y)] = (cst + cstg< thr) ? cst  + cstg: thr; }
		 }

}



void EDP::ROTH_lb(int n_lb, float * buf_lb, unsigned char * I_b)
{

	int rule_cnst =3; int r = 32;  int div = 256/r; int size = r*r*r;  POINT p_lb0 = {n_lb%rule_cnst, n_lb/rule_cnst};
	 POINT_F p_l = {float(p_lb0.x*(m_width/rule_cnst) + m_width/rule_cnst/2), float(p_lb0.y*(m_height/rule_cnst) + m_height/rule_cnst/2)};
     int r_sig = m_width/rule_cnst/2;
	
	 double * exp_sig = new double [m_width+ m_height];
	 for(int i =0; i <m_width+ m_height; i++ ) exp_sig[i] = exp(-(double)i*i/(2*r_sig*r_sig));
	  double logmin =  1./((double)m_nPixels*m_nPixels);
	   double  logmax =  -log(logmin);
	 double max = 0;
	double  * buf2srt = new double[size];  for(int p = 0; p < size; p++)buf2srt[p] = 0;
         for(int y = 0; y < m_height; y++)for(int x = 0; x < m_width; x++)
		 {

            int  dist =(int)sqrt((float) (p_l.y-y)*(p_l.y-y) + (p_l.x-x)*(p_l.x-x)); double w = exp_sig[dist];
			 int vl[3] = {CUTFL((float)I_b[IND_IC(x,y,0)]/div, r ), CUTFL((float)I_b[IND_IC(x,y,1)]/div, r ),CUTFL((float)I_b[IND_IC(x,y,2)]/div, r )} ;
             buf2srt[IND_HC(vl[0],vl[1],vl[2], r)]    +=  w ;
		 }
		 float sum =0; for(int p = 0; p < size; p++) {sum += buf2srt[p]; max = (max> buf2srt[p]) ? max : buf2srt[p];}
		 buf2srt[0] = -log(buf2srt[0]/sum + logmin); float min = buf2srt[0];   max =0;
		 for(int p = 1; p < size; p++) {buf2srt[p] = -log(buf2srt[p]/sum + logmin);  min = (min < buf2srt[p] ) ? min: buf2srt[p]; }
		/* for(int p = 0; p < size; p++) {buf2srt[p] -= min; max = (max > buf2srt[p] ) ? max : buf2srt[p];}*/
		 sum = 0;
		          for(int y = 0; y < m_height; y++)for(int x = 0; x < m_width; x++)
		 {

			 int vl[3] = {CUTFL((float)I_b[IND_IC(x,y,0)]/div, r ), CUTFL((float)I_b[IND_IC(x,y,1)]/div, r ),CUTFL((float)I_b[IND_IC(x,y,2)]/div, r )} ;
             buf_lb[IND(x,y)] = buf2srt[IND_HC(vl[0],vl[1],vl[2], r)] ; max = (max >  buf_lb[IND(x,y)] ) ? max : buf_lb[IND(x,y)]; sum += buf_lb[IND(x,y)];
		 }
delete [] exp_sig,
delete [] buf2srt;
}


void EDP::K_means(int cl_Q, COLOR_F * K_cls, int r, int * buf_to_sort, short int * clr_cls)
{
	int size = r*r*r; int dev = 256/r;
     int * ind_out = new int [size];
	 float * dist = new float [cl_Q];
 ////////////////////////////////////////////
SetSort(size,  buf_to_sort, ind_out);

COLOR_F a[8]; /////debug
COLOR_F * k_cls[2];
k_cls[0]  = new COLOR_F [cl_Q];
k_cls[1]  = new COLOR_F [cl_Q];
for(int i =0; i < cl_Q; i ++)for (int c =0; c< 3; c++) k_cls[0][i].cl[c] = IND_HC_B(c, ind_out[i], r);
int * k_cls_i =  new int [cl_Q];
delete [] ind_out;
short int *  cl_c  = new short int [size];

int It =15;

///////////////////////////////////////////
for(int it =1; it <= It; it++){ int it0 = (it-1)%2, it1 = (it)%2;  int cc[3]; float  max;
for(cc[0] = 0; cc[0]< r; cc[0]++)for(cc[1] = 0; cc[1]< r; cc[1]++) for(cc[2] = 0; cc[2]< r; cc[2]++)
{
	float min; int  ind_min;
	for(int d = 0; d < cl_Q; d++)
	{
	dist[d] =0;
	for(int c = 0; c <3; c++){float vl =  k_cls[it0][d].cl[c]- cc[c]; dist[d] += vl*vl; }
	if(!d){min = dist[d]; ind_min =d;}
	else  {if(min > dist[d]){min = dist[d]; ind_min =d;}}
	}
    cl_c[IND_HC(cc[0],cc[1], cc[2],r)] = ind_min;
}
for(int ci =0; ci < cl_Q; ci++){k_cls_i[ci] =0; k_cls[it1][ci].cl[0]=0; k_cls[it1][ci].cl[1]=0; k_cls[it1][ci].cl[2]=0;}
for(cc[0] = 0; cc[0]< r; cc[0]++)for(cc[1] = 0; cc[1]< r; cc[1]++) for(cc[2] = 0; cc[2]< r; cc[2]++)
{

int d = cl_c[IND_HC(cc[0],cc[1], cc[2],r)];
int add_v = buf_to_sort[IND_HC(cc[0],cc[1], cc[2],r)];
k_cls_i[d] += add_v;
for(int c =0; c <3; c++)k_cls[it1][d].cl[c] += cc[c]*add_v;
}
for(int d = 0; d < cl_Q; d++) for(int c =0; c <3; c++)k_cls[it1][d].cl[c] /= k_cls_i[d];
	for(int d = 0; d < cl_Q; d++)
	{
	dist[d] =0;
	for(int c = 0; c <3; c++){float vl =  k_cls[it0][d].cl[c]- k_cls[it1][d].cl[c]; dist[d] += vl*vl; }
	if(!d) max = dist[d];
	else  {if(max < dist[d]){max = dist[d];}}
	}
	for(int ii =0; ii < 8; ii ++) a[ii] = k_cls[it1][ii];
	if(it == It || max < 1./dev)
	{
		for(int d = 0; d < cl_Q; d++) for(int c =0; c <3; c++){K_cls[d].cl[c] = round_fl(k_cls[it1][d].cl[c]*dev); if(K_cls[d].cl[c] >255) K_cls[d].cl[c] =255;}
		for(int i = 0; i < size; i++) clr_cls[i] = cl_c [i];
		it =It+1;
	}
}
///////////////////////////////////////////
delete [] k_cls_i;
delete [] dist;
delete [] k_cls[0];
delete [] k_cls[1];
delete [] cl_c;
}


void EDP::SetSort(int n_p, int * buf_to_sort, int * ind_out)
{
	int *i_exp2 = new int [31]; i_exp2[0] = 1; 
	for(int i = 1; i< 31; i++ ) i_exp2[i] = 2*i_exp2[i-1];
	
	int st =1; int N  = 0; 
	while(st){ if(n_p <= i_exp2[N]) st = 0; else N++; }
	
	int size= i_exp2[N];
	int *buf_max[2];
	int * buf_ind[2];
	
	buf_max[0] = new int  [ size];
	buf_max[1] = new int  [ size];
	buf_ind[0] =   new int  [ size];
	buf_ind[1] =   new int  [ size];
	
	for(int i = 0 ; i < size; i++){
		buf_ind[0][i] = i; 
		if(i<n_p) buf_max[0][i] = buf_to_sort[i]; 
		else buf_max[0][i] = 0;
	}
	////////////////////////////////////////////////////////////

	st=0;
	for(int in=1; in <= N ; in ++ ) {
		SetSort_n(i_exp2, in, N,  buf_max[(st+1)%2],  buf_ind[(st+1)%2],  buf_max[st],  buf_ind[st]); 
		st=(st+1)%2;  
	}
	for(int i = 0 ; i < n_p; i++){
		ind_out[i] = buf_ind[st][i]; 
	}
	///////////////////////////////////////////
	delete [] buf_max[0];
	delete [] buf_max[1];
	delete [] buf_ind[0];
	delete [] buf_ind[1];
	delete [] i_exp2;
}


inline void EDP::SetSort_n(int * i_exp2, int n,  int N, int* buf_out_max, int* buf_out_ind, int * buf_in_max, int* buf_in_ind) {
	int J_max=i_exp2[N-n];
	int I_max = i_exp2[n];
	int I_max_2 = I_max/2;
	for(int j=0; j<J_max; j++) { 
		int k=0, l=0;

		for(int i=0;i<I_max;i++) { 
			//------------------------------------
			if(k!=-1&&l!=-1) {	
				if(buf_in_max[j*I_max+k] > buf_in_max[j*I_max +I_max_2  +l]) {
					buf_out_max[j*I_max+i] = buf_in_max[j*I_max+k];
					buf_out_ind[j*I_max+i] = buf_in_ind[j*I_max+k];
					k++; if(k>=I_max_2) k=-1;   
				} else {
					buf_out_max[j*I_max+i] = buf_in_max [j*I_max +I_max_2  +l];
					buf_out_ind[j*I_max+i] = buf_in_ind [j*I_max +I_max_2  +l];
					l++; if(l>=I_max_2) l=-1;   
				}
			} else {
				if(k!=-1){
					buf_out_max[j*I_max+i] = buf_in_max[j*I_max+k];
					buf_out_ind[j*I_max+i] = buf_in_ind[j*I_max+k];
					k++;  
				}
				if(l!=-1){
					buf_out_max[j*I_max+i] = buf_in_max [j*I_max +I_max_2  +l];
					buf_out_ind[j*I_max+i] = buf_in_ind [j*I_max +I_max_2  +l];
					l++;
				}
			}
		}//-------------------------------------
	}
}

inline float  EDP::Get_mul_gr( float *buf, EDP::BYTE* bufb) {
	int size =m_nPixels;
	float mn[3]={0,0,0};
	for( int i= 0; i<size; i++) {	
		for( int ci= 0; ci<3; ci++)
			mn[ci]+=bufb[i+ci*m_nPixels];
	}
	for( int ci= 0; ci<3; ci++)
		mn[ci]/=m_nPixels;
	
	float sig_b=0;
	for( int i= 0; i<size; i++) {
		for( int ci= 0; ci<3; ci++)
			sig_b+=(mn[ci]-bufb[i+ci*m_nPixels])*(mn[ci]-bufb[i+ci*m_nPixels]);
	}
	
	
	float sig_f=0;
	for( int i= 0; i<size; i++) {
		{
			for( int ci= 0; ci<6; ci++)
				sig_f+=(buf[i+ci*m_nPixels])*(buf[i+ci*m_nPixels]);
		}
	}
	
	sig_b=sqrt(sig_b/m_nPixels);
	sig_f=sqrt(sig_f/m_nPixels); 
	
	float ret = sig_b/sig_f;
	sigma_img = 1./sigma_img; 
	sigma_img *= sig_f;
	return ret;
}
void EDP::GaussCosConv2DFst(float sigma,int x, int y, float* buf, float * mask)
{
	//------------------------------
    alpha_gauss_cos_opt=1.3867895353;
    double pisq = asin(1.)*asin(1.)*4;
	double alphasq2 = alpha_gauss_cos_opt*alpha_gauss_cos_opt*2;
    beta_gauss_cos = exp(-pisq/alphasq2);
    A_gauss_cos=(1.+beta_gauss_cos)/2;;
    B_gauss_cos=(1.-beta_gauss_cos)/2;
    gamma_gauss_cos=asin(1.)*2/alpha_gauss_cos_opt;
	iW_gauss_cos=(int)(gamma_gauss_cos*sigma);
	double thr = 1./ (x*y*x);
	cmplxCosSin= new cmplx [iW_gauss_cos+2];
	teta_gauss_cos=alpha_gauss_cos_opt/sigma;
	for(int i=0;i<iW_gauss_cos+2;i++)
	{cmplxCosSin[i].re=cos(teta_gauss_cos*i); cmplxCosSin[i].im=sin(teta_gauss_cos*i);}
	//------------------------------
	double * buf_out = new double [x*y];
	double * wgh_out = new double [x*y];
	for(int i = 0; i < x*y; i++){buf_out[i] = buf[i]*mask[i];  wgh_out[i] = mask[i];}
	/*float mul = (float)x*y/cnt;*/
     dirBuf = new double [x+y];
	 invBuf = new double [x+y];
	 cmplxBuf=new cmplx [x+y];
	 cmplxBuf2=new cmplx [x+y];
	 double * yBuf=new  double  [x+y];
	 double * yWgh=new  double  [x+y];
	//==============
	
	for(int  yy = 0; yy <  y; yy++)
	{
	GaussCosConvFst(x,  &wgh_out[yy*x]);
    GaussCosConvFst(x,  &buf_out[yy*x]);
	}

	//===============
for(int xx = 0; xx < x; xx++)
{
	for(int yy = 0; yy < y ; yy++){ yBuf[yy]   = buf_out[yy*x+ xx]; yWgh[yy] =  wgh_out[yy*x+ xx];}
   GaussCosConvFst(y, yWgh);
   GaussCosConvFst(y, yBuf);
   for(int yy = 0; yy < y; yy++){ buf[yy*x+ xx] = yBuf[yy];/* (yWgh[yy]) ? yBuf[yy]/yWgh[yy]: 0;*/ mask[yy*x+ xx] = yWgh[yy]; }
}
//---------------------------------------
	//for(int i = 0; i < x*y; i++){buf[i] *= mask[i]; }
//---------------------------------------
	delete [] buf_out;
	delete [] wgh_out;
	delete [] dirBuf;
	delete [] invBuf;
	delete [] yBuf;
	delete [] cmplxBuf;
	delete [] cmplxBuf2;
	delete [] yWgh;
	delete [] cmplxCosSin;

}
void EDP::GaussConv2D(float sigma,int x, int y, float* buf)
{
	//------------------------------
	
   int iW = 4*sigma;
	float *wt_r = new float  [iW+1];
	teta_gauss_cos=alpha_gauss_cos_opt/sigma;
	for(int i=0;i<iW+1;i++)
	{wt_r[i] = exp(-(float)i*i/2/(sigma*sigma));}
	//------------------------------
	double * buf_out = new double [x*y];
	for(int i = 0; i < x*y; i++){buf_out[i] = buf[i];}
	 double * yBuf=new  double  [x+y];

	for(int  yy = 0; yy <  y; yy++)
	{
    GaussConv( x,  &buf_out[yy*x], wt_r, iW);
	}

	//===============
for(int xx = 0; xx < x; xx++)
{
	for(int yy = 0; yy < y ; yy++){ yBuf[yy]   = buf_out[yy*x+ xx]; }
    GaussConv( y,  yBuf, wt_r, iW);
   for(int yy = 0; yy < y; yy++){ buf[yy*x+ xx] =  yBuf[yy];}
}
//---------------------------------------
//---------------------------------------
	delete [] buf_out;
	delete [] yBuf;
	delete [] wt_r;

}
void EDP::GaussCosConv2DFst(float sigma,int x, int y, float* buf)
{
	//------------------------------
	
	//double cosalph = acos((exp(-0.5) - 0.538)/0.462);
	//double cosalphr = acos((exp(-0.5) - 0.5)/0.5);
    alpha_gauss_cos_opt=1.25;//1.3867895353;
    double pisq = asin(1.)*asin(1.)*4;
	double alphasq2 = alpha_gauss_cos_opt*alpha_gauss_cos_opt*2;
    beta_gauss_cos = exp(-pisq/alphasq2);
    A_gauss_cos=(1.+beta_gauss_cos)/2;
    B_gauss_cos=(1.-beta_gauss_cos)/2;
    gamma_gauss_cos=asin(1.)*2/alpha_gauss_cos_opt;
	iW_gauss_cos=(int)(gamma_gauss_cos*sigma);
	double thr = 1./(iW_gauss_cos*iW_gauss_cos*50);
	cmplxCosSin= new cmplx [iW_gauss_cos+2];
	teta_gauss_cos=alpha_gauss_cos_opt/sigma;
	for(int i=0;i<iW_gauss_cos+2;i++)
	{cmplxCosSin[i].re=cos(teta_gauss_cos*i); cmplxCosSin[i].im=sin(teta_gauss_cos*i);}
	//------------------------------
	double * buf_out = new double [x*y];
	for(int i = 0; i < x*y; i++){buf_out[i] = buf[i];}
	 double * yBuf=new  double  [x+y];
	//==============
	
	for(int  yy = 0; yy <  y; yy++)
	{
    GaussCosConvFst_(x,  &buf_out[yy*x]);
	}

	//===============
for(int xx = 0; xx < x; xx++)
{
	for(int yy = 0; yy < y ; yy++){ yBuf[yy]   = buf_out[yy*x+ xx]; }
   GaussCosConvFst_(y, yBuf);
   for(int yy = 0; yy < y; yy++){ buf[yy*x+ xx] =(yBuf[yy]>= thr)? yBuf[yy] : 0;}
}
//---------------------------------------
//---------------------------------------
	delete [] buf_out;
	delete [] yBuf;
	delete [] cmplxCosSin;

}
void EDP::GaussCosConv2DFst(float sigma,  int iq, short int * map,  float* lmbK,  int x, int y, float* buf)
{
	//------------------------------
	
	//double cosalph = acos((exp(-0.5) - 0.538)/0.462);
	//double cosalphr = acos((exp(-0.5) - 0.5)/0.5);
    alpha_gauss_cos_opt=1.25;//1.3867895353;
    double pisq = asin(1.)*asin(1.)*4;
	double alphasq2 = alpha_gauss_cos_opt*alpha_gauss_cos_opt*2;
    beta_gauss_cos = exp(-pisq/alphasq2);
    A_gauss_cos=(1.+beta_gauss_cos)/2;
    B_gauss_cos=(1.-beta_gauss_cos)/2;
    gamma_gauss_cos=asin(1.)*2/alpha_gauss_cos_opt;
	iW_gauss_cos=(int)(gamma_gauss_cos*sigma);
	double thr = 1./(iW_gauss_cos*iW_gauss_cos*50);
	cmplxCosSin= new cmplx [iW_gauss_cos+2];
	teta_gauss_cos=alpha_gauss_cos_opt/sigma;
	for(int i=0;i<iW_gauss_cos+2;i++)
	{cmplxCosSin[i].re=cos(teta_gauss_cos*i); cmplxCosSin[i].im=sin(teta_gauss_cos*i);}
	//------------------------------
	double * buf_out = new double [x*y];
	for(int i = 0; i < x*y; i++){buf_out[i] = buf[i];}
	 double * yBuf=new  double  [x+y];
	//==============
	
	for(int  yy = 0; yy <  y; yy++)
	{
    GaussCosConvFst_(x,  &buf_out[yy*x]);
	}

	//===============
for(int xx = 0; xx < x; xx++)
{
	for(int yy = 0; yy < y ; yy++){ yBuf[yy]   = buf_out[yy*x+ xx]; }
   GaussCosConvFst_(y, yBuf);
   for(int yy = 0; yy < y; yy++){ float a =(yBuf[yy]>= thr)? yBuf[yy] : 0; buf[yy*x+ xx] =  (iq==map[yy*x+ xx])? buf[yy*x+ xx]*lmbK[yy*x+ xx]+a : a; }
}
//---------------------------------------
//---------------------------------------
	delete [] buf_out;
	delete [] yBuf;
	delete [] cmplxCosSin;

}
void EDP::GaussCosConv2D_Clr(float sigma, float sig_cl, int x, int y, float* buf, unsigned char* I_b)
{
	   int n_c = 1000;
	  int cc[3];
	  float * gss_c_wt = new float [n_c];
     get_std( sig_cl, I_b,  gss_c_wt,  n_c );
	//------------------------------
    alpha_gauss_cos_opt=1.25;
    double pisq = asin(1.)*asin(1.)*4;
	double alphasq2 = alpha_gauss_cos_opt*alpha_gauss_cos_opt*2;
    beta_gauss_cos = exp(-pisq/alphasq2);
    A_gauss_cos=(1.+beta_gauss_cos)/2;;
    B_gauss_cos=(1.-beta_gauss_cos)/2;
    gamma_gauss_cos=asin(1.)*2/alpha_gauss_cos_opt;
	iW_gauss_cos=(int)(gamma_gauss_cos*sigma);
	double thr = 1./ (x*y*x);
	float *CosSin= new float [iW_gauss_cos+2];
	teta_gauss_cos=alpha_gauss_cos_opt/sigma;
	for(int i=0;i<iW_gauss_cos+2;i++)
	{CosSin[i] = (A_gauss_cos + B_gauss_cos*cos(teta_gauss_cos*i));}
	//------------------------------
	float * buf_out = new float [x*y];
	for(int i = 0; i < x*y; i++){buf_out[i] = buf[i];}
	for(int i = 0; i < x*y; i++){buf[i] =gss_wnd_c( i,    I_b,   iW_gauss_cos, gss_c_wt, CosSin, buf_out);}

//---------------------------------------

    delete [] gss_c_wt;
	delete [] buf_out;
	delete [] CosSin;

}
inline void EDP:: GaussCosConvFst(int N,  double *InBuf)
{

		for(int i = 0; i < N; i++) {invBuf[i] = cmplxBuf2[i].re =InBuf[i];cmplxBuf2[i].im =0;  }
       CmplxExpFst(  teta_gauss_cos,  N, cmplxBuf2);
	   MuFst(  N, invBuf);
	   for(int i = 0; i < N; i++) InBuf[i] = cmplxBuf2[i].re*B_gauss_cos + A_gauss_cos*invBuf[i];


}
inline double EDP::CmplxExpFst(double teta, int iN, EDP:: cmplx * inbuf)
{
	for(int i=0;i<iN;i++) {cmplxBuf[i]=inbuf[i]; }
	cmplx sum={0,0};
	//for(int i=0;i<=iW_gauss_cos;i++) sum=addCmplx(sum, mulCmplx(fromBufN( i,  iN, cmplxBuf), expTet(i)));
	int fn = (iW_gauss_cos+1 < iN) ?  iW_gauss_cos+1 : iN;
	//for(int i=0;i < fn ;i++){ sum = addCmplx(sum, mulCmplx (cmplxBuf[i], expTet(i)) );}
	for(int i=0;i < fn ;i++){ if(inbuf[i].re){sum.re += inbuf[i].re * cmplxCosSin[i].re; sum.im += inbuf[i].re * cmplxCosSin[i].im; }}
	inbuf[0]=sum;
    cmplx et_1 = expTet(-1);
	cmplx et_m = expTet(-iW_gauss_cos-1);
	cmplx et_p  = expTet(iW_gauss_cos);
	for(int i=1, ii = -iW_gauss_cos, iii =iW_gauss_cos+1; i<iN; i++, ii++, iii++)
	{
		//cmplx s1 = mulCmplx(inbuf[i-1], expTet(-1));
	 //   cmplx s2 = mulCmplx(fromBufN( i-iW_gauss_cos-1,  iN, cmplxBuf), expTet(-iW_gauss_cos-1));
	 //   s2.re = -s2.re; s2.im = -s2.im;
	 //   cmplx s3 = mulCmplx(fromBufN( i+iW_gauss_cos,  iN, cmplxBuf), expTet(iW_gauss_cos));
	 //   s1=addCmplx(s1, s2);
	 //   inbuf[i]=addCmplx(s1, s3);
		 inbuf[i] = mulCmplx(inbuf[i-1], et_1);
		 if(ii >= 0){inbuf[i].re -= cmplxBuf[ii].re * et_m.re; inbuf[i].im -= cmplxBuf[ii].re * et_m.im; }
	    //cmplx s2 = mulCmplx(fromBufN( ii,  iN, cmplxBuf), et_m);
	    //s2.re = -s2.re; s2.im = -s2.im;
         if(iii< iN){inbuf[i].re += cmplxBuf[iii].re * et_p.re;  inbuf[i].im += cmplxBuf[iii].re * et_p.im; }
	}
return 0;
}
inline void EDP::GaussCosConvFst_( int N, double * inbuf)
{
    double * ib = new double [N];
    for(int i=0;i< N;i++) { ib[i] = inbuf[i]; }
	cmplx sum={0,0}; double sm = 0;
	double a = A_gauss_cos; double b = B_gauss_cos;
	int fn = (iW_gauss_cos+1 < N) ?  iW_gauss_cos+1 : N;
	for(int  i = 0; i < fn ; i++) { sum.re += ib [i] * cmplxCosSin[i].re; sum.im += ib [i] * cmplxCosSin[i].im;  sm += ib [i]; }
	inbuf[0]=sum.re * b + sm * a;
    cmplx et_1 = expTet(-1);
	cmplx et_m = expTet(-iW_gauss_cos-1);
	cmplx et_p  = expTet(iW_gauss_cos);
	for(int i=1, ii = -iW_gauss_cos, iii =iW_gauss_cos+1; i<N; i++, ii++, iii++)
	{
		 sum = mulCmplx(sum, et_1);
		 if(ii >= 0){ sum.re -= ib[ii] * et_m.re; sum.im -= ib[ii] * et_m.im;  sm -= ib[ii]; }
         if(iii <  N){sum.re += ib[iii] * et_p.re;  sum.im += ib[iii] * et_p.im;  sm += ib[iii]; }
		 inbuf[i] = sum.re *b + sm * a;
	}
	 delete [] ib;
}
inline void EDP::GaussConv( int N, double * inbuf, float * wt, int iW)
{
    double * ib = new double [N];
    for(int i=0;i< N;i++) { ib[i] = inbuf[i]; }
	for(int i=0;i< N;i++)
	{ double sum = 0;
	for(int ii=-iW; ii<=iW; ii++)
	{    int ind = i + ii;
		 if(ind>=0&& ind<N)sum += wt[abs(ii)]*ib[ind] ;
	}
	inbuf[i] = sum;
	}
	 delete [] ib;
}
inline EDP::cmplx EDP::mulCmplx( EDP::cmplx a, EDP::cmplx b)
{

	cmplx ret = {a.re*b.re-a.im*b.im,  a.re*b.im+a.im*b.re} ;
	return ret ;
}
inline EDP::cmplx  EDP::addCmplx( EDP::cmplx a, EDP::cmplx b)
{
	cmplx ret= {a.re+b.re,  a. im+b.im} ;
	return ret;
}

inline EDP::cmplx  EDP::expTet( int i)
{
cmplx ret;
if(i>=0){ret.re=cmplxCosSin[i].re; ret.im=cmplxCosSin[i].im; return ret;}
else  {ret.re=cmplxCosSin[abs(i)].re; ret.im=- cmplxCosSin[abs(i)].im; return ret;}

}
inline void EDP::MuFst( int iN, double * inbuf)
{
	for(int i=0;i<iN;i++) {dirBuf[i]=inbuf[i]; }
	double sum=0;
	for(int i=0;i<=iW_gauss_cos;i++) if(i<iN) sum += dirBuf[i]; //fromBufN( i,  iN, dirBuf);
	inbuf[0] = sum;
	for(int i=1, ii = -iW_gauss_cos, iii =iW_gauss_cos+1; i<iN; i++, ii++, iii++)
	{
	inbuf[i] = inbuf[i-1];
	if(ii>=0 && dirBuf[ii])inbuf[i] -= dirBuf[ii];
	if(iii<iN && dirBuf[iii])  inbuf[i] +=  dirBuf[iii];/*fromBufN( i-iW_gauss_cos-1,  iN, dirBuf)*/ /* fromBufN( i+iW_gauss_cos,  iN, dirBuf)*/
	}
}
inline float EDP::fromBufN(int i, int iN, double* buf)
{
	return ((i>=0&&i<iN)? buf[i]:0);

}
inline EDP::cmplx EDP::fromBufN(int i, int iN, EDP::cmplx* buf)
{
	cmplx zr={0,0};  return ((i>=0&&i<iN)? buf[i]:zr);
}

inline void EDP:: SampledFlt(float sigma_c_mul, float sigma, int cl_Q, unsigned char * I_b, float *buf_to_flt, short int * clr_cls,  EDP::COLOR_F * K_cls, int r, int nb_ch )
{ ///////////////////////////

	  float ** buf = new float* [cl_Q];
	  float ** bufs = new float* [cl_Q];
	  float ** mask = new float * [cl_Q];
	  float ** maskc = new float * [cl_Q];
	  double sigc = sigma_c_mul*255;

	  int div = 256/r;
	  int n_c = 1000;
	  int cc[1000] ;
	  float ss[1000];
	  float ccf[1000];
	  float * gss_c_wt = new float [n_c];
	  //float * ijK = new float [cl_Q*cl_Q];
	  //double *sigmaK = new double  [cl_Q];
	  //float * sigmaSK = new float [cl_Q];
	  //int *sigma_cnt = new int [cl_Q];
	  short int * q_map = new short int [m_nPixels];
	  for(int i = 0; i < cl_Q; i++) {mask[i] = new float [m_nPixels]; maskc[i] = new float [m_nPixels]; }
	  get_std( sigma_c_mul, I_b,  gss_c_wt,  n_c );
	  //------------------------
						  for(int p = 0; p < m_nPixels; p++)
						  {
							  for(int c =0; c < 3; c ++){cc[c] = CUTFL((float)I_b[p + c*m_nPixels]/div, r); }
							     int cls = q_map[p] =  clr_cls[IND_HC(cc[0], cc[1], cc[2], r)];
								 for(int i = 0; i < cl_Q; i++) maskc[i][p]=mask[i][p] = (cls == i) ? 1 : 0;
								 //float cwt =0;
								 //for(int i = 0; i < cl_Q; i++)
								 //{double dlt =0; for(int c =0; c < 3; c ++)
								 //{float dlt0 =  K_cls[i].cl[c] -  I_b[p + c*m_nPixels]; dlt +=dlt0*dlt0;}
								 //dlt = sqrt(dlt); if(dlt>=n_c) dlt = n_c-1; cwt += mask[i][p] = gss_c_wt[(int)dlt]; }
								 // for(int i = 0; i < cl_Q; i++) maskc[i][p] = mask[i][p] /= cwt;
						  }


	 //------------------------------------------------------------
	for(int i = 0; i < cl_Q; i++)buf[i] = new float [m_nPixels*3];
    for(int i = 0; i < cl_Q; i++)bufs[i] = new float [m_nPixels];
	   //---------------------------------------------------------
/*	  	                  for(int ic = 0; ic < nb_ch; ic++)
						  for(int i = 0; i < cl_Q; i++)
	                       {
						    for(int p = 0; p < m_nPixels; p++) bufs[i][p+ic*m_nPixels] = buf_to_flt[p+ic*m_nPixels]*mask[i][p];
						   GaussCosConv2DFst( sigma,  m_width, m_height,  &bufs[i][ic*m_nPixels] );
							}   */
						  for(int ic = 0; ic < 3; ic++)
						  for(int i = 0; i < cl_Q; i++)
	                       {
						    for(int p = 0; p < m_nPixels; p++) buf[i][p+ic*m_nPixels] = (float)I_b[p+ic*m_nPixels]*maskc[i][p];
						   GaussCosConv2DFst( sigma,  m_width, m_height,  &buf[i][ic*m_nPixels] );
							}
	  for(int i = 0; i < cl_Q; i++) GaussCosConv2DFst( sigma, m_width, m_height, mask[i]);
  /*    float * ksi_b = new float [m_nPixels];
	  float * ksi_m = new float [m_nPixels]; for(int p = 0; p < m_nPixels; p++) ksi_m[p] =1;  GaussCosConv2DFst( sigma,  m_width, m_height, ksi_m);*/

/////////////////////////////////////////////////////////////////////////////
	                          //double ksi_mean = 0; double ml =  mul_tab(cl_Q);
							  //---------------


	  					//	 for(int p = 0; p < m_nPixels; p++)
						  //  {
							 // int iq = q_map[p]; float sum =0;  float ds, dst = 0;
							 // for(int c =0; c < 3; c ++){ ds = buf[iq][p +c*m_nPixels]/mask[iq][p] - I_b[p + c*m_nPixels]; ds *= ds; sum += ds; }
        //                     ksi_b[p] = sum;
							 //}
							 //GaussCosConv2DFst( sigma,  m_width, m_height, ksi_b);
							 //for(int p = 0; p < m_nPixels; p++){ ksi_b[p] /= ksi_m[p]; ksi_b[p] = sqrt(ksi_b[p])*ml;}
							 //for(int p = 0; p < m_nPixels; p++)
							 //{
							 //ksi_b[p] = fi_erf(ksi_b[p] , sigc);
							 //}
/////////////////////////////////////////////////////////////////////////////
       for(int ic = 0; ic < nb_ch; ic++)
	  {
      	  for(int i = 0; i < cl_Q; i++)
	                       {
							   for(int p = 0; p < m_nPixels; p++){ bufs[i][p] = buf_to_flt[p+ic*m_nPixels]*maskc[i][p];  }
						      GaussCosConv2DFst( sigma,  m_width, m_height, bufs[i] );
							}
         for(int p = 0; p < m_nPixels; p++){ buf_to_flt[p+ic*m_nPixels] = 0; }
		  for(int p = 0; p < m_nPixels; p++)
						  {

							  for(int c =0; c < 3; c ++){ cc[c] = I_b[p + c*m_nPixels]; }
	/*						  for(int c =0; c < nb_ch; c ++){ ss[c] = buf_to_flt[p + c*m_nPixels];  buf_to_flt[p+ c*m_nPixels] =0;} */
							   float w_sum =0;   							   							
								for(int i = 0; i <cl_Q; i++)
								{
                                    float sv, cv, w,  dist = 0;  int iq = q_map[p];

									if(mask[i][p])
									{

                                     sv = bufs[i][p]/mask[i][p];
	/*								 if(iq==i){ cv = cv*ksi_b[p] + (1. - ksi_b[p])*ss[c];}		*/	
									//ccf[c] = cv;
						
									for(int c =0; c < 3; c ++)
									{
                                     cv =  buf[i][p +c*m_nPixels]/mask[i][p];
/*									 if(iq==i){ cv = cv*ksi_b[p] + (1. - ksi_b[p])*cc[c];}		*/	
									dist +=(cv - cc[c])*(cv - cc[c]);
									}
									dist = sqrt(dist); w =  gss_c_wt[round_fl(dist)]*mask[i][p];
									
									}
									else {w = 0;}
									
									w_sum += w;
									buf_to_flt[p +ic*m_nPixels] += (w)? sv*w :0;
								}
								    buf_to_flt[p +ic*m_nPixels] /= (w_sum)? w_sum:1;
								//if(buf_to_flt[p +c*m_nPixels]<0)buf_to_flt[p +c*m_nPixels]=0; if(buf_to_flt[p +c*m_nPixels]>255)buf_to_flt[p +c*m_nPixels]=255;
								
						  }
	  }
     //------------------------------------
	 delete [] gss_c_wt;
	 for(int i = 0; i < cl_Q; i++){ delete [] buf[i]; delete [] bufs[i]; delete [] mask[i]; delete [] maskc[i];}
	 delete [] buf;
	 delete [] mask;
	 delete [] maskc;
	 delete [] bufs;
	 delete [] q_map;
	 //delete [] sigma_cnt;
	 //delete [] sigmaSK;
	 //delete [] sigmaK;
	 //delete [] ksi_b;
	 //delete [] ksi_m;

	}////////////////////////////////////
	inline double EDP:: mul_tab(int q)
	{                        double ret;
							  if(q ==1) ret = 2.45;
							  if(q ==2) ret = 2.4;
							  if(q>=3&&q<6) ret = 2.30 -  (q-3)*0.15/3;
							  if(q>=6&&q<12) ret = 2.15 -  (q-6)*0.1/6;
							  if(q>=12&&q<24) ret = 2.05 -  (q-12)*0.25/12;
							  if(q>=24&&q<48) ret = 1.8 -  (q-24)*0.2/24;
							  if(q>=48&&q<96) ret = 1.6 -  (q-48)*0.05/48;
							  if(q>=96) ret = 1.55;
							  return ret;
	}
inline void EDP:: SampledFlt__(float sigma_c_mul, float sigma, int cl_Q, unsigned char * I_b, float *buf_to_flt, short int * clr_cls,  EDP::COLOR_F * K_cls, int r, int nb_ch )
{ ///////////////////////////
	  int div = 256/r;
	  int n_c = 1000;
	  int cc[3] ;
	  float ** mask = new float * [cl_Q];
	  float * gss_c_wt = new float [n_c];
	  short int * q_map = new short int [m_nPixels];
	  for(int i = 0; i < cl_Q; i++) {mask[i] = new float [m_nPixels]; }
	  //------------------------
						  for(int p = 0; p < m_nPixels; p++)
						  {
							     for(int c =0; c < 3; c ++)cc[c] = CUTFL((float)I_b[p + c*m_nPixels]/div, r);
							     int cls = q_map[p] =  clr_cls[IND_HC(cc[0], cc[1], cc[2], r)];
								 for(int i = 0; i < cl_Q; i++)mask[i][p] = (cls == i) ? 1 : 0;
						  }
	  //------------------------
     get_std( sigma_c_mul, I_b,  gss_c_wt,  n_c );
	   //---------------------------------------------------------
	  for(int i = 0; i < cl_Q; i++)
	  {
	  GaussCosConv2DFst( sigma, m_width, m_height, mask[i]);
	  }

	
////////////////---------------------------------------------------------------------------------
	 	                      for(int ic = 0; ic < nb_ch; ic++)
						  {

						  for(int p = 0; p < m_nPixels; p++)
						  {

							    for(int c =0; c < 3; c ++)cc[c] = I_b[p + c*m_nPixels];
							   buf_to_flt[p+ic*m_nPixels] =0; float w_sum =0;   							   							
								for(int i = 0; i <cl_Q; i++)
								{
                                    float w,  dist = 0;
									for(int c =0; c < 3; c ++)
									{
									dist +=(cc[c] - K_cls[i].cl[c])*(cc[c] - K_cls[i].cl[c]);
									}
									dist = sqrt(dist); w = gss_c_wt[round_fl(dist)];
									
									
									w_sum += w*mask[i][p];
                                     buf_to_flt[p +ic*m_nPixels] += K_cls[i].cl[ic]*w*mask[i][p];
								}
								buf_to_flt[p +ic*m_nPixels] /= (w_sum)? w_sum:1;
						  }
						  }

     //------------------------------------
	 delete [] gss_c_wt;
	 for(int i = 0; i < cl_Q; i++){   delete [] mask[i]; }
	 delete [] mask;
	 delete [] q_map;
	}////////////////////////////////////
inline int EDP:: Gss_rand( float sk, int r_w)
{         int i, ret; float sum;
           for(i=0, sum=0;i<12;i++) {sum+=(double)rand()/RAND_MAX;}
          sum=( sum - 6 ) * sk;
		  ret = round_fl(sum);
          if(abs(ret)>r_w)ret = r_w*ret/abs(ret);
		  return ret;
}
void EDP::MakeMK_mask(int nsp, int nmw, float sk, int r_w, float ** crd_w, POINT ** crd)
{
srand( 27644437 );
int i, dm = r_w*2+1;
float sum;
BYTE * mask = new BYTE [dm*dm]; for (int i =0; i < dm*dm; i++) mask[i] = 0;

for (int  n=0; n < nmw; n++)
{
	for(int p = 0 ; p < nsp; p++)
	{
		int st = 1; while(st){
		   int x = Gss_rand( sk, r_w);
		   int y =  Gss_rand( sk, r_w);
		   int ind = (x+r_w) + (y + r_w)*dm;
		   if(!mask[ind])
		   { mask[ind] =1; st=0; crd_w[n][p] = exp(-(x*x+y*y)/2/sk/sk);  crd[n][p].x = x;  crd[n][p].y = y; }
		}

	}
/////////////////////////////////////
for (int i =0; i < dm*dm; i++) mask[i] = 0;
}


	
//=======================
delete [] mask;
}


float  EDP::dev_map( BYTE * l_b,  float sigma_c_mul, float sigma, short int * map,  float * mul_buf)
{
    #define NMBS 6
	int nmbs = NMBS*4+1;
	POINT ind_w[NMBS*4+1];
	 float  ind_w_w[NMBS*4+1 ];
    ind_w[nmbs-1].x = 0;
	ind_w[nmbs-1].y = 0;
	ind_w_w[nmbs-1] = 1;
	//      int nsp =64, nmw = 23; float sk = 10; int r_w = 25, dm = 51;
 //    float ** crd_w = new float* [nmw]; POINT ** crd = new POINT* [nmw];
	// for(int i = 0; i < nmw; i ++){crd_w[i] = new float [nsp]; crd[i] = new POINT [nsp];}
	//MakeMK_mask( nsp,  nmw,  sk,   r_w,   crd_w,  crd);

 //   for(int p =0; p < m_nPixels*3; p++ )I_ims[1][p] =0;
 //   for(int i =r_w; i < m_width - r_w ; i += dm )
	//for(int j =r_w; j < m_height - r_w ; j += dm )
	//for(int p =0; p < nsp; p ++ )
	//{ int x = i + crd[(i+j)%nmw][p].x;
	//   int y = j + crd[(i+j)%nmw][p].y;
 //     for(int c =0; c < 3; c ++ ) I_ims[1][x + y*m_width +c*m_nPixels ] =  I_ims[0][x + y*m_width +c*m_nPixels ];
	//}
	//for(int i = 0; i < nmw; i ++){delete [] crd_w[i]; delete [] crd[i];}
	//delete [] crd_w ; delete [] crd;


	 float rf = 0.67*sigma; int ind_cnt =0;
  for(int p= 0; p< m_nPixels*3; p++){ mul_buf[p]  =1;}
	 for(int i = -2; i<=2; i++) for(int j = -2; j<=2; j++)
	 { if(/*!(abs(i)==2&&abs(j)==2)&&*/!(((i)==0&&(j)==0)))
	 {
		 float x = ind_w[ind_cnt].x = round_fl(rf*i);
	     float y = ind_w[ind_cnt].y = round_fl(rf*j);
		 if((abs(i)==2&&abs(j)==2))
		 {
		 float x = ind_w[ind_cnt].x = i/2;
	     float y = ind_w[ind_cnt].y = j/2;
		 }
		 ind_w_w[ind_cnt++] = exp(-(x*x+y*y)/2/sigma/sigma);
	 }
	 }
     int cc[3], cc_[3]; float sumr[3], sums[3], vr[3], vs[3];
	 int q;  float wr,  wrsm, wssm,  dist,  ww;
	 float sgc = 255*sigma_c_mul; sgc *= 2*sgc; sgc = 1./sgc;
	 for(int p = 0 ; p <m_nPixels; p++)
	 {
					  int x = p%m_width, y = p/m_width;
					  for(int c = 0; c< 3; c++)cc[c] = l_b[p+c*m_nPixels];
					 wrsm = wssm = 0; q = map[p];
					 for(int c = 0; c< 3; c++){sumr[c] =0; sums[c] =0;}
					 for(int i =0; i<nmbs; i++)
					 {
					 int  xx =  ind_w[i].x + x, yy = ind_w[i].y + y;
					 if(xx>=0&&xx< m_width && yy>=0 && yy < m_height)
					 { int pp = xx+ yy*m_width;
					 dist =0;
					 if(q == map[pp])  {
					for(int c = 0; c< 3; c++){float aa = cc[c] - (cc_[c] = l_b[pp+c*m_nPixels]); dist += aa*aa;}
					 dist *= sgc; wr = exp(-dist)*ind_w_w[i];
					 for(int c = 0; c< 3; c++){sumr[c] += wr*cc_[c]; sums[c] +=  ind_w_w[i]*cc_[c];}
					 wrsm += wr; wssm += ind_w_w[i];
					 }}}
					 if(wrsm){
						 for(int c = 0; c< 3; c++){
							 mul_buf[p + c*m_nPixels]= sumr[c]/wrsm;
							  //ww = wssm/wrsm;
         //                  mul_buf[p + c*m_nPixels]= ww*sumr[c]/sums[c];
						 //float vl = -(float)cc[c] + (sums[c]/wssm); if(vl) mul_buf[p + c*m_nPixels]  =((-(float)cc[c]+ sumr[c]/wrsm)/vl);
	      //               if(mul_buf[p + c*m_nPixels] > 1) mul_buf[p + c*m_nPixels]=1; if(mul_buf[p + c*m_nPixels] <0) mul_buf[p + c*m_nPixels]=0;
					 }}

	 }
	 //-------------------------------------------------
 //dev_map3x3( l_b,  sigma_c_mul, sigma, map,  mul_buf) ;
   return 0;// sum/m_nPixels;
	 }
void  EDP::dev_map3x3( BYTE * l_b,  float sigma_c_mul, float sigma, short int * map,  float * mul_buf)
{
	int nmbs = 9;
	POINT ind_w[9];
	 float  ind_w_w[9];
	 float *mul_buf_c = new float [ m_nPixels*3];
	   for(int p= 0; p< m_nPixels*3; p++){ mul_buf_c[p]  =mul_buf[p] ;}


      int ind_cnt =0;
	 for(int i = -1; i<=1; i++) for(int j = -1; j<=1; j++)
	 {

		 float x = ind_w[ind_cnt].x;// = round_fl(i);
	     float y = ind_w[ind_cnt].y ;//= round_fl(j);
		 ind_w_w[ind_cnt++] = exp(-(x*x+y*y)/2/sigma/sigma);

	 }
     int cc[3]; float sumr[3];
	 float wr,  wrsm,  dist,  ww;
	 float sgc = 255*sigma_c_mul; sgc *= 2*sgc; sgc = 1./sgc;
	 for(int p = 0 ; p <m_nPixels; p++)
	 {
					  int x = p%m_width, y = p/m_width;
					  for(int c = 0; c< 3; c++)cc[c] = l_b[p+c*m_nPixels];
					 wrsm  = 0;
					 for(int c = 0; c< 3; c++){sumr[c] =0;}
					 for(int i =0; i<nmbs; i++)
					 {
					 int  xx =  ind_w[i].x + x, yy = ind_w[i].y + y;
					 if(xx>=0&&xx< m_width && yy>=0 && yy < m_height)
					 { int pp = xx+ yy*m_width;
					 dist =0;
					  {
					for(int c = 0; c< 3; c++){float aa = cc[c] - ( l_b[pp+c*m_nPixels]); dist += aa*aa;}
					 dist *= sgc; wr = exp(-dist)*ind_w_w[i];
					 for(int c = 0; c< 3; c++){sumr[c] += wr*mul_buf_c[pp+c*m_nPixels]; }
					 wrsm += wr;
					 }}}
					 {
						 for(int c = 0; c< 3; c++){mul_buf[p + c*m_nPixels] = sumr[c]/wrsm;
					 }}

	 }
delete [] mul_buf_c;
	 }
void EDP::K_mean_Flt( EDP::BYTE * I_b,  float sigma_c_mul, float sigma, int cl_Q, float * buf_to_flt, int n_bufs)
{///------------------------ K_mean

	

	 COLOR_F * K_cls = new COLOR_F [cl_Q];
     int r = 32;  int div = 256/r; int size = r*r*r;
	 short int * clr_cls  = new short int [size];
	 int * buf_to_sort  = new int [size]; memset(buf_to_sort, 0, sizeof(int)*size);
         for(int y = 0; y < m_height; y++)for(int x = 0; x < m_width; x++)
		 {
			 int vl[3] = {CUTFL((float)I_b[IND_IC(x,y,0)]/div, r ), CUTFL((float)I_b[IND_IC(x,y,1)]/div, r ),CUTFL((float)I_b[IND_IC(x,y,2)]/div, r )} ;
             buf_to_sort[IND_HC(vl[0],vl[1],vl[2], r)]++;
		 }
	 K_means(cl_Q,  K_cls, r, buf_to_sort, clr_cls);
	 //	 for(int p =0; p < m_nPixels; p++ )
	 //{  int cc[3];
		// for(int c =0; c < 3; c ++)cc[c] = CUTFL((float)I_b[p + c*m_nPixels]/div, r);
		// int cls = clr_cls[IND_HC(cc[0], cc[1], cc[2], r)];
		// for(int c =0; c < 3; c ++)I_ims[1][p + c*m_nPixels ] = (cls==3)?/* K_cls[cls].cl[c]*/I_ims[0][p + c*m_nPixels ]:255;
	 //}
    SampledFlt(sigma_c_mul, sigma,  cl_Q,  I_b, buf_to_flt,  clr_cls,  K_cls, r,n_bufs );
	 delete [] buf_to_sort;
	 delete [] clr_cls;
	 delete [] K_cls;

	} ///------------------ end K means
void EDP::BLF_spr( EDP::BYTE * l_b,  float sigma_c_mul, float sigma,  float * buf_to_flt, int n_bufs)
{///------------------------ K_mean


    float * bfw = new float [m_nPixels*n_bufs];   for(int p= 0; p< m_nPixels*3; p++){ bfw[p]  =buf_to_flt[p];}
	float * cc_ = new float [n_bufs];
     #define NMBS 6
	int nmbs = NMBS*4+1;
	POINT ind_w[NMBS*4+1];
	 float  ind_w_w[NMBS*4+1 ];
    ind_w[nmbs-1].x = 0;
	ind_w[nmbs-1].y = 0;
	ind_w_w[nmbs-1] = 1;
	 float rf = 0.67*sigma; int ind_cnt =0;

	 for(int i = -2; i<=2; i++) for(int j = -2; j<=2; j++)
	 { if(/*!(abs(i)==2&&abs(j)==2)&&*/!(((i)==0&&(j)==0)))
	 {
		 float x = ind_w[ind_cnt].x = round_fl(rf*i);
	     float y = ind_w[ind_cnt].y = round_fl(rf*j);
		 if((abs(i)==2&&abs(j)==2))
		 {
		 float x = ind_w[ind_cnt].x = i/2;
	     float y = ind_w[ind_cnt].y = j/2;
		 }
		 ind_w_w[ind_cnt++] = exp(-(x*x+y*y)/2/sigma/sigma);
	 }
	 }
     int cc[3]; float sumr[3];
	 int q;  float wr,  wrsm,  dist,  ww;
	 float sgc = 255*sigma_c_mul; sgc *= 2*sgc; sgc = 1./sgc;
	 for(int p = 0 ; p <m_nPixels; p++)
	 {
					  int x = p%m_width, y = p/m_width;
					  for(int c = 0; c< 3; c++)cc[c] = l_b[p+c*m_nPixels];
					 wrsm  = 0;
					 for(int c = 0; c< 3; c++){sumr[c] =0; }
					 for(int i =0; i<nmbs; i++)
					 {
					 int  xx =  ind_w[i].x + x, yy = ind_w[i].y + y;
					 if(xx>=0&&xx< m_width && yy>=0 && yy < m_height)
					 { int pp = xx+ yy*m_width;
					 dist =0;
					  {
					for(int c = 0; c< 3; c++){float aa = cc[c] -  l_b[pp+c*m_nPixels]; dist += aa*aa;}
					 dist *= sgc; wr = exp(-dist)*ind_w_w[i];
					 for(int c = 0; c< n_bufs; c++){sumr[c] += wr*bfw[pp+c*m_nPixels];}
					 wrsm += wr;
					 }}}
					 if(wrsm)  for(int c = 0; c< n_bufs; c++) buf_to_flt[p + c*m_nPixels]= sumr[c]/wrsm;
					 else          for(int c = 0; c< n_bufs; c++)buf_to_flt[p + c*m_nPixels]= 0;

	 }
	 //-------------------------------------------------
delete [] bfw;
delete [] cc_;

	} ///------------------ end K means
void EDP::K_mean_Flt_cl( EDP::BYTE * I_b,  float sigma_c_mul, float sigma, int cl_Q, float * buf_to_flt)
{///------------------------ K_mean

	
     int n_bufs =3;
	 COLOR_F * K_cls = new COLOR_F [cl_Q];
     int r = 32;  int div = 256/r; int size = r*r*r;
	 short int * clr_cls  = new short int [size];
	 int * buf_to_sort  = new int [size]; memset(buf_to_sort, 0, sizeof(int)*size);
         for(int y = 0; y < m_height; y++)for(int x = 0; x < m_width; x++)
		 {
			 int vl[3] = {CUTFL((float)I_b[IND_IC(x,y,0)]/div, r ), CUTFL((float)I_b[IND_IC(x,y,1)]/div, r ),CUTFL((float)I_b[IND_IC(x,y,2)]/div, r )} ;
             buf_to_sort[IND_HC(vl[0],vl[1],vl[2], r)]++;
		 }
	 K_means(cl_Q,  K_cls, r, buf_to_sort, clr_cls);

	 //for(int p =0; p < m_nPixels; p++ )
	 //{  int cc[3];
		// for(int c =0; c < 3; c ++)cc[c] = CUTFL((float)I_b[p + c*m_nPixels]/div, r);
		// int cls = clr_cls[IND_HC(cc[0], cc[1], cc[2], r)];
		//  for(int c =0; c < 3; c ++)I_ims[2][p + c*m_nPixels ] = K_cls[cls].cl[c];
	 //}
      SampledFlt(sigma_c_mul, sigma,  cl_Q,  I_b, buf_to_flt,  clr_cls,  K_cls, r,n_bufs);
	 for(int p =0; p < m_nPixels; p++ )
	 {  int cc[3];
		 for(int c =0; c < 3; c ++)cc[c] = CUTFL(buf_to_flt[p + c*m_nPixels]/div, r);
		 int cls = clr_cls[IND_HC(cc[0], cc[1], cc[2], r)];
		  for(int c =0; c < 3; c ++)buf_to_flt[p + c*m_nPixels ] = K_cls[cls].cl[c];
	 }
	 delete [] buf_to_sort;
	 delete [] clr_cls;
	 delete [] K_cls;
	}


void   EDP:: Gr_I(unsigned char *in, POINT_F * b_0, POINT_F* Histo, int g_Ql, double   * mean, double * dsp )
{
//-------------------------------------------------
	POINT_F  H[256];
     for(int i = 0; i < 256; i++) {H[i].x =0; H[i].y =0;}
	 POINT_F * Vl = new POINT_F [g_Ql];  for(int i = 0; i < g_Ql; i++) {Vl[i].x = 512; Vl[i].y = 512; }

//-------------------------------------------------

		for(int y = 0; y< m_height ; y++)
		for(int x = 0; x< m_width ; x++)
		{
		m_G0[IND(x,y)].m =0;
         if(x == m_width -1)m_G0[IND(x,y)].x =0;
		 else                           m_G0[IND(x,y)].x  = (in[IND(x+1,y)] - in[IND(x,y)]) ;
		 if(y == m_height -1)m_G0[IND(x,y)].y =0;
		 else                           m_G0[IND(x,y)].y  = (in[IND(x,y+1)] - in[IND(x,y)]);
		}
        mean[0] =0; for(int p = 0; p < m_nPixels; p++) mean[0] +=  in[p ]; mean[0] /= m_nPixels;
		dsp[0] =0;    for(int p = 0; p < m_nPixels; p++) dsp[0] += (mean[0] - in[p ]) *( mean[0] - in[p]); dsp[0] = sqrt(dsp[0]/m_nPixels);
		POINT_F mn = {0,0}; for(int p = 0; p < m_nPixels; p++){ mn.x +=  m_G0[p].x; mn.y +=  m_G0[p].y;  } mn.x /= m_nPixels; mn.y /= m_nPixels;
		/////
		float mul = 1./(g_Ql - 1);
		 for(int p = 0; p < m_nPixels; p++) {H[(int)fabs(m_G0[p].x)].x += 1; H[(int)fabs(m_G0[p].y)].y += 1;}
         for(int p = 1; p < 256; p++) {H[p].x += H[p-1].x; H[p].y += H[p-1].y;}
		 for(int p = 0; p < 256; p++) {H[p].x /= m_nPixels*mul; H[p].y /= m_nPixels*mul;}
		 POINT st = {0,0};
		 for(int p = 0; p < 256; p++)
		 {
			 POINT fn = {(int)H[p].x, (int)H[p].y};  if(p == g_Ql-1) {fn.x =  g_Ql-1; fn.y =  g_Ql-1;}
			 for(int i = st.x; i<= fn.x; i++)if(Vl[i].x > p) Vl[i].x  = p;
			 for(int i = st.x; i<= fn.x; i++)if(Vl[i].y > p) Vl[i].y  = p;
		     st = fn;
		 }
		 for(int y = 0; y< m_height ; y++)
		for(int x = 0; x< m_width ; x++)
		{
		if(m_G0[IND(x,y)].x > 0) m_G0[IND(x,y)].x = Vl[(int)Histo[(int)b_0[IND(x,y)].x].x].x;
		if(m_G0[IND(x,y)].x < 0) m_G0[IND(x,y)].x = - Vl[(int)Histo[(int)b_0[IND(x,y)].x].x].x;
		if(m_G0[IND(x,y)].y > 0) m_G0[IND(x,y)].y = Vl[(int)Histo[(int)b_0[IND(x,y)].y].y].y;
		if(m_G0[IND(x,y)].y < 0) m_G0[IND(x,y)].y = - Vl[(int)Histo[(int)b_0[IND(x,y)].y].y].y;
		}
		POINT_F mn_ = {0,0}; for(int p = 0; p < m_nPixels; p++){ mn_.x +=  m_G0[p].x; mn_.y +=  m_G0[p].y;  } mn_.x /= m_nPixels; mn_.y /= m_nPixels;
		for(int p = 0; p < m_nPixels; p++){ m_G0[p].x += mn.x - mn_.x; m_G0[p].y += mn.y - mn_.y; }
//---------------------------------------------------------------------------------------------------------------------------------


delete [] Vl;
}
void   EDP:: Gr_I( float *in,  double   * mean, double * dsp )
{


//-------------------------------------------------

		for(int y = 0; y< m_height ; y++)
		for(int x = 0; x< m_width ; x++)
		{
		m_G0[IND(x,y)].m =0;
         if(x == m_width -1)m_G0[IND(x,y)].x =0;
		 else                           m_G0[IND(x,y)].x  = (in[IND(x+1,y)] - in[IND(x,y)]) ;
		 if(y == m_height -1)m_G0[IND(x,y)].y =0;
		 else                           m_G0[IND(x,y)].y  = (in[IND(x,y+1)] - in[IND(x,y)]);
		}
        mean[0] =0; for(int p = 0; p < m_nPixels; p++) mean[0] +=  in[p ]; mean[0] /= m_nPixels;
		dsp[0] =0;    for(int p = 0; p < m_nPixels; p++) dsp[0] += (mean[0] - in[p ]) *( mean[0] - in[p]); dsp[0] = sqrt(dsp[0]/m_nPixels);
		
}
void   EDP:: Gr_I_BF( float *in,  double   * mean, double * dsp, float sigma, float sigc )
{


//-------------------------------------------------

		for(int y = 0; y< m_height ; y++)
		for(int x = 0; x< m_width ; x++)
		{
		m_G0[IND(x,y)].m =0;
         if(x == m_width -1)m_G0[IND(x,y)].x =0;
		 else                           m_G0[IND(x,y)].x  =Gr_I_BF( in[IND(x+1,y)] ,  in[IND(x,y)],   sigma,  sigc );
		 if(y == m_height -1)m_G0[IND(x,y)].y =0;
		 else                           m_G0[IND(x,y)].y  = Gr_I_BF( in[IND(x,y + 1)] ,  in[IND(x,y)],   sigma,  sigc );
		}
        mean[0] =0; for(int p = 0; p < m_nPixels; p++) mean[0] +=  in[p ]; mean[0] /= m_nPixels;
		dsp[0] =0;    for(int p = 0; p < m_nPixels; p++) dsp[0] += (mean[0] - in[p ]) *( mean[0] - in[p]); dsp[0] = sqrt(dsp[0]/m_nPixels);
		
}
void   EDP:: Gr_I_BF__(unsigned char * I_r,  float *in,  double   * mean, double * dsp, float sigma, float sigc_ , int nch)
{
int wn = sigma*2.5;
int n_wr = 3000;
int thrr = 6;
float str = (float)n_wr/thrr;
float sigc = sigc_*255;
float sigc_c =1./(sigc*sigc*2);
float * rfl_x = new float [m_nPixels*nch];
float * rfl_y = new float [m_nPixels*nch];
float * ws = new float [wn+1];
float * wr = new float [n_wr+1];
for(int gi = 0; gi <=wn; gi++ )ws[gi] = exp(- (float)gi*gi/2/sigma/sigma);
for(int gi = 0; gi <=n_wr; gi++ )wr[gi] = exp(- (float)thrr*gi/n_wr);
        float * sum_x  =  new float [nch];
		float  *sum_y  =new float   [nch];
		for(int y = 0; y< m_height ; y++)
		for(int x = 0; x< m_width ; x++)
		{
		int ind = IND(x,y);
        float sum_w_x =0;
		float sum_w_y =0;
		for(int c =0; c< nch; c++){sum_x[c]  =0; sum_y[c]  =0;}
		for(int gi = - wn; gi <=wn; gi++ )
		{
			int ind_x =(gi+x>=0&&gi+x<m_width) ? IND(x+gi,y) : -1;
			int ind_y =(gi+y>=0&&gi+y<m_height) ? IND(x,y+gi) : -1;

			float dist_x  = 0, dist_y =0; for(int c = 0; c < 3; c++ )
			{if(ind_x+1) dist_x +=( I_r[c*m_nPixels + ind_x]-I_r[c*m_nPixels + ind])*( I_r[c*m_nPixels + ind_x]-I_r[c*m_nPixels + ind]);
			if(ind_y+1) dist_y +=( I_r[c*m_nPixels + ind_y]-I_r[c*m_nPixels + ind])*( I_r[c*m_nPixels + ind_y]-I_r[c*m_nPixels + ind]);}
	        dist_x *=sigc_c; dist_y *=sigc_c; if(dist_x>thrr)dist_x = thrr;  if(dist_y>thrr)dist_y = thrr;
         float w_x =  ws[abs(gi)]*wr[(int)(dist_x*str)];//exp(- (float)gi*gi/2/sigma/sigma)*exp( -dist_x/2/sigc/sigc);
		 float w_y = ws[abs(gi)]*wr[(int)(dist_y*str)];//exp(- (float)gi*gi/2/sigma/sigma)*exp( -dist_y/2/sigc/sigc);
		 if(ind_x+1) { for(int c =0; c< nch; c++){ sum_x[c] += in[ind_x + c*m_nPixels]*w_x;}  sum_w_x += w_x; }
		 if(ind_y+1) { for(int c =0; c< nch; c++) { sum_y[c] += in[ind_y + c*m_nPixels]*w_y;} sum_w_y += w_y; }
		}

		 for(int c =0; c< nch; c++){rfl_x[ind + c*m_nPixels] =  sum_x[c]/sum_w_x;
		 rfl_y[ind + c*m_nPixels] = sum_y[c]/sum_w_y;}
		}


//-------------------------------------------------
		float mll = 0.5;
        for(int c =0; c< nch; c++)
		for(int y = 0; y< m_height ; y++)
		for(int x = 0; x< m_width ; x++)
		{
		m_G0[IND(x,y)+ c*m_nPixels].m =0;
         if(x == m_width -1)m_G0[IND(x,y)+ c*m_nPixels].x =0;
		 else                           m_G0[IND(x,y)+ c*m_nPixels].x  =rfl_x[IND(x+1,y)+ c*m_nPixels] - rfl_x[IND(x,y)+ c*m_nPixels] ;
		 if(y == m_height -1)m_G0[IND(x,y)+ c*m_nPixels].y =0;
		 else                           m_G0[IND(x,y)+ c*m_nPixels].y  = rfl_y[IND(x,y+1)+ c*m_nPixels] - rfl_y[IND(x,y)+ c*m_nPixels] ;
		}
		for(int c =0; c< nch; c++){mean[c] =0; for(int p = 0; p < m_nPixels; p++) mean[c] += (rfl_x[p + c*m_nPixels]+rfl_y[p + c*m_nPixels])*mll; mean[c] /= m_nPixels; dsp[c] =0;}
		//for(int c =0; c< nch; c++){mean[c] =0; for(int p = 0; p < m_nPixels; p++) mean[c] += (rfl_x[p + c*m_nPixels]+rfl_y[p + c*m_nPixels])*mll; mean[c] /= m_nPixels; dsp[c] =0;}
		   /* for(int p = 0; p < m_nPixels; p++) dsp[0] += (mean[0] - in[p ]) *( mean[0] - in[p]); dsp[0] = sqrt(dsp[0]/m_nPixels);*/
		//for(int c =0; c< nch; c++){for(int p = 0; p < m_nPixels; p++) in[p + c*m_nPixels] = (rfl_x[p + c*m_nPixels]+rfl_x[p + c*m_nPixels])/2; }
delete [] rfl_x;		
delete [] rfl_y;	
delete [] ws;
delete [] wr;
delete [] sum_x;
delete [] sum_y;
}
inline float   EDP:: Gr_I_BF( float v,  float v0,  float sigma, float sigc )
{
	
//float w =  v - v0; w *=w;  w = exp(-(w/2/sigc/sigc) -((float)1./2/sigma/sigma) );
////return (((v + v0*w)-(v*w+v0))/(1+w));
//if(v>v0)return 1;
//else return -1;
float sg =-1;	if(v>v0) sg =1;
return log(fabs(v-v0)+1)*sg;
}
void   EDP:: FilterGr( unsigned char *out,  unsigned char *in)
{
	int r_w =2; float c_w = 1; int g_Ql = 2048;
	POINT_F* Histo = new POINT_F [g_Ql];
    m_G0 = new POINT_G [m_nPixels];
	POINT_F * b_0 = new POINT_F [m_nPixels];
	//--------------------------------------------------------------
    Gr(  r_w, c_w,  in,  b_0,  Histo, g_Ql );

	for(int c = 0; c < 3; c++)
    {
	double   mean,  dsp;
     Gr_I(&in[c*m_nPixels], b_0, Histo,  g_Ql,  &mean, &dsp );
    //---------------------------------------------------------------------------------------------------------------------------------

        MkRecFromG( &out[c*m_nPixels], mean, dsp);
		//for(int p = 0; p < m_nPixels; p++) {out[p + c*m_nPixels]  = abs(out[p + c*m_nPixels] - in[p + c*m_nPixels])*16;}
   }
//-----------------///
delete [] m_G0;
delete [] b_0;
delete [] Histo;
}
void   EDP:: BFilterGr( unsigned char * r_I, float *out, int nb_ch, float sigma, float sigc)
{

 m_G0 = new POINT_G [m_nPixels*nb_ch];
	//--------------------------------------------------------------
   double   * mean = new double [nb_ch]; double   * dsp = new double [nb_ch];
    Gr_I_BF__(r_I, out,  mean, dsp, sigma,sigc, nb_ch );
	for(int c = 0; c < nb_ch; c++) { if(c)    for(int p =0; p< m_nPixels; p++ )m_G0[p] = m_G0[p + c*m_nPixels];
		MkRecFromG( &out[c*m_nPixels], mean[c], 0);
	}
//-----------------///
delete [] m_G0;
delete [] mean;
delete [] dsp;
}
void EDP:: Gr( int r_w, float c_w,  unsigned char *in, POINT_F * b_0, POINT_F* Histo, int g_Ql )
{
//-------------------------------------------------
     for(int i = 0; i < g_Ql; i++) {Histo[i].x =0; Histo[i].y =0;}
	POINT_F * b_2 = new POINT_F [m_nPixels];
	Make_Inflc_Ptn(b_0,  in, r_w,  c_w, 0);
    Make_Inflc_Ptn(b_2,  in, r_w,  c_w, 2);
	float mm, pp;
    for(int x =0; x< m_width; x++ )
	for(int y =0; y< m_height; y++ )
	{
     if(x < m_width-1)
	 {
		 pp = b_0[IND(x,y)].x; mm = b_2[IND(x+1,y)].x;
	     b_0[IND(x,y)].x  = (pp>mm) ?(1. - pp)*(g_Ql-1): (1. - mm)*(g_Ql-1);
	 }
	 else b_0[IND(x,y)].x = 0 ;
	      if(y < m_height-1)
	 {
		 pp = b_0[IND(x,y)].y; mm = b_2[IND(x+1,y)].y;
	     b_0[IND(x,y)].y  = (pp>mm) ?(1. - pp)*(g_Ql-1): (1. - mm)*(g_Ql-1);
	 }
	else b_0[IND(x,y)].y = 0 ;
	}
	float mul = 1./(g_Ql-1);
	for(int p = 0; p < m_nPixels; p++) {Histo[(int)b_0[p].x].x += 1; Histo[(int)b_0[p].y].y += 1;}
    for(int p = 1; p < g_Ql; p++) {Histo[p].x += Histo[p-1].x; Histo[p].y += Histo[p-1].y;}
	for(int p = 0; p < g_Ql; p++) {Histo[p].x /= m_nPixels*mul; Histo[p].y /= m_nPixels*mul;}

//-----------------------

delete [] b_2;
}
void EDP::MkRecFromG(unsigned char *out_b, double out_m, double out_d)
{
	int minxy = (m_width < m_height) ? m_width : m_height;
	int sz = minxy; m_nRnk =0 ; while (sz){sz = (sz!=1)? sz/2 + sz%2 :0; m_nRnk++; }
	//-------------------------------------------------------------------
	im_x = new int [m_nRnk];
	im_y = new int [m_nRnk];
	im_sh = new int [m_nRnk];
    im_x[0] = m_width;
	im_y[0] = m_height;
	im_sh[0] = im_x[0]*im_y[0];
	for( int i = 1;  i < m_nRnk;  i++)
			{
			im_x[i]   =  im_x[i-1]/2 + im_x[i-1]%2 ;
			im_y[i]   =  im_y[i-1]/2 + im_y[i-1]%2 ;
			im_sh[i] =  im_x[i]*im_y[i];
			}
///////////////////////////////////////////////////////////////
	m_MD = new float [m_nPixels];
	m_H    = new float  [m_nPixels];
	m_G =   new POINT_G *[m_nRnk];
	for( int i = 0;  i < m_nRnk;  i++){m_G[i] = new POINT_G [im_sh[i]];}
/////////////////////////////////////////////////////////////
		for(int p = 0; p< m_nPixels; p++)m_G[0][p] = m_G0[p];
	int dx[3][2];
	int dy[2][3];
	for( int i = 1;  i < m_nRnk;  i++)
	{
		for(int y = 0; y< im_y[i] ; y++)
		for(int x = 0; x< im_x[i] ; x++)
		{
         int pnt = x + y*im_x[i];
		 int xx[3] = {x*2, x*2 + 1, x*2 +2};
		 int yy[3] = {y*2, y*2 + 1, y*2 +2};
         for(int p = 0; p < 3 ; p++)for(int q = 0; q < 2 ; q++) dx[p][q] = xx[p] + yy[q]*im_x[i-1];
		 for(int p = 0; p < 2 ; p++)for(int q = 0; q < 3 ; q++) dy[p][q] = xx[p] + yy[q]*im_x[i-1];
		 if(xx[2] >= im_x[i-1])
		 {
          dx[2][0] = -1;
		  dx[2][1] = -1;
		 }
		 if(yy[2] >= im_y[i-1])
		 {
          dy[0][2] = -1;
		  dy[1][2] = -1;
		 }
		 if(xx[1] >= im_x[i-1])
		 {
          dx[1][0] = -1;
		  dx[1][1] = -1;
		  dy[1][0] = -1;
		  dy[1][1] = -1;
		  dy[1][2] = -1;
		 }
		 if(yy[1] >= im_y[i-1])
		 {
          dy[0][1] = -1;
		  dy[1][1] = -1;
		  dx[0][1] = -1;
		  dx[1][1] = -1;
		  dx[2][1] = -1;
		 }
		 if(xx[1] >= im_x[i-1])
		 {
		  dy[1][0] =  dy[0][0];
		  dy[1][1] =  dy[0][1];
		  dy[1][2] =  dy[0][2];
		 }
		 if(yy[1] >= im_y[i-1])
		 {
		  dx[0][1] = dx[0][0];
		  dx[1][1] = dx[1][0];
		  dx[2][1] = dx[2][0];
		 }
///////////////////////////////////////
	if( yy[1]< im_y[i-1] && xx[1] < im_x[i-1] ) m_G[i][pnt].m  = ((m_G[i-1][dx[0][0]].x +  m_G[i-1][dx[0][0]].y) *3 +  m_G[i-1][dy[1][0]].y +  m_G[i-1][dx[0][1]].x )/8;
	else
	{
		  if( yy[1]>= im_y[i-1] && xx[1] < im_x[i-1] )  m_G[i][pnt].m  =m_G[i-1][dx[0][0]].x/2;
		  if( yy[1]< im_y[i-1] && xx[1] >= im_x[i-1] )  m_G[i][pnt].m  =m_G[i-1][dy[0][0]].y/2;
		  if( yy[1]>= im_y[i-1] && xx[1] >= im_x[i-1]) m_G[i][pnt].m = 0;
	}
////////////////////////////////////

			if(x == im_x[i]-1) m_G[i][pnt].x = 0;
			else
			{
			float sum = 0;
			for(int p = 0; p < 3 ; p++)
		    for(int q = 0; q < 2 ; q++)
			if(dx[p][q] != -1) {sum += (p == 1) ?  m_G[i-1][dx[p][q]].x*2 : m_G[i-1][dx[p][q]].x;}
            m_G[i][pnt].x = sum/4;
			}
         	if(y == im_y[i]-1) m_G[i][pnt].y = 0;
			else
			{
			float sum = 0;
			for(int p = 0; p < 2 ; p++)
		    for(int q = 0; q < 3 ; q++)
			if(dy[p][q] != -1){sum += (q == 1) ?  m_G[i-1][dy[p][q]].y*2 : m_G[i-1][dy[p][q]].y; }
            m_G[i][pnt].y = sum/4;
			}
		}
	}
//---------------------------------------------------------------
    if(im_x[m_nRnk-1] == im_y[m_nRnk-1]) m_H[0] = 0;
	if(im_x[m_nRnk-1] > im_y[m_nRnk-1]){ m_H[0] = 0; for(int p = 1; p< im_x[m_nRnk-1]; p++ )m_H[p] = m_H[p-1] + m_G[m_nRnk-1][p-1].x; }
	if(im_x[m_nRnk-1] < im_y[m_nRnk-1]){ m_H[0] = 0; for(int p= 1; p< im_y[m_nRnk-1]; p++ )m_H[p] = m_H[p-1] + m_G[m_nRnk-1][p-1].y; }


	for( int i = m_nRnk-2;  i >= 0;  i--)
	{
        for(int y = 0; y< im_y[i] ; y++)
		for(int x = 0; x< im_x[i] ; x++)
        m_MD[x + y *im_x[i]] = m_H[(x/2) + (y/2)*im_x[i+1]] - m_G[i+1][(x/2) + (y/2)*im_x[i+1]].m;

		for(int y = 0; y< im_y[i] ; y++)
		for(int x = 0; x< im_x[i] ; x++)
		{
	    int pnt = x + y*im_x[i];
        int  x__  =   x%2;
	    int  y__  =   y%2;
    //---------------------------------- 0 0
	if(!x__ && !y__){
            int cnt = 1;
			m_H[pnt] =  m_MD[pnt];
			if(x)
				{
				m_H[pnt]  += m_MD[pnt - 1] + m_G[i][pnt-2].x + m_G[i][pnt -1].x;
				cnt++;
				}
			if(y)
				{
				m_H[pnt]  += m_MD[pnt - im_x[i]] + m_G[i][pnt - 2*im_x[i]].y + m_G[i][pnt - im_x[i]].y;
				cnt++;
				}
             m_H[pnt] /= cnt;
		}
	//---------------------------------------10
		if( x__ && !y__){
			int cnt = 1;
			m_H[pnt] =  m_MD[pnt]  +m_G[i][pnt - 1].x;
		  if(x < im_x[i]-1)
				{
                m_H[pnt]  += m_MD[pnt + 1] - m_G[i][pnt].x;
				cnt ++;
				}
			if(y)
				{
				m_H[pnt]  += m_MD[pnt - im_x[i]] + m_G[i][(x-1) + (y-2)*im_x[i]].y + m_G[i][(x-1) + (y-1)*im_x[i]].x + m_G[i][pnt - im_x[i]].y;
				cnt ++;
				}
              m_H[pnt] /= cnt;
		}
		//------------------------------------01
		if( !x__ &&  y__){
			int cnt = 1;
			m_H[pnt] =  m_MD[pnt]  + m_G[i][pnt  - im_x[i]].y;
			if(x)
				{
				m_H[pnt]  += m_MD[pnt - 1] + m_G[i][(x-2) + (y-1)*im_x[i]].y + m_G[i][(x-2) + (y)*im_x[i]].x + m_G[i][pnt - 1].x;
				cnt++;
				}
			if(y<im_y[i]-1)
				{
				m_H[pnt]  += m_MD[pnt + im_x[i]] - m_G[i][pnt].y;
				cnt++;
				}
                m_H[pnt] /= cnt;
		}
		//---------------------------------------11
		if( x__ &&  y__){
			int cnt = 1;
			m_H[pnt] =  m_MD[pnt]  + ( m_G[i][(x-1) + (y-1)*im_x[i]].y  + m_G[i][(x-1) + (y)*im_x[i]].x);
			if(x < im_x[i] - 1)
				{
				m_H[pnt]  += m_MD[pnt +1] + m_G[i][(x+1) + (y-1)*im_x[i]].y - m_G[i][pnt].x;
				cnt ++;
				}
			if(y < im_y[i] - 1)
				{
				m_H[pnt]  += m_MD[pnt +im_x[i]] + m_G[i][(x-1) + (y+1)*im_x[i]].x - m_G[i][pnt].y;
				cnt ++;
				}
			    m_H[pnt] /= cnt;
		}
		}

	}
//---------------------------------------------------------------
    double in_mean =0, in_dsp = 0;
	for(int p = 0; p< im_sh[0] ; p++) in_mean += m_H[p];
	in_mean /= im_sh[0];
	for(int p = 0; p< im_sh[0] ; p++) in_dsp += (m_H[p]- in_mean) * (m_H[p]- in_mean);
    in_dsp = sqrt(in_dsp/ im_sh[0]);

	double  Bmd = (out_d)? out_d/in_dsp: 1.;
    double  Cmd = (out_m - in_mean*Bmd);


for( int p = 0; p < im_sh[0]; p++)
{
	int vl = round_fl((float)(Bmd*m_H[p] + Cmd));
	out_b[p] = (vl < 0 )? 0: ((vl >255) ? 255: vl);
}
//---------------------------------------------------------------
	 delete [] im_x ;
	 delete [] im_y ;
	 delete [] im_sh ;
	 delete [] m_MD;
	for( int i = 0;  i < m_nRnk;  i++) delete [] m_G[i];
	delete [] m_G;
}

void EDP::MkRecFromG(float *out_b, double out_m, double out_d)
{
	int minxy = (m_width < m_height) ? m_width : m_height;
	int sz = minxy; m_nRnk =0 ; while (sz){sz = (sz!=1)? sz/2 + sz%2 :0; m_nRnk++; }
	//-------------------------------------------------------------------
	im_x = new int [m_nRnk];
	im_y = new int [m_nRnk];
	im_sh = new int [m_nRnk];
    im_x[0] = m_width;
	im_y[0] = m_height;
	im_sh[0] = im_x[0]*im_y[0];
	for( int i = 1;  i < m_nRnk;  i++)
			{
			im_x[i]   =  im_x[i-1]/2 + im_x[i-1]%2 ;
			im_y[i]   =  im_y[i-1]/2 + im_y[i-1]%2 ;
			im_sh[i] =  im_x[i]*im_y[i];
			}
///////////////////////////////////////////////////////////////
	m_MD = new float [m_nPixels];
	m_H    = new float  [m_nPixels];
	m_G =   new POINT_G *[m_nRnk];
	for( int i = 0;  i < m_nRnk;  i++){m_G[i] = new POINT_G [im_sh[i]];}
/////////////////////////////////////////////////////////////
		for(int p = 0; p< m_nPixels; p++)m_G[0][p] = m_G0[p];
	int dx[3][2];
	int dy[2][3];
	for( int i = 1;  i < m_nRnk;  i++)
	{
		for(int y = 0; y< im_y[i] ; y++)
		for(int x = 0; x< im_x[i] ; x++)
		{
         int pnt = x + y*im_x[i];
		 int xx[3] = {x*2, x*2 + 1, x*2 +2};
		 int yy[3] = {y*2, y*2 + 1, y*2 +2};
         for(int p = 0; p < 3 ; p++)for(int q = 0; q < 2 ; q++) dx[p][q] = xx[p] + yy[q]*im_x[i-1];
		 for(int p = 0; p < 2 ; p++)for(int q = 0; q < 3 ; q++) dy[p][q] = xx[p] + yy[q]*im_x[i-1];
		 if(xx[2] >= im_x[i-1])
		 {
          dx[2][0] = -1;
		  dx[2][1] = -1;
		 }
		 if(yy[2] >= im_y[i-1])
		 {
          dy[0][2] = -1;
		  dy[1][2] = -1;
		 }
		 if(xx[1] >= im_x[i-1])
		 {
          dx[1][0] = -1;
		  dx[1][1] = -1;
		  dy[1][0] = -1;
		  dy[1][1] = -1;
		  dy[1][2] = -1;
		 }
		 if(yy[1] >= im_y[i-1])
		 {
          dy[0][1] = -1;
		  dy[1][1] = -1;
		  dx[0][1] = -1;
		  dx[1][1] = -1;
		  dx[2][1] = -1;
		 }
		 if(xx[1] >= im_x[i-1])
		 {
		  dy[1][0] =  dy[0][0];
		  dy[1][1] =  dy[0][1];
		  dy[1][2] =  dy[0][2];
		 }
		 if(yy[1] >= im_y[i-1])
		 {
		  dx[0][1] = dx[0][0];
		  dx[1][1] = dx[1][0];
		  dx[2][1] = dx[2][0];
		 }
///////////////////////////////////////
	if( yy[1]< im_y[i-1] && xx[1] < im_x[i-1] ) m_G[i][pnt].m  = ((m_G[i-1][dx[0][0]].x +  m_G[i-1][dx[0][0]].y) *3 +  m_G[i-1][dy[1][0]].y +  m_G[i-1][dx[0][1]].x )/8;
	else
	{
		  if( yy[1]>= im_y[i-1] && xx[1] < im_x[i-1] )  m_G[i][pnt].m  =m_G[i-1][dx[0][0]].x/2;
		  if( yy[1]< im_y[i-1] && xx[1] >= im_x[i-1] )  m_G[i][pnt].m  =m_G[i-1][dy[0][0]].y/2;
		  if( yy[1]>= im_y[i-1] && xx[1] >= im_x[i-1]) m_G[i][pnt].m = 0;
	}
////////////////////////////////////

			if(x == im_x[i]-1) m_G[i][pnt].x = 0;
			else
			{
			float sum = 0;
			for(int p = 0; p < 3 ; p++)
		    for(int q = 0; q < 2 ; q++)
			if(dx[p][q] != -1) {sum += (p == 1) ?  m_G[i-1][dx[p][q]].x*2 : m_G[i-1][dx[p][q]].x;}
            m_G[i][pnt].x = sum/4;
			}
         	if(y == im_y[i]-1) m_G[i][pnt].y = 0;
			else
			{
			float sum = 0;
			for(int p = 0; p < 2 ; p++)
		    for(int q = 0; q < 3 ; q++)
			if(dy[p][q] != -1){sum += (q == 1) ?  m_G[i-1][dy[p][q]].y*2 : m_G[i-1][dy[p][q]].y; }
            m_G[i][pnt].y = sum/4;
			}
		}
	}
//---------------------------------------------------------------
    if(im_x[m_nRnk-1] == im_y[m_nRnk-1]) m_H[0] = 0;
	if(im_x[m_nRnk-1] > im_y[m_nRnk-1]){ m_H[0] = 0; for(int p = 1; p< im_x[m_nRnk-1]; p++ )m_H[p] = m_H[p-1] + m_G[m_nRnk-1][p-1].x; }
	if(im_x[m_nRnk-1] < im_y[m_nRnk-1]){ m_H[0] = 0; for(int p= 1; p< im_y[m_nRnk-1]; p++ )m_H[p] = m_H[p-1] + m_G[m_nRnk-1][p-1].y; }


	for( int i = m_nRnk-2;  i >= 0;  i--)
	{
        for(int y = 0; y< im_y[i] ; y++)
		for(int x = 0; x< im_x[i] ; x++)
        m_MD[x + y *im_x[i]] = m_H[(x/2) + (y/2)*im_x[i+1]] - m_G[i+1][(x/2) + (y/2)*im_x[i+1]].m;

		for(int y = 0; y< im_y[i] ; y++)
		for(int x = 0; x< im_x[i] ; x++)
		{
	    int pnt = x + y*im_x[i];
        int  x__  =   x%2;
	    int  y__  =   y%2;
    //---------------------------------- 0 0
	if(!x__ && !y__){
            int cnt = 1;
			m_H[pnt] =  m_MD[pnt];
			if(x)
				{
				m_H[pnt]  += m_MD[pnt - 1] + m_G[i][pnt-2].x + m_G[i][pnt -1].x;
				cnt++;
				}
			if(y)
				{
				m_H[pnt]  += m_MD[pnt - im_x[i]] + m_G[i][pnt - 2*im_x[i]].y + m_G[i][pnt - im_x[i]].y;
				cnt++;
				}
             m_H[pnt] /= cnt;
		}
	//---------------------------------------10
		if( x__ && !y__){
			int cnt = 1;
			m_H[pnt] =  m_MD[pnt]  +m_G[i][pnt - 1].x;
		  if(x < im_x[i]-1)
				{
                m_H[pnt]  += m_MD[pnt + 1] - m_G[i][pnt].x;
				cnt ++;
				}
			if(y)
				{
				m_H[pnt]  += m_MD[pnt - im_x[i]] + m_G[i][(x-1) + (y-2)*im_x[i]].y + m_G[i][(x-1) + (y-1)*im_x[i]].x + m_G[i][pnt - im_x[i]].y;
				cnt ++;
				}
              m_H[pnt] /= cnt;
		}
		//------------------------------------01
		if( !x__ &&  y__){
			int cnt = 1;
			m_H[pnt] =  m_MD[pnt]  + m_G[i][pnt  - im_x[i]].y;
			if(x)
				{
				m_H[pnt]  += m_MD[pnt - 1] + m_G[i][(x-2) + (y-1)*im_x[i]].y + m_G[i][(x-2) + (y)*im_x[i]].x + m_G[i][pnt - 1].x;
				cnt++;
				}
			if(y<im_y[i]-1)
				{
				m_H[pnt]  += m_MD[pnt + im_x[i]] - m_G[i][pnt].y;
				cnt++;
				}
                m_H[pnt] /= cnt;
		}
		//---------------------------------------11
		if( x__ &&  y__){
			int cnt = 1;
			m_H[pnt] =  m_MD[pnt]  + ( m_G[i][(x-1) + (y-1)*im_x[i]].y  + m_G[i][(x-1) + (y)*im_x[i]].x);
			if(x < im_x[i] - 1)
				{
				m_H[pnt]  += m_MD[pnt +1] + m_G[i][(x+1) + (y-1)*im_x[i]].y - m_G[i][pnt].x;
				cnt ++;
				}
			if(y < im_y[i] - 1)
				{
				m_H[pnt]  += m_MD[pnt +im_x[i]] + m_G[i][(x-1) + (y+1)*im_x[i]].x - m_G[i][pnt].y;
				cnt ++;
				}
			    m_H[pnt] /= cnt;
		}
		}

	}
//---------------------------------------------------------------
    double in_mean =0, in_dsp = 0;
	for(int p = 0; p< im_sh[0] ; p++) in_mean += m_H[p];
	in_mean /= im_sh[0];
	for(int p = 0; p< im_sh[0] ; p++) in_dsp += (m_H[p]- in_mean) * (m_H[p]- in_mean);
    in_dsp = sqrt(in_dsp/ im_sh[0]);

	double  Bmd = (out_d)? out_d/in_dsp: 1.;
    double  Cmd = (out_m - in_mean*Bmd);
    double ddd = out_m - in_mean;

for( int p = 0; p < im_sh[0]; p++)
{
	//int vl = round_fl((float)(Bmd*m_H[p] + Cmd));
	//out_b[p] =  (vl < 0 )? 0: ((vl >255) ? 255: vl);
	out_b[p] = (float)(Bmd*m_H[p] + Cmd);
}
//---------------------------------------------------------------
	 delete [] im_x ;
	 delete [] im_y ;
	 delete [] im_sh ;
	 delete [] m_MD;
	for( int i = 0;  i < m_nRnk;  i++) delete [] m_G[i];
	delete [] m_G;
}

/*
	PS: I use "variance" freely.
	Function to calculate variance of disparities
	
	Also sets global variable sigma_img.

	@params IN dL - disparity per pixel. Dimensions= m_nPixels*1
	@returns ret - variance of disparity.


*/
double EDP::sigma_img_get( int *dL) {
	float mean =0;
	double ret;
	for(int p =0; p < m_nPixels; p++) {
		mean += dL[p];
	}
	mean /= m_nPixels;
    sigma_img = 0;
    for(int p =0; p < m_nPixels; p++) {
    	sigma_img += (mean-dL[p])*(mean - dL[p]);
    }
    ret = sigma_img /= m_nPixels;
    sigma_img =9.* (sigma_img)/m_nLabels/m_nLabels;
	return ret;
}
int EDP::FillStepClr(int p_n, int thr,  POINT* cr_out,  POINT* cr_in, BYTE* FillBf, POINT ij_img) {
	int ImgSh = m_nPixels;
	int ret_p_n=0;
	////////////////////////
	POINT  dxdy[4];
	dxdy[0].x = 1;  dxdy[0].y = 0;
	dxdy[1].x = 0;  dxdy[1].y =  1;
	dxdy[2].x = -1; dxdy[2].y = 0;
	dxdy[3].x = 0;  dxdy[3].y = -1;
	/////////////////////////
	for(int p = 0; p<p_n; p++)
	for(int k = 0; k<4;k++) {
		int i = cr_in[p].x + dxdy[k].x; int j = cr_in[p].y + dxdy[k].y;
		POINT ij = {i,j};
		if(i>=0&&i<ij_img.x &&j>=0&&j<ij_img.y) {//chk_nbrs
			int i_sh = i + j*ij_img.x ;
			int i_0  = cr_in[p].x + cr_in[p].y*ij_img.x;
			int st_msk = (WMask[i_sh])?  0: 1;
			//-----------------------------------------------
			int t_vl =0; /**/
			if(st_msk) {
				for(int c=0; c<3; c++, i_0 += ImgSh, i_sh += ImgSh) {
					t_vl += abs(FillBf[i_0]-FillBf[i_sh]);
				}
			}
			//---------------------------------------------
			if(t_vl<=thr&&st_msk) {
				WMask[i + j*ij_img.x]=IndCount;
				cr_out[ret_p_n++]=ij;
				Count++;
			}
		} /*chk_nbrs*/
	}
	return ret_p_n;
}

int EDP::FillFrPoinClr(int thr, POINT p, POINT**cr_in_out,  BYTE* FillBf, POINT ij_img) {
	Count=1;
	int s =0;
	WMask[p.y*ij_img.x+p.x]=IndCount;
	cr_in_out[s][0]= p;
	int prFl=1;
	//===============
	while((prFl=FillStepClr(prFl, thr, cr_in_out[(s+1)%2],cr_in_out[s], FillBf, ij_img))){
		s=(s+1)%2;
	}
	//===============
	ResInd[IndCount].x=Count;
	IndCount++;


	return Count;
}
void EDP::GetClrMask(int thr, BYTE* FillBf) {
	POINT ij_img = {m_width, m_height};
	int size = ij_img.x* ij_img.y;
	POINT *cr_in_out[2];
	///////////////////////////////
	for(int i=0;i<2;i++){
		cr_in_out[i] =new POINT [size];
	}


	memset(WMask, 0, size*sizeof(int));
	//----------------------------------------------------------------------------------------
	IndCount=1;
	POINT p;
	////////////////////////////////
	for(int i=0;i<size;i++) {
		p.x=i%ij_img.x;
		p.y=i/ij_img.x;
		if((!WMask[i]))
			FillFrPoinClr(thr, p, cr_in_out, FillBf, ij_img);
	}
	IndCount--;
	///////////////////////////////////////
	for(int i=0;i<2;i++){
		delete []  cr_in_out[i];
	}

	FOR_PX_p {
		WMask[p]--;
	}
}
void EDP::Make_Gr_fl( int thr, int thr2, int lr) {
	WMask = new int [m_nPixels];
	ResInd   = new POINT [m_nPixels];
	BYTE* FillBf = (!lr) ? I_ims[0] : I_ims[1];
	
	POINT_CH * msk_grp  = (!lr) ? B_L_buf_pp_ : B_R_buf_pp_;
	POINT_CH * msk_grm = (!lr) ? B_L_buf_mm_ : B_R_buf_mm_;
	POINT * dgr   = new POINT    [m_nPixels];
	//POINT * dgr2 = new POINT  [m_nPixels];
	//--------------------------------
	GetClrMask(thr, FillBf);
	for(int i =0; i < m_nPixels; i++ ) {
		int x = i%m_width; int y = i/m_width;
		int xpp = (x+1<m_width)? x+1:x; int ypp = (y+1<m_height)? y+1:y;
		int iypp = ypp*m_width + x;
		int ixpp = y*m_width + xpp;
		dgr[i].x = (WMask[i] - WMask[ixpp]) ? 1 : 0; dgr[i].y = (WMask[i] - WMask[iypp]) ? 1 : 0;
	}

	for(int i =0; i < m_nPixels; i++ ) {
		int x = i%m_width; int y = i/m_width;
		int xpp = (x+1<m_width)? x+1:x; int ypp = (y+1<m_height)? y+1:y;
		int xmm = (x > 0)? x-1 : x; int ymm = (y > 0)? y-1 : y;
		int iypp = ypp*m_width + x;
		int ixpp = y*m_width + xpp;
		int iymm = ymm*m_width + x;
		int ixmm = y*m_width + xmm;
		if(msk_grp[i].x && dgr[i].x && ResInd[WMask[i]].x > thr2) { 
			msk_grp[i].x = 0; msk_grm[ixpp].x = 0;
		}
		if(msk_grp[i].y && dgr[i].y && ResInd[WMask[i]].x > thr2) { 
			msk_grp[i].y = 0; msk_grm[ixpp].y = 0;
		}
	}
	delete [] dgr;
	//delete [] dgr2;
	delete [] WMask;
	delete [] ResInd;
}
void EDP::Make_Gr_fl_buf( int thr, int lr, double *cost, int L, float alp )
{
	WMask = new int [m_nPixels];
	ResInd   = new POINT [m_nPixels];
	BYTE* FillBf = (!lr) ? I_ims[0] : I_ims[1];
	GetClrMask(thr, FillBf);
	
    double * mvl = new double [IndCount];
	int  *cnt 	 = new int [IndCount];
	double * test= new double [m_nPixels];
	for(int l = 0; l < L; l++) {
		//------------------
		FOR_I(IndCount) mvl[i] = 0;
		FOR_I(IndCount) cnt[i] = 0;

		double * tb = &cost[l*N_PX];
		FOR_PX_p {test[p] = tb[p];}
		FOR_PX_p {
			int ind = WMask[p];
			mvl[ind] += test[p]; 
			cnt[ind]++;
		}
		FOR_I(IndCount) mvl[i] /= cnt[i];
		FOR_PX_p {
			tb[p] = tb[p]*(1-alp) + mvl[WMask[p]]*alp;
		}
	}//-----------------

	delete [] mvl;
	delete [] test;
	delete [] cnt;
	delete [] WMask;
	delete [] ResInd;
}

inline void  EDP::UpdMsgC_Z(int x, int y, float * Di_hat, float* M) {

	int K = m_nLabels;

	REAL  lambda = 1;// Cnst_1;//Smthnss[0].semi*5;
	REAL  smoothMax = 3; //Smthnss[0].max*5;
	{
		int k;
		REAL vlm,  vlp,   delta =   Di_hat[0]; 
		for (k=1; k<K; k++)
			TRUNCATE(delta,  Di_hat[k]);
		for (k=0; k<K; k++) {
			M[k] = Di_hat[k] - delta;
			if(k){
				vlm = Di_hat[k-1] - delta + lambda; 
				TRUNCATE(M[k],  vlm);
			}
			if(k<K-1){
				vlp = Di_hat[k+1] - delta + lambda; 
				TRUNCATE(M[k],  vlp);
			}
		}
		for (k=0; k<K; k++) 
			TRUNCATE(M[k],  smoothMax);
	}
//----------------------------------------

}

inline void  EDP::UpdMsgC_Z( double * Di_hat, double* M)
{

  int K = m_nLabels;

  REAL  lambda =  0.5/6;//Smthnss[0].semi*5;
  REAL  smoothMax = 0.5; //Smthnss[0].max*5;
{
     int k;
    REAL vlm,  vlp,   delta =   Di_hat[0]; for (k=1; k<K; k++)TRUNCATE(delta,  Di_hat[k]);
    for (k=0; k<K; k++)
	{
      M[k] = Di_hat[k] - delta;
	  if(k){vlm = Di_hat[k-1] - delta + lambda; TRUNCATE(M[k],  vlm);}
	  if(k<K-1){vlp = Di_hat[k+1] - delta + lambda; TRUNCATE(M[k],  vlp);}
	}
     for (k=0; k<K; k++) TRUNCATE(M[k],  smoothMax);

}
//----------------------------------------

}

inline void  EDP::UpdMsgC_Z_full( double * C,  double* M)
{
double  *CC = new double [m_nPixels*m_nLabels];
double  *MM = new double [m_nPixels*m_nLabels];
for(int p= 0; p<m_nPixels; p++)
for(int l= 0; l<m_nLabels; l++)
CC[p*m_nLabels+l]= C[p+l*m_nPixels];
//-------------------------------------------------
for(int p= 0; p<m_nPixels; p++)UpdMsgC_Z( &CC[p*m_nLabels], &MM[p*m_nLabels]);
//-------------------------------------------------
for(int p= 0; p<m_nPixels; p++)
for(int l= 0; l<m_nLabels; l++)
M[p+l*m_nPixels] = MM[p*m_nLabels+l];
  delete [] MM;
 delete [] CC;
}

inline void  EDP::UpdMsgC_F( float a, float * D, float* M, float*cM, float* W)
{
	
	int k, K =m_nLabels;
	cM[0] = D[0]; float m = exp(-a);M[K-1]= D[K-1];
 for (k=1; k<K; k++) cM[k] = D[k] + m*cM[k-1];
 for (k=K-2; k>=0; k--)M[k] = D[k] + m*M[k+1];
 for (k=0; k<K; k++)M[k] = (M[k]+cM[k])/W[k];

}


void EDP:: TRWS_CST_Z( float *m_cst )
{
	float a =1;
float* E = new float [m_nLabels];for(int l= 0; l<m_nLabels; l++)E[l]=1;
float* W = new float [m_nLabels];float* WW = new float [m_nLabels];
UpdMsgC_F(a, E,W,WW,E); delete [] WW;
////////////////////////////////////////////7
float * Mrg = new float [m_nLabels*m_nPixels];
float * cMrg = new float [m_nLabels*m_nPixels];
for(int p= 0; p<m_nPixels; p++)
for(int l= 0; l<m_nLabels; l++)
Mrg[p*m_nLabels+l]= m_cst[p+l*m_nPixels];

for(int y =0; y<m_height; y++)
for(int x =0; x<m_width; x ++)
{
UpdMsgC_F( a,
		  &Mrg[IND(x,y)*m_nLabels],
		  &cMrg[IND(x,y)*m_nLabels], E,W);
//UpdMsgC_Z( x, y,&Mrg[IND(x,y)*m_nLabels], &cMrg[IND(x,y)*m_nLabels]);
}// end dir 0 ////////////////
for(int p= 0; p<m_nPixels; p++)
for(int l= 0; l<m_nLabels; l++)
m_cst[p+l*m_nPixels]  = cMrg[p*m_nLabels+l];
//for(int p= m_width*100; p<m_width*101; p++)
//for(int l= 0; l<m_nLabels; l++)
//I_ims[2][p%m_width+l*m_width +m_nPixels]  = I_ims[2][p%m_width+l*m_width]=
//I_ims[2][p%m_width+l*m_width+ 2*m_nPixels]= m_cst[p+l*m_nPixels]*255;
//////////////////////
delete [] Mrg;
delete [] cMrg;
delete [] E;
delete [] W;

}
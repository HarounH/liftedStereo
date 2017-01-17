#ifndef __EDP_H__
#define __EDP_H__

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <png++/png.hpp>

typedef struct {
	int tau;
	short int d;
	unsigned char rgb[4];
	unsigned char lb; // 0 -
}
SLTN_CRV;


typedef struct {
	double  semi;
    double max;
}
DELTA_D;

typedef struct {
	float x;
    float y;
}
POINT_F;
typedef struct {
	float x;
    float y;
	float m;
}
POINT_G;
typedef struct {
	char x;
    char y;
}
POINT_CH;


class EDP {
public:
	typedef struct {
		int x;
		int y;
	} POINT;
	typedef struct { float cl[3]; } COLOR_F;
    typedef double REAL;
	typedef unsigned char BYTE;
	typedef int Label;
    typedef float CostVal;
	typedef struct {double  re; double im;} cmplx;
	EDP(int nmIm, unsigned char ** i_ims, int width, int height, int nLabels );
    ~EDP();



    // For general smoothness functions, this code tries to cache all function values in an array
    // for efficiency.  To prevent this, call the following function before calling initialize():
	void findSlt(float * dsp);
		void GR_filter( int thr,  unsigned char *in, double * fltb);
	//void  findSlt(float * dsp, char * c1, char * c2, float scl);
	int Size_fw;
	int Lt;
	int Sc_out;
protected:
	int nm_Ims;
	unsigned char **I_ims;
	int m_width;
	int m_height;
	int m_nLabels;
	int m_nPixels;
	
///////////////////////
	 inline double erf(double x);
      inline void drw();
	  inline double   UpdateMessageL1(REAL* M, REAL* Di_hat, int K, REAL gamma, CostVal lambda, CostVal smoothMax);
	   inline void CopyVector(REAL* to, CostVal* from, int K);
	   inline void  AddVector(REAL* to, REAL* from, int K);
	   inline double   SubtractMin(REAL *D, int K);
	  inline double    fi_erf(double ksi,  double sigc);
	  inline void    fi_erf_ini();
	  inline double    fi_erf_tab_f(double ksi, double vi, double sigc);
	  	inline double   mul_tab(int q);
	  inline double    fi_bf(double ksi, double vi, double sigc);
	  inline double    fi_erf(double ksi, double vi, double sigc);
	 inline void GaussConv( int N, double * inbuf, float * wt, int iW);
	 void  GaussConv2D(float sigma,int x, int y, float* buf);
	 inline void   SampledFlt_Add(float sigma_c_mul, float sigma, int cl_Q, unsigned char * I_b, float *buf_to_flt, short int * clr_cls,   COLOR_F * K_cls, int r, int nb_ch );
	 void  K_mean_Flt_Add(  BYTE * I_b,  float sigma_c_mul, float sigma, int cl_Q, float * buf_to_flt, int n_bufs) ;
	 //////////////////////////////////////
	   inline void    copy_m_D( float * m_cost, int dr);
	      inline void    copy_m_D( double * m_cost, int dr);
	      inline void    min_cost( float * m_cost, int * dsp);
		     inline void    med_cost( float * m_cost, int * dsp);
			   inline void    max_cost( float * m_cost, int * dsp);
			    inline void   prf_cost( float * m_cost);
				inline void    max_cost_do( float * m_cost);
	 //////////////////////////////////////
	void   Make_Gr_inf(int Ql, int thr1, int thr2);
    void Allocate();
    void initializeAlg();
	void  grBuf_( float* buf_out,  BYTE* buf_in);
	float * bufcc;

	void view_D(int st, double * cst);
		void view_D(int st, float * cst);
	void transform_D_(float thr);
	inline void MSG_UP_(float dsp_st, float dsp_fs,  REAL * outb,  REAL * outb_c,  REAL mul);
	int inline gss_wnd_ML(int ind, int min, int max,  int * im, int * buf0, int r_w, float * gss_c_wt,float * gss_r_wt, float *cst, float* Histo);
	void cost_subpix_(int r, float stp, int iRL);
	void Gss_wnd_ML( int r_w, float cl_prc, int * sup,  int * DspM, int o_s, int rg );
	inline void min_cost( double * m_cost, int * dsp, int lbs, int sc);
	
	void MRF__z(int itr, int * dL, int * dR);
	void modified_MRF__z(int itr, int * dL, int * dR, png::image< png::ga_pixel >& timg);
	void hybrid_MRF__z(int itr, int * dL, int * dR, png::image< png::ga_pixel >& truthimg);

	void CorrectSub_new(int *outL, int *outR, /*double * m_cost,*/ int r,  int sc);
	float inline cst_sub_pix(  float dsp0, int ind0,  BYTE * buf0,  BYTE * buf1);
	void Make_Gr_inf_pp(float thr_ar, int thr,  int thr1, int thr2);
	void Outlr_reg_filter(float thr_ar, int thr, int * rez, int lr, int r_w, float sig_c);
	//short int *K_ind;
  void  K_mean_Flt_Add_new(  BYTE * I_b,  float sigma_c_mul, float sigmax,float sigmay, int cl_Q, double * buf_to_flt, int n_bufs) ;
void     FilterBF_C( unsigned char *outc,float sig, float sigc, int cl_Q);
  inline void   SampledFlt_Add_new(float sigma_c_mul, float sigmaxy, int cl_Q, unsigned char * I_b, double *buf_to_flt, short int * K_ind,   COLOR_F * K_cls, int nb_ch );
void   Make_Gr_fl_buf( int thr, int lr, double *cost, int L, float alp )  ;	
  //---------------------------
  inline int   round_fl(double vl);
	inline void   UpdMsgC_F(float a, float * D, float* M, float*cM, float* W);
	int *        WMask;
	POINT *  ResInd;
	int    IndCount;
	int    Count;
	float sigma_img;
	void   comb_cst( double *out_cst, double ** in_cst, int lr);
	void   mk_BF_cost( float sigma_c_mul, float sigmax,float sigmay, int cl_Q, double **intr_cst, int n_bufs);
	double      sigma_img_get( int *dL);
	int    FillFrPoinClr(int thr, POINT p, POINT**cr_in_out,  BYTE* FillBf, POINT ij_img);
	void GetClrMask(int thr, BYTE* FillBf);
	int    FillStepClr(int p_n, int thr,  POINT* cr_out,  POINT* cr_in, BYTE* FillBf, POINT ij_img);
	void Make_Gr_fl( int thr, int thr2, int lr)  ;
	void BFfilter(int * inout, float sigma, float sig_cl, int m_q, int lr);
	//---------------------------
	inline void   UpdMsgC_Z( double * Di_hat, double* M);
	inline void   UpdMsgC_Z_full( double * C,  double* M);
	void  ITR_BF(int iti,   BYTE * I_b,  float sigma_c_mul, float sigmax,float sigmay, int cl_Q, double * C, int n_bufs);
   void  K_mean_Flt_Add_gr(  BYTE * I_b,  float sigma_c_mul, float sigmax,float sigmay, int cl_Q, double * buf_to_flt, int n_bufs) ;
   void     OneD_filter(unsigned char * RI, double * cst, int r, float sig_c);
	  void   BFfilter(unsigned char * inout, float sigma, float sig_cl, int m_q);
     inline void   UpdMsgC_Z(int x, int y, float * Di_hat, float* M);
	   void   TRWS_CST_Z( float *m_cst );
	void  Intermed(int *outL, int *outR);
	float  min_mean();
	 inline void   prf_cost_( float * m_cost);
	 void    med_fl( int r_w,  int r_c, unsigned char * inb);
	 void   transform_D(float thr,  float pw);
	void    BF_1D( int r_w, float cl_prc,int iLR, float * fb );
	float inline  BF_1D(int ind,   float * buf0, int r_w, float * gss_c_wt,float * gss_r_wt, int iRL);
	  inline void   prf_cost_inv_( float * m_cost);
void     FilterGrDsp( unsigned char * sup, int *outc);////GRADIENT
void     FilterGrC( unsigned char * sup, unsigned char *outc, int thr); ///grfilter;
void     FilterGrDsp( unsigned char * sup, int *outc, int thr);
void     FilterGrC( unsigned char * sup, double *outc, int thr);
	 inline double    mean_bf( float * m_cost, int pixs);
	 	  	  inline void    cross_chkL( float * m_cost, int * dL, int * dR);
			   inline void    cross_chkL_cst( float * m_cost, int * dL);
			   void  K_mean_Flt_Add_new_( int * I_b,  float sigma_c_mul, float sigmax,float sigmay, int cl_Q, double * buf_to_flt, int n_bufs) ;
	
	int inline  mask_BF_med(int ind,  unsigned char * msk, int * buf0, int r_w, float * gss_c_wt,  float * gss_r_wt, int iRL);
	 void    Make_Gr_inf ();
	float *adapX;
	float *adapY;
	float adapx, adapy;
	inline void    min_cost( double * m_cost, int * dsp);
	    void   adap_d( double *m_cst, int *DspM, int lr);
	void    get_std( double cl_r, unsigned char * buf,  double * gss_c_wt, int n_c );
	inline void   SampledFlt_Add(float sigma_c_mul, float sigma, int cl_Q, unsigned char * I_b, double *buf_to_flt, short int * clr_cls,   COLOR_F * K_cls, int r, int nb_ch );
	void  K_mean_Flt_Add(  BYTE * I_b,  float sigma_c_mul, float sigmax,float sigmay, int cl_Q, double * buf_to_flt, int n_bufs);
	void     OneD_filterX(unsigned char * RI, float * cst, int r, float sig_c);
	void     OneD_filter(unsigned char * RI, float * cst, int r, float sig_c);
	 inline float   mean_min_sm( float * m_cost);
	float  L_occ(int *outL );
	 int  Gss_wnd_Z( int r_w, float cl_prc, int * DspM, int iLR, int o_s, float div );
	 void inline  Histo_min_max_Z(int ind, int *min, int *max,  int * buf0, int r_w, float* Histo, int Z);
	int inline  gss_wnd_Z(int ind, int min, int max,  int * buf0, int r_w, float * gss_c_wt,float * gss_r_wt, int iRL, float* Histo, float thrw);
	void   TRWS_CST(int itr,   REAL * Mrg, double *m_cst, int *DspM, int lr);
	void   TRWS_CST(int itr,   REAL * Mrg, float *m_cst, int *DspM, int lr);
	 inline void  addMsgC( REAL* M,  REAL* Mrg, int K);
	 inline void   MkMrgC(int x, int y,  REAL* Msg,  REAL* Mrg);
	 inline void   UpdMsgC(int x, int y, REAL* Di,  REAL* Msg,  REAL* Mrg, int dir);
	 inline void   subMsgC( REAL* Di,  REAL* M,  REAL* Mrg,int K,  REAL gamma);
	 inline  REAL    UpdMsgCL1(  REAL* M,  REAL* Di_hat, int K,  REAL lambda,  REAL smoothMax);
	//--------------------------
	 inline void   prf_cost_inv( float * m_cost);
	inline double    copy_m_D_b( float * m_cost, int dr);
	void  BLF_spr(  BYTE * I_b,  float sigma_c_mul, float sigma,  float * buf_to_flt, int n_bufs) ;
	inline int   Gss_rand( float sk, int r_w);
	void   dev_map3x3( BYTE * l_b,  float sigma_c_mul, float sigma, short int * map,  float * mul_buf) ;
	void  MakeMK_mask(int nsp, int nmw, float sk, int r_w, float ** crd_w, POINT ** crd);
	inline void   SampledFlt__(float sigma_c_mul, float sigma, int cl_Q, unsigned char * I_b, float *buf_to_flt, short int * clr_cls,   COLOR_F * K_cls, int r, int nb_ch );
	float   dev_map(  BYTE * I_b,  float sigma_c_mul, float sigma, short int * map,  float * mul_buf);
	void  GaussCosConv2DFst(float sigma, int iq, short int * map,  float* lmbK, int x, int y, float* buf);
	void  K_mean_Flt_cl(  BYTE * I_b,  float sigma_c_mul, float sigma, int cl_Q, float * buf_to_flt) ;
	void  GaussCosConv2DFst(float sigma,int x, int y, float* buf);
	REAL getErrGTI(float *in, float * out);
	void PrjectSlt(int st_intrp, int *outL, int *outR);
	void MkDspCorr_RL( int *inL, int *inR);
	void  Gss_wnd_( int r_w, float cl_prc, int * DspM, int iLR, int o_s );
	int inline gss_wnd_(int ind,  int min, int max, int * buf0, int r_w, float * gss_c_wt,float * gss_r_wt, int iRL, float* Histo);
  void  inline color_med(int ind, int r_w, unsigned char * out_vl, unsigned char * im );
  void  inline color_med_( int r_w,unsigned char * im );
  void inline  Histo_min_max(int ind, int *min, int *max,  int * buf0, int r_w, float* Histo);
  void Mk_DRL(int * DspL,  int *DspR, int thr);
  void  MskRL_Corr( int *in,int iRL);
  void  filterDSI_tau(unsigned char *im1,    int rw,     float cl_fc);
  void  SubpixelAcc(int itr, float strt_dsp, float pw2mul, int depth, int iRL);
  void  get_dsp_subpix(int r, float st, int itr, int iRL);
  void  cost_subpix(int r, float stp, int iRL);
  void  inline cst_subpix_pix(  float *out_buf, int ind0, BYTE * buf0, BYTE * buf1);
  void  get_dsp_sbpx_mrf(  int iti);
  void dsp_sp_mrf(int itr, REAL **ms);
  inline void   MSG_UP(float dsp_st, float dsp_fs,  REAL * outb,  REAL * outb_c,  REAL mul);
  inline int   Get_min( REAL *buf);
  void    Add_DspFL( REAL *tmp_vl,  REAL ** ms);
  void   CorrectSub(int *outL, int *outR, int dpth, int sc);
  void  grBuf( float* buf_out,  BYTE* buf_in);
  inline float   Get_mul_gr( float *buf, BYTE* bufb);
  void  inline   cst_subpix_pixG(  float *out_buf, int ind0,  BYTE * buf0,  BYTE * buf1);
  void GrBuf();
  void  Mk_DG(int * DspL);
  void GT_gnr(int *outL, int *outR);
  void MkPnlCst_G(int *LMp);
  void PrjectSltG( int *outL);
  void ProjectImCrvNew(int yy, unsigned char *img, float dir);
  void ProjectImCrvNew_(int yy, unsigned char *img, float dir);
  void Make_G_lw_hi_( unsigned char * G_buf,  unsigned char * I_buf, int dir, float thr);
  void Make_Inflc_Ptn( POINT_F * fltb,  unsigned char * I_buf, int r_w, float cl_prc, int dir);
  float inline  BiF_wnd_( int ind, int dir,   unsigned char * im, int r_w, float * gss_c_wt, float * gss_r_wt);
  double  Cost_mean_vl();
  double  Mean_tr(double * max);
  void MRF__(int itr, int * dL, int * dR);

  void MkDspCorr_pp(int *inL, int *inR, int iLR);
  void   Gss_wnd_pp( int r_w, float cl_prc, int * DspM, int iLR, int o_s );
  void inline  Histo_min_max_pp(int ind, int *min, int *max,  int * buf0, int r_w, float* Histo, int *map);
  int inline  gss_wnd_pp(int ind, int min, int max,  int * buf0, int r_w, float * gss_c_wt,float * gss_r_wt, int iRL, float* Histo, int * map);
  ///-------------------------------------------------------------
  void GaussCosConv2D_Clr(float sigma, float sig_cl, int x, int y, float* buf, unsigned char* I_b);
  int  get_lb(int p, float * buf_lb, int N);
  void  getCst_lft(int n_lb, float * buf_lb, float thr);
  float inline gss_wnd_c(int ind,    unsigned char * im, int r_w, float * gss_c_wt,float * gss_r_wt, float *buf_f);
 private:
	  double sq2;
	  double mPi;
	  double * fi_erf_tab;
	  int fi_erf_q;
	 float INFL_TR;

    CostVal m_smoothMax; // used only if
    CostVal m_lambda;    // m_type == L1 or m_type == L2
	/////////////////////////////////////////////
	double *Bi_ws;
	double *TBi_ws;
	short int *TBi_p0;
	int sigmaix;
	int sigmai2;
	int sigmaiy;
	int sigmai;
	int sigmaX;
	int sigmaY;
	double Lmb_Y;
	int Col_outl_thr;
	float Lmbd_trn;
	int ST_Cub;
	int chTr;
	int subpix_cnst;
	float subpix_cst_tune;
	float * fine_BC[4];
	int * Dst_B;
	float * DspFL;
	int Lbs_S;
	float Lbs_st;
	float * Lbs_S_map;
	float * spx_D;
    short int  * XL_tau;
	short int  * XR_tau;
	int *tau_XL_d;
	int *tau_XR_d;
	int *V_tau_bck;
	int *R_tau_bck;
	int *L_tau_bck;
	int *V_tau_frw;
	int *R_tau_frw;
	int *L_tau_frw;
	int Tau;
	DELTA_D *Smthnss;
	float  * Cost_map;
	int WpH;

	int R_INF_W;
	float CL_INF_R;
    Label *m_answer;
    Label *gtL_answer;
	Label *gtR_answer;
	CostVal *m_D;
	float * m_D_pR;
	float * m_D_pL;
	float PCst_mlt;
	int Itr_num;
	char Trn_n;
	REAL Lmbd_y;
	REAL Pnlt_y_max;
	float C_app_1;
	float C_app_2;
	float C_dsp_1;
	float C_dsp_2;
	float C_gr_1;
	float C_gr_2;
	REAL Cost_mean;
	int Cost_max;
	float Gr_prc;
	unsigned char * Mask_RL[2];
	unsigned char *IG_ims[2];
	POINT_F * B_L_buf_pp;
	POINT_F * B_L_buf_mm;
	POINT_F * B_R_buf_pp;
	POINT_F * B_R_buf_mm;
	POINT_CH * B_L_buf_pp_;
	POINT_CH * B_L_buf_mm_;
	POINT_CH * B_R_buf_pp_;
	POINT_CH * B_R_buf_mm_;
	REAL * Msg_X_pp;
	REAL * Msg_X_mm;
	REAL * Mrg_bf;
	int * T_max_y;
	SLTN_CRV  * Sltn_crv;
	SLTN_CRV  * Sltn_crv_tst;
	///////////////////////////

    REAL* m_messages;
	REAL   computeDSI_tau(unsigned char *im1,   unsigned char *im2,   int birchfield,  int squaredDiffs, int truncDiffs);
	REAL   computeDSI_tauG(unsigned char *im1,   unsigned char *im2,   int birchfield,  int squaredDiffs, int truncDiffs, float alph);
	///////////////////////////////
	REAL  findS(int y);
	REAL  findSb(int y);
	inline  REAL  getErrGT();
	void  XL_XR_Tau();
	void   drwCost(int yy);
	void  PaintCrv(int yy );
	void  Make_G_lw_hi( unsigned char * G_buf,  unsigned char * I_buf, int dir, float thr);
	/* void computeDSI_Pnt_tau();     */
	 //REAL findMrg(int y,  int dir);
	 inline int   round_fl(float vl);
	 void ProjectImCrv(int yy, unsigned char *img, float dir);
	 void  ProjectImCrvB(int yy, unsigned char *img, float dir);
	 void  ProjectDCrv(int yy, int dir,int *out );
	 /*REAL   computeDSI_tauGT( int truncDiffs, int st);  */
	 REAL getErrGTI();
	 void  MkPnlCst_RL(int *LMp, int *RMp);
	 REAL   TRWS_R(int itr,   REAL * Mrgnl, int * DspM);
	 inline void  MkMrgR(int x, int y, REAL* Msg, REAL* Msg_add, REAL* Mrg);
	 inline void  UpdMsgYR_mm(int x, int y, REAL* Di,REAL* Msg, REAL* Msg_add,  REAL* Mrg);
	 inline void  UpdMsgYR_pp(int x, int y, REAL* Di,REAL* Msg, REAL* Msg_add,  REAL* Mrg);
	 inline void  UpdMsgXR_pp(int x, int y, REAL* Di,REAL* Msg, REAL* Msg_add,  REAL* Mrg);
	 inline void  UpdMsgXR_mm(int x, int y, REAL* Di,REAL* Msg, REAL* Msg_add,  REAL* Mrg);
	 ///////////////////////////////////////////////////////////////////
	 REAL   TRWS_L(int itr,   REAL * Mrgnl, int * DspM);
	 inline void  MkMrgL(int x, int y, REAL* Msg, REAL* Msg_add, REAL* Mrg);
	 inline void  UpdMsgYL_mm(int x, int y, REAL* Di,REAL* Msg, REAL* Msg_add,  REAL* Mrg);
	 inline void  UpdMsgYL_pp(int x, int y, REAL* Di,REAL* Msg, REAL* Msg_add,  REAL* Mrg);
	 inline void  UpdMsgXL_pp(int x, int y, REAL* Di,REAL* Msg, REAL* Msg_add,  REAL* Mrg);
	 inline void  UpdMsgXL_mm(int x, int y, REAL* Di,REAL* Msg, REAL* Msg_add,  REAL* Mrg);
	 inline void  subMsg(int * xr,int * xsh, REAL* Di,REAL* M, REAL* Mrg,int K, REAL gamma);
	 inline void  subMsg_(int * xr, int * xsh, REAL* Di,  REAL* M,  REAL* Mrg,int Km,REAL add, REAL gamma);
	 inline void  addMsg(int * xr, int * xsh,REAL* M, REAL* Mrg, int K);
	 inline void  addMsg_(int * xr, int * xsh,  REAL* M,  REAL* Mrg,int Km, REAL add);
	 inline REAL  UpdMsgL1(int * xr,  REAL* M,  REAL* Di_hat, int K, REAL lambda, REAL smoothMax);
	 inline REAL  UpdMsgL1_(int * xr,  REAL* M,  REAL* Di_hat, REAL* add,int K, REAL lambda, REAL smoothMax);
	 ///////////////////////////////////
	 void  Make_G_fromBiF( float * fltb,  unsigned char * I_buf, int r, float cl_r, int dir);
	void inline  BiF_wnd( float *ret, int ind, int dir,   unsigned char * im, int r_w, float * gss_c_wt, float * gss_r_wt);
	//////////////////////// GR //////////////////////////////
	int * im_sh;
	int * im_x;
	int * im_y;
	int    m_nRnk;
	float * m_H;
	float * m_MD;
	POINT_G ** m_G;
	POINT_G *m_G0;
	//---------------------------------
	void ROTH_lb(int n_lb, float * buf_lb, unsigned char * I_b, int * answer);
    void  ROTH_lb(int n_lb, float * buf_lb, unsigned char * I_b);
	void K_means(int cl_Q, COLOR_F * K_cls, int r, int * buf_to_sort, short int * clr_cls);
	void SetSort(int n_p, int * buf_to_sort, int * ind_out);
	inline void SetSort_n(int * i_exp2, int n,  int N, int* vl_buf_out, int* ind_buf_out, int * vl_buf_in, int* ind_buf_in);
	////////////////////////////GAUSS FAST /////////////////////////////////////////
	double * gss_vl_buf;
	double alpha_gauss_cos_opt;
	double A_gauss_cos;
	double B_gauss_cos;
	double beta_gauss_cos;
	double gamma_gauss_cos;
	double teta_gauss_cos;
	cmplx * cmplxBuf;
	cmplx * cmplxBuf2;
	double *   wght_gauss_cos_buf;
	cmplx * cmplxCosSin;
	int         iW_gauss_cos;
	double *  dirBuf;
	double *  invBuf;
  //------------------------------
	void K_mean_Flt( BYTE * I_b,  float sigma_c_mul, float sigma, int cl_Q, float * buf_to_flt, int n_bufs) ;
	inline void  SampledFlt(float sigma_c_mul, float sigma, int cl_Q, unsigned char * I_b, float *buf_to_flt, short int * clr_cls,   COLOR_F * K_cls, int r, int nb_ch );
	void  get_std( float cl_r, unsigned char * buf,  float * gss_c_wt, int n_c );
	inline cmplx fromBufN(int i, int iN, cmplx* buf);
	inline float fromBufN(int i, int iN, double * buf);
	inline void MuFst( int iN, double * inbuf);
	inline cmplx mulCmplx( cmplx a, cmplx b);
	inline cmplx  addCmplx( cmplx a, cmplx b);
	inline   cmplx  expTet( int i);
	inline double  CmplxExpFst( double teta, int N,  cmplx * inbuf);
	inline void       GaussCosConvFst(int N,  double *InBuf);
	inline void       GaussCosConvFst_( int N, double * inbuf);
	void                 GaussCosConv2DFst(float sigma,int x, int y, float* buf, float *mask);
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
	void   Gr_I(unsigned char *in, POINT_F * b_0, POINT_F* Histo, int g_Ql, double   * mean, double * dsp );
    void  Gr( int r_w, float c_w,  unsigned char *in, POINT_F * b_0, POINT_F* Histo, int g_Ql );
	void   FilterGr( unsigned char *out,  unsigned char *in);
	void  MkRecFromG(unsigned char *out_b, double  out_m, double out_d);
	void  MkRecFromG(float *out_b, double out_m, double out_d);
	void    Gr_I( float *in,  double   * mean, double * dsp );
	void    BFilterGr( unsigned char * r_I, float *out, int nb_ch, float sigma, float sigc);
	void    Gr_I_BF( float *in,  double   * mean, double * dsp, float sigma, float sigc );
	inline float    Gr_I_BF( float v,  float v0,  float sigma, float sigc );
	void    Gr_I_BF__(unsigned char * I_r,  float *in,  double   * mean, double * dsp, float sigma, float sigc, int nch );
};

#endif /*  __EDP_H__ */

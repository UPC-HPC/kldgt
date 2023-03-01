/***************************************************************************
   Copyright (c) CGGVERITAS SINGAPORE
   Created by Yang Kunlun on Nov. 30, 2014.
* **************************************************************************/
#include "msdginterp.h"
#include <math.h>
#include <fftw3.h>
#include <PFL_C.h>
#include <complex.h>
#include <mkl.h>
#include "invmatrix.h"
#include "utility.h"

/**
 * @brief compute dot product of solo and dual vector
 * */
static inline double dotp2d_solo(int m, int n, float ** restrict a)
{
    double sum = 0.0f;
    for (int i=0; i<m; i++)
      for (int j=0; j<n; j++)
        sum += a[i][j]*a[i][j];
    return sum;
}
static inline double dotp2d_dual(int m, int n, float ** restrict a, float ** restrict b)
{
    double sum = 0.0f;
    for (int i=0; i<m; i++)
      for (int j=0; j<n; j++)
        sum += a[i][j]*b[i][j];
    return sum;
}

#define MSDGINTERP_WARP_DISPATCH __declspec(cpu_specific(core_2nd_gen_avx))
#include "l1inv_lsqr.h"
#undef MSDGINTERP_WARP_DISPATCH

#define MSDGINTERP_WARP_DISPATCH __declspec(cpu_specific(core_i7_sse4_2))
#include "l1inv_lsqr.h"
#undef MSDGINTERP_WARP_DISPATCH

#define MSDGINTERP_WARP_DISPATCH __declspec(cpu_dispatch(core_i7_sse4_2, core_2nd_gen_avx))
MSDGINTERP_WARP_DISPATCH inline static void 
taup_warp_dual_p(const int nx, const int np, const int hfid, float complex **pfpT, 
        float complex **AP, float complex *xPrimary, float *preflectivity, 
        float complex **AG, float complex *xGhost, float complex **pfx) {};

MSDGINTERP_WARP_DISPATCH
static void instantaneous_weight(void *ptask, const float * restrict hw, int nt, int nr, int np,
                                 float ** restrict inp, float ** restrict weight,
                                 float spower, float * restrict conv) {}; 
MSDGINTERP_WARP_DISPATCH
static void c_l1inv_TAUP_FW_DG(int nsgtrace, int hfid, int lfid, int npx, 
				  float **TX_p, float **TX_az, float **TX_ay, float ** taup, l1inv_t *l1para) {};
MSDGINTERP_WARP_DISPATCH
static void c_l1inv_TAUP_FW_CMSTOS_STOS_p(int imode, int nsgtrace, int hfid, int lfid, int npx, 
				 float **TX, float ** taup, l1inv_t *l1para) {};
MSDGINTERP_WARP_DISPATCH
static void c_l1inv_TAUP_FW_CMSTOS_STOS_z(int imode, int nsgtrace, int hfid, int lfid, int npx, 
				 float **TX, float ** taup, l1inv_t *l1para) {};
MSDGINTERP_WARP_DISPATCH
static void c_l1inv_matrix_TAUP_FW_RV_CMSTOS_STOS_p(int imode, int nsgtrace, int hfid, int npx, 
					   float **taup, float **taupou, l1inv_t *l1para) {};
MSDGINTERP_WARP_DISPATCH
static void c_l1inv_matrix_TAUP_FW_RV_CMSTOS_STOS_z(int imode, int nsgtrace, int hfid, int npx, 
					   float **taup, float **taupou, l1inv_t *l1para) {};
MSDGINTERP_WARP_DISPATCH
static void c_l1inv_matrix_TAUP_FW_RV_DG(int nsgtrace, int hfid, int npx, float **taup, 
                                         float **taupou, l1inv_t *l1para) {};
MSDGINTERP_WARP_DISPATCH
void c_l1inv_TAUP_RV_STOS_MS(int nsgtrace, int hfid, int lfid, int npx, float ** taup, 
			  float **RFX,  l1inv_t *l1para, int datamode, int datum, int outtype) {};
MSDGINTERP_WARP_DISPATCH
void c_l1inv_TAUP_RV_STOS(int nsgtrace, int hfid, int lfid, int npx, float ** taup, 
			  float **RFX,  l1inv_t *l1para, int datamode, int datum) {};
MSDGINTERP_WARP_DISPATCH
void c_l1inv_TAUP_RV_CTOS(int nsgtrace, int hfid, int lfid, int npx, float ** taup, 
			  float **RFX,  l1inv_t *l1para, int datamode, int datum) {};
MSDGINTERP_WARP_DISPATCH
void c_l1inv_TAUP_RV_BASIC(int nsgtrace, int hfid, int lfid, int npx, float ** taup, 
			   float **RFX,  float ** TPP, l1inv_t *l1para, int datum) {};
MSDGINTERP_WARP_DISPATCH
void c_l1inv_TAUP_RV_MTOS(int nsgtrace, int hfid, int lfid, int npx, float ** taup, 
			     float **RFX,  l1inv_t *l1para, int datamode, int datum) {};
MSDGINTERP_WARP_DISPATCH
void c_l1inv_lsqr(int imode, float repstop, int niter, 
		       int nsgtrace, int hfid, int lfid, int npx, float sppow, 
		       float ** pweight, float ** dataxt_p, float ** dataxt_az, float ** dataxt_ay,
		       float **taup, float **tweight, l1inv_t *l1para) {};
MSDGINTERP_WARP_DISPATCH
void c_l1inv_lsqr_cg(int imode, float repstop, int niter, 
		       int nsgtrace, int hfid, int lfid, int npx, float sppow, 
		       float ** pweight, float ** dataxt_p, float ** dataxt_az, float ** dataxt_ay,
		       float **taup, float **tweight, l1inv_t *l1para) {};
MSDGINTERP_WARP_DISPATCH
void c_l1inv_lsqr_bcg(int imode, float repstop, int niter, 
		       int nsgtrace, int hfid, int lfid, int npx, float sppow, 
		       float ** pweight, float ** dataxt_p, float ** dataxt_az, float ** dataxt_ay,
		       float **taup, float **tweight, l1inv_t *l1para) {};
MSDGINTERP_WARP_DISPATCH
void ms_lowcutf(float **input, int nsgtrace, l1inv_t *l1para, int datamode) {};
MSDGINTERP_WARP_DISPATCH
void l1inv3d_taup_merge(l1inv_t *l1para, int hfid, int npx, float **taup, float **taup1) {};

MSDGINTERP_WARP_DISPATCH
void taup_warp_solo_dsdg3d(int nsgtrace, int npx, int hfid, float **pfpT0in, float **pfpT1in, 
			   float **TG, float *xGhost, float *obliqcorr, float *recz, float **pfx, l1inv_t *l1para) {};
MSDGINTERP_WARP_DISPATCH
void c_dsdg3d_TAUP_RV_CTOS(int nsgtrace, int hfid, int lfid, int npx, float *recz, float ** taup_p, float ** taup_az,
			   float **RFX,  l1inv_t *l1para, int datamode, int datum) {};
MSDGINTERP_WARP_DISPATCH
void c_dsdg3d_TAUP_RV_MTOS(int nsgtrace, int hfid, int lfid, int npx, float *recz, float ** taup_p, float ** taup_az,
			  float **RFX,  l1inv_t *l1para, int datamode) {};
MSDGINTERP_WARP_DISPATCH
float compute_scalar(int ip, int ifreq, float recz, l1inv_t *l1para) {};

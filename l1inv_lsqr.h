//#include "test1.h"
/***************************************************************************
  Copyright (c) CGGVERITAS SINGAPORE
  Created by Yang Kunlun on Nov. 30, 2014.
 * **************************************************************************/
/*
 * @brief kernel of taup transform, forward and reverse all use it.
 * @param nx            [R] number of spatial trace
 * @param np            [R] number of p trace
 * @param hfid          [R] high bound of frequency component
 * @param reflectivity  [R] usually -1.0f for free surface
 * @param pfpT          [R] Pointer of Frequency Plan wave trace Transposed
 * @param TP            [R] phase table for primary
 * @param yPrimary      [S] temporal primary phase spin variable
 * @param TG            [R] phase table for ghost
 * @param yGhost        [S] temporal ghost phase spin variable
 * @param ztop          [R] ?
 * @param pfx           [W] pointer to frequency spatial matrix
 *
 * Warp core for all kinds of taup transform.
 * The deghost is actually a dual redatum from surface to cable and mirror cable.
 * solo is a total subset of dual, let's keep exactly that way.
 */

#define taup_warp_solo(nx, np, hfid, pfpT, TP, yPrimary, ztop, pfx)                 \
{ for(int ix=0; ix<nx; ix++) {                                                      \
	float complex tmp;                                                              \
	for(int ip=0; ip<np; ip++)                                                      \
	yPrimary[ip] = 1.0f*ztop[ip];                                                 \
	for(int ifreq=0; ifreq<hfid; ifreq++) {                                       \
		tmp = 0.0f;                                                               \
		for(int ip=0; ip<np; ip++) {                                              \
			tmp += yPrimary[ip]*pfpT[ifreq][ip];                                  \
			yPrimary[ip] *= TP[ix][ip];                                           \
		} pfx[ix][ifreq] = tmp;                                                   \
	} } }
 
    MSDGINTERP_WARP_DISPATCH 
__attribute__((always_inline))   
//        float *restrict preflectivity, float complex **restrict  AG, 
//        float complex *restrict xGhost, float complex **restrict pfx)
//{
 inline static void 
taup_warp_dual_p(const int nx, const int np, const int hfid, float complex **restrict pfpT, 
        float complex **restrict   AP, float complex *restrict xPrimary, 
        float *restrict  preflectivity, float complex **restrict   AG, 
        float complex *restrict xGhost, float complex **restrict  pfx)
{ 
    __assume_aligned(pfpT,32);
    __assume_aligned(AP,32);
    __assume_aligned(xPrimary,32);
    __assume_aligned(preflectivity,32);
    __assume_aligned(AG,32);
    __assume_aligned(xGhost,32);
    __assume_aligned(pfx,32);

    for(int ix=0; ix<nx; ix++) { 
        float complex tmp;

#pragma vector aligned         
        for(int ip=0; ip<np; ip++) {
            xPrimary[ip] = 1.0f;
            xGhost[ip] = 1.0f; 
        }
	for(int ifreq=0; ifreq<hfid; ifreq++) {
            tmp = 0.0f;
#pragma vector aligned             
            for(int ip=0; ip<np; ip++) {
                tmp+=pfpT[ifreq][ip]*(xPrimary[ip]+preflectivity[ifreq]*xGhost[ip]); 
                xPrimary[ip] *= AP[ix][ip];
                xGhost[ip] *= AG[ix][ip];   
            }
            pfx[ix][ifreq] = tmp;                                         
	}
	tmp = 0.0f;                                                              
#pragma simd reduction (+:tmp)
      	for(int ip=0;ip<np;ip++)                                                 
	    tmp += pfpT[0][ip];

	pfx[ix][0] = tmp*(1.0+preflectivity[0]); 
    } 
}

#define taup_warp_solo_z1(npx, nsgtrace, hfid, pfxT, AP, yPrimary, method, zytop, integration, pfp) \
{ for(int ip=0; ip<npx; ip++) {                                                     \
	float complex tmp;                                                              \
	for(int itrace=0; itrace<nsgtrace; itrace++){                                   \
		if(method==1)	yPrimary[itrace] = 1.0f;				    \
		else yPrimary[itrace] = 1.0f*zytop[ip];			            \
	}									            \
	for(int ifreq=0; ifreq<hfid; ifreq++) {                                       \
		tmp = 0.0f;                                                               \
		for(int itrace=0; itrace<nsgtrace; itrace++) {                            \
			tmp += yPrimary[itrace]*pfxT[ifreq][itrace];                          \
			yPrimary[itrace] *= AP[ip][itrace];                                   \
		} pfp[ip][ifreq] = tmp*integration[ifreq];                                \
	} } }                                                                                                                                                  

#define taup_warp_solo_z2(nsgtrace, np, hfid, pfpT, TP, yPrimary, method, zytop, integration, pfx) \
{ for(int itrace=0; itrace<nsgtrace; itrace++) {                                    \
	float complex tmp;                                                              \
	for(int ip=0; ip<np; ip++){                                                     \
		if(method==1)	yPrimary[ip] = 1.0f;				            \
		else yPrimary[ip] = 1.0f*zytop[ip];				            \
	}									            \
	for(int ifreq=0; ifreq<hfid; ifreq++) {                                       \
		tmp = 0.0f;                                                               \
		for(int ip=0; ip<np; ip++) {                                              \
			tmp += yPrimary[ip]*pfpT[ifreq][ip];                                  \
			yPrimary[ip] *= TP[itrace][ip];                                       \
		} pfx[itrace][ifreq] = tmp*integration[ifreq];                            \
	} } }                                                                         
// zxue
#define taup_warp_quad_p(nx, np, hfid, pfpT, TP, yPrimary, reflectivity, reflectivity2, TG, \
		yGhost, TGS, yGhostS, TGR, yGhostR, pfx)                                      \
{ for(int ix=0; ix<nx; ix++) {                                                      \
	float complex tmp;                                                              \
	for(int ip=0; ip<np; ip++) {                                                    \
		yPrimary[ip] = 1.0f;                                                          \
		yGhost[ip]   = 1.0f;                                   \
		yGhostS[ip]  = 1.0f;                                                  \
		yGhostR[ip]  = 1.0f; }                                                \
	for(int ifreq=0; ifreq<hfid; ifreq++) {                                       \
		tmp = 0.0f;                                                               \
		for(int ip=0; ip<np; ip++) {                                              \
			tmp += yPrimary[ip]*pfpT[ifreq][ip];                                  \
			yPrimary[ip] *= TP[ix][ip];                                           \
			tmp += reflectivity[ifreq] *reflectivity2[ifreq]*yGhost[ip]*pfpT[ifreq][ip];                                    \
			yGhost[ip] *= TG[ix][ip];                                             \
			tmp += reflectivity2[ifreq]*yGhostS[ip]*pfpT[ifreq][ip];                                   \
			yGhostS[ip] *= TGS[ix][ip];                                           \
			tmp += reflectivity[ifreq] *yGhostR[ip]*pfpT[ifreq][ip];                                   \
			yGhostR[ip] *= TGR[ix][ip];                                           \
		} pfx[ix][ifreq] = tmp;                                                   \
	}                                                                             \
							 } }

// zxue
#define taup_warp_dual_azay1(npx, nsgtrace, hfid, pfxT, AP, xPrimary, reflectivity, AG, xGhost, method, zytop, integration, pfp)	\
{ for(int ip=0; ip<npx; ip++) {						\
	float complex tmp;                    			        \
	for(int itrace=0; itrace<nsgtrace; itrace++) {                        \
		if(method==1)                                                       \
		{									\
			xPrimary[itrace] = 1.0f;					\
			xGhost[itrace] = 1.0f;					\
			\
		}									\
		else								\
		{									\
			xPrimary[itrace] = 1.0f*zytop[ip];				\
			xGhost[itrace] = zytop[ip];			\
		}									\
	}									\
	for(int ifreq=0; ifreq<hfid; ifreq++) {				\
		tmp = 0.0f;								\
		for(int itrace=0; itrace<nsgtrace; itrace++) {			\
			tmp += xPrimary[itrace]*pfxT[ifreq][itrace];		        \
			tmp += reflectivity[ifreq]*xGhost[itrace]*pfxT[ifreq][itrace];			\
			xPrimary[itrace] *= AP[ip][itrace];				\
			xGhost[itrace] *= AG[ip][itrace];	                    		\
		}pfp[ip][ifreq] = tmp*integration[ifreq];                           \
	}                                                                     \
							  }}									\

// zxue
#define taup_warp_dual_azay2(nsgtrace, npx, hfid, pfpT, TP, xPrimary, reflectivity, TG, xGhost, method, zytop, integration, pfx)	\
{ for(int itrace=0; itrace<nsgtrace; itrace++) {			\
	float complex tmp;                    			        \
	for(int ip=0; ip<npx; ip++) {                                         \
		if(method==1)                                                       \
		{									\
			xPrimary[ip] = 1.0f;						\
			xGhost[ip] = 1.0f;					\
			\
		}									\
		else								\
		{									\
			xPrimary[ip] = 1.0f*zytop[ip];					\
			xGhost[ip] = zytop[ip];				\
		}									\
	}									\
	for(int ifreq=1; ifreq<hfid; ifreq++) {				\
		tmp = 0.0f;								\
		for(int ip=0; ip<npx; ip++) {                                       \
			xPrimary[ip] *= TP[itrace][ip];					\
			xGhost[ip] *= TG[itrace][ip];	                    	        \
			tmp += xPrimary[ip]*pfpT[ifreq][ip];				\
			tmp += reflectivity[ifreq]*xGhost[ip]*pfpT[ifreq][ip];				\
		}pfx[itrace][ifreq] = tmp*integration[ifreq];                       \
	}                                                                     \
	tmp = 0.0f;								\
	for(int ip=0;ip<npx;ip++)						\
	tmp += pfpT[0][ip];							\
	pfx[itrace][0] = tmp*(1.0+reflectivity[0]); } }                          \


/*
 *   There is no need for special treat for the loops, icc is really good at 
 * vectorize and unroll them.
 * 
 */

#define FFT_R2C_PP(plan, t, m, f, nf2)                                                            \
	for(int i=0; i<m; i++) {                                                                        \
		PFL_RCFFT(plan, (nf2-2),(float*)(t)[i],(PFLComplex*)(f)[i],1,-1,(nf2),(nf2),1.0f/((nf2)-2));  \
	}

#define FFT_C2R_PP(plan, f, m, t, nf2)                                                            \
	for(int i=0; i<m; i++)                                                                          \
PFL_CRFFT(plan, (nf2-2),(PFLComplex*)(f)[i],(float*)(t)[i],1,+1,(nf2),(nf2),1.0f);            \

#define MAT_TRANSPOSE(a, t, m, n)                                                   \
	for(int i=0; i<m; i++) { for(int j=0; j<n; j++) t[i][j] = a[j][i]; }

/**
 * @brief Sparsity weighting calculated from instantaneous_amplitude^spower
 * @param rfplan    [R] r2c forward plan of length nr
 * @param crplan    [R] c2c reverse plan of length nr
 * @param inp       [S] input data, same as output when exit
 * @param weight    [S] output data
 * @param nr        [R] fft logical length, and fast dimention of inp buffer
 * @param np        [R] number of plan wave trace
 * @param spower    [R] the power to rise
 * @param tmp1      [S] scratch, size of nr float complex
 * @param tmp2      [S] scratch, size of nr float complex
 *
 * Hopefully this will be faster than hilbert transform in time domain. (verified result correct)
 *    The accurate instantaneous amplitude is not good for deghost, and create more ringing,
 * the time domain truncated filter version has better time domain resolution.
 * */
	MSDGINTERP_WARP_DISPATCH
static void instantaneous_weight(void *ptask, const float * restrict hw, int nt, int nr, int np,
		float ** restrict inp, float ** restrict weight,
		float spower, float * restrict conv)
{
	float * restrict hilbert = conv+nt;
	for(int i=0; i<np; i++) {
		int status = vslsConvExec1D(ptask, hw, 1, inp[i], 1, conv, 1);
		assert(status==VSL_STATUS_OK);
		for(int j=0; j<nr; j++)
			inp[i][j] = powf(inp[i][j]*inp[i][j]+hilbert[j]*hilbert[j], spower);
		memcpy(weight[i], inp[i], nr*sizeof(float));
	}
}

	MSDGINTERP_WARP_DISPATCH
static void c_l1inv_TAUP_FW_DG(int nsgtrace, int hfid, int lfid, int npx, 
		float **TX_p, float **TX_az, float **TX_ay, float ** taup, l1inv_t *l1para)
{
	int ip, ix, ifreq, idata, nsgtrace4d, nsgtracep, nsgtrace_az, nsgtrace_ay; 

	float * restrict ztop = l1para->ztop;
	float * restrict ytop = l1para->ytop;

	float complex ** restrict AP       = (float complex **)l1para->AP_p;
	float complex ** restrict AG       = (float complex **)l1para->AG_p;
	float complex ** restrict AS       = (float complex **)l1para->AS_p;

	float complex ** restrict pfp      = (float complex **)l1para->pfp;
	float complex ** restrict pfx      = (float complex **)l1para->pfx;
	float complex ** restrict pfxT     = (float complex **)l1para->pfxT;
	float complex  * restrict yGhost   = (float complex  *)l1para->xGhost;
	float complex  * restrict yPrimary = (float complex  *)l1para->xPrimary;

	float complex ** restrict pfpadd    = (float complex **)l1para->pfpadd;

	float complex * restrict integration = (float complex  *)l1para->integration;
	float **TX_p_4d;

	float complex ** restrict AGS     = (float complex **)l1para->AGS_p;
	float complex ** restrict AGR     = (float complex **)l1para->AGR_p;
	float complex  * restrict yGhostS = (float complex  *)l1para->xGhostS;
	float complex  * restrict yGhostR = (float complex  *)l1para->xGhostR;

	float complex ** restrict pSrcFilter = (float complex **) l1para->pSrcFilter;

        const int fftnr=l1para->fftnr;

	for (ifreq=0;ifreq<l1para->fftnc;ifreq++) integration[ifreq] = 1.0f;

	for(ip=0;ip<npx;ip++)
		memset(pfpadd[ip], 0.0f, l1para->fftNC*2*sizeof(float));

	if(l1para->choosemethod==CHOOSEMETHOD_MSDGI && l1para->zerotraceattr == 1)
	{
		if (l1para->choosemethod_msdgi_pseudo == 0)
		{
			nsgtracep   = l1para->nsgtracearr4d[0];
			nsgtrace_az = l1para->nsgtracearr4d[1];
			nsgtrace_ay = l1para->nsgtracearr4d[2];
		}
		else
		{
			nsgtracep = nsgtrace;
			nsgtrace_az = l1para->nsgtracearr4d[l1para->ndata];
		}
	}
	else
	{
		nsgtracep   = nsgtrace;
		nsgtrace_az = nsgtrace;
		nsgtrace_ay = nsgtrace;
	}

	if (l1para->ms_pyes==1 && (l1para->choosemethod==CHOOSEMETHOD_MSDGI||l1para->choosemethod==CHOOSEMETHOD_BSDG3D
				||l1para->choosemethod==CHOOSEMETHOD_BSDG2D))
	{

		FFT_R2C_PP(l1para->fftplan[0], TX_p, nsgtracep, pfx, fftnr);

		MAT_TRANSPOSE(pfx, pfxT, hfid, nsgtracep);

		taup_warp_dual_p(npx, nsgtracep, hfid, pfxT, AP, yPrimary, l1para->prfl2, AG, yGhost, pfp);

		if(l1para->choosemethod==CHOOSEMETHOD_MSDGI && l1para->jointmc == YESYES)
		{
			for (ip = 0; ip < npx; ip++)
				for (ifreq = 0; ifreq < hfid; ifreq++)
					pfp[ip][ifreq] *= conj (pSrcFilter[ip][ifreq]);
		}

		for(ip=0;ip<npx;ip++) {
			for(ifreq=0;ifreq<hfid;ifreq++)
			{
				pfpadd[ip][ifreq] = pfp[ip][ifreq];
				pfp[ip][ifreq] = 0.0f;
			}
		}
	}
	else if (l1para->ms_pyes==1 && (l1para->choosemethod==CHOOSEMETHOD_JOINTSR3D || l1para->choosemethod==CHOOSEMETHOD_JOINTSR2D))
	{
		FFT_R2C_PP(l1para->fftplan[0], TX_p, nsgtrace, pfx, fftnr);

		MAT_TRANSPOSE(pfx, pfxT, hfid, nsgtrace);

		if (l1para->computemode == 0)
		{
			taup_warp_quad_p(npx, nsgtrace, hfid, pfxT, AP, yPrimary, l1para->prfl2, l1para->srcrfl,AG, yGhost,
					AGS, yGhostS, AGR, yGhostR, pfp);
		}
		else
		{
			if ((l1para->choosemethod==CHOOSEMETHOD_JOINTSR3D||l1para->choosemethod==CHOOSEMETHOD_JOINTSR2D) 
					&& (l1para->jointmethod == JOINTDESIGNATURE || l1para->jointmethod == JOINTDESIG_SRC_DEGHOST)
					|| l1para->jointmethod == JOINT_SRC_DEGHOST)
			{
				float *obliqcorr = (float *) calloc (MAX (nsgtrace, npx), sizeof (float));
				for (ip = 0; ip < MAX (nsgtrace, npx);ip++) obliqcorr[ip] = 1.0f;
				taup_warp_solo (npx, nsgtrace, hfid, pfxT, AS, yPrimary, obliqcorr, pfp);
				free (obliqcorr);
			}
			else
			{
				taup_warp_dual_p(npx, nsgtrace, hfid, pfxT, AP, yPrimary, l1para->prfl2, AG, yGhost, pfp);
			}

			for (ip = 0; ip < npx; ip++)
				for (ifreq = 0; ifreq < hfid; ifreq++)
					pfp[ip][ifreq] *= conj (pSrcFilter[ip][ifreq]);
		}

		for(ip=0;ip<npx;ip++) {
			for(ifreq=0;ifreq<hfid;ifreq++)
			{
				pfpadd[ip][ifreq] = pfp[ip][ifreq];
				pfp[ip][ifreq] = 0.0f;
			}
		}
	}
	else if (l1para->ms_pyes==1 && l1para->choosemethod==CHOOSEMETHOD_FP4D)
	{

		for(idata=0;idata<l1para->ndata;idata++)
		{

			TX_p_4d = &TX_p[idata*nsgtrace];

			l1para->AP_p = (float **)l1para->AP_parray[idata];
			l1para->AG_p = (float **)l1para->AG_parray[idata];

			AP       = (float complex **)l1para->AP_p;
			AG       = (float complex **)l1para->AG_p;

			nsgtrace4d = l1para->nsgtracearr4d[idata]; //PS: check check

			FFT_R2C_PP(l1para->fftplan[0], TX_p_4d, nsgtrace4d, pfx, fftnr);

			MAT_TRANSPOSE(pfx, pfxT, hfid, nsgtrace4d);

			taup_warp_dual_p(npx, nsgtrace4d, hfid, pfxT, AP, yPrimary, l1para->prfl2, AG, yGhost, pfp);

			if(l1para->jointmc == YESYES)
			{

				l1para->pSrcFilter = (float **)l1para->pSrcFilter_array[idata];
				pSrcFilter = (float complex **) l1para->pSrcFilter;

				for (ip = 0; ip < npx; ip++)
					for (ifreq = 0; ifreq < hfid; ifreq++)
						pfp[ip][ifreq] *= conj (pSrcFilter[ip][ifreq]);   //// BUG FIX : NEW
			}


			for(ip=0;ip<npx;ip++) {
				for(ifreq=0;ifreq<hfid;ifreq++)
				{
					pfpadd[ip][ifreq] = pfp[ip][ifreq]+pfpadd[ip][ifreq];
					pfp[ip][ifreq] = 0.0f;
				}
			}	

		}

		if (l1para->ms_zyes==1) nsgtrace_az = nsgtrace4d = l1para->nsgtracearr4d[l1para->ndata];
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////

	if (l1para->ms_zyes==1)
	{

		FFT_R2C_PP(l1para->fftplan[0], TX_az, nsgtrace_az, pfx, fftnr);

		MAT_TRANSPOSE(pfx, pfxT, hfid, nsgtrace_az);

		for(ifreq=0; ifreq<hfid; ifreq++) {        
			if (l1para->zmethod==1||l1para->zmethod==2)
				integration[ifreq] = l1para->zlowf[ifreq];	  
			else
			{
				integration[ifreq] = (-I)*l1para->zwtilt[ifreq];
			}
		}

		if (l1para->choosemethod==CHOOSEMETHOD_FP4D)
		{
			l1para->AP_p = (float **)l1para->AP_parray[0];
			l1para->AG_p = (float **)l1para->AG_parray[0];
			AP       = (float complex **)l1para->AP_p;
			AG       = (float complex **)l1para->AG_p;

			if(l1para->jointmc == YESYES)
			{
				l1para->pSrcFilter = (float **)l1para->pSrcFilter_array[0];
				pSrcFilter = (float complex **) l1para->pSrcFilter;
			}

		}

		// zxue
		taup_warp_dual_azay1(npx, nsgtrace_az, hfid, pfxT, AP, yPrimary, l1para->zrfl, AG, yGhost, 
				l1para->zmethod, ztop, integration, pfp);

		if(l1para->jointmc == YESYES)
		{
			for (ip = 0; ip < npx; ip++)
				for (ifreq = 0; ifreq < hfid; ifreq++)
					pfp[ip][ifreq] *= conj (pSrcFilter[ip][ifreq]); 
		}

		for(ip=0;ip<npx;ip++) {
			for(ifreq=0;ifreq<hfid;ifreq++)
			{
				pfpadd[ip][ifreq] = pfp[ip][ifreq]+pfpadd[ip][ifreq];
				pfp[ip][ifreq] = 0.0f;  //// <---
			}
		}	
	}

	if (l1para->ms_yyes==1)
	{

		FFT_R2C_PP(l1para->fftplan[0], TX_ay, nsgtrace_ay, pfx, fftnr);

		MAT_TRANSPOSE(pfx, pfxT, hfid, nsgtrace_ay);

		for(ifreq=0; ifreq<hfid; ifreq++) {        
			if (l1para->ymethod==1||l1para->ymethod==2)
				integration[ifreq] = l1para->ylowf[ifreq];	  
			else
			{
				integration[ifreq] = (-I)*l1para->ywtilt[ifreq];
			}
		}

		// zxue
		taup_warp_dual_azay1(npx, nsgtrace_ay, hfid, pfxT, AP, yPrimary, l1para->yrfl, AG, yGhost, 
				l1para->ymethod, ytop, integration, pfp);

		if(l1para->jointmc == YESYES)
		{
			for (ip = 0; ip < npx; ip++)
				for (ifreq = 0; ifreq < hfid; ifreq++)
					pfp[ip][ifreq] *= conj (pSrcFilter[ip][ifreq]);  
		}

		for(ip=0;ip<npx;ip++) {
			for(ifreq=0;ifreq<hfid;ifreq++)
				pfpadd[ip][ifreq] = pfp[ip][ifreq]+pfpadd[ip][ifreq];
			pfp[ip][ifreq] = 0.0f;  //// <---
		}	
	}


	for(ip=0;ip<npx;ip++) {
		for(ifreq=hfid;ifreq<l1para->fftnc;ifreq++)
			pfpadd[ip][ifreq] = 0.0f;

		pfpadd[ip][0] = pfpadd[ip][1] = 0.0f;  //// <---                    
		for (ifreq=2;ifreq<lfid;ifreq++) 
			pfpadd[ip][ifreq]=pfpadd[ip][ifreq]*((float)ifreq-1.0f)/((float)lfid-2.0f); //// <---
	}

	FFT_C2R_PP(l1para->fftplan[1], pfpadd, npx, taup, fftnr);

	for(ip=0;ip<npx;ip++) 
		for(ifreq=fftnr-2;ifreq<=fftnr-1;ifreq++) 
			taup[ip][ifreq] = 0.0f;
}

	MSDGINTERP_WARP_DISPATCH
static void c_l1inv_matrix_TAUP_FW_RV_DG(int nsgtrace, int hfid, int npx, float **taup, 
		float **taupou, l1inv_t *l1para)
{
	//    perf_start(MSDGINTERP_PERF_CTAUP);
	int ip, ix, ifreq, itrace, i, j, idata, nsgtrace4d, nsgtracep, nsgtrace_az, nsgtrace_ay; 

	float * restrict ztop = l1para->ztop;
	float * restrict ytop = l1para->ytop;

	float complex ** restrict AG       = (float complex**)l1para->AG_p;
	float complex ** restrict TG       = (float complex**)l1para->TG_p;
	float complex ** restrict AP       = (float complex**)l1para->AP_p;
	float complex ** restrict TP       = (float complex**)l1para->TP_p;
	float complex ** restrict AS       = (float complex**)l1para->AS_p;
	float complex ** restrict TS       = (float complex**)l1para->TS_p;

	float complex ** restrict pfp      = (float complex**)l1para->pfp;
	float complex ** restrict pfpinit  = (float complex**)l1para->pfpinit; // LKL
	float complex ** restrict pfpin    = (float complex**)l1para->pfpin;
	float complex ** restrict pfpadd   = (float complex**)l1para->pfpadd;
	float complex ** restrict pfx      = (float complex**)l1para->pfx;
	float complex ** restrict pfpT     = (float complex**)l1para->pfpT;
	float complex ** restrict pfxT     = (float complex**)l1para->pfxT;
	float complex  * restrict yGhost   = (float complex *)l1para->xGhost;
	float complex  * restrict yPrimary = (float complex *)l1para->xPrimary;

	float complex * restrict integration = (float complex  *)l1para->integration;

	float complex tmp;

	float complex ** restrict AGS     = (float complex **)l1para->AGS_p;
	float complex ** restrict AGR     = (float complex **)l1para->AGR_p;
	float complex ** restrict TGS     = (float complex **)l1para->TGS_p;
	float complex ** restrict TGR     = (float complex **)l1para->TGR_p;
	float complex  * restrict yGhostS = (float complex  *)l1para->xGhostS;
	float complex  * restrict yGhostR = (float complex  *)l1para->xGhostR;

	float complex ** restrict pSrcFilter = (float complex **) l1para->pSrcFilter;
        const int fftnr=l1para->fftnr;

	for (ifreq=0;ifreq<l1para->fftnc;ifreq++) integration[ifreq] = 1.0f;

	for(ip=0;ip<npx;ip++)
		memset(pfpadd[ip], 0.0f, l1para->fftNC*2*sizeof(float));

	if(l1para->choosemethod==CHOOSEMETHOD_MSDGI && l1para->zerotraceattr == 1)
	{
		if (l1para->choosemethod_msdgi_pseudo == 0)
		{
			nsgtracep   = l1para->nsgtracearr4d[0];
			nsgtrace_az = l1para->nsgtracearr4d[1];
			nsgtrace_ay = l1para->nsgtracearr4d[2];
		}
		else
		{
			nsgtracep = nsgtrace;
			nsgtrace_az = l1para->nsgtracearr4d[l1para->ndata];
		}
	}
	else
	{
		nsgtracep   = nsgtrace;
		nsgtrace_az = nsgtrace;
		nsgtrace_ay = nsgtrace;
	}

	FFT_R2C_PP(l1para->fftplan[0], taup, npx, pfpin, fftnr);  

	//LKL: keep the original taup model for later use in regularization
	if(l1para->regparam > 0.0f||l1para->choosemethod==CHOOSEMETHOD_FP4D)
	{
		for (ip = 0; ip < npx; ip++)
			for (ifreq = 0; ifreq < hfid; ifreq++)
				pfpinit[ip][ifreq] = pfpin[ip][ifreq];
	}

	if (l1para->ms_pyes==1 && (l1para->choosemethod==CHOOSEMETHOD_MSDGI||l1para->choosemethod==CHOOSEMETHOD_BSDG3D
				||l1para->choosemethod==CHOOSEMETHOD_BSDG2D))  //// BUG FIX
	{

		if(l1para->choosemethod==CHOOSEMETHOD_MSDGI && l1para->jointmc == YESYES)
		{
			for (ip = 0; ip < npx; ip++)
				for (ifreq = 0; ifreq < hfid; ifreq++)
					pfpin[ip][ifreq] *= pSrcFilter[ip][ifreq];
		}

		MAT_TRANSPOSE(pfpin, pfpT, hfid, npx);

		taup_warp_dual_p(nsgtracep, npx, hfid, pfpT, TP, yPrimary, l1para->prfl2, TG, yGhost, pfx);

		////////// above is reverse taup, and below is forward taup =======================
		MAT_TRANSPOSE(pfx, pfxT, hfid, nsgtracep);

		taup_warp_dual_p(npx, nsgtracep, hfid, pfxT, AP, yPrimary, l1para->prfl2, AG, yGhost, pfp);

		if(l1para->choosemethod==CHOOSEMETHOD_MSDGI && l1para->jointmc == YESYES)
		{
			for (ip = 0; ip < npx; ip++)
				for (ifreq = 0; ifreq < hfid; ifreq++)
					pfp[ip][ifreq] *= conj (pSrcFilter[ip][ifreq]);
		}

		for (i=0;i<npx;i++) {   //mute high frequency, may not be necessary in the cg loop!
			for (j=hfid;j<l1para->fftnc;j++) 
				pfp[i][j] = 0.0f;
		}

		for (i=0;i<npx;i++) {   //mute high frequency, may not be necessary in the cg loop!
			for (j=0;j<l1para->fftnc;j++) 
				pfpadd[i][j] = pfp[i][j];
		}

	}
	else if (l1para->ms_pyes==1 && (l1para->choosemethod==CHOOSEMETHOD_JOINTSR3D || l1para->choosemethod==CHOOSEMETHOD_JOINTSR2D))
	{
		if (l1para->computemode == 1)
			for (ip = 0; ip < npx; ip++)
				for (ifreq = 0; ifreq < hfid; ifreq++)
					pfpin[ip][ifreq] *= pSrcFilter[ip][ifreq];

		MAT_TRANSPOSE(pfpin, pfpT, hfid, npx);

		if (l1para->computemode == 0)
		{
			taup_warp_quad_p(nsgtrace, npx, hfid, pfpT, TP, yPrimary, l1para->prfl2, l1para->srcrfl, TG, yGhost,
					TGS, yGhostS, TGR, yGhostR, pfx);
		}
		else
		{
			if ((l1para->choosemethod==CHOOSEMETHOD_JOINTSR3D||l1para->choosemethod==CHOOSEMETHOD_JOINTSR2D) 
					&& (l1para->jointmethod == JOINTDESIGNATURE || l1para->jointmethod == JOINTDESIG_SRC_DEGHOST
						||l1para->jointmethod == JOINT_SRC_DEGHOST))
			{
				float *obliqcorr = (float *) calloc (MAX (nsgtrace, npx), sizeof (float));
				for (ip = 0; ip < MAX (nsgtrace, npx);ip++) obliqcorr[ip] = 1.0f;
				taup_warp_solo (nsgtrace, npx, hfid, pfpT, TS, yPrimary, obliqcorr, pfx);
				free (obliqcorr);
			}
			else
			{
				taup_warp_dual_p(nsgtrace, npx, hfid, pfpT, TP, yPrimary, l1para->prfl2, TG, yGhost, pfx);
			}
		}

		////////// above is reverse taup, and below is forward taup =======================
		MAT_TRANSPOSE(pfx, pfxT, hfid, nsgtrace);

		if (l1para->computemode == 0)
		{
			taup_warp_quad_p(npx, nsgtrace, hfid, pfxT, AP, yPrimary, l1para->prfl2, l1para->srcrfl, AG, yGhost,
					AGS, yGhostS, AGR, yGhostR, pfp);
		}
		else
		{
			if ((l1para->choosemethod==CHOOSEMETHOD_JOINTSR3D||l1para->choosemethod==CHOOSEMETHOD_JOINTSR2D) 
					&& (l1para->jointmethod == JOINTDESIGNATURE || l1para->jointmethod == JOINTDESIG_SRC_DEGHOST
						||l1para->jointmethod == JOINT_SRC_DEGHOST))
			{
				float *obliqcorr = (float *) calloc (MAX (nsgtrace, npx), sizeof (float));
				for (ip = 0; ip < MAX (nsgtrace, npx);ip++) obliqcorr[ip] = 1.0f;
				taup_warp_solo (npx, nsgtrace, hfid, pfxT, AS, yPrimary, obliqcorr, pfp);
				free (obliqcorr);
			}
			else
			{
				taup_warp_dual_p(npx, nsgtrace, hfid, pfxT, AP, yPrimary, l1para->prfl2, AG, yGhost, pfp);
			}
			for (ip = 0; ip < npx; ip++)
				for (ifreq = 0; ifreq < hfid; ifreq++)
					pfp[ip][ifreq] *= conj (pSrcFilter[ip][ifreq]);
		}

		for (i=0;i<npx;i++) {   //mute high frequency, may not be necessary in the cg loop!
			for (j=hfid;j<l1para->fftnc;j++) 
				pfp[i][j] = 0.0f;
		}

		for (i=0;i<npx;i++) {   //mute high frequency, may not be necessary in the cg loop!
			for (j=0;j<l1para->fftnc;j++) 
				pfpadd[i][j] = pfp[i][j];
		}

	}
	else if (l1para->ms_pyes==1 && l1para->choosemethod==CHOOSEMETHOD_FP4D)
	{

		for(idata=0;idata<l1para->ndata;idata++) 
		{

			for (ip = 0; ip < npx; ip++)
				for (ifreq = 0; ifreq < hfid; ifreq++)
					pfpin[ip][ifreq] = pfpinit[ip][ifreq];

			if(l1para->jointmc == YESYES)
			{
				l1para->pSrcFilter = (float **)l1para->pSrcFilter_array[idata];
				pSrcFilter = (float complex **) l1para->pSrcFilter;

				for (ip = 0; ip < npx; ip++)
					for (ifreq = 0; ifreq < hfid; ifreq++)
						pfpin[ip][ifreq] *= pSrcFilter[ip][ifreq];  //// BUG FIX : NEW
			}

			nsgtrace4d = l1para->nsgtracearr4d[idata];

			MAT_TRANSPOSE(pfpin, pfpT, hfid, npx);

			l1para->TP_p = (float **)l1para->TP_parray[idata];
			l1para->TG_p = (float **)l1para->TG_parray[idata];

			TP       = (float complex **)l1para->TP_p;
			TG       = (float complex **)l1para->TG_p;


			taup_warp_dual_p(nsgtrace4d, npx, hfid, pfpT, TP, yPrimary, l1para->prfl2, TG, yGhost, pfx);

			////////// above is reverse taup, and below is forward taup =======================
			MAT_TRANSPOSE(pfx, pfxT, hfid, nsgtrace4d);

			l1para->AP_p = (float **)l1para->AP_parray[idata];
			l1para->AG_p = (float **)l1para->AG_parray[idata];

			AP       = (float complex **)l1para->AP_p;
			AG       = (float complex **)l1para->AG_p;

			taup_warp_dual_p(npx, nsgtrace4d, hfid, pfxT, AP, yPrimary, l1para->prfl2, AG, yGhost, pfp);

			if(l1para->jointmc == YESYES)
			{
				for (ip = 0; ip < npx; ip++)
					for (ifreq = 0; ifreq < hfid; ifreq++)
						pfp[ip][ifreq] *= conj (pSrcFilter[ip][ifreq]);   //// BUG FIX : NEW
			}

			for (i=0;i<npx;i++) {   //mute high frequency, may not be necessary in the cg loop!
				for (j=hfid;j<l1para->fftnc;j++) 
					pfp[i][j] = 0.0f;
			}

			for (i=0;i<npx;i++) {   //mute high frequency, may not be necessary in the cg loop!
				for (j=0;j<l1para->fftnc;j++) 
					pfpadd[i][j] = pfp[i][j] + pfpadd[i][j];
			}

		}

		if (l1para->ms_zyes==1) nsgtrace_az = nsgtrace4d = l1para->nsgtracearr4d[l1para->ndata];
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////

	if (l1para->ms_zyes==1)
	{
		if (l1para->ms_pyes == 0 && l1para->jointmc == YESYES)
			for (ip = 0; ip < npx; ip++)
				for (ifreq = 0; ifreq < hfid; ifreq++)
					pfpin[ip][ifreq] *= pSrcFilter[ip][ifreq];  

		if (l1para->choosemethod==CHOOSEMETHOD_FP4D)
		{

			for (ip = 0; ip < npx; ip++)
				for (ifreq = 0; ifreq < hfid; ifreq++)
					pfpin[ip][ifreq] = pfpinit[ip][ifreq];

			l1para->TP_p = (float **)l1para->TP_parray[0];
			l1para->TG_p = (float **)l1para->TG_parray[0];
			TP       = (float complex **)l1para->TP_p;
			TG       = (float complex **)l1para->TG_p;

			if(l1para->jointmc == YESYES)
			{
				l1para->pSrcFilter = (float **)l1para->pSrcFilter_array[0];
				pSrcFilter = (float complex **) l1para->pSrcFilter;

				for (ip = 0; ip < npx; ip++)
					for (ifreq = 0; ifreq < hfid; ifreq++)
						pfpin[ip][ifreq] *= pSrcFilter[ip][ifreq];  

			}

		}

		MAT_TRANSPOSE(pfpin, pfpT, hfid, npx);         //// <--- DO WE NEED TO DO THE MATRIX TRANSPOSE AGAIN ??

		for (ifreq=0;ifreq<hfid;ifreq++)
		{
			if (l1para->zmethod==0)
			{
				for (ip=0;ip<npx;ip++)
				{
					pfpT[ifreq][ip] = (+I)*pfpT[ifreq][ip]*l1para->zwtilt[ifreq];
				}
			}
		}

		// zxue
		taup_warp_dual_azay2(nsgtrace_az, npx, hfid, pfpT, TP, yPrimary, l1para->zrfl, TG, yGhost, 
				l1para->zmethod, ztop, integration, pfx);

		////////// above is reverse taup, and below is forward taup =======================
		MAT_TRANSPOSE(pfx, pfxT, hfid, nsgtrace_az);

		for(int ifreq=0; ifreq<hfid; ifreq++) {        
			if (l1para->zmethod==0)
				integration[ifreq] = (-I)*l1para->zwtilt[ifreq];

		}

		if (l1para->choosemethod==CHOOSEMETHOD_FP4D)
		{
			l1para->AP_p = (float **)l1para->AP_parray[0];
			l1para->AG_p = (float **)l1para->AG_parray[0];
			AP       = (float complex **)l1para->AP_p;
			AG       = (float complex **)l1para->AG_p;
		}

		// zxue
		taup_warp_dual_azay1(npx, nsgtrace_az, hfid, pfxT, AP, yPrimary, l1para->zrfl, AG, yGhost, 
				l1para->zmethod, ztop, integration, pfp);

		if(l1para->jointmc == YESYES)
		{
			for (ip = 0; ip < npx; ip++)
				for (ifreq = 0; ifreq < hfid; ifreq++)
					pfp[ip][ifreq] *= conj (pSrcFilter[ip][ifreq]); 
		}

		// zxue
		for (ip=0;ip<npx;ip++) {
			tmp = 0.0f;							
			for(itrace=0;itrace<nsgtrace_az;itrace++)			
				tmp += pfxT[0][itrace];					
			pfp[ip][0] = tmp*(1.0+l1para->zrfl[0]);				
		}

		for (i=0;i<npx;i++) {   //mute high frequency, may not be necessary in the cg loop!
			for (j=hfid;j<l1para->fftnc;j++) 
				pfp[i][j] = 0.0f;
		}

		for (i=0;i<npx;i++) {   //mute high frequency, may not be necessary in the cg loop!
			for (j=0;j<l1para->fftnc;j++) 
				pfpadd[i][j] = pfpadd[i][j] + pfp[i][j];
		}

	}

	////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////

	if (l1para->ms_yyes==1)
	{

		if (l1para->ms_pyes == 0 && l1para->ms_zyes == 0 && l1para->jointmc == YESYES)
			for (ip = 0; ip < npx; ip++)
				for (ifreq = 0; ifreq < hfid; ifreq++)
					pfpin[ip][ifreq] *= pSrcFilter[ip][ifreq];  

		MAT_TRANSPOSE(pfpin, pfpT, hfid, npx);         //// <--- DO WE NEED TO DO THE MATRIX TRANSPOSE AGAIN ??

		for (ifreq=0;ifreq<hfid;ifreq++)
		{
			if (l1para->ymethod==0)
			{
				for (ip=0;ip<npx;ip++)
				{
					pfpT[ifreq][ip] = (+I)*pfpT[ifreq][ip]*l1para->ywtilt[ifreq];
				}
			}
		}

		// zxue
		taup_warp_dual_azay2(nsgtrace_ay, npx, hfid, pfpT, TP, yPrimary, l1para->yrfl, TG, yGhost, 
				l1para->ymethod, ytop, integration, pfx);

		////////// above is reverse taup, and below is forward taup =======================
		MAT_TRANSPOSE(pfx, pfxT, hfid, nsgtrace_ay);

		for(int ifreq=0; ifreq<hfid; ifreq++) {        
			if (l1para->ymethod==0)
				integration[ifreq] = (-I)*l1para->ywtilt[ifreq];

		}

		// zxue
		taup_warp_dual_azay1(npx, nsgtrace_ay, hfid, pfxT, AP, yPrimary, l1para->yrfl, AG, yGhost, 
				l1para->ymethod, ytop, integration, pfp);

		if(l1para->jointmc == YESYES)
		{
			for (ip = 0; ip < npx; ip++)
				for (ifreq = 0; ifreq < hfid; ifreq++)
					pfp[ip][ifreq] *= conj (pSrcFilter[ip][ifreq]);   
		}

		// zxue
		for (ip=0;ip<npx;ip++) {
			tmp = 0.0f;							
			for(itrace=0;itrace<nsgtrace_ay;itrace++)			
				tmp += pfxT[0][itrace];					
			pfp[ip][0] = tmp*(1.0+l1para->yrfl[0]);			
		}

		for (i=0;i<npx;i++) {   //mute high frequency, may not be necessary in the cg loop!
			for (j=hfid;j<l1para->fftnc;j++) 
				pfp[i][j] = 0.0f;
		}

		for (i=0;i<npx;i++) {   //mute high frequency, may not be necessary in the cg loop!
			for (j=0;j<l1para->fftnc;j++) 
				pfpadd[i][j] = pfpadd[i][j] + pfp[i][j];
		}

	}

	////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////

	//LKL: Add regparam * taup_initial to regularize the CG solutions
	float lambda = l1para->regparam * npx;
	if(lambda > 0.0f) { 
		for (i=0;i<npx;i++)
			for (j=0;j<hfid;j++) // high-frequency supposed to be muted
				pfpadd[i][j] += lambda * pfpinit[i][j];
	}

	FFT_C2R_PP(l1para->fftplan[1], pfpadd, npx, taupou, fftnr);

	for (i=0; i<npx; i++) {
		taupou[i][fftnr-2] = 0.0f;
		taupou[i][fftnr-1] = 0.0f;
	}

}


	MSDGINTERP_WARP_DISPATCH
static void c_l1inv_TAUP_FW_CMSTOS_STOS_p(int imode, int nsgtrace, int hfid, int lfid, int npx, 
		float **TX, float ** taup, l1inv_t *l1para)
{
	int ip, ifreq;

	float complex ** restrict AA;

	float complex ** restrict pfp    = (float complex **)l1para->pfp;
	float complex ** restrict pfx    = (float complex **)l1para->pfx;
	float complex ** restrict pfxT   = (float complex **)l1para->pfxT;
	float complex  * restrict xPrimary = (float complex  *)l1para->xPrimary;

	float * restrict obliqcorr = l1para->obliqcorr;
        const int fftnr=l1para->fftnr;

	if (imode == 0)
		AA = (float complex **)l1para->AP_p;
	else if (imode == 1)
		AA = (float complex **)l1para->AG_p;
	else if (imode == 2)
		AA = (float complex **)l1para->AS_p;

	if (imode == 0 || imode == 2)
		for (ip=0;ip<MAX(nsgtrace,npx);ip++) obliqcorr[ip] = 1.0f;
	else if (imode == 1)
		for (ip=0;ip<MAX(nsgtrace,npx);ip++) obliqcorr[ip] = -1.0f; // zxue

	FFT_R2C_PP(l1para->fftplan[0], TX, nsgtrace, pfx, fftnr)

	MAT_TRANSPOSE(pfx, pfxT, hfid, nsgtrace);

	taup_warp_solo(npx, nsgtrace, hfid, pfxT, AA, xPrimary, obliqcorr, pfp);

	// zxue
	if(imode == 1){
		for(ip=0;ip<npx;ip++) {
			for(ifreq=0;ifreq<hfid;ifreq++){
				pfp[ip][ifreq] *= l1para->prfl2[ifreq];
			}
		}
	}

	for(ip=0;ip<npx;ip++) {
		for(ifreq=hfid;ifreq<l1para->fftnc;ifreq++)
			pfp[ip][ifreq] = 0.0f;

		pfp[ip][0] = pfp[ip][1] = 0.0f;

		for (ifreq=2;ifreq<lfid;ifreq++) 
			pfp[ip][ifreq]=pfp[ip][ifreq]*((float)ifreq-1.0f)/((float)lfid-2.0f);
	}

	FFT_C2R_PP(l1para->fftplan[1], pfp, npx, taup, fftnr)

		for(ip=0;ip<npx;ip++) 
			for(ifreq=fftnr-2;ifreq<=fftnr-1;ifreq++) 
				taup[ip][ifreq] = 0.0f;

}

	MSDGINTERP_WARP_DISPATCH
static void c_l1inv_matrix_TAUP_FW_RV_CMSTOS_STOS_p(int imode, int nsgtrace, int hfid, int npx, 
		float **taup, float **taupou, l1inv_t *l1para)
{

	int ip, ifreq;

	float complex ** restrict AA;      
	float complex ** restrict TT;      

	float complex ** restrict pfp    = (float complex**)l1para->pfp; 
	float complex ** restrict pfx    = (float complex**)l1para->pfx;
	float complex ** restrict pfpT   = (float complex**)l1para->pfpT;
	float complex ** restrict pfxT   = (float complex**)l1para->pfxT;

	float complex  * restrict xPrimary = (float complex *)l1para->xPrimary;
        const int fftnr=l1para->fftnr;

	float * restrict obliqcorr = l1para->obliqcorr;

	if (imode == 0)
	{
		AA      = (float complex **)l1para->AP_p;
		TT      = (float complex **)l1para->TP_p;
	}
	else if (imode == 1)
	{
		AA      = (float complex **)l1para->AG_p;
		TT      = (float complex **)l1para->TG_p;

	}
	else if (imode == 2)
	{
		AA      = (float complex **)l1para->AS_p;
		TT      = (float complex **)l1para->TS_p;
	}

	if (imode == 0 || imode == 2)      for (ip=0;ip<MAX(nsgtrace,npx);ip++) obliqcorr[ip] = 1.0f;
	else if (imode == 1)               for (ip=0;ip<MAX(nsgtrace,npx);ip++) obliqcorr[ip] = -1.f; // zxue

	FFT_R2C_PP(l1para->fftplan[0], taup, npx, pfp, fftnr)

	MAT_TRANSPOSE(pfp, pfpT, hfid, npx);

	taup_warp_solo(nsgtrace, npx, hfid, pfpT, TT, xPrimary, obliqcorr, pfx);

	// zxue
	if(imode == 1){
		for(ip=0;ip<nsgtrace;ip++) {
			for(ifreq=0;ifreq<hfid;ifreq++){
				pfx[ip][ifreq] *= l1para->prfl2[ifreq];
			}
		}
	}

	/// above is reverse taup, and below is forward taup ==
	MAT_TRANSPOSE(pfx, pfxT, hfid, nsgtrace);

	taup_warp_solo(npx, nsgtrace, hfid, pfxT, AA, xPrimary, obliqcorr, pfp);

	// zxue
	if(imode == 1){
		for(ip=0;ip<npx;ip++) {
			for(ifreq=0;ifreq<hfid;ifreq++){
				pfp[ip][ifreq] *= l1para->prfl2[ifreq];
			}
		}
	}

	for (ip=0;ip<npx;ip++) {   //mute high frequency, may not be necessary in the cg loop!
		for (ifreq=hfid;ifreq<l1para->fftnc;ifreq++) 
			pfp[ip][ifreq] = 0.0f;
	}

	FFT_C2R_PP(l1para->fftplan[1], pfp, npx, taupou, fftnr)

		for (ip=0; ip<npx; ip++) {
			taupou[ip][fftnr-2] = 0.0f;
			taupou[ip][fftnr-1] = 0.0f;
		}

}

	MSDGINTERP_WARP_DISPATCH
static void c_l1inv_TAUP_FW_CMSTOS_STOS_z(int imode, int nsgtrace, int hfid, int lfid, int npx, 
		float **TX, float ** taup, l1inv_t *l1para)
{
	int ip, ifreq;

	float complex ** restrict AA;

	float complex ** restrict pfp    = (float complex **)l1para->pfp;
	float complex ** restrict pfx    = (float complex **)l1para->pfx;
	float complex ** restrict pfxT   = (float complex **)l1para->pfxT;
	float complex  * restrict xPrimary = (float complex  *)l1para->xPrimary;

	float * restrict obliqcorr = l1para->obliqcorr;
        const int fftnr=l1para->fftnr;

	float complex * restrict integration = (float complex  *)l1para->integration;

	if (imode == 0)
		AA = (float complex **)l1para->AP_p;
	else if (imode == 1)
		AA = (float complex **)l1para->AG_p;
	else if (imode == 2)
		AA = (float complex **)l1para->AS_p;

	if (imode == 0 || imode == 2)
		for (ip=0;ip<MAX(nsgtrace,npx);ip++) obliqcorr[ip] = 1.0f;
	else if (imode == 1)
		for (ip=0;ip<MAX(nsgtrace,npx);ip++) obliqcorr[ip] = -l1para->prfl2[0]; // zxue

	for (ifreq=0;ifreq<l1para->fftnc;ifreq++) integration[ifreq] = 1.0f;

	if (l1para->zmethod == 0||l1para->zmethod == 2)
		for (ip=0;ip<npx;ip++) 
			obliqcorr[ip] = l1para->ztop[ip];

	for(int ifreq=0; ifreq<hfid; ifreq++)   
		if (l1para->zmethod==0) integration[ifreq] = (-I)*l1para->zwtilt[ifreq];

	FFT_R2C_PP(l1para->fftplan[0], TX, nsgtrace, pfx, fftnr)

		MAT_TRANSPOSE(pfx, pfxT, hfid, nsgtrace);

	taup_warp_solo_z1(npx, nsgtrace, hfid, pfxT, AA, xPrimary, l1para->zmethod, obliqcorr, integration, pfp);
	////taup_warp_solo(npx, nsgtrace, hfid, pfxT, AA, xPrimary, obliqcorr, pfp);

	for(ip=0;ip<npx;ip++) {
		for(ifreq=hfid;ifreq<l1para->fftnc;ifreq++)
			pfp[ip][ifreq] = 0.0f;

		pfp[ip][0] = pfp[ip][1] = 0.0f;

		for (ifreq=2;ifreq<lfid;ifreq++) 
			pfp[ip][ifreq]=pfp[ip][ifreq]*((float)ifreq-1.0f)/((float)lfid-2.0f);
	}

	FFT_C2R_PP(l1para->fftplan[1], pfp, npx, taup, fftnr)

		for(ip=0;ip<npx;ip++) 
			for(ifreq=fftnr-2;ifreq<=fftnr-1;ifreq++) 
				taup[ip][ifreq] = 0.0f;

}

	MSDGINTERP_WARP_DISPATCH
static void c_l1inv_matrix_TAUP_FW_RV_CMSTOS_STOS_z(int imode, int nsgtrace, int hfid, int npx, 
		float **taup, float **taupou, l1inv_t *l1para)
{

	int ip, ifreq;

	float complex ** restrict AA;      
	float complex ** restrict TT;      

	float complex ** restrict pfp    = (float complex**)l1para->pfp; 
	float complex ** restrict pfx    = (float complex**)l1para->pfx;
	float complex ** restrict pfpT   = (float complex**)l1para->pfpT;
	float complex ** restrict pfxT   = (float complex**)l1para->pfxT;

	float complex  * restrict xPrimary = (float complex *)l1para->xPrimary;

	float * restrict obliqcorr = l1para->obliqcorr;
        const int fftnr=l1para->fftnr;

	float complex * restrict integration = (float complex  *)l1para->integration;

	if (imode == 0)
	{
		AA      = (float complex **)l1para->AP_p;
		TT      = (float complex **)l1para->TP_p;
	}
	else if (imode == 1)
	{
		AA      = (float complex **)l1para->AG_p;
		TT      = (float complex **)l1para->TG_p;

	}
	else if (imode == 2)
	{
		AA      = (float complex **)l1para->AS_p;
		TT      = (float complex **)l1para->TS_p;
	}

	if (imode == 0 || imode == 2)      for (ip=0;ip<MAX(nsgtrace,npx);ip++) obliqcorr[ip] = 1.0f;
	else if (imode == 1)               for (ip=0;ip<MAX(nsgtrace,npx);ip++) obliqcorr[ip] = -l1para->prfl2[0];

	if (l1para->zmethod == 0||l1para->zmethod == 2)
		for (ip=0;ip<npx;ip++) obliqcorr[ip] = l1para->ztop[ip];

	for (ifreq=0;ifreq<l1para->fftnc;ifreq++) integration[ifreq] = 1.0f;

	for(int ifreq=0; ifreq<hfid; ifreq++)
		if (l1para->zmethod==0) integration[ifreq] = (-I)*l1para->zwtilt[ifreq];

	FFT_R2C_PP(l1para->fftplan[0], taup, npx, pfp, fftnr)

		MAT_TRANSPOSE(pfp, pfpT, hfid, npx);

	taup_warp_solo_z2(nsgtrace, npx, hfid, pfpT, TT, xPrimary, l1para->zmethod, obliqcorr, integration, pfx);
	////taup_warp_solo(nsgtrace, npx, hfid, pfpT, TT, xPrimary, obliqcorr, pfx);

	/// above is reverse taup, and below is forward taup ==
	MAT_TRANSPOSE(pfx, pfxT, hfid, nsgtrace);

	taup_warp_solo_z1(npx, nsgtrace, hfid, pfxT, AA, xPrimary, l1para->zmethod, obliqcorr, integration, pfp);
	////taup_warp_solo(npx, nsgtrace, hfid, pfxT, AA, xPrimary, obliqcorr, pfp);

	for (ip=0;ip<npx;ip++) {   //mute high frequency, may not be necessary in the cg loop!
#pragma ivdep
		for (ifreq=hfid;ifreq<l1para->fftnc;ifreq++) 
			pfp[ip][ifreq] = 0.0f;
	}

	FFT_C2R_PP(l1para->fftplan[1], pfp, npx, taupou, fftnr)

		for (ip=0; ip<npx; ip++) {
			taupou[ip][fftnr-2] = 0.0f;
			taupou[ip][fftnr-1] = 0.0f;
		}

}

// This is for joint source receiver deghosting purpose 
// May change the name later 
	MSDGINTERP_WARP_DISPATCH
void c_l1inv_TAUP_RV_STOS(int nsgtrace, int hfid, int lfid, int npx, float ** taup, 
		float **RFX,  l1inv_t *l1para, int datamode, int datum)
{
	float complex ** restrict TP; 
	float complex ** restrict TG; 
	float complex ** restrict TS; 
	float complex ** restrict TGS; 
	float complex ** restrict TGR; 
	float complex ** restrict pfp;
	float complex ** restrict pfx;
	float complex ** restrict pfpT;
	float complex *  restrict xPrimary;
	float complex *  restrict xGhost;
	float complex *  restrict xGhostS;
	float complex *  restrict xGhostR;

	float complex **  FX;   //// <---
	float complex ** restrict pSrcFilter = (float complex **) l1para->pSrcFilter;
	float complex ** restrict pTgtFilter = (float complex **) l1para->pTgtFilter;
        const int fftnr=l1para->fftnr;

	pfp      = (float complex **)l1para->pfp;
	pfx      = (float complex **)l1para->pfx;
	pfpT     = (float complex **)l1para->pfpT;
	xPrimary = (float complex *)l1para->xPrimary;
	xGhost   = (float complex *)l1para->xGhost;
	xGhostS  = (float complex *)l1para->xGhostS;
	xGhostR  = (float complex *)l1para->xGhostR;

	FX      = (float complex **)RFX;

	int i,ip,ifreq,j; float * restrict obliqcorr = l1para->obliqcorr;

	for (int ip=0;ip<npx;ip++) obliqcorr[ip] = 1.0f;

	FFT_R2C_PP(l1para->fftplan[0], taup, npx, pfp, fftnr);


	if (datamode == 0 && l1para->computemode == 1)
		for (ip = 0; ip < npx; ip++)
			for (ifreq = 0; ifreq < hfid; ifreq++)
				pfp[ip][ifreq] *= pSrcFilter[ip][ifreq];

	MAT_TRANSPOSE(pfp, pfpT, hfid, npx);

	// datamode = 0 no redatum (for residual)
	// datamode = 1 redatum (for output)
	if (datamode == 0)
	{
		TP  = (float complex **)l1para->TP_p; 
		TG  = (float complex **)l1para->TG_p; 
		TS  = (float complex **)l1para->TS_p; 
		TGS = (float complex **)l1para->TGS_p; 
		TGR = (float complex **)l1para->TGR_p; 

		if (l1para->computemode == 0)
		{
			taup_warp_quad_p(nsgtrace, npx, hfid, pfpT, TP, xPrimary, l1para->prfl2, l1para->srcrfl, TG, xGhost,
					TGS, xGhostS, TGR, xGhostR, pfx);
		}
		else
		{
			if ((l1para->choosemethod==CHOOSEMETHOD_JOINTSR3D||l1para->choosemethod==CHOOSEMETHOD_JOINTSR2D) 
					&& (l1para->jointmethod == JOINTDESIGNATURE || l1para->jointmethod == JOINTDESIG_SRC_DEGHOST
						||l1para->jointmethod == JOINT_SRC_DEGHOST))
			{
				float *obliqcorr = (float *) calloc (MAX (nsgtrace, npx), sizeof (float));
				for (ip = 0; ip < MAX (nsgtrace, npx);ip++) obliqcorr[ip] = 1.0f;
				taup_warp_solo (nsgtrace, npx, hfid, pfpT, TS, xPrimary, obliqcorr, pfx);
				free (obliqcorr);
			}
			else
			{
				taup_warp_dual_p (nsgtrace, npx, hfid, pfpT, TP, xPrimary, l1para->prfl2, TG, xGhost, pfx);
			}
		}
	}
	else if (datamode == 1)
	{
		if (datum == CABLE90)
			TP = (float complex **)l1para->TP_p; 
		else if (datum == SURFACE0 || datum == SURFACE90)
			TP = (float complex **)l1para->TS_p; 

		if ((l1para->choosemethod==CHOOSEMETHOD_JOINTSR3D||l1para->choosemethod==CHOOSEMETHOD_JOINTSR2D) 
				&& (l1para->jointmethod == JOINTDESIGNATURE || l1para->jointmethod == JOINTDESIG_SRC_DEGHOST
					||l1para->jointmethod == JOINT_SRC_DEGHOST))
			TP = (float complex **)l1para->TS_p; 

		float complex phaseshift = -1.0f;
		if ((l1para->choosemethod==CHOOSEMETHOD_JOINTSR3D||l1para->choosemethod==CHOOSEMETHOD_JOINTSR2D) 
				&& (l1para->jointmethod == JOINTDESIG_SRC_DEGHOST || l1para->jointmethod == JOINTDESIG_RCV_DEGHOST
					||l1para->jointmethod == JOINT_SRC_DEGHOST))
			phaseshift = I;

		if (datum == SURFACE0)
			for (ifreq=0;ifreq<hfid;ifreq++)
				for (ip=0;ip<npx;ip++)
					pfpT[ifreq][ip] *= phaseshift;


		///////////////////////
		if (l1para->target_type == TARGETTYPE_GAPDECON)
		{
			for (ifreq = 0; ifreq < hfid; ifreq++)
				for (ip = 0; ip < npx; ip++)
					pfpT[ifreq][ip] *= pTgtFilter[ip][ifreq];
		}
		///////////////////////


		taup_warp_solo(nsgtrace, npx, hfid, pfpT, TP, xPrimary, obliqcorr, pfx);
	}

	for (i=0;i<nsgtrace;i++) //mute high frequency!
		for (j=hfid;j<l1para->fftnc;j++) 
			pfx[i][j] = 0.0f;

	if (l1para->lownear>0) {
		for (i=0;i<nsgtrace;i++) {  
			pfx[i][0] = 0.0f;
			pfx[i][1] = 0.0f;

			for (j=2;j<lfid;j++) 
				pfx[i][j] *= ((float)j-1.0f)/((float)lfid-2.0f);
		}
	}

	FFT_C2R_PP(l1para->fftplan[1], pfx, nsgtrace, RFX, fftnr);

	for (i=0;i<nsgtrace;i++) 
		FX[i][l1para->fftnc-1] = 0.0f;  //// <---
}

// This is for joint source receiver deghosting of multi sensor
// May change the name later 
	MSDGINTERP_WARP_DISPATCH
void c_l1inv_TAUP_RV_STOS_MS(int nsgtrace, int hfid, int lfid, int npx, float ** taup, 
		float **RFX,  l1inv_t *l1para, int datamode, int datum, int outtype)
{
	float complex ** restrict TP; 
	float complex ** restrict TG; 
	float complex ** restrict TGS; 
	float complex ** restrict TGR; 
	float complex ** restrict pfp;
	float complex ** restrict pfx;
	float complex ** restrict pfpT;
	float complex *  restrict xPrimary;
	float complex *  restrict xGhost;
	float complex *  restrict xGhostS;
	float complex *  restrict xGhostR;
        const int fftnr=l1para->fftnr;

	float complex **  FX;   //// <---
	float complex ** restrict pSrcFilter = (float complex **) l1para->pSrcFilter;
	float complex * restrict integration = (float complex  *)l1para->integration;
	float complex ** restrict pTgtFilter = (float complex **) l1para->pTgtFilter;

	pfp      = (float complex **)l1para->pfp;
	pfx      = (float complex **)l1para->pfx;
	pfpT     = (float complex **)l1para->pfpT;
	xPrimary = (float complex *)l1para->xPrimary;
	xGhost   = (float complex *)l1para->xGhost;
	xGhostS  = (float complex *)l1para->xGhostS;
	xGhostR  = (float complex *)l1para->xGhostR;

	FX      = (float complex **)RFX;

	int i,ip,ifreq,j,method; float * restrict obliqcorr = l1para->obliqcorr;
	float * restrict wtilt;
	float * restrict reflect; // zxue

	for (int ip=0;ip<npx;ip++) obliqcorr[ip] = 1.0f;
	for (ifreq=0;ifreq<l1para->fftnc;ifreq++) integration[ifreq] = 1.0f;

	if (datamode == 1) // az
	{
		wtilt   = (float *)l1para->zwtilt;
		method  = l1para->zmethod;	
		reflect = (float *)l1para->zrfl; // zxue
		if (l1para->zmethod == 0||l1para->zmethod == 2)
			for (ip=0;ip<npx;ip++) obliqcorr[ip] = l1para->ztop[ip];
	}
	else if (datamode == 2) // ay
	{
		wtilt   = (float *)l1para->ywtilt;
		method  = l1para->ymethod;	
		reflect = (float *)l1para->yrfl; // zxue
		if (l1para->ymethod == 0||l1para->ymethod == 2)
			for (ip=0;ip<npx;ip++) obliqcorr[ip] = l1para->ytop[ip];
	}

	FFT_R2C_PP(l1para->fftplan[0], taup, npx, pfp, fftnr);

	// outtype = 0 no redatum (for residual)
	// outtype = 1 redatum (for output)
	if (outtype == 0)
		for (ip = 0; ip < npx; ip++)
			for (ifreq = 0; ifreq < hfid; ifreq++)
				pfp[ip][ifreq] *= pSrcFilter[ip][ifreq];

	MAT_TRANSPOSE(pfp, pfpT, hfid, npx);

	if (outtype == 0)
	{
		TP  = (float complex **)l1para->TP_p; 
		TG  = (float complex **)l1para->TG_p; 

		if (datamode == 0) // p
		{
			taup_warp_dual_p (nsgtrace, npx, hfid, pfpT, TP, xPrimary, l1para->prfl2, TG, xGhost, pfx);
		}
		else if (datamode == 1 || datamode == 2) // az || ay
		{
			if (method==0)
				for (ifreq=0;ifreq<hfid;ifreq++)
					for (ip=0;ip<npx;ip++)
						pfpT[ifreq][ip] = (+I)*pfpT[ifreq][ip]*wtilt[ifreq];

			taup_warp_dual_azay2(nsgtrace, npx, hfid, pfpT, TP, xPrimary, reflect, TG, xGhost, 
					method, obliqcorr, integration, pfx);
		}
	}
	else if (outtype == 1)
	{
		if (datum == CABLE90)
			TP = (float complex **)l1para->TP_p; 
		else if (datum == SURFACE0 || datum == SURFACE90)
			TP = (float complex **)l1para->TS_p; 

		if (datum == SURFACE0)
			for (ifreq=0;ifreq<hfid;ifreq++)
				for (ip=0;ip<npx;ip++)
					pfpT[ifreq][ip] *= -1.0f;

		if ((datamode == 1)&&(l1para->zmethod==0)||((datamode == 2)&&(l1para->ymethod==0)))
			for (ifreq=0;ifreq<hfid;ifreq++)
				for (ip=0;ip<npx;ip++)
					pfpT[ifreq][ip] = (I)*pfpT[ifreq][ip]*wtilt[ifreq];

		///////////////////////
		if (l1para->target_type == TARGETTYPE_GAPDECON) 
		{
			for (ifreq = 0; ifreq < hfid; ifreq++)
				for (ip = 0; ip < npx; ip++)
					pfpT[ifreq][ip] *= pTgtFilter[ip][ifreq];
		}
		///////////////////////

		taup_warp_solo(nsgtrace, npx, hfid, pfpT, TP, xPrimary, obliqcorr, pfx);
	}

	for (i=0;i<nsgtrace;i++) //mute high frequency!
		for (j=hfid;j<l1para->fftnc;j++) 
			pfx[i][j] = 0.0f;

	if (l1para->lownear>0) {
		for (i=0;i<nsgtrace;i++) {  
			pfx[i][0] = 0.0f;
			pfx[i][1] = 0.0f;

			for (j=2;j<lfid;j++) 
				pfx[i][j] *= ((float)j-1.0f)/((float)lfid-2.0f);
		}
	}

	FFT_C2R_PP(l1para->fftplan[1], pfx, nsgtrace, RFX, fftnr);

	for (i=0;i<nsgtrace;i++) 
		FX[i][l1para->fftnc-1] = 0.0f;  //// <---
}

	MSDGINTERP_WARP_DISPATCH
void c_l1inv_TAUP_RV_CTOS(int nsgtrace, int hfid, int lfid, int npx, float ** taup, 
		float **RFX,  l1inv_t *l1para, int datamode, int datum)
{

	float complex ** restrict TP; 
	float complex ** restrict pfp;
	float complex ** restrict pfx;
	float complex ** restrict pfpT;
	float complex *  restrict xPrimary;
        const int fftnr=l1para->fftnr;

	float complex **  FX;   //// <---

	pfp      = (float complex **)l1para->pfp;
	pfx      = (float complex **)l1para->pfx;
	pfpT     = (float complex **)l1para->pfpT;
	xPrimary = (float complex *)l1para->xPrimary;

	FX      = (float complex **)RFX;

	if (datum == CABLE90)
		TP       = (float complex **)l1para->TP_p; 
	else if (datum == SURFACE0 || datum == SURFACE90)
		TP       = (float complex **)l1para->TS_p; 
	else if (datum == DGRGCABLE90)
		TP       = (float complex **)l1para->TPR_p; 

	int i,ip,ifreq,j; float * restrict obliqcorr = l1para->obliqcorr;

	float * restrict wtilt;

	for (int ip=0;ip<npx;ip++) obliqcorr[ip] = 1.0f;

	if (datamode == 1)
	{
		wtilt = (float *)l1para->zwtilt;

		if (l1para->zmethod == 0||l1para->zmethod == 2)
			for (ip=0;ip<npx;ip++) obliqcorr[ip] = l1para->ztop[ip];
	}
	else if (datamode == 2)
	{
		wtilt = (float *)l1para->ywtilt;

		if (l1para->ymethod == 0||l1para->ymethod == 2)
			for (ip=0;ip<npx;ip++) obliqcorr[ip] = l1para->ytop[ip];
	}

	FFT_R2C_PP(l1para->fftplan[0], taup, npx, pfp, fftnr);

	MAT_TRANSPOSE(pfp, pfpT, hfid, npx);

	if (datum == SURFACE0)
		for (ifreq=0;ifreq<hfid;ifreq++)
			for (ip=0;ip<npx;ip++)
				pfpT[ifreq][ip] = (I)*pfpT[ifreq][ip];

	if ((datamode == 1)&&(l1para->zmethod==0)||((datamode == 2)&&(l1para->ymethod==0)))
	{
		for (ifreq=0;ifreq<hfid;ifreq++)
		{
			for (ip=0;ip<npx;ip++)
			{
				pfpT[ifreq][ip] = (I)*pfpT[ifreq][ip]*wtilt[ifreq];
			}
		}
	}

	taup_warp_solo(nsgtrace, npx, hfid, pfpT, TP, xPrimary, obliqcorr, pfx);

	if (l1para->gundelay>0.0)
	{
		float complex z;

		float fnyq = 500000.0f/l1para->srate;
		float deltaf = fnyq/(float)(l1para->fftnc-1);

		for (i=0;i<nsgtrace;i++)
		{
			for (j=0;j<l1para->fftnc;j++)
			{
				z = (float complex) (cosf(2*PI*j*deltaf*l1para->gundelay),-sinf(2*PI*j*deltaf*l1para->gundelay));
				pfx[i][j] = pfx[i][j]*z;
			}
		}
	}

	for (i=0;i<nsgtrace;i++) //mute high frequency!
		for (j=hfid;j<l1para->fftnc;j++) 
			pfx[i][j] = 0.0f;

	if (l1para->lownear>0) {
		for (i=0;i<nsgtrace;i++) {  
			pfx[i][0] = 0.0f;
			pfx[i][1] = 0.0f;

			for (j=2;j<lfid;j++) 
				pfx[i][j] *= ((float)j-1.0f)/((float)lfid-2.0f);
		}
	}


	FFT_C2R_PP(l1para->fftplan[1], pfx, nsgtrace, RFX, fftnr);

	for (i=0;i<nsgtrace;i++) 
		FX[i][l1para->fftnc-1] = 0.0f;  //// <---

}

	MSDGINTERP_WARP_DISPATCH
void c_l1inv_TAUP_RV_MTOS(int nsgtrace, int hfid, int lfid, int npx, float ** taup, 
		float **RFX,  l1inv_t *l1para, int datamode, int datum)
{
	int ip,ifreq,i,j;  

	float complex ** restrict TG; 
	float complex ** restrict pfp;
	float complex ** restrict pfx;
	float complex ** restrict pfpT;
	float complex *  restrict xGhost;
        const int fftnr=l1para->fftnr;

	float complex **  FX;  //// <---

	pfp      = (float complex **)l1para->pfp;
	pfx      = (float complex **)l1para->pfx;
	pfpT     = (float complex **)l1para->pfpT;
	xGhost = (float complex *)l1para->xGhost;

	FX      = (float complex **)RFX;

	if (datum == DGRGCABLE90)
		TG       = (float complex **)l1para->TGR_p; 
	else
		TG       = (float complex **)l1para->TG_p; 

	float * restrict obliqcorr = l1para->obliqcorr;

	float * restrict wtilt;

	for (ip=0;ip<npx;ip++) obliqcorr[ip] = 1.0f;

	if (datamode == 1)
	{
		wtilt = (float *)l1para->zwtilt;
		if (l1para->zmethod==1) 
		{
			for (ip=0;ip<npx;ip++) obliqcorr[ip] = -1.f; // zxue     
		}
		else 
		{
			for (ip=0;ip<npx;ip++) obliqcorr[ip] = -l1para->ztop[ip]; // zxue
		}
	}
	if (datamode == 2)
	{
		wtilt = (float *)l1para->ywtilt;
		if (l1para->ymethod==1) 
		{
			for (ip=0;ip<npx;ip++) obliqcorr[ip] = -1.f; // zxue // bug fix with missing for loop 
		}
		else 
		{
			for (ip=0;ip<npx;ip++) obliqcorr[ip] = -l1para->ytop[ip];
		}
	}
	else if (datamode == 0)
	{
		for (ip=0;ip<npx;ip++) obliqcorr[ip] = -1.f; // zxue
	}

	FFT_R2C_PP(l1para->fftplan[0], taup, npx, pfp, fftnr);

	MAT_TRANSPOSE(pfp, pfpT, hfid, npx);

	if ((datamode == 1&&l1para->zmethod==0)||(datamode == 2&&l1para->ymethod==0))
	{
		for (ifreq=0;ifreq<hfid;ifreq++)
		{
			for (ip=0;ip<npx;ip++)
			{
				pfpT[ifreq][ip] = (I)*pfpT[ifreq][ip]*wtilt[ifreq];
			}
		}
	}

	taup_warp_solo(nsgtrace, npx, hfid, pfpT, TG, xGhost, obliqcorr, pfx);

	// zxue
	if(datamode == 1){
		for (ip=0;ip<nsgtrace;ip++)
			for (ifreq=0; ifreq<hfid; ifreq++)
				pfx[ip][ifreq] *= l1para->zrfl[ifreq];
	}
	if(datamode == 2){
		for (ip=0;ip<nsgtrace;ip++)
			for (ifreq=0; ifreq<hfid; ifreq++)
				pfx[ip][ifreq] *= l1para->yrfl[ifreq];
	}
	if(datamode == 0){
		for (ip=0;ip<nsgtrace;ip++)
			for (ifreq=0; ifreq<hfid; ifreq++)
				pfx[ip][ifreq] *= l1para->prfl2[ifreq];
	}

	if (l1para->gundelay>0.0)
	{
		float complex z;

		float fnyq = 500000.0f/l1para->srate;
		float deltaf = fnyq/(float)(l1para->fftnc-1);

		for (i=0;i<nsgtrace;i++)
		{
			for (j=0;j<l1para->fftnc;j++)
			{
				z = (float complex) (cosf(2*PI*j*deltaf*l1para->gundelay),-sinf(2*PI*j*deltaf*l1para->gundelay));
				pfx[i][j] = pfx[i][j]*z;
			}
		}
	}

	for (i=0;i<nsgtrace;i++) //mute high frequency!
		for (j=hfid;j<l1para->fftnc;j++) 
			pfx[i][j] = 0.0f;

	if (l1para->lownear>0) {
		for (i=0;i<nsgtrace;i++) {  
			pfx[i][0] = 0.0f;
			pfx[i][1] = 0.0f;

			for (j=2;j<lfid;j++) 
				pfx[i][j] *= ((float)j-1.0f)/((float)lfid-2.0f);
		}
	}


	FFT_C2R_PP(l1para->fftplan[1], pfx, nsgtrace, RFX, fftnr);

	for (i=0;i<nsgtrace;i++) 
		FX[i][l1para->fftnc-1] = 0.0f;  //// <---

}

	MSDGINTERP_WARP_DISPATCH
void c_l1inv_TAUP_RV_BASIC(int nsgtrace, int hfid, int lfid, int npx, float ** taup, 
		float **RFX,  float ** TPP, l1inv_t *l1para, int datum)
{
	int ifreq,ip; 
	float complex ** restrict TP;  
	float complex ** restrict pfp;
	float complex ** restrict pfx;
	float complex ** restrict pfpT;
	float complex *  restrict xPrimary;
        const int fftnr=l1para->fftnr;

	pfp      = (float complex **)l1para->pfp;
	pfx      = (float complex **)l1para->pfx;
	pfpT     = (float complex **)l1para->pfpT;
	xPrimary = (float complex *)l1para->xPrimary;

	TP       = (float complex **)TPP;

	int i,j; float * restrict obliqcorr = l1para->obliqcorr;

	for (ip=0;ip<npx;ip++) obliqcorr[ip] = 1.0f;

	FFT_R2C_PP(l1para->fftplan[0], taup, npx, pfp, fftnr);

	MAT_TRANSPOSE(pfp, pfpT, hfid, npx);

	if (datum == SURFACE0)
		for (ifreq=0;ifreq<hfid;ifreq++)
			for (ip=0;ip<npx;ip++)
				pfpT[ifreq][ip] = (I)*pfpT[ifreq][ip];

	taup_warp_solo(nsgtrace, npx, hfid, pfpT, TP, xPrimary, obliqcorr, pfx);

	for (i=0;i<nsgtrace;i++) //mute high frequency!
		for (j=hfid;j<l1para->fftnc;j++) 
			pfx[i][j] = 0.0f;

	if (l1para->lownear>0) {
		for (i=0;i<nsgtrace;i++) {  
			pfx[i][0] = 0.0f;
			pfx[i][1] = 0.0f;

			for (j=2;j<lfid;j++) 
				pfx[i][j] *= ((float)j-1.0f)/((float)lfid-2.0f);
		}
	}

	FFT_C2R_PP(l1para->fftplan[1], pfx, nsgtrace, RFX, fftnr);

}

	MSDGINTERP_WARP_DISPATCH
void c_l1inv_lsqr(int imode, float repstop, int niter, 
		int nsgtrace, int hfid, int lfid, int npx, float sppow, 
		float ** pweight, float ** dataxt_p, float ** dataxt_az, float ** dataxt_ay,
		float **taup, float **tweight, l1inv_t *l1para)
{   

	if (l1para->cgmethod == 0)
		c_l1inv_lsqr_cg(imode, repstop, niter, nsgtrace, hfid, lfid, npx, sppow, 
				pweight, dataxt_p, dataxt_az, dataxt_ay,
				taup, tweight, l1para);
	else
		c_l1inv_lsqr_bcg(imode, repstop, niter, nsgtrace, hfid, lfid, npx, sppow, 
				pweight, dataxt_p, dataxt_az, dataxt_ay,
				taup, tweight, l1para);

}

	MSDGINTERP_WARP_DISPATCH
void c_l1inv_lsqr_bcg(int imode, float repstop, int niter, 
		int nsgtrace, int hfid, int lfid, int npx, float sppow, 
		float ** pweight, float ** dataxt_p, float ** dataxt_az, float ** dataxt_ay,
		float **taup, float **tweight, l1inv_t *l1para)
{   

	//// t,As,nu,ws can be potentially removed

	int i, j, iter;
	float pTAp, gammar, alpha, beta, gammar_last, gammar0;
        const int fftnr=l1para->fftnr;

	float ** restrict p  = (float **)l1para->cg_d;      // bbb->d
	float ** restrict r  = (float **)l1para->cg_r;      // ddd->r
	float ** restrict s  = (float **)l1para->cg_s;      
	float ** restrict Ap = (float **)l1para->cg_Ad;     // ccc->Ad
	float ** restrict As = (float **)l1para->cg_As; 
	float ** restrict rhat0 = (float **)l1para->cg_rhat0;     
	float ** wp = (imode!=MSDGINTERP_CG_WEIGHTED)?(p):(l1para->cg_wd);    //an alias if not weighted!
	float ** ws = (imode!=MSDGINTERP_CG_WEIGHTED)?(s):(l1para->cg_wd);    //an alias if not weighted!

	float rho,rho_last,omega;

	if (l1para->fpredatum==FPREDATUM_NO||l1para->fpredatum==FPREDATUM_DGRG) 
		c_l1inv_TAUP_FW_DG  (nsgtrace, hfid, lfid, npx, dataxt_p, dataxt_az, dataxt_ay, taup, l1para);
	else if (l1para->fpredatum==FPREDATUM_CTOS||l1para->fpredatum==FPREDATUM_CIN)
		c_l1inv_TAUP_FW_CMSTOS_STOS_p(0, nsgtrace, hfid, lfid, npx, dataxt_p, taup, l1para);
	else if (l1para->fpredatum==FPREDATUM_XG_Z)
		c_l1inv_TAUP_FW_CMSTOS_STOS_z(0, nsgtrace, hfid, lfid, npx, dataxt_az, taup, l1para);
	else if (l1para->fpredatum==FPREDATUM_MTOS)
		c_l1inv_TAUP_FW_CMSTOS_STOS_p(1, nsgtrace, hfid, lfid, npx, dataxt_p, taup, l1para);
	////else if (l1para->fpredatum==FPREDATUM_SIN)
	else if (l1para->fpredatum==FPREDATUM_STOM||l1para->fpredatum==FPREDATUM_STOC||l1para->fpredatum==FPREDATUM_SIN) 
		c_l1inv_TAUP_FW_CMSTOS_STOS_p(2, nsgtrace, hfid, lfid, npx, dataxt_p, taup, l1para);

	float ** b = (float **)taup;   //// Solving A x = b

	if (imode==MSDGINTERP_CG_WEIGHTED) {
		for(i=0; i<npx; i++) 
			for(j=0; j<fftnr; j++)
				b[i][j] *= pweight[i][j];
	}

	for(i=0; i<npx; i++)
		for(j=0; j<fftnr; j++)
			r[i][j] = b[i][j];   

	for(i=0; i<npx; i++)
		for(j=0; j<fftnr; j++)
			rhat0[i][j] = r[i][j];

	rho = rho_last = alpha = omega = 1.0f;

	for(i=0; i<npx; i++)
		for(j=0; j<fftnr; j++)
			Ap[i][j] = p[i][j] = 0.0f;   //// Ap = nu

	float ** restrict x = (imode==MSDGINTERP_CG_WEIGHTED)?(taup):(pweight);
	for(i=0; i<npx; i++)
		for(j=0; j<fftnr; j++)
			x[i][j] = 0.0f;

	gammar = gammar0 = gammar_last = dotp2d_solo(npx, fftnr, r);
	if (sqrt(gammar)<1.0e-15) 
		return;

	for (iter=1;iter<=niter;iter++) {

		rho = dotp2d_dual(npx, fftnr, rhat0, r);

		beta = (rho/rho_last)*(alpha/omega);

		for (i=0;i<npx;i++) 
			for (j=0;j<fftnr;j++) {
				p[i][j] = r[i][j] + beta*(p[i][j]-omega*Ap[i][j]);
			}

		if (imode==MSDGINTERP_CG_WEIGHTED) {
			for (i=0;i<npx;i++) 
				for (j=0;j<fftnr;j++) 
					wp[i][j] = p[i][j] * pweight[i][j]; 
		}

		if (l1para->fpredatum==FPREDATUM_NO||l1para->fpredatum==FPREDATUM_DGRG) 
			c_l1inv_matrix_TAUP_FW_RV_DG    (nsgtrace, hfid, npx, wp, Ap, l1para);//// p -> A p
		else if (l1para->fpredatum==FPREDATUM_CTOS||l1para->fpredatum==FPREDATUM_CIN)
			c_l1inv_matrix_TAUP_FW_RV_CMSTOS_STOS_p(0, nsgtrace, hfid, npx, wp, Ap, l1para);
		else if (l1para->fpredatum==FPREDATUM_XG_Z)
			c_l1inv_matrix_TAUP_FW_RV_CMSTOS_STOS_z(0, nsgtrace, hfid, npx, wp, Ap, l1para);
		else if (l1para->fpredatum==FPREDATUM_MTOS)
			c_l1inv_matrix_TAUP_FW_RV_CMSTOS_STOS_p(1, nsgtrace, hfid, npx, wp, Ap, l1para);
		////else if (l1para->fpredatum==FPREDATUM_SIN)
		else if (l1para->fpredatum==FPREDATUM_STOM||l1para->fpredatum==FPREDATUM_STOC||l1para->fpredatum==FPREDATUM_SIN) 
			c_l1inv_matrix_TAUP_FW_RV_CMSTOS_STOS_p(2, nsgtrace, hfid, npx, wp, Ap, l1para);


		if (imode==MSDGINTERP_CG_WEIGHTED) {
			for (i=0;i<npx;i++) 
				for (j=0;j<fftnr;j++) {
					Ap[i][j] = Ap[i][j]*pweight[i][j];
				}
		}

		alpha = rho/(dotp2d_dual(npx, fftnr, rhat0, Ap));

		rho_last =  rho;

		for(i=0; i<npx; i++)
			for(j=0; j<fftnr; j++) {
				s[i][j] = r[i][j] -  alpha * Ap[i][j];
			}

		if (imode==MSDGINTERP_CG_WEIGHTED) {
			for (i=0;i<npx;i++) 
				for (j=0;j<fftnr;j++) 
					ws[i][j] = s[i][j] * pweight[i][j]; 
		}

		if (l1para->fpredatum==FPREDATUM_NO||l1para->fpredatum==FPREDATUM_DGRG) 
			c_l1inv_matrix_TAUP_FW_RV_DG    (nsgtrace, hfid, npx, ws, As, l1para);//// s -> A s
		else if (l1para->fpredatum==FPREDATUM_CTOS||l1para->fpredatum==FPREDATUM_CIN)
			c_l1inv_matrix_TAUP_FW_RV_CMSTOS_STOS_p(0, nsgtrace, hfid, npx, ws, As, l1para);
		else if (l1para->fpredatum==FPREDATUM_XG_Z)
			c_l1inv_matrix_TAUP_FW_RV_CMSTOS_STOS_z(0, nsgtrace, hfid, npx, ws, As, l1para);
		else if (l1para->fpredatum==FPREDATUM_MTOS)
			c_l1inv_matrix_TAUP_FW_RV_CMSTOS_STOS_p(1, nsgtrace, hfid, npx, ws, As, l1para);
		////else if (l1para->fpredatum==FPREDATUM_SIN)
		else if (l1para->fpredatum==FPREDATUM_STOM||l1para->fpredatum==FPREDATUM_STOC||l1para->fpredatum==FPREDATUM_SIN) 
			c_l1inv_matrix_TAUP_FW_RV_CMSTOS_STOS_p(2, nsgtrace, hfid, npx, ws, As, l1para);

		if (imode==MSDGINTERP_CG_WEIGHTED) {
			for (i=0;i<npx;i++) 
				for (j=0;j<fftnr;j++) {
					As[i][j] = As[i][j]*pweight[i][j]; //// t = As
				}
		}

		omega = dotp2d_dual(npx, fftnr, As, s)/dotp2d_solo(npx, fftnr, As);

		for(i=0; i<npx; i++)
			for(j=0; j<fftnr; j++) {
				x[i][j] += alpha*p[i][j]+omega*s[i][j];
				r[i][j] = s[i][j] - omega*As[i][j];
			}

		if (sqrt(gammar/gammar0) > 10.0f)
			G2_logCritical("bsdg","Unstable iteration: %d %f %f %f",iter,sqrt(gammar/gammar0),gammar,gammar0);

		if (sqrt(gammar/gammar0) <= repstop) {
			break;     
		}
		if (iter!=niter) {   //no update for last iteration
			gammar = dotp2d_solo(npx, fftnr, r);
		}
	} 

	if (imode==MSDGINTERP_CG_WEIGHTED) {
		for(i=0; i<npx; i++)
			for(j=0; j<fftnr; j++)
				x[i][j] *= pweight[i][j];
	} else {
		instantaneous_weight(l1para->hilbert_task, l1para->hilbert_hw, l1para->hilbert_nt,
				fftnr, npx, pweight, tweight, sppow, 
				l1para->hilbert_conv);  //TODO: move out of cg
	}
}

	MSDGINTERP_WARP_DISPATCH
void c_l1inv_lsqr_cg(int imode, float repstop, int niter, 
		int nsgtrace, int hfid, int lfid, int npx, float sppow, 
		float ** pweight, float ** dataxt_p, float ** dataxt_az, float ** dataxt_ay,
		float **taup, float **tweight, l1inv_t *l1para)
{   
	int i, j, iter;
	float gammar, alpha, beta, gammar_last, gammar0,tolcheck;
	double xcov;
        const int fftnr=l1para->fftnr;

	FILE *fp_rmsqc;
	if (l1para->rmsqc) fp_rmsqc = fopen(l1para->rmsqcfile,"a")
		;
	float ** restrict d  = (float **)l1para->cg_d;      // bbb->d
	float ** restrict r  = (float **)l1para->cg_r;      // ddd->r
	float ** restrict Ad = (float **)l1para->cg_Ad;     // ccc->Ad
	float ** wd = (imode!=MSDGINTERP_CG_WEIGHTED)?(d):(l1para->cg_wd);    //an alias if not weighted!

	if (l1para->fpredatum==FPREDATUM_NO||l1para->fpredatum==FPREDATUM_DGRG) 
		c_l1inv_TAUP_FW_DG  (nsgtrace, hfid, lfid, npx, dataxt_p, dataxt_az, dataxt_ay, taup, l1para);
	else if (l1para->fpredatum==FPREDATUM_CTOS||l1para->fpredatum==FPREDATUM_CIN)
		c_l1inv_TAUP_FW_CMSTOS_STOS_p(0, nsgtrace, hfid, lfid, npx, dataxt_p, taup, l1para);
	else if (l1para->fpredatum==FPREDATUM_XG_Z)
		c_l1inv_TAUP_FW_CMSTOS_STOS_z(0, nsgtrace, hfid, lfid, npx, dataxt_az, taup, l1para);
	else if (l1para->fpredatum==FPREDATUM_MTOS)
		c_l1inv_TAUP_FW_CMSTOS_STOS_p(1, nsgtrace, hfid, lfid, npx, dataxt_p, taup, l1para);
	////else if (l1para->fpredatum==FPREDATUM_SIN)
	else if (l1para->fpredatum==FPREDATUM_STOM||l1para->fpredatum==FPREDATUM_STOC||l1para->fpredatum==FPREDATUM_SIN) 
		c_l1inv_TAUP_FW_CMSTOS_STOS_p(2, nsgtrace, hfid, lfid, npx, dataxt_p, taup, l1para);

	float ** b = (float **)taup;

	if (imode==MSDGINTERP_CG_WEIGHTED) {
		for(i=0; i<npx; i++) 
			for(j=0; j<fftnr; j++)
				b[i][j] *= pweight[i][j];
	}
	for(i=0; i<npx; i++)
#pragma vector aligned
		for(j=0; j<fftnr; j++)
			d[i][j] = r[i][j] = b[i][j];
        //memcpy(d[0],b[0],npx*fftnr*sizeof(float));
        //memcpy(r[0],b[0],npx*fftnr*sizeof(float));

	float ** restrict x = (imode==MSDGINTERP_CG_WEIGHTED)?(taup):(pweight);

	for(i=0; i<npx; i++)
#pragma vector aligned  	
            for(j=0; j<fftnr; j++)
			x[i][j] = 0.0f;
        //memset(x[0],0,npx*fftnr);		

	gammar = gammar0 = gammar_last = dotp2d_solo(npx, fftnr, r);
	if (sqrt(gammar)<1.0e-15) 
		return;

	for (iter=1;iter<=niter;iter++) {

		if (imode==MSDGINTERP_CG_WEIGHTED) {
			for (i=0;i<npx;i++) 
				for (j=0;j<fftnr;j++) 
					wd[i][j] = d[i][j] * pweight[i][j]; 
		}

		if (l1para->fpredatum==FPREDATUM_NO||l1para->fpredatum==FPREDATUM_DGRG) 
			c_l1inv_matrix_TAUP_FW_RV_DG    (nsgtrace, hfid, npx, wd, Ad, l1para);
		else if (l1para->fpredatum==FPREDATUM_CTOS||l1para->fpredatum==FPREDATUM_CIN)
			c_l1inv_matrix_TAUP_FW_RV_CMSTOS_STOS_p(0, nsgtrace, hfid, npx, wd, Ad, l1para);
		else if (l1para->fpredatum==FPREDATUM_XG_Z)
			c_l1inv_matrix_TAUP_FW_RV_CMSTOS_STOS_z(0, nsgtrace, hfid, npx, wd, Ad, l1para);
		else if (l1para->fpredatum==FPREDATUM_MTOS)
			c_l1inv_matrix_TAUP_FW_RV_CMSTOS_STOS_p(1, nsgtrace, hfid, npx, wd, Ad, l1para);
		////else if (l1para->fpredatum==FPREDATUM_SIN)
		else if (l1para->fpredatum==FPREDATUM_STOM||l1para->fpredatum==FPREDATUM_STOC||l1para->fpredatum==FPREDATUM_SIN) 
			c_l1inv_matrix_TAUP_FW_RV_CMSTOS_STOS_p(2, nsgtrace, hfid, npx, wd, Ad, l1para);

		if (imode==MSDGINTERP_CG_WEIGHTED) {
			xcov = 0.0;
			for (i=0;i<npx;i++) 
				for (j=0;j<fftnr;j++) {
					Ad[i][j] *= pweight[i][j];
					xcov += d[i][j] * Ad[i][j];
				} 
			alpha = gammar / xcov; 
		} else {
			alpha = gammar / dotp2d_dual(npx, fftnr, d, Ad);
		}

		for(i=0; i<npx; i++){
#pragma ivdep
			for(j=0; j<fftnr; j++) {
				x[i][j] += alpha*d[i][j];
				r[i][j] -= alpha*Ad[i][j];
			}
                }

		tolcheck = sqrtf(gammar/gammar0);
		if(l1para->rmsqc) fprintf(fp_rmsqc,"%3d %.6f\n",iter,tolcheck);

		if (tolcheck > 10.0f)
			G2_logCritical("bsdg","Unstable iteration: %d %f %f %f",iter,sqrt(gammar/gammar0),gammar,gammar0);

		if (sqrt(gammar/gammar0) <= repstop) {
			break;     
		}
		if (iter!=niter) {   //no update for last iteration

			gammar = dotp2d_solo(npx, fftnr, r);

			beta = gammar/gammar_last;
			for (i=0;i<npx;i++) 
				for (j=0;j<fftnr;j++) 
					d[i][j]=r[i][j]+beta*d[i][j];
			gammar_last = gammar;

			if (imode==MSDGINTERP_CG_WEIGHTED) {
				for(i=0; i<npx; i++)
					for(j=0; j<fftnr; j++)
						wd[i][j] = d[i][j]*pweight[i][j];
			}
		}
	} 

	if (imode==MSDGINTERP_CG_WEIGHTED) {
		for(i=0; i<npx; i++)
			for(j=0; j<fftnr; j++)
				x[i][j] *= pweight[i][j];
	} else {
		instantaneous_weight(l1para->hilbert_task, l1para->hilbert_hw, l1para->hilbert_nt,
				fftnr, npx, pweight, tweight, sppow, 
				l1para->hilbert_conv);  //TODO: move out of cg
	}
	if (l1para->rmsqc) fclose(fp_rmsqc);
}

	MSDGINTERP_WARP_DISPATCH
void ms_lowcutf(float **input, int nsgtrace, l1inv_t *l1para, int datamode)
{

	int itrc,ifreq;
        const int fftnr=l1para->fftnr;

	float complex ** restrict pfx      = (float complex **)l1para->pfx;

	float * restrict xzlowf;

	if (datamode == 1)
		xzlowf = (float *)l1para->zlowf;
	else if (datamode == 2)
		xzlowf = (float *)l1para->ylowf;

	FFT_R2C_PP(l1para->fftplan[0], input, nsgtrace, pfx, fftnr);

	for (itrc=0;itrc<nsgtrace;itrc++)
		for (ifreq=0;ifreq<l1para->fftnc;ifreq++)
			pfx[itrc][ifreq] = pfx[itrc][ifreq]*xzlowf[ifreq];


	FFT_C2R_PP(l1para->fftplan[1], pfx, nsgtrace, input, fftnr);

}

	MSDGINTERP_WARP_DISPATCH
void l1inv3d_taup_merge(l1inv_t *l1para, int hfid, int npx, float **taup, float **taup1)
{
	int i,j,ip,ifreq;
	float x;
        const int fftnr=l1para->fftnr;

	float complex ** restrict pfp      = (float complex**)l1para->pfp;
	float complex ** restrict pfpT     = (float complex**)l1para->pfpT;

	FFT_R2C_PP(l1para->fftplan[0], taup, npx, pfp, fftnr);

	MAT_TRANSPOSE(pfp, pfpT, hfid, npx);

	FFT_R2C_PP(l1para->fftplan[0], taup1, npx, pfp, fftnr);

	for (ip=0;ip<npx;ip++)
	{
		for (ifreq=0;ifreq<l1para->fftnc;ifreq++)
		{
			x=l1para->zlowf[ifreq];
			pfp[ip][ifreq] =  pfp[ip][ifreq]*(1.0f-x)+pfpT[ifreq][ip]*x;
			////pfp[ip][ifreq] =  pfp[ip][ifreq]*(1.0f-x)+pfpT[ifreq][ip]*x;
		}
	}

	FFT_C2R_PP(l1para->fftplan[1], pfp, npx, taup, fftnr);

	for (i=0;i<npx;i++)
		for (j=l1para->nsamp;j<=fftnr-1;j++)
			taup[i][j] = 0.0f;

	return;
}

	MSDGINTERP_WARP_DISPATCH
float compute_scalar(int ip, int ifreq, float recz, l1inv_t *l1para)
{

	float px,py,pxy,tg,fsin,fcos,x,tmpr,tmpi,fmid,fmid_n,fmid_p;

	float invwater = 1.0f/l1para->vwater;

	float fnyq = 500000.0f/l1para->srate;
	float deltaf = fnyq/(float)(l1para->fftnc-1);
	float twopif=deltaf*2.0f*3.1415927f;

	int ftap = 5.0f/deltaf;
	float ftaphalf = 0.5f*l1para->ftap/deltaf;

	int inotch0,inotch1;
	float fnotch0,fnotch1;
	float fpeak0,fpeak1;
	float tapbeg,tapend,frac;

	if (l1para->ipn_type != -1)
	{
		if (l1para->ipn_type == 1)      //// DSDGOUT_Puppeak
		{
			inotch1 = (l1para->ipn + 1)/2;
			inotch0 = MAX(0,inotch1 - 1);
		}
		else                           //// DSDGOUT_Pupnotch
		{
			inotch1 = (l1para->ipn)/2;
			inotch0 = MAX(0,inotch1 - 1);
		}
	}

	px = l1para->pxsave[ip];
	py = l1para->pysave[ip];
	pxy = sqrtf(px*px+py*py);

	fsin=MAX(-0.99999f,MIN(0.99999f,pxy*l1para->vwater));
	fcos=sqrtf(1.0f-fsin*fsin);
	tg = MAX(l1para->ztminscl*0.001f,2.0f*fcos*recz*invwater);

	if (l1para->ipn_type != -1)
	{

		fnotch0 = (float)(inotch0)*(1.0f/tg/deltaf);
		fnotch1 = (float)(inotch1)*(1.0f/tg/deltaf);

		fmid_n = 0.5f*(fnotch0+fnotch1);

		if (inotch1 == 0)   //// should be relevant only for Pupnotch and Zuppeak
			fpeak0 = 0.0f;
		else
			fpeak0 = fnotch0 + (0.5f/tg/deltaf);

		fpeak1 = fnotch1 + (0.5f/tg/deltaf);

		fmid_p = 0.5f*(fpeak0+fpeak1);

		if (l1para->ipn_type == 1)
		{

			tapbeg = fmid_n-ftaphalf;
			tapend = fmid_n+ftaphalf;

			if (ifreq >= tapbeg && ifreq <= tapend)  ////Puppeak; also Zupnotch ?
			{
				frac = ((float)ifreq - tapbeg)/ftaphalf;
				x=sinf(0.5*PI*frac);
				x = MAX(0.0f,x);
				return 1.0f-x;
			}
			else
			{
				x = 0.0f;
				return 1.0f-x;
			}	  

		}
		else
		{
			tapbeg = fmid_p-ftaphalf;
			tapend = fmid_p+ftaphalf;

			if (ifreq >= tapbeg && ifreq <= tapend)  ////Pupnotch; also Zuppeak ?
			{
				frac = ((float)ifreq - tapbeg)/ftaphalf;
				x=sinf(0.5*PI*frac);
				x = MAX(0.0f,x);
				return x;
			}
			else
			{
				x = 0.0f;
				return x;
			}

		}
	}

	if (l1para->zlowmerge<1.0e-15)
		fmid = 0.5f/tg/deltaf;
	else if (l1para->zlowmerge<1.0e-2)
		fmid = 0.0f;
	else
		fmid = l1para->zlowmerge/deltaf;

	x=0.5f*(1.0f+cosf(twopif*ifreq*tg));
	x=powf(x,l1para->zpower);
	if (fmid>0)
	{
		if (ifreq<fmid) x=0.0f;
		else if (ifreq<fmid+ftap)
			x=x*0.5f*(1.0f-cosf(PI*(ifreq-fmid)*1.0f/ftap));
	}

	return x;

}

	MSDGINTERP_WARP_DISPATCH
void taup_warp_solo_dsdg3d(int nsgtrace, int npx, int hfid, float **pfpT0in, float **pfpT1in, 
		float **TGG, float *xxGhost, float *obliqcorr, float *recz, float **pfxou, l1inv_t *l1para)
{

	int ix,ip,ifreq;

	float complex * restrict xGhost;

	float complex ** restrict TG; 

	float complex ** restrict pfpT0;
	float complex ** restrict pfpT1;

	float complex ** restrict pfx;

	float complex tmp;
	float complex tmp1;

	float x;

	TG       = (float complex **)TGG; 

	xGhost   = (float complex *)xxGhost;

	pfpT0    = (float complex **)pfpT0in;
	pfpT1    = (float complex **)pfpT1in;

	pfx      = (float complex **)pfxou;

	for(ix=0; ix<nsgtrace; ix++) 
	{					
		for(ip=0; ip<npx; ip++)            
			xGhost[ip] = 1.0f*obliqcorr[ip];                                                 

		for(ifreq=0; ifreq<hfid; ifreq++) 
		{                                       
			tmp = 0.0f;                                                               
			for(ip=0; ip<npx; ip++) 
			{            
				if (l1para->dsdgout==DSDGOUT_P||l1para->dsdgout==DSDGOUT_Z)
					x = 1.0f;
				else
					x = compute_scalar(ip, ifreq, recz[ix], l1para);

				if (l1para->dsdgout==DSDGOUT_PZ)
					tmp1=pfpT0[ifreq][ip]*(1.0f-x)+pfpT1[ifreq][ip]*x;
				else if (l1para->dsdgout==DSDGOUT_P||l1para->dsdgout==DSDGOUT_Pnotch||l1para->dsdgout==DSDGOUT_Pupnotch)
					tmp1=pfpT0[ifreq][ip]*x;
				else if (l1para->dsdgout==DSDGOUT_Z)
					tmp1=pfpT0[ifreq][ip]*x;
				else if (l1para->dsdgout==DSDGOUT_Ppeak||l1para->dsdgout==DSDGOUT_Puppeak)
					tmp1=pfpT0[ifreq][ip]*(1.0f-x);
				else if (l1para->dsdgout==DSDGOUT_Znotch||l1para->dsdgout==DSDGOUT_Zupnotch)
					tmp1=pfpT0[ifreq][ip]*(1.0f-x);
				else if (l1para->dsdgout==DSDGOUT_Zpeak||l1para->dsdgout==DSDGOUT_Zuppeak)
					tmp1=pfpT0[ifreq][ip]*x;

				tmp += xGhost[ip]*tmp1;  

				xGhost[ip] *= TG[ix][ip];                                           
			} 
			pfx[ix][ifreq] = tmp;                                                   
		} 
	}

}

	MSDGINTERP_WARP_DISPATCH
void c_dsdg3d_TAUP_RV_CTOS(int nsgtrace, int hfid, int lfid, int npx, float *recz, float **taup_p, float **taup_az,
		float **RFX,  l1inv_t *l1para, int datamode, int datum)
{

	float complex ** restrict pfp0;
	float complex ** restrict pfp1;

	float complex ** restrict pfx;

	float complex ** restrict pfpT0;
	float complex ** restrict pfpT1;
        const int fftnr=l1para->fftnr;

	float complex **  FX;   //// <---

	pfp0     = (float complex **)l1para->pfp;
	pfp1     = (float complex **)l1para->pfpin;

	pfx      = (float complex **)l1para->pfx;

	pfpT0    = (float complex **)l1para->pfpT;
	pfpT1    = (float complex **)l1para->pfpT2;

	FX      = (float complex **)RFX;

	int i,ip,ifreq,j; 
	float * restrict obliqcorr = l1para->obliqcorr;

	float * restrict wtilt;

	for (int ip=0;ip<npx;ip++) obliqcorr[ip] = 1.0f;

	if (datamode == 0) 
	{

		FFT_R2C_PP(l1para->fftplan[0], taup_p, npx, pfp0, fftnr);
		MAT_TRANSPOSE(pfp0, pfpT0, hfid, npx);

		if (datum == SURFACE0)
			for (ifreq=0;ifreq<hfid;ifreq++)
				for (ip=0;ip<npx;ip++)
					pfpT0[ifreq][ip] = (I)*pfpT0[ifreq][ip];

	}

	if (datamode == 1)
	{
		wtilt = (float *)l1para->zwtilt;

		if (l1para->zmethod == 0||l1para->zmethod == 2)
			for (ip=0;ip<npx;ip++) obliqcorr[ip] = l1para->ztop[ip];

		FFT_R2C_PP(l1para->fftplan[0], taup_az, npx, pfp1, fftnr);
		MAT_TRANSPOSE(pfp1, pfpT1, hfid, npx);

		if (datum == SURFACE0)
			for (ifreq=0;ifreq<hfid;ifreq++)
				for (ip=0;ip<npx;ip++)
					pfpT1[ifreq][ip] = (I)*pfpT1[ifreq][ip];

		if (l1para->zmethod==0)
		{
			for (ifreq=0;ifreq<hfid;ifreq++)
			{
				for (ip=0;ip<npx;ip++)
				{
					pfpT1[ifreq][ip] = (I)*pfpT1[ifreq][ip]*wtilt[ifreq];
				}
			}
		}

	}

	if (datamode == 10) 
	{

		FFT_R2C_PP(l1para->fftplan[0], taup_p, npx, pfp0, fftnr);
		MAT_TRANSPOSE(pfp0, pfpT0, hfid, npx);

		if (datum == SURFACE0)
			for (ifreq=0;ifreq<hfid;ifreq++)
				for (ip=0;ip<npx;ip++)
					pfpT0[ifreq][ip] = (I)*pfpT0[ifreq][ip];

		FFT_R2C_PP(l1para->fftplan[0], taup_az, npx, pfp1, fftnr);
		MAT_TRANSPOSE(pfp1, pfpT1, hfid, npx);

		if (datum == SURFACE0)
			for (ifreq=0;ifreq<hfid;ifreq++)
				for (ip=0;ip<npx;ip++)
					pfpT1[ifreq][ip] = (I)*pfpT1[ifreq][ip];

	}

	if (datum == CABLE90)
	{
		taup_warp_solo_dsdg3d(nsgtrace, npx, hfid, l1para->pfpT, l1para->pfpT2, 
				l1para->TP_p, l1para->xPrimary, obliqcorr, recz, l1para->pfx, l1para);
	}
	else if (datum == SURFACE0 || datum == SURFACE90)
	{
		taup_warp_solo_dsdg3d(nsgtrace, npx, hfid, l1para->pfpT, l1para->pfpT2, 
				l1para->TS_p, l1para->xPrimary, obliqcorr, recz, l1para->pfx, l1para);
	}

	for (i=0;i<nsgtrace;i++) //mute high frequency!
		for (j=hfid;j<l1para->fftnc;j++) 
			pfx[i][j] = 0.0f;

	if (l1para->lownear>0) {
		for (i=0;i<nsgtrace;i++) {  
			pfx[i][0] = 0.0f;
			pfx[i][1] = 0.0f;

			for (j=2;j<lfid;j++) 
				pfx[i][j] *= ((float)j-1.0f)/((float)lfid-2.0f);
		}
	}

	FFT_C2R_PP(l1para->fftplan[1], pfx, nsgtrace, RFX, fftnr);

	for (i=0;i<nsgtrace;i++) 
		FX[i][l1para->fftnc-1] = 0.0f;  //// <---

}

	MSDGINTERP_WARP_DISPATCH
void c_dsdg3d_TAUP_RV_MTOS(int nsgtrace, int hfid, int lfid, int npx, float *recz, float ** taup_p, float ** taup_az,
		float **RFX,  l1inv_t *l1para, int datamode)
{
	int ip,ifreq,i,j;  

	float complex ** restrict pfp0;
	float complex ** restrict pfp1;

	float complex ** restrict pfx;

	float complex ** restrict pfpT0;
	float complex ** restrict pfpT1;
        const int fftnr=l1para->fftnr;

	float complex **  FX;  //// <---

	pfp0      = (float complex **)l1para->pfp;
	pfp1      = (float complex **)l1para->pfpin;

	pfx      = (float complex **)l1para->pfx;

	pfpT0     = (float complex **)l1para->pfpT;
	pfpT1     = (float complex **)l1para->pfpT2;

	FX      = (float complex **)RFX;

	float * restrict obliqcorr = l1para->obliqcorr;

	float * restrict wtilt;

	for (ip=0;ip<npx;ip++) obliqcorr[ip] = 1.0f;

	if (datamode == 0) 
	{

		for (ip=0;ip<npx;ip++) obliqcorr[ip] = -1.f; // zxue

		FFT_R2C_PP(l1para->fftplan[0], taup_p, npx, pfp0, fftnr);

		MAT_TRANSPOSE(pfp0, pfpT0, hfid, npx);

	}

	if (datamode == 1)
	{
		wtilt = (float *)l1para->zwtilt;
		if (l1para->zmethod==1) 
			for (ip=0;ip<npx;ip++) obliqcorr[ip] = -1.0f; // zxue    
		else 
			for (ip=0;ip<npx;ip++) obliqcorr[ip] = -l1para->ztop[ip]; // zxue

		FFT_R2C_PP(l1para->fftplan[0], taup_az, npx, pfp1, fftnr);

		MAT_TRANSPOSE(pfp1, pfpT1, hfid, npx);

		if (l1para->zmethod==0)
		{
			for (ifreq=0;ifreq<hfid;ifreq++)
			{
				for (ip=0;ip<npx;ip++)
				{
					pfpT1[ifreq][ip] = (I)*pfpT1[ifreq][ip]*wtilt[ifreq];
				}
			}
		}
	}

	if (datamode == 10) 
	{

		for (ip=0;ip<npx;ip++) obliqcorr[ip] = -1.f; // zxue

		FFT_R2C_PP(l1para->fftplan[0], taup_p, npx, pfp0, fftnr);

		MAT_TRANSPOSE(pfp0, pfpT0, hfid, npx);

		FFT_R2C_PP(l1para->fftplan[0], taup_az, npx, pfp1, fftnr);

		MAT_TRANSPOSE(pfp1, pfpT1, hfid, npx);

	}

	taup_warp_solo_dsdg3d(nsgtrace, npx, hfid, l1para->pfpT, l1para->pfpT2, 
			l1para->TG_p, l1para->xGhost, obliqcorr, recz, l1para->pfx, l1para);

	// zxue
	if(datamode == 0 || datamode == 10){
		for (ip=0;ip<nsgtrace;ip++)
			for (ifreq=0;ifreq<hfid;ifreq++)
				l1para->pfx[ip][ifreq] *= l1para->prfl2[ifreq];
	}
	if(datamode == 1){
		for (ip=0;ip<nsgtrace;ip++)
			for (ifreq=0;ifreq<hfid;ifreq++)
				l1para->pfx[ip][ifreq] *= l1para->zrfl[ifreq];
	}

	for (i=0;i<nsgtrace;i++) //mute high frequency!
		for (j=hfid;j<l1para->fftnc;j++) 
			pfx[i][j] = 0.0f;

	if (l1para->lownear>0) {
		for (i=0;i<nsgtrace;i++) {  
			pfx[i][0] = 0.0f;
			pfx[i][1] = 0.0f;

			for (j=2;j<lfid;j++) 
				pfx[i][j] *= ((float)j-1.0f)/((float)lfid-2.0f);
		}
	}

	FFT_C2R_PP(l1para->fftplan[1], pfx, nsgtrace, RFX, fftnr);

	for (i=0;i<nsgtrace;i++) 
		FX[i][l1para->fftnc-1] = 0.0f;  //// <---

}

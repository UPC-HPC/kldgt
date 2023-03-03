#include "msdginterp.h"
#include "hilbert.h"
#include "invmatrix.h"
#include <fftw3.h>
#include <PFL_C.h>
#include "utility.h"
#include <complex>
#include <fenv.h>

using namespace std;

#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )

#define FFT_R2C_PP(plan, t, m, f, nf2)                                                    \
	for(int i=0; i<m; i++)                                                                  \
PFL_RCFFT(plan, (nf2-2),(float*)(t)[i],(PFLComplex*)(f)[i],1,-1,(nf2),(nf2),1.0f/((nf2)-2));

#define FFT_C2R_PP(plan, f, m, t, nf2)                                              \
	for(int i=0; i<m; i++)                                                            \
PFL_CRFFT(plan, (nf2-2),(PFLComplex*)(f)[i],(float*)(t)[i],1,+1,(nf2),(nf2),1.0f);

extern "C" {
	int reorder(int npx, l1inv_t *l1para, int *indx, float **taup, float pscale);
}

void get_lo_and_hi_f (l1inv_t *l1para, float rpx, float &flocut, float &fhicut);
void butterworth_bandpass (int flt_len, float *flt, int lfid, int hfid, int flt_order);
void apply_dipfilter3d (int npx, int *indx, float **taup, l1inv_t *l1para);
void convolve_dataxt_with_target (int nsgtrace, int fftnr, float **dataxt_p, l1inv_t *l1para);
void merge_lowfreq_filter (int nsgtrace, l1inv_t *l1para, float **input_p, float **output, float *recz);

void bsdg3d_deghosting(l1inv_t *l1para, int nsgtrace, float *offsetx, float *offsety, float *recz,
		float **input_p, 
		float **output, float sppow, int lfid)
{
	int npx, npy, nq;
	int ip, itrc, isamp;
	float fnyq, deltaf;

	float maxamp=-1.0f;

	float resscale = l1para->resscale;

	float **taup;
	float **taup_filter;
	float **pweight,**tweight;

	float **dataxt_p, **residual1_p;
	float **dataxt_az, **residual1_az;
	float **dataxt_ay, **residual1_ay;

	int np_orig,np_reduce;
	int *indx;

	int hfid,lfid0,datamode,datum;

	l1para->ms_pyes=1;
	l1para->ms_zyes=0;
	l1para->ms_yyes=0;

	int pyes=l1para->ms_pyes;
	int zyes=l1para->ms_zyes;
	int yyes=l1para->ms_yyes;

	for (itrc=0;itrc<nsgtrace;itrc++)
		for (isamp=0;isamp<l1para->nsamp;isamp++)
			if (fabsf(input_p[itrc][isamp])>maxamp) maxamp=fabsf(input_p[itrc][isamp]);
	test_amp(maxamp);
	if (maxamp<1e-5)
	{
		for (itrc=0;itrc<nsgtrace;itrc++)
			for (isamp=0;isamp<l1para->nsamp;isamp++) output[itrc][isamp]=input_p[itrc][isamp];
		return;
	}

        if (nsgtrace < l1para->iprocessing_nthres)
        {
                for (itrc=0;itrc<nsgtrace;itrc++)
                   for (isamp=0;isamp<l1para->nsamp;isamp++) output[itrc][isamp]=input_p[itrc][isamp];
                return;
        }

	calc_matrix_A_T_p (l1para, nsgtrace, offsetx, offsety, recz, lfid);

	PFL_GET(&l1para->fftplan[0], PFL_PORTABLE_DEF); //kn
	PFL_GET(&l1para->fftplan[1], PFL_PORTABLE_DEF); //kn


	fnyq = 500000.0/l1para->srate;
	deltaf = fnyq/(float)(l1para->fftnc-1);

	lfid0 = 2.0f/deltaf; //// 2 Hz;

	////index of the maximum/minimum frequency
	hfid = min(floor((float)(l1para->fftnc)*l1para->highfreq/fnyq),l1para->fftnc);

	npx = l1para->npxsave;
	npy = l1para->npysave;
	nq  = l1para->nqsave;

	npx = npx*npy*nq; 
	np_orig = npx;

	dataxt_p    = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));
	residual1_p = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));

	mycopy(input_p,dataxt_p,nsgtrace,l1para->nsamp);
	for (itrc=0;itrc<nsgtrace;itrc++) 
		for (isamp=l1para->nsamp;isamp<l1para->fftnr;isamp++) 
			dataxt_p[itrc][isamp] = 0.0;

	l1para->nsgtracepad   = ceil((float)(nsgtrace  /4.))*4;
	l1para->npxpad = ceil((float)(npx/4.))*4;

	taup          = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));
	taup_filter   = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));

	pweight = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));
	tweight = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));

	indx = (int *)calloc(np_orig,sizeof(int));

	myinit(pweight,npx,l1para->fftnr,1.0f);
	myinit(tweight,npx,l1para->fftnr,1.0f);

	myinit(output,nsgtrace,l1para->nsamp,0.0f);

	myinit(taup,  npx, l1para->fftnr, 0.0f); 

	for (ip = 0; ip < np_orig; ip++) indx[ip] = ip;


	c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(5,l1para->niter*0.35), nsgtrace, MIN(l1para->lowstart/deltaf,hfid),
			lfid0, npx, MIN(1.0f,sppow*0.5f), pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight, l1para);  

	c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter*0.35), nsgtrace, MIN((l1para->lowstart+5.0f)/deltaf,hfid),
			lfid0, npx, MIN(1.0f,sppow*0.5f), pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight, l1para);

	np_reduce = reorder(np_orig, l1para, indx, taup, 0.5f);

	recalc_matrix_A_T_p(l1para, nsgtrace, offsetx, offsety, recz, indx, np_reduce, lfid);

	////G2_logInfo("bsdg","nptotal=%d, np_reduce=%d",npx,np_reduce);fflush(0);
	npx = np_reduce;

	c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(5,l1para->niter*0.4), nsgtrace, MIN(l1para->highfreq*0.5f/deltaf,hfid),
			lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight, l1para); 

	c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter), nsgtrace, MIN(l1para->highfreq/deltaf,hfid),    
			lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight,l1para);  

	mycopy (taup, taup_filter, npx, l1para->fftnr);
	if (l1para->ldipf > 0)
		apply_dipfilter3d (npx, indx, taup_filter, l1para);

	if (l1para->fpredatum==FPREDATUM_NO)
	{
		datamode = 0;
		c_l1inv_TAUP_RV_MTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, 0);
		ARRAY_MATH(residual1_p, input_p, +, dataxt_p, nsgtrace, l1para->nsamp);

		datamode = 0;
		datum = 0;
		c_l1inv_TAUP_RV_CTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
		ARRAY_MATH(residual1_p, residual1_p, -, dataxt_p, nsgtrace, l1para->nsamp);

		datamode = 0;
		datum = l1para->datum;
		c_l1inv_TAUP_RV_CTOS(nsgtrace, hfid, lfid, npx, taup_filter, dataxt_p, l1para, datamode, datum);
		ARRAY_MATH(output, output, +, dataxt_p, nsgtrace, l1para->nsamp);
	}
	else if (l1para->fpredatum==FPREDATUM_CTOS)
	{
		datamode = 0;
		datum = 0;
		c_l1inv_TAUP_RV_CTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
		ARRAY_MATH(residual1_p, input_p, -, dataxt_p, nsgtrace, l1para->nsamp);

		datamode = 0;
		datum = l1para->datum;
		c_l1inv_TAUP_RV_CTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
		ARRAY_MATH(output, output, +, dataxt_p, nsgtrace, l1para->nsamp);

	}
	else if (l1para->fpredatum==FPREDATUM_MTOS)
	{
		datamode = 0;
		datum = 0;
		c_l1inv_TAUP_RV_MTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
		ARRAY_MATH(residual1_p, input_p, -, dataxt_p, nsgtrace, l1para->nsamp);

		datamode = 0;
		datum = l1para->datum;
		c_l1inv_TAUP_RV_CTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
		ARRAY_MATH(output, output, +, dataxt_p, nsgtrace, l1para->nsamp);

	}
	else if (l1para->fpredatum==FPREDATUM_SIN) 
	{

		datamode = 0;
		datum = SURFACE90;
		c_l1inv_TAUP_RV_CTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
		ARRAY_MATH(residual1_p, input_p, -, dataxt_p, nsgtrace, l1para->nsamp);

		datamode = 0;
		datum = 0;
		c_l1inv_TAUP_RV_CTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
		ARRAY_MATH(output, output, +, dataxt_p, nsgtrace, l1para->nsamp);

		datamode = 0;
		datum = 0;
		c_l1inv_TAUP_RV_MTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
		ARRAY_MATH(output, output, -, dataxt_p, nsgtrace, l1para->nsamp);

	}
	else if (l1para->fpredatum==FPREDATUM_CIN) 
	{

		datamode = 0;
		datum = 0;
		c_l1inv_TAUP_RV_CTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
		ARRAY_MATH(residual1_p, input_p, -, dataxt_p, nsgtrace, l1para->nsamp);
		ARRAY_MATH(output, output, +, dataxt_p, nsgtrace, l1para->nsamp);

		datamode = 0;
		datum = 0;
		c_l1inv_TAUP_RV_MTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
		ARRAY_MATH(output, output, -, dataxt_p, nsgtrace, l1para->nsamp);

	}

	if ( !(l1para->resscale < 0.0f) ) // user choose to skip the fitting of weak energy if <0
	{
		npx = np_orig;

		myinit(pweight,npx,l1para->fftnr,1.0f);
		myinit(tweight,npx,l1para->fftnr,1.0f);

		myinit(taup,  npx, l1para->fftnr, 0.0f); 

		for (ip = 0; ip < np_orig; ip++) indx[ip] = ip;

		calc_matrix_A_T_p(l1para, nsgtrace, offsetx, offsety, recz, lfid);

		c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(3,l1para->niter*0.2*resscale), nsgtrace, MIN(l1para->lowstart/deltaf,hfid),
				lfid0, npx, MIN(0.8f,sppow*0.5f), pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para); 

		c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(3,l1para->niter*0.2*resscale), nsgtrace, MIN((l1para->lowstart+5.0f)/deltaf,hfid),
				lfid0, npx, MIN(0.8f,sppow*0.5f), pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para); 


		np_reduce = reorder(np_orig, l1para, indx, taup, 1.0f);
		npx = np_reduce;

		recalc_matrix_A_T_p(l1para, nsgtrace, offsetx, offsety, recz, indx, np_reduce, lfid);

		c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(4,l1para->niter*0.2*resscale), nsgtrace, MIN(l1para->highfreq*0.5f/deltaf,hfid),
				lfid, npx, sppow*0.5f, pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para); 

		c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter*0.4*resscale), nsgtrace, MIN(l1para->highfreq/deltaf,hfid),    
				lfid, npx, sppow*0.5f, pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para);  

		mycopy (taup, taup_filter, npx, l1para->fftnr);
		if (l1para->ldipf > 0)
			apply_dipfilter3d (npx, indx, taup_filter, l1para);

		if (l1para->fpredatum==FPREDATUM_NO)
		{
			datamode = 0;
			datum = 0;
			c_l1inv_TAUP_RV_CTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
			ARRAY_MATH(residual1_p, residual1_p, -, dataxt_p, nsgtrace, l1para->nsamp);

			datamode = 0;
			c_l1inv_TAUP_RV_MTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, 0);
			ARRAY_MATH(residual1_p, residual1_p, +, dataxt_p, nsgtrace, l1para->nsamp); 

			datamode = 0;
			datum = l1para->datum;
			c_l1inv_TAUP_RV_CTOS(nsgtrace, hfid, lfid, npx, taup_filter, dataxt_p, l1para, datamode, datum);
			ARRAY_MATH(output, output, +, dataxt_p, nsgtrace, l1para->nsamp); //=> output = output + dataxt;
		}
		else if (l1para->fpredatum==FPREDATUM_CTOS)
		{
			datamode = 0;
			datum = 0;
			c_l1inv_TAUP_RV_CTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
			ARRAY_MATH(residual1_p, residual1_p, -, dataxt_p, nsgtrace, l1para->nsamp);

			datamode = 0;
			datum = l1para->datum;
			c_l1inv_TAUP_RV_CTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
			ARRAY_MATH(output, output, +, dataxt_p, nsgtrace, l1para->nsamp);
		}
		else if (l1para->fpredatum==FPREDATUM_MTOS)
		{
			datamode = 0;
			datum = 0;
			c_l1inv_TAUP_RV_MTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
			ARRAY_MATH(residual1_p, residual1_p, -, dataxt_p, nsgtrace, l1para->nsamp);

			datamode = 0;
			datum = l1para->datum;
			c_l1inv_TAUP_RV_CTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
			ARRAY_MATH(output, output, +, dataxt_p, nsgtrace, l1para->nsamp);
		}
		else if (l1para->fpredatum==FPREDATUM_CIN||l1para->fpredatum==FPREDATUM_SIN) //// DIFFERENT from BSDG (there is possibly a bug in BSDG)
		{
			datamode = 0;
			datum = 0;
			c_l1inv_TAUP_RV_CTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
			ARRAY_MATH(output, output, +, dataxt_p, nsgtrace, l1para->nsamp);

			datamode = 0;
			datum = 0;
			c_l1inv_TAUP_RV_MTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
			ARRAY_MATH(output, output, -, dataxt_p, nsgtrace, l1para->nsamp);
		}
	} // this ends the check on resscale

	if (l1para->fpredatum==FPREDATUM_NO)
	{
		//=> output = output + 0.5*l1para->residual*residual1;
		ARRAY_MATH(output, 0.5*l1para->residual*residual1_p, +, output, nsgtrace, l1para->nsamp); 
	}
	else if (l1para->fpredatum==FPREDATUM_CTOS)
	{
		ARRAY_MATH(output, residual1_p, +, output, nsgtrace, l1para->nsamp); 
	}
	else if (l1para->fpredatum==FPREDATUM_MTOS)
	{
		ARRAY_MATH(output, residual1_p, +, output, nsgtrace, l1para->nsamp); 
	}

	if (l1para->taupQC==1) 
	{
		for (itrc=0;itrc<nsgtrace;itrc++){
			for (isamp=0;isamp<l1para->nsamp;isamp++){
				output[itrc][isamp] = residual1_p[itrc][isamp];
			}
		}
	}

	free(indx);

	flexible_free_array2d(dataxt_p);        //myfree (dataxt,nsgtrace);
	flexible_free_array2d(residual1_p);     //myfree (residual1,nsgtrace);

	flexible_free_array2d(taup);          //myfree (taup,   np_orig);
	flexible_free_array2d(taup_filter);   //myfree (taup,   np_orig);
	flexible_free_array2d(pweight);       //myfree (pweight,np_orig);
	flexible_free_array2d(tweight);       //myfree (tweight,np_orig);

	PFL_FREE(l1para->fftplan[0]); //kn
	PFL_FREE(l1para->fftplan[1]); //kn
}

void jointsr3d_deghosting(l1inv_t *l1para, int nsgtrace, float *offsetx, float *offsety, float *recz,
		float **input_p, 
		float **output, float sppow, int lfid)
{
	int npx, npy, nq;
	int ip, itrc, isamp;
	float fnyq, deltaf;

	float maxamp=-1.0f;

	float **taup;
	float **taup_filter;
	float **pweight,**tweight;

	float **dataxt_p, **residual1_p, **input_p_save;
	float **dataxt_az, **residual1_az;
	float **dataxt_ay, **residual1_ay;

	int np_orig,np_reduce;
	int *indx;

	int hfid,lfid0,datum,operation;

	l1para->ms_pyes=1;
	l1para->ms_zyes=0;
	l1para->ms_yyes=0;

	int pyes=l1para->ms_pyes;
	int zyes=l1para->ms_zyes;
	int yyes=l1para->ms_yyes;

	for (itrc=0;itrc<nsgtrace;itrc++)
		for (isamp=0;isamp<l1para->nsamp;isamp++)
			if (fabsf(input_p[itrc][isamp])>maxamp) maxamp=fabsf(input_p[itrc][isamp]);
	test_amp(maxamp);
	if (maxamp<1e-5)
	{
		for (itrc=0;itrc<nsgtrace;itrc++)
			for (isamp=0;isamp<l1para->nsamp;isamp++) output[itrc][isamp]=input_p[itrc][isamp];
		return;
	}

	if (l1para->jointmethod != JOINTDEGHOST)
		calc_matrix_A_T_p_joint (l1para, nsgtrace, offsetx, offsety, recz, lfid);
	else
		calc_matrix_A_T_p (l1para, nsgtrace, offsetx, offsety, recz, lfid);

	PFL_GET(&l1para->fftplan[0], PFL_PORTABLE_DEF); //kn
	PFL_GET(&l1para->fftplan[1], PFL_PORTABLE_DEF); //kn


	fnyq = 500000.0/l1para->srate;
	deltaf = fnyq/(float)(l1para->fftnc-1);

	lfid0 = 2.0f/deltaf; //// 2 Hz;

	////index of the maximum/minimum frequency
	hfid = min(floor((float)(l1para->fftnc)*l1para->highfreq/fnyq),l1para->fftnc);

	npx = l1para->npxsave;
	npy = l1para->npysave;
	nq = l1para->nqsave;

	npx = npx*npy*nq; 
	np_orig = npx;

	dataxt_p    = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));
	input_p_save    = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));
	residual1_p = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));

	mycopy(input_p,dataxt_p,nsgtrace,l1para->nsamp);
	mycopy(input_p,input_p_save,nsgtrace,l1para->nsamp);
	for (itrc=0;itrc<nsgtrace;itrc++) 
		for (isamp=l1para->nsamp;isamp<l1para->fftnr;isamp++) 
			dataxt_p[itrc][isamp] = 0.0;

	if (l1para->jointmethod != JOINTDEGHOST && l1para->targetconvopt == 1)
		convolve_dataxt_with_target (nsgtrace, l1para->fftnr, dataxt_p, l1para);

	mycopy(dataxt_p,input_p,nsgtrace,l1para->nsamp);

	l1para->nsgtracepad   = ceil((float)(nsgtrace  /4.))*4;
	l1para->npxpad = ceil((float)(npx/4.))*4;

	taup          = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));
	taup_filter   = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));

	pweight = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));
	tweight = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));

	indx = (int *)calloc(np_orig,sizeof(int));

	myinit(pweight,npx,l1para->fftnr,1.0f);
	myinit(tweight,npx,l1para->fftnr,1.0f);

	myinit(output,nsgtrace,l1para->nsamp,0.0f);

	myinit(taup,  npx, l1para->fftnr, 0.0f); 

	for (ip = 0; ip < np_orig; ip++) indx[ip] = ip;

	c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(5,l1para->niter*0.35), nsgtrace, MIN(l1para->lowstart/deltaf,hfid),
			lfid0, npx, MIN(1.0f,sppow*0.5f), pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight, l1para);  

	c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter*0.35), nsgtrace, MIN((l1para->lowstart+5.0f)/deltaf,hfid),
			lfid0, npx, MIN(1.0f,sppow*0.5f), pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight, l1para);

	np_reduce = reorder(np_orig, l1para, indx, taup, 0.5f);

	if (l1para->jointmethod != JOINTDEGHOST)
		recalc_matrix_A_T_p_joint(l1para, nsgtrace, offsetx, offsety, recz, indx, np_reduce, lfid);
	else
		recalc_matrix_A_T_p(l1para, nsgtrace, offsetx, offsety, recz, indx, np_reduce, lfid);

	////G2_logInfo("bsdg","nptotal=%d, np_reduce=%d",npx,np_reduce);fflush(0);
	npx = np_reduce;

	c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(5,l1para->niter*0.4), nsgtrace, MIN(l1para->highfreq*0.5f/deltaf,hfid),
			lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight, l1para); 

	c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter), nsgtrace, MIN(l1para->highfreq/deltaf,hfid),    
			lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight,l1para);  

	mycopy (taup, taup_filter, npx, l1para->fftnr);
	if (l1para->ldipf > 0)
		apply_dipfilter3d (npx, indx, taup_filter, l1para);

	operation = 0;
	c_l1inv_TAUP_RV_STOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, operation, 0);
	ARRAY_MATH(residual1_p, input_p, -, dataxt_p, nsgtrace, l1para->nsamp);

	operation = 1;
	datum = l1para->datum;
	c_l1inv_TAUP_RV_STOS(nsgtrace, hfid, lfid, npx, taup_filter, dataxt_p, l1para, operation, datum);

	if (l1para->jointmethod != JOINTDEGHOST && l1para->targetconvopt == 0)
		convolve_dataxt_with_target (nsgtrace, l1para->fftnr, dataxt_p, l1para);

	ARRAY_MATH(output, output, +, dataxt_p, nsgtrace, l1para->nsamp);


	if( !(l1para->resscale < 0.0f) ){ // user choose to skip the fitting of weak energy if < 0
		npx = np_orig;

		myinit(pweight,npx,l1para->fftnr,1.0f);
		myinit(tweight,npx,l1para->fftnr,1.0f);

		myinit(taup,  npx, l1para->fftnr, 0.0f); 

		for (ip = 0; ip < np_orig; ip++) indx[ip] = ip;

		if (l1para->jointmethod != JOINTDEGHOST)
			calc_matrix_A_T_p_joint(l1para, nsgtrace, offsetx, offsety, recz, lfid);
		else
			calc_matrix_A_T_p(l1para, nsgtrace, offsetx, offsety, recz, lfid);

		c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(3,l1para->niter*0.2*l1para->resscale), nsgtrace, MIN(l1para->lowstart/deltaf,hfid),
				lfid0, npx, MIN(0.8f,sppow*0.5f), pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para); 

		c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(3,l1para->niter*0.2*l1para->resscale), nsgtrace, MIN((l1para->lowstart+5.0f)/deltaf,hfid),
				lfid0, npx, MIN(0.8f,sppow*0.5f), pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para); 


		np_reduce = reorder(np_orig, l1para, indx, taup, 1.0f);
		npx = np_reduce;

		if (l1para->jointmethod != JOINTDEGHOST)
			recalc_matrix_A_T_p_joint(l1para, nsgtrace, offsetx, offsety, recz, indx, np_reduce, lfid);
		else
			recalc_matrix_A_T_p(l1para, nsgtrace, offsetx, offsety, recz, indx, np_reduce, lfid);


		c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(4,l1para->niter*0.2*l1para->resscale), nsgtrace, MIN(l1para->highfreq*0.5f/deltaf,hfid),
				lfid, npx, sppow*0.5f, pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para); 

		c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter*0.4*l1para->resscale), nsgtrace, MIN(l1para->highfreq/deltaf,hfid),    
				lfid, npx, sppow*0.5f, pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para);  

		mycopy (taup, taup_filter, npx, l1para->fftnr);
		if (l1para->ldipf > 0)
			apply_dipfilter3d (npx, indx, taup_filter, l1para);

		operation = 0;
		c_l1inv_TAUP_RV_STOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, operation, 0);
		ARRAY_MATH(residual1_p, residual1_p, -, dataxt_p, nsgtrace, l1para->nsamp);

		operation = 1;
		datum = l1para->datum;
		c_l1inv_TAUP_RV_STOS(nsgtrace, hfid, lfid, npx, taup_filter, dataxt_p, l1para, operation, datum);

		if (l1para->jointmethod != JOINTDEGHOST && l1para->targetconvopt == 0)
			convolve_dataxt_with_target (nsgtrace, l1para->fftnr, dataxt_p, l1para);

		ARRAY_MATH(output, output, +, dataxt_p, nsgtrace, l1para->nsamp);
	} // this ends the check on resscale

	// merge 1d filter for low frequency //
	if (l1para->lmerge1d == YESYES)
		merge_lowfreq_filter (nsgtrace, l1para, input_p_save, output, recz);

	float factor = 0.25f;
	if (l1para->jointmethod == JOINTDEGHOST || l1para->jointmethod == JOINTDESIG_SRCRCV_DEGHOST)
	{
		if(l1para->datum == SURFACE90)
			factor = -1.00f; 
		else
			factor = 1.00f;
	}
	else if (l1para->jointmethod == JOINTDESIG_SRC_DEGHOST ||
			l1para->jointmethod == JOINTDESIG_RCV_DEGHOST)    factor = 0.50f; 
	else if (l1para->jointmethod == JOINTDESIGNATURE)          factor = 1.00f; 

	if (l1para->jointmethod == JOINTDEGHOST) l1para->designorm_factor = 1.0f;
	factor = factor * l1para->designorm_factor;

	ARRAY_MATH(output, factor*l1para->residual*residual1_p, +, output, nsgtrace, l1para->nsamp); 

	if (l1para->taupQC==1) 
	{
		for (itrc=0;itrc<nsgtrace;itrc++){
			for (isamp=0;isamp<l1para->nsamp;isamp++){
				output[itrc][isamp] = l1para->designorm_factor * residual1_p[itrc][isamp];
			}
		}
	}

	free(indx);

	flexible_free_array2d(dataxt_p);        //myfree (dataxt,nsgtrace);
	flexible_free_array2d(input_p_save);
	flexible_free_array2d(residual1_p);     //myfree (residual1,nsgtrace);

	flexible_free_array2d(taup);          //myfree (taup,   np_orig);
	flexible_free_array2d(taup_filter);   //myfree (taup,   np_orig);
	flexible_free_array2d(pweight);       //myfree (pweight,np_orig);
	flexible_free_array2d(tweight);       //myfree (tweight,np_orig);

	PFL_FREE(l1para->fftplan[0]); //kn
	PFL_FREE(l1para->fftplan[1]); //kn
}


//void apply_dipfilter3d (int **hpfl, int npx, int *ip_sorted, float **taup, l1inv_t *l1para){
void apply_dipfilter3d (int npx, int *ip_sorted, float **taup, l1inv_t *l1para){
	int ip, ipx, lfid, hfid, ifreq;
	float rpx, px, py, flocut, fhicut;
	complex<float>** pfp  = (complex<float>**)(l1para->pfp);
	float *bwflt = (float *) calloc (l1para->fftnr, sizeof (float));
	//  float **pfp_r, **pfp_i;
	//  pfp_r = l1para->pfp_r;
	//  pfp_i = l1para->pfp_i;

	//  MY_RCFFTNEW(hpfl[0][0], taup, pfp_r, pfp_i, npx, l1para->fftnr);
	FFT_R2C_PP(l1para->fftplan[0], taup, npx, pfp, l1para->fftnr)

		for (ipx = 0; ipx < npx; ipx++){
			ip = ip_sorted[ipx];
			px = l1para->pxsave[ip];
			py = l1para->pysave[ip];
			rpx=sqrtf(px*px+py*py);
			////rpx = l1para->pxmin + (float) ipx * deltapx;

			get_lo_and_hi_f (l1para, rpx, flocut, fhicut);
			// apply band pass filter
			flocut = MAX (0.0f, flocut);
			fhicut = MIN (fhicut, 500000.0f / l1para->srate);
			lfid = (int) flocut * l1para->srate * l1para->fftnc / 500000.0;
			hfid = (int) fhicut * l1para->srate * l1para->fftnc / 500000.0;

			//G2_logInfo("bsdg","%s: rpx %f pxmin %f pxmax %f lfid %d hfid %d fftnc %d", __func__,
			//    rpx, l1para->pxmin, l1para->pxmax, lfid, hfid, l1para->fftnc);
			butterworth_bandpass (l1para->fftnc, bwflt, lfid, hfid, 4);
			for (ifreq = 0; ifreq < l1para->fftnc; ifreq++){
				//        pfp_r[ipx][ifreq] *= bwflt[2 * ifreq];
				//        pfp_i[ipx][ifreq] *= bwflt[2 * ifreq];
				pfp[ipx][ifreq] *= bwflt[2 * ifreq];
			}
		}
	//  MY_CRFFTNEW(hpfl[1][1], pfp_r, pfp_i, taup, npx, l1para->fftnr);
	FFT_C2R_PP(l1para->fftplan[1], pfp, npx, taup, l1para->fftnr)
		free (bwflt);
	return;
}


/////////////////////////////////////////////////////////////////////////////
// 3. butterworth band-pass filter
/////////////////////////////////////////////////////////////////////////////
void butterworth_bandpass (int flt_len, float *flt, int lfid, int hfid, int flt_order){
	// Note: flt_len is for complex float
	float w[2]; w[0] = 0.0f; w[1] = 0.0f;
	int i;

	// Vectorization of the nested loop below involves (masked) divisions by zero, and computations with the
	// resulting infinities. Even though these results are thrown away, if FPEs are enabled, the program will halt.
	// Therefore, we temporarily disable FPEs, and restore the FPE settings after the loop.

	int old_exceptions = fegetexcept();
	fedisableexcept(old_exceptions);

	flt[0] = 0.0f; flt[1] = 0.0f;
	for (i = 1; i < flt_len; i++){
		w[1] = (float) lfid / (float) i;
		flt[2 * i]     = 1.0f / sqrt (1.0f + powf (w[1], flt_order * 2.0f));
		w[1] = (float) i / (float) hfid;
		flt[2 * i]    *= 1.0f / sqrt (1.0f + powf (w[1], flt_order * 2.0f));
		flt[2 * i + 1] = 0.0f;
	}

	feenableexcept(old_exceptions);
}


void get_lo_and_hi_f (l1inv_t *l1para, float rpx, float &flocut, float &fhicut){
	int i;
	float rtap;
	if (l1para->pxlo[0] <= rpx && rpx <= l1para->pxhi[0])
	{
		flocut = 0.0f;
		fhicut = 9999.9f;
	}
	else if (rpx < l1para->pxlo[l1para->ldipf-1] || rpx > l1para->pxhi[l1para->ldipf-1])
	{
		flocut = l1para->flocut[l1para->ldipf-1];
		fhicut = l1para->fhicut[l1para->ldipf-1];
	}
	else
	{
		for (i=l1para->ldipf-2; i>=0; i--)
		{
			if (rpx < l1para->pxlo[i] || rpx > l1para->pxhi[i])
			{
				if (rpx < l1para->pxlo[i])
					rtap = (rpx - l1para->pxlo[i]) / (l1para->pxlo[i+1] - l1para->pxlo[i]);
				else
					rtap = (rpx - l1para->pxhi[i]) / (l1para->pxhi[i+1] - l1para->pxhi[i]);
				//        flocut = l1para->flocut[i] + rtap * l1para->flocut[i+1];
				//        fhicut = l1para->fhicut[i] + rtap * l1para->fhicut[i+1];
				flocut = l1para->flocut[i] + rtap * (l1para->flocut[i+1] - l1para->flocut[i]);
				fhicut = l1para->fhicut[i] + rtap * (l1para->fhicut[i+1] - l1para->fhicut[i]);
				break;
			}
		}
	}
	return;
}

void convolve_dataxt_with_target (int nsgtrace, int fftnr, float **dataxt_p, l1inv_t *l1para)
{
	if (l1para->target_type == TARGETTYPE_GAPDECON) 
	  {
	    l1para->designorm_factor = 1.0f; 
	    return;
	  }

	int itrc, isamp, ifreq, trlen = NHEAD + l1para->nsampsrc;
	float fnyq, deltaf, twopif, tdelay = l1para->tzerotgt, maxtgtamp = 0.0f;
	complex<float> **pfx = (complex<float> **) flexible_array2d (nsgtrace, l1para->fftnc, sizeof (complex<float>));
	complex<float> **pfxtarget = (complex<float> **) flexible_array2d (1, l1para->fftnc, sizeof (complex<float>));
	float **ptarget = (float **) flexible_array2d (1, l1para->fftnr, sizeof (float)); 
	complex<float> **p1dTgt = (complex<float> **) l1para->p1dTgt;

	fnyq = 500000.0f / l1para->srate;
	deltaf = fnyq / (float) (l1para->fftnc - 1);
	twopif = deltaf * 2.0f * 3.1415927f;

	if (l1para->target_type == TARGETTYPE_WAVELET)
		memcpy (ptarget[0], l1para->psig + l1para->ntrcsrc * trlen + NHEAD, l1para->nsampsrc * sizeof (float));
	else if (l1para->target_type == TARGETTYPE_SPIKE)
		ptarget[0][0] = 1.0f;

	//for (isamp = 0; isamp < l1para->nsampsrc; isamp++)
	//  maxtgtamp = MAX (maxtgtamp, fabsf (ptarget[0][isamp]));
	//for (isamp = 0; isamp < l1para->nsampsrc; isamp++)
	//  ptarget[0][isamp] /= maxtgtamp;
	if (l1para->designorm == DESIGNORM_MAX||l1para->designorm == DESIGNORM_NONE)
	{
		for (isamp = 0; isamp < l1para->nsampsrc; isamp++)
			maxtgtamp = MAX (maxtgtamp, fabsf (ptarget[0][isamp]));
	}
	else if (l1para->designorm == DESIGNORM_RMS)
	{
		for (isamp = 0; isamp < l1para->nsampsrc; isamp++)
			maxtgtamp += ptarget[0][isamp] * ptarget[0][isamp];
		maxtgtamp = sqrtf (maxtgtamp / l1para->nsampsrc);
	}

	if(l1para->designorm == DESIGNORM_NONE)
	  {
	    if (l1para->targetconvopt == 0) l1para->designorm_factor = maxtgtamp/l1para->total_gunvol;
	    else l1para->designorm_factor = 1.0f/l1para->total_gunvol;
	  }
	else l1para->designorm_factor = 1.0f;

	if (l1para->designorm != DESIGNORM_NONE)
		for (isamp = 0; isamp < l1para->nsampsrc; isamp++)
			ptarget[0][isamp] /= maxtgtamp;



	FFT_R2C_PP (l1para->fftplan[0], dataxt_p, nsgtrace, pfx, l1para->fftnr);

	FFT_R2C_PP (l1para->fftplan[0], ptarget, 1, pfxtarget, l1para->fftnr);
	for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
		pfxtarget[0][ifreq] *= (l1para->fftnr - 2);

	if (l1para->lmerge1d == YESYES)
		for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
			p1dTgt[0][ifreq] = pfxtarget[0][ifreq] * exp (complex<float> (0.0f, twopif * ifreq * tdelay));

	for (itrc = 0; itrc < nsgtrace; itrc++)
		for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
			pfx[itrc][ifreq] *= pfxtarget[0][ifreq] * exp (complex<float> (0.0f, twopif * ifreq * tdelay));

	FFT_C2R_PP (l1para->fftplan[1], pfx, nsgtrace, dataxt_p, l1para->fftnr);
	flexible_free_array2d (pfx); flexible_free_array2d (pfxtarget);
	flexible_free_array2d (ptarget);
}


void merge_lowfreq_filter (int nsgtrace, l1inv_t *l1para, float **input_p, float **output, float *recz)
{
	// input_p is original input and output is output from 3D inversion //
	int itrc, isamp, ifreq, lfid, hpfl[2] ;
	int trlensig = NHEAD + l1para->nsampsrc;
	float phase0, phase1, phase2, invwater = 1.0f / l1para->vwater;
	float fnyq, deltaf, total_gunvol;
	float **inputxt_p  = (float**) flexible_array2d (nsgtrace, l1para->fftnr, sizeof (float));
	float **outputxt_p = (float**) flexible_array2d (nsgtrace, l1para->fftnr, sizeof (float));
	complex<float> **inputxf_p    = (complex<float> **) flexible_array2d (nsgtrace, l1para->fftnr, sizeof (float));
	complex<float> **outputxf_p   = (complex<float> **) flexible_array2d (nsgtrace, l1para->fftnr, sizeof (float));
	complex<float> **pSrcFilterAz = (complex<float> **) flexible_array2d (1, l1para->fftnr, sizeof (float));
	complex<float> **pFilter      = (complex<float> **) flexible_array2d (1, l1para->fftnr, sizeof (float));
	head_t *src;

	mycopy (input_p, inputxt_p,  nsgtrace, l1para->nsamp);
	mycopy (output,  outputxt_p, nsgtrace, l1para->nsamp);
	for (itrc = 0; itrc < nsgtrace; itrc++) 
	{
		for (isamp = l1para->nsamp; isamp <l1para->fftnr; isamp++) 
		{
			inputxt_p[itrc][isamp]  = 0.0f;
			outputxt_p[itrc][isamp] = 0.0f;
		}
	}

	PFL_GET (&hpfl[0], PFL_PORTABLE_DEF);
	PFL_GET (&hpfl[1], PFL_PORTABLE_DEF);
	fnyq = 500000.0f / l1para->srate;
	deltaf = fnyq / (float) (l1para->fftnc - 1);

	float twopif = deltaf * 2.0f * 3.1415927f;

	FFT_R2C_PP (hpfl[0], inputxt_p, nsgtrace,        inputxf_p, l1para->fftnr); 
	FFT_R2C_PP (hpfl[0], outputxt_p,nsgtrace,        outputxf_p,l1para->fftnr); 

	complex<float> ** p1dSrc = (complex<float> **) l1para->p1dSrc;
	complex<float> ** p1dTgt = (complex<float> **) l1para->p1dTgt;
	float rFilter = 0.0f;
	pFilter[0][0] = pFilter[0][l1para->fftnc - 1] = 0.0f; // no change at 0 Hz
	for (ifreq = 1; ifreq < l1para->fftnc - 1; ifreq++)
	{
		pFilter[0][ifreq] = p1dTgt[0][ifreq] / p1dSrc[0][ifreq];
		rFilter = abs (pFilter[0][ifreq]);

		if (rFilter > l1para->mrg1dmaxg)
			pFilter[0][ifreq] *= l1para->mrg1dmaxg / rFilter;
	}

	// merge with output from inversion at l1para->mrg1dfreq Hz //
	lfid = (int) (l1para->mrg1dfreq / deltaf);

	if (l1para->jointmethod == JOINTDESIG_SRCRCV_DEGHOST || l1para->jointmethod == JOINTDESIG_RCV_DEGHOST)
	{
		float *l_f = (float *) calloc (l1para->fftnc, sizeof (float));
		for (ifreq = 0; ifreq < l1para->fftnc; ifreq++) l_f[ifreq] = twopif * ifreq; 

		complex<float> **p1dRcv = (complex<float> **) flexible_array2d (1, l1para->fftnc, sizeof (complex<float>));
		for (itrc = 0; itrc < nsgtrace; itrc++)
		{
			float trcvghost = recz[itrc] * invwater;
			for (ifreq = 1; ifreq < l1para->fftnc; ifreq++)
			{
				if (l1para->datum == CABLE90)
					p1dRcv[0][ifreq] = (1.0f + l1para->prfl2[ifreq] * exp (complex<float> (0.0f, - (l_f[ifreq] * 2.0f * trcvghost))));
				else
					p1dRcv[0][ifreq] = exp (complex<float> (0.0f,  l_f[ifreq] * trcvghost)) + 
						l1para->prfl2[ifreq] * exp (complex<float> (0.0f, -l_f[ifreq] * trcvghost));

				pFilter[0][ifreq] = p1dTgt[0][ifreq] / p1dSrc[0][ifreq] / p1dRcv[0][ifreq];
				rFilter = abs (pFilter[0][ifreq]);

				if (rFilter > l1para->mrg1dmaxg)
					pFilter[0][ifreq] *= l1para->mrg1dmaxg / rFilter;

				float rtap = sqrtf (1.0f / (1.0f + powf ((float) ifreq / lfid, 6.0f)));
				inputxf_p[itrc][ifreq] *= pFilter[0][ifreq];
				outputxf_p[itrc][ifreq] = (1.0f - rtap) * outputxf_p[itrc][ifreq] +
					rtap  * inputxf_p[itrc][ifreq]; 
			}
		}
		free (l_f);
		flexible_free_array2d (p1dRcv);
	}
	else
	{
		for (itrc = 0; itrc < nsgtrace; itrc++)
		{
			for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
			{
				float rtap = sqrtf (1.0f / (1.0f + powf ((float) ifreq / lfid, 6.0f)));
				inputxf_p[itrc][ifreq] *= pFilter[0][ifreq];
				outputxf_p[itrc][ifreq] = (1.0f - rtap) * outputxf_p[itrc][ifreq] +
					rtap  * inputxf_p[itrc][ifreq]; 
			}
		}
	}

	FFT_C2R_PP (hpfl[1], outputxf_p, nsgtrace, outputxt_p, l1para->fftnr);
	mycopy (outputxt_p, output,  nsgtrace, l1para->nsamp);

	PFL_FREE (hpfl[0]);
	PFL_FREE (hpfl[1]);

	flexible_free_array2d (inputxt_p);
	flexible_free_array2d (inputxf_p);
	flexible_free_array2d (outputxt_p);
	flexible_free_array2d (outputxf_p);
	flexible_free_array2d (pSrcFilterAz);
	flexible_free_array2d (pFilter);
}

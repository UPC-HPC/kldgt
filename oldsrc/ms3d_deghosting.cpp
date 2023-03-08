#include "msdginterp.h"
#include "hilbert.h"
#include "invmatrix.h"
#include <fftw3.h>
#include <PFL_C.h>
#include "utility.h"
#include <complex>

using namespace std;

#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )

extern "C" {
  int reorder(int npx, l1inv_t *l1para, int *indx, float **taup, float pscale);
}

void convolve_dataxt_with_target (int nsgtrace, int fftnr, float **dataxt_p, l1inv_t *l1para);
void apply_dipfilter3d (int npx, int *indx, float **taup, l1inv_t *l1para);

 void ms3d_deghosting(l1inv_t *l1para, int nsgtrace, float *offsetx, float *offsety, float *recz,
		       float **input_p, float **input_az, float **input_ay, 
		       float **output, float sppow, int lfid)
{
  int npx, npy, nsgtracep, nsgtrace_az, nsgtrace_ay;
  int ip, itrc, isamp;
  float fnyq, deltaf;

  float maxamp=-1.0f;
  
  float **taup,**taup1;
  float **taup_filter;
  float **pweight,**tweight;

  float **dataxt_p, **residual1_p, **residual1_ponly;
  float **dataxt_az, **residual1_az;
  float **dataxt_ay, **residual1_ay;

  float **dataxt_tmp;

  int np_orig,np_reduce;
  int *indx;

  int hfid,lfid0,datamode,datum;

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

  if(l1para->zerotraceattr == 1)
    {
      nsgtracep   = l1para->nsgtracearr4d[0];
      nsgtrace_az = l1para->nsgtracearr4d[1];
      nsgtrace_ay = l1para->nsgtracearr4d[2];
    }
  else
    {
      nsgtracep   = nsgtrace;
      nsgtrace_az = nsgtrace;
      nsgtrace_ay = nsgtrace;
    }

  calc_matrix_A_T_p (l1para, nsgtracep, offsetx, offsety, recz, lfid);

  PFL_GET(&l1para->fftplan[0], PFL_PORTABLE_DEF); //kn
  PFL_GET(&l1para->fftplan[1], PFL_PORTABLE_DEF); //kn

  
  fnyq = 500000.0/l1para->srate;
  deltaf = fnyq/(float)(l1para->fftnc-1);

  lfid0 = 2.0f/deltaf; //// 2 Hz;

  ////index of the maximum/minimum frequency
  hfid = min(floor((float)(l1para->fftnc)*l1para->highfreq/fnyq),l1para->fftnc);

  npx = l1para->npxsave;
  npy = l1para->npysave;
  
  npx = npx*npy; 
  np_orig = npx;
  
  dataxt_p    = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));
  residual1_p = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));
  residual1_ponly = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));

  dataxt_tmp    = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));

 
  mycopy(input_p,dataxt_p,nsgtrace,l1para->nsamp);
  for (itrc=0;itrc<nsgtrace;itrc++) 
    for (isamp=l1para->nsamp;isamp<l1para->fftnr;isamp++) 
      dataxt_p[itrc][isamp] = 0.0;

  if ( l1para->msdataz == YESYES )
    {	
      dataxt_az    = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));
      residual1_az = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));

      mycopy(input_az,dataxt_az,nsgtrace,l1para->nsamp);
      for (itrc=0;itrc<nsgtrace;itrc++) 
	for (isamp=l1para->nsamp;isamp<l1para->fftnr;isamp++) 
	  dataxt_az[itrc][isamp] = 0.0;
      
      datamode = 1;
      if (l1para->iexternlc == 0) ms_lowcutf(dataxt_az, nsgtrace, l1para, datamode);

      mycopy(dataxt_az,input_az,nsgtrace,l1para->nsamp);

    }

  if ( l1para->msdatay == YESYES )
    {	
      dataxt_ay    = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));
      residual1_ay = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));
      
      mycopy(input_ay,dataxt_ay,nsgtrace,l1para->nsamp);
      for (itrc=0;itrc<nsgtrace;itrc++) 
	for (isamp=l1para->nsamp;isamp<l1para->fftnr;isamp++) 
	  dataxt_ay[itrc][isamp] = 0.0;
      
      datamode = 2;
      ms_lowcutf(dataxt_ay, nsgtrace, l1para, datamode);//PS: check check
      
      mycopy(dataxt_ay,input_ay,nsgtrace,l1para->nsamp);
      
    }

  l1para->nsgtracepad   = ceil((float)(nsgtrace  /4.))*4;
  l1para->npxpad = ceil((float)(npx/4.))*4;
  
  taup   = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));
  taup1 = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));
  taup_filter   = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));

  pweight = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));
  tweight = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));

  indx = (int *)calloc(np_orig,sizeof(int));

  myinit(pweight,npx,l1para->fftnr,1.0f);
  myinit(tweight,npx,l1para->fftnr,1.0f);

  myinit(output,nsgtrace,l1para->nsamp,0.0f);

  myinit(taup,  npx, l1para->fftnr, 0.0f); 
  myinit(taup1, npx, l1para->fftnr, 0.0f); 

  for (ip = 0; ip < np_orig; ip++) indx[ip] = ip;

  if (l1para->lowguide==0)
    {

      l1para->ms_pyes=1;
      l1para->ms_zyes=0;
      l1para->ms_yyes=0;
      
      c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(5,l1para->niterlow*0.35), nsgtrace, MIN(l1para->lowstart/deltaf,hfid),
		   lfid0, npx, MIN(1.0f,sppow*0.5f), pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight, l1para);  
      
      c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niterlow*0.35), nsgtrace, MIN((l1para->lowstart+5.0f)/deltaf,hfid),
		   lfid0, npx, MIN(1.0f,sppow*0.5f), pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight, l1para);
      
      l1para->ms_pyes=pyes;
      l1para->ms_zyes=zyes;
      l1para->ms_yyes=yyes;
      
    }
  else
    {
      
      c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(5,l1para->niter*0.35), nsgtrace, MIN(l1para->lowstart/deltaf,hfid),
		   lfid0, npx, MIN(1.0f,sppow*0.5f), pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight, l1para);  
      
      c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter*0.35), nsgtrace, MIN((l1para->lowstart+5.0f)/deltaf,hfid),
		   lfid0, npx, MIN(1.0f,sppow*0.5f), pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight, l1para);
    }

  np_reduce = reorder(np_orig, l1para, indx, taup, 0.5f);

  recalc_matrix_A_T_p(l1para, nsgtracep, offsetx, offsety, recz, indx, np_reduce, lfid);

  ////G2_logInfo("bsdg","nptotal=%d, np_reduce=%d",npx,np_reduce);fflush(0);
  npx = np_reduce;

  if ((l1para->outtype==3||l1para->outtype==5)&&(l1para->zlowcut>=1.0f||l1para->lzfilter>0))
    {
      
      l1para->ms_pyes=1;
      l1para->ms_zyes=0;
      l1para->ms_yyes=0;
      
      c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(5,l1para->niterlow*0.4), nsgtrace, MIN(l1para->highfreq_p*0.5f/deltaf,hfid),
		   lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup1, tweight, l1para); 
      
      c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niterlow), nsgtrace, MIN(l1para->highfreq_p/deltaf,hfid),    
		   lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup1, tweight,l1para);  
      
      if (l1para->lowguide==0)
	{
	  
	  datamode = 0;
	  c_l1inv_TAUP_RV_MTOS(nsgtracep, hfid, lfid, npx, taup1, dataxt_tmp, l1para, datamode, 0);
	  ARRAY_MATH(residual1_ponly, input_p, +, dataxt_tmp, nsgtracep, l1para->nsamp);
	  	  
	  datamode = 0;
	  datum = 0;
	  c_l1inv_TAUP_RV_CTOS(nsgtracep, hfid, lfid, npx, taup1, dataxt_tmp, l1para, datamode, datum);
	  ARRAY_MATH(residual1_ponly, residual1_ponly, -, dataxt_tmp, nsgtracep, l1para->nsamp);
            
	}

      l1para->ms_pyes=pyes;
      l1para->ms_zyes=zyes;
      l1para->ms_yyes=yyes;

      c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(5,l1para->niter*0.4), nsgtrace, MIN(l1para->highfreq*0.5f/deltaf,hfid),
		   lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight, l1para); 
      
      c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter), nsgtrace, MIN(l1para->highfreq/deltaf,hfid),    
		   lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight,l1para);  
      
      l1inv3d_taup_merge(l1para, hfid, npx, taup, taup1);    //// output is taup
      
    }
  else
    {

      if (l1para->lowguide==0)
	{

	  l1para->ms_pyes=1;
	  l1para->ms_zyes=0;
	  l1para->ms_yyes=0;
	  
	  c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(5,l1para->niter*0.4), nsgtrace, MIN(l1para->highfreq*0.5f/deltaf,hfid),
		       lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup1, tweight, l1para); 
	  
	  c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter), nsgtrace, MIN(l1para->highfreq/deltaf,hfid),    
		       lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup1, tweight,l1para);  
	  
	  
	  datamode = 0;
	  c_l1inv_TAUP_RV_MTOS(nsgtracep, hfid, lfid, npx, taup1, dataxt_tmp, l1para, datamode, 0);
	  ARRAY_MATH(residual1_ponly, input_p, +, dataxt_tmp, nsgtracep, l1para->nsamp);
	  
	  
	  datamode = 0;
	  datum = 0;
	  c_l1inv_TAUP_RV_CTOS(nsgtracep, hfid, lfid, npx, taup1, dataxt_tmp, l1para, datamode, datum);
	  ARRAY_MATH(residual1_ponly, residual1_ponly, -, dataxt_tmp, nsgtracep, l1para->nsamp);

	  l1para->ms_pyes=pyes;
	  l1para->ms_zyes=zyes;
	  l1para->ms_yyes=yyes;
	
	}

      c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(5,l1para->niter*0.4), nsgtrace, MIN(l1para->highfreq*0.5f/deltaf,hfid),
		   lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight, l1para); 
      
      c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter), nsgtrace, MIN(l1para->highfreq/deltaf,hfid),    
		   lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight,l1para);  

    }
 
  mycopy (taup, taup_filter, npx, l1para->fftnr);
  if (l1para->ldipf > 0)
    apply_dipfilter3d (npx, indx, taup_filter, l1para);
 
  if (l1para->fpredatum==FPREDATUM_NO)
    {

      datamode = 0;
      c_l1inv_TAUP_RV_MTOS(nsgtracep, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, 0);
      ARRAY_MATH(residual1_p, input_p, +, dataxt_p, nsgtracep, l1para->nsamp);
  
      datamode = 0;
      datum = 0;
      c_l1inv_TAUP_RV_CTOS(nsgtracep, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
      ARRAY_MATH(residual1_p, residual1_p, -, dataxt_p, nsgtracep, l1para->nsamp);
      
      datamode = 0;
      datum = l1para->datum;
      c_l1inv_TAUP_RV_CTOS(nsgtracep, hfid, lfid, npx, taup_filter, dataxt_p, l1para, datamode, datum);
      ARRAY_MATH(output, output, +, dataxt_p, nsgtracep, l1para->nsamp);
    }
  else if (l1para->fpredatum==FPREDATUM_CTOS)
    {
      datamode = 0;
      datum = 0;
      c_l1inv_TAUP_RV_CTOS(nsgtracep, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
      ARRAY_MATH(residual1_p, input_p, -, dataxt_p, nsgtracep, l1para->nsamp);
      
      datamode = 0;
      datum = l1para->datum;
      c_l1inv_TAUP_RV_CTOS(nsgtracep, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
      ARRAY_MATH(output, output, +, dataxt_p, nsgtracep, l1para->nsamp);

    }

  if (l1para->ms_zyes==1)
    {
      datamode = 1;
      c_l1inv_TAUP_RV_MTOS(nsgtrace_az, hfid, lfid, npx, taup, dataxt_az, l1para, datamode, 0);
      ARRAY_MATH(residual1_az, input_az, +, dataxt_az, nsgtrace_az, l1para->nsamp);
      
      datamode = 1;
      datum = 0;
      c_l1inv_TAUP_RV_CTOS(nsgtrace_az, hfid, lfid, npx, taup, dataxt_az, l1para, datamode, datum);
      ARRAY_MATH(residual1_az, residual1_az, -, dataxt_az, nsgtrace_az, l1para->nsamp);

    }

  if (l1para->ms_yyes==1)
    {
      datamode = 2;
      c_l1inv_TAUP_RV_MTOS(nsgtrace_ay, hfid, lfid, npx, taup, dataxt_ay, l1para, datamode, 0);
      ARRAY_MATH(residual1_ay, input_ay, +, dataxt_ay, nsgtrace_ay, l1para->nsamp);
      
      datamode = 2;
      datum = 0;
      c_l1inv_TAUP_RV_CTOS(nsgtrace_ay, hfid, lfid, npx, taup, dataxt_ay, l1para, datamode, datum);
      ARRAY_MATH(residual1_ay, residual1_ay, -, dataxt_ay, nsgtrace_ay, l1para->nsamp);
      
    }
  
  if( !(l1para->resscale<0.0f) ) // user choose to skip the fitting of weak energy if <0
    {
      npx = np_orig;
      
      myinit(pweight,npx,l1para->fftnr,1.0f);
      myinit(tweight,npx,l1para->fftnr,1.0f);
      
      myinit(taup,  npx, l1para->fftnr, 0.0f); 
      myinit(taup1, npx, l1para->fftnr, 0.0f); 

      for (ip = 0; ip < np_orig; ip++) indx[ip] = ip;
      
      calc_matrix_A_T_p(l1para, nsgtracep, offsetx, offsety, recz, lfid);
      
      if (l1para->lowguide==0)
	{
	  
	  l1para->ms_pyes=1;
	  l1para->ms_zyes=0;
	  l1para->ms_yyes=0;
	  
	  c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(3,l1para->niter*0.2*l1para->resscale), nsgtrace, MIN(l1para->lowstart/deltaf,hfid),
		   lfid0, npx, MIN(0.8f,sppow*0.5f), pweight, residual1_ponly, residual1_az, residual1_ay, taup, tweight, l1para); 
	  
	  c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(3,l1para->niter*0.2*l1para->resscale), nsgtrace, MIN((l1para->lowstart+5.0f)/deltaf,hfid),
		       lfid0, npx, MIN(0.8f,sppow*0.5f), pweight, residual1_ponly, residual1_az, residual1_ay, taup, tweight, l1para);
	  
	  l1para->ms_pyes=pyes;
	  l1para->ms_zyes=zyes;
	  l1para->ms_yyes=yyes;
	  
	}
      else
	{
	  c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(3,l1para->niter*0.2*l1para->resscale), nsgtrace, MIN(l1para->lowstart/deltaf,hfid),
		       lfid0, npx, MIN(0.8f,sppow*0.5f), pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para); 
	  
	  c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(3,l1para->niter*0.2*l1para->resscale), nsgtrace, MIN((l1para->lowstart+5.0f)/deltaf,hfid),
		       lfid0, npx, MIN(0.8f,sppow*0.5f), pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para); 
	}
      
      np_reduce = reorder(np_orig, l1para, indx, taup, 1.0f);
      npx = np_reduce;
      
      recalc_matrix_A_T_p(l1para, nsgtracep, offsetx, offsety, recz, indx, np_reduce, lfid);
      
      if ((l1para->outtype==3||l1para->outtype==5)&&(l1para->zlowcut>=1.0f||l1para->lzfilter>0)) //// <-- condition possibly different in production version v203.03
	{
	  
	  l1para->ms_pyes=1;
	  l1para->ms_zyes=0;
	  l1para->ms_yyes=0;
      
	  c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(4,l1para->niter*0.2*l1para->resscale), nsgtrace, MIN(l1para->highfreq_p*0.5f/deltaf,hfid),
		       lfid, npx, sppow*0.5f, pweight, residual1_ponly, residual1_az, residual1_ay, taup1, tweight, l1para); 
      
	  c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter*0.4*l1para->resscale), nsgtrace, MIN(l1para->highfreq_p/deltaf,hfid),    
		       lfid, npx, sppow*0.5f, pweight, residual1_ponly, residual1_az, residual1_ay, taup1, tweight, l1para);  

	  l1para->ms_pyes=pyes;
	  l1para->ms_zyes=zyes;
	  l1para->ms_yyes=yyes;

	  c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(4,l1para->niter*0.2*l1para->resscale), nsgtrace, MIN(l1para->highfreq*0.5f/deltaf,hfid),
		       lfid, npx, sppow*0.5f, pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para); 
      
	  c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter*0.4*l1para->resscale), nsgtrace, MIN(l1para->highfreq/deltaf,hfid),    
		       lfid, npx, sppow*0.5f, pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para);  

	  l1inv3d_taup_merge(l1para, hfid, npx, taup, taup1);    

	}

      else
	{

	  c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(4,l1para->niter*0.2*l1para->resscale), nsgtrace, MIN(l1para->highfreq*0.5f/deltaf,hfid),
		       lfid, npx, sppow*0.5f, pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para); 
      
	  c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter*0.4*l1para->resscale), nsgtrace, MIN(l1para->highfreq/deltaf,hfid),    
		       lfid, npx, sppow*0.5f, pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para);  

	}
 
      mycopy (taup, taup_filter, npx, l1para->fftnr);
      if (l1para->ldipf > 0)
	apply_dipfilter3d (npx, indx, taup_filter, l1para);

      if (l1para->fpredatum==FPREDATUM_NO)
	{
	  datamode = 0;
	  datum = 0;
	  c_l1inv_TAUP_RV_CTOS(nsgtracep, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
	  ARRAY_MATH(residual1_p, residual1_p, -, dataxt_p, nsgtracep, l1para->nsamp);
      
	  datamode = 0;
	  c_l1inv_TAUP_RV_MTOS(nsgtracep, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, 0);
	  ARRAY_MATH(residual1_p, residual1_p, +, dataxt_p, nsgtracep, l1para->nsamp); 
      
	  datamode = 0;
	  datum = l1para->datum;
	  c_l1inv_TAUP_RV_CTOS(nsgtracep, hfid, lfid, npx, taup_filter, dataxt_p, l1para, datamode, datum);
	  ARRAY_MATH(output, output, +, dataxt_p, nsgtracep, l1para->nsamp); //=> output = output + dataxt;
	}
      else if (l1para->fpredatum==FPREDATUM_CTOS)
	{
	  datamode = 0;
	  datum = 0;
	  c_l1inv_TAUP_RV_CTOS(nsgtracep, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
	  ARRAY_MATH(residual1_p, residual1_p, -, dataxt_p, nsgtracep, l1para->nsamp);
      
	  datamode = 0;
	  datum = l1para->datum;
	  c_l1inv_TAUP_RV_CTOS(nsgtracep, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
	  ARRAY_MATH(output, output, +, dataxt_p, nsgtracep, l1para->nsamp);
	}
    } // this ends the check on resscale

  if ((l1para->ms_pyes==0&&(l1para->zmethod==0))||(l1para->ms_pyes==0&&(l1para->ymethod==0))) //// <--- ** CHECK **
    {
      ARRAY_MATH(output, 0.0*residual1_p, +, output, nsgtracep, l1para->nsamp); 
    }
  else if (l1para->fpredatum==FPREDATUM_CTOS)
    {
      ARRAY_MATH(output, residual1_p, +, output, nsgtracep, l1para->nsamp); 
    }
  else
    {
      //=> output = output + 0.5*l1para->residual*residual1;
      ARRAY_MATH(output, 0.5*l1para->residual*residual1_p, +, output, nsgtracep, l1para->nsamp); 
    }

  if (l1para->taupQC==2)
    {
      datamode = 1;
      c_l1inv_TAUP_RV_MTOS(nsgtrace_az, hfid, lfid, npx, taup, dataxt_az, l1para, datamode, 0);
      ARRAY_MATH(residual1_az, residual1_az, +, dataxt_az, nsgtrace_az, l1para->nsamp);
      
      datamode = 1;
      c_l1inv_TAUP_RV_CTOS(nsgtrace_az, hfid, lfid, npx, taup, dataxt_az, l1para, datamode, 0);
      ARRAY_MATH(residual1_az, residual1_az, -, dataxt_az, nsgtrace_az, l1para->nsamp);
    }

  if (l1para->taupQC==1) 
    {
      for (itrc=0;itrc<nsgtracep;itrc++){
	for (isamp=0;isamp<l1para->nsamp;isamp++){
	  output[itrc][isamp] = residual1_p[itrc][isamp];
	}
      }
    }
  if (l1para->taupQC==2) 
    {
      for (itrc=0;itrc<nsgtrace_az;itrc++){
	for (isamp=0;isamp<l1para->nsamp;isamp++){
	  output[itrc][isamp] = residual1_az[itrc][isamp];
	}
      }
    }

  free(indx);

  flexible_free_array2d(dataxt_p);        //myfree (dataxt,nsgtrace);
  flexible_free_array2d(residual1_p);     //myfree (residual1,nsgtrace);
  flexible_free_array2d(residual1_ponly); //myfree (residual1,nsgtrace);

  flexible_free_array2d(dataxt_tmp); 

  if ( l1para->msdataz == YESYES ) 
    {
      flexible_free_array2d(dataxt_az);       //myfree (dataxt,nsgtrace);
      flexible_free_array2d(residual1_az);    //myfree (residual1,nsgtrace);
    }

  if ( l1para->msdatay == YESYES ) 
    {
      flexible_free_array2d(dataxt_ay);       //myfree (dataxt,nsgtrace);
      flexible_free_array2d(residual1_ay);    //myfree (residual1,nsgtrace);
    }

  flexible_free_array2d(taup);          //myfree (taup,   np_orig);
  flexible_free_array2d(taup1);        //myfree (taup,   np_orig);
  flexible_free_array2d(taup_filter);   //myfree (taup,   np_orig);
  flexible_free_array2d(pweight);       //myfree (pweight,np_orig);
  flexible_free_array2d(tweight);       //myfree (tweight,np_orig);

  PFL_FREE(l1para->fftplan[0]); //kn
  PFL_FREE(l1para->fftplan[1]); //kn
}

 /// Multi-sensor deghost-reghost
void ms3d_dgrg(l1inv_t *l1para, int nsgtrace, float *offsetx, float *offsety, float *recz,
	       float **input_p, float **input_az, float **input_ay, 
	       float **output, float sppow, int lfid)
{
  int npx, npy, nsgtracep, nsgtrace_az, nsgtrace_ay;
  int ip, itrc, isamp;
  float fnyq, deltaf;

  float maxamp=-1.0f;
  
  float **taup,**taup1;
  float **pweight,**tweight;

  float **dataxt_p, **residual1_p, **residual1_ponly;
  float **dataxt_az, **residual1_az;
  float **dataxt_ay, **residual1_ay;

  float **dataxt_tmp;

  int np_orig,np_reduce;
  int *indx;

  int hfid,lfid0,datamode,datum;

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

  if(l1para->zerotraceattr == 1)
    {
      nsgtracep   = l1para->nsgtracearr4d[0];
      nsgtrace_az = l1para->nsgtracearr4d[1];
      nsgtrace_ay = l1para->nsgtracearr4d[2];
    }
  else
    {
      nsgtracep   = nsgtrace;
      nsgtrace_az = nsgtrace;
      nsgtrace_ay = nsgtrace;
    }

  calc_matrix_A_T_p (l1para, nsgtracep, offsetx, offsety, recz, lfid);

  PFL_GET(&l1para->fftplan[0], PFL_PORTABLE_DEF); //kn
  PFL_GET(&l1para->fftplan[1], PFL_PORTABLE_DEF); //kn

  
  fnyq = 500000.0/l1para->srate;
  deltaf = fnyq/(float)(l1para->fftnc-1);

  lfid0 = 2.0f/deltaf; //// 2 Hz;

  ////index of the maximum/minimum frequency
  hfid = min(floor((float)(l1para->fftnc)*l1para->highfreq/fnyq),l1para->fftnc);

  npx = l1para->npxsave;
  npy = l1para->npysave;
  
  npx = npx*npy; 
  np_orig = npx;
  
  dataxt_p    = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));
  residual1_p = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));
  residual1_ponly = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));

  dataxt_tmp    = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));

 
  mycopy(input_p,dataxt_p,nsgtrace,l1para->nsamp);
  for (itrc=0;itrc<nsgtrace;itrc++) 
    for (isamp=l1para->nsamp;isamp<l1para->fftnr;isamp++) 
      dataxt_p[itrc][isamp] = 0.0;

  if ( l1para->msdataz == YESYES )
    {	
      dataxt_az    = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));
      residual1_az = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));

      mycopy(input_az,dataxt_az,nsgtrace,l1para->nsamp);
      for (itrc=0;itrc<nsgtrace;itrc++) 
	for (isamp=l1para->nsamp;isamp<l1para->fftnr;isamp++) 
	  dataxt_az[itrc][isamp] = 0.0;
      
      datamode = 1;
      ms_lowcutf(dataxt_az, nsgtrace, l1para, datamode);

      mycopy(dataxt_az,input_az,nsgtrace,l1para->nsamp);

    }

  if ( l1para->msdatay == YESYES )
    {	
      dataxt_ay    = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));
      residual1_ay = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));
      
      mycopy(input_ay,dataxt_ay,nsgtrace,l1para->nsamp);
      for (itrc=0;itrc<nsgtrace;itrc++) 
	for (isamp=l1para->nsamp;isamp<l1para->fftnr;isamp++) 
	  dataxt_ay[itrc][isamp] = 0.0;
      
      datamode = 2;
      ms_lowcutf(dataxt_ay, nsgtrace, l1para, datamode);
      
      mycopy(dataxt_ay,input_ay,nsgtrace,l1para->nsamp);
      
    }

  l1para->nsgtracepad   = ceil((float)(nsgtrace  /4.))*4;
  l1para->npxpad = ceil((float)(npx/4.))*4;
  
  taup   = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));
  taup1 = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));
  
  pweight = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));
  tweight = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));

  indx = (int *)calloc(np_orig,sizeof(int));

  myinit(pweight,npx,l1para->fftnr,1.0f);
  myinit(tweight,npx,l1para->fftnr,1.0f);

  myinit(output,nsgtrace,l1para->nsamp,0.0f);

  myinit(taup,  npx, l1para->fftnr, 0.0f); 
  myinit(taup1, npx, l1para->fftnr, 0.0f); 

  for (ip = 0; ip < np_orig; ip++) indx[ip] = ip;

  if (l1para->lowguide==0)
    {

      l1para->ms_pyes=1;
      l1para->ms_zyes=0;
      l1para->ms_yyes=0;
      
      c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(5,l1para->niter*0.35), nsgtrace, MIN(l1para->lowstart/deltaf,hfid),
		   lfid0, npx, MIN(1.0f,sppow*0.5f), pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight, l1para);  
      
      c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter*0.35), nsgtrace, MIN((l1para->lowstart+5.0f)/deltaf,hfid),
		   lfid0, npx, MIN(1.0f,sppow*0.5f), pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight, l1para);
      
      l1para->ms_pyes=pyes;
      l1para->ms_zyes=zyes;
      l1para->ms_yyes=yyes;
      
    }
  else
    {
      
      c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(5,l1para->niter*0.35), nsgtrace, MIN(l1para->lowstart/deltaf,hfid),
		   lfid0, npx, MIN(1.0f,sppow*0.5f), pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight, l1para);  
      
      c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter*0.35), nsgtrace, MIN((l1para->lowstart+5.0f)/deltaf,hfid),
		   lfid0, npx, MIN(1.0f,sppow*0.5f), pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight, l1para);
    }

  np_reduce = reorder(np_orig, l1para, indx, taup, 0.5f);

  recalc_matrix_A_T_p(l1para, nsgtracep, offsetx, offsety, recz, indx, np_reduce, lfid);

  ////G2_logInfo("bsdg","nptotal=%d, np_reduce=%d",npx,np_reduce);fflush(0);
  npx = np_reduce;

  if ((l1para->outtype==3||l1para->outtype==5)&&l1para->zlowcut>=1.0f)
    {
      
      l1para->ms_pyes=1;
      l1para->ms_zyes=0;
      l1para->ms_yyes=0;
      
      c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(5,l1para->niter*0.4), nsgtrace, MIN(l1para->highfreq*0.5f/deltaf,hfid),
		   lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup1, tweight, l1para); 
      
      c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter), nsgtrace, MIN(l1para->highfreq/deltaf,hfid),    
		   lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup1, tweight,l1para);  
      
      if (l1para->lowguide==0)
	{
	  
	  datamode = 0;
	  c_l1inv_TAUP_RV_MTOS(nsgtracep, hfid, lfid, npx, taup1, dataxt_tmp, l1para, datamode, 0);
	  ARRAY_MATH(residual1_ponly, input_p, +, dataxt_tmp, nsgtracep, l1para->nsamp);
	  	  
	  datamode = 0;
	  datum = 0;
	  c_l1inv_TAUP_RV_CTOS(nsgtracep, hfid, lfid, npx, taup1, dataxt_tmp, l1para, datamode, datum);
	  ARRAY_MATH(residual1_ponly, residual1_ponly, -, dataxt_tmp, nsgtracep, l1para->nsamp);
            
	}

      l1para->ms_pyes=pyes;
      l1para->ms_zyes=zyes;
      l1para->ms_yyes=yyes;

      c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(5,l1para->niter*0.4), nsgtrace, MIN(l1para->highfreq*0.5f/deltaf,hfid),
		   lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight, l1para); 
      
      c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter), nsgtrace, MIN(l1para->highfreq/deltaf,hfid),    
		   lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight,l1para);  
      
      l1inv3d_taup_merge(l1para, hfid, npx, taup, taup1);    //// output is taup
      
    }
  else
    {

      if (l1para->lowguide==0)
	{

	  l1para->ms_pyes=1;
	  l1para->ms_zyes=0;
	  l1para->ms_yyes=0;
	  
	  c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(5,l1para->niter*0.4), nsgtrace, MIN(l1para->highfreq*0.5f/deltaf,hfid),
		       lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup1, tweight, l1para); 
	  
	  c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter), nsgtrace, MIN(l1para->highfreq/deltaf,hfid),    
		       lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup1, tweight,l1para);  
	  
	  
	  datamode = 0;
	  c_l1inv_TAUP_RV_MTOS(nsgtracep, hfid, lfid, npx, taup1, dataxt_tmp, l1para, datamode, 0);
	  ARRAY_MATH(residual1_ponly, input_p, +, dataxt_tmp, nsgtracep, l1para->nsamp);
	  
	  
	  datamode = 0;
	  datum = 0;
	  c_l1inv_TAUP_RV_CTOS(nsgtracep, hfid, lfid, npx, taup1, dataxt_tmp, l1para, datamode, datum);
	  ARRAY_MATH(residual1_ponly, residual1_ponly, -, dataxt_tmp, nsgtracep, l1para->nsamp);

	  l1para->ms_pyes=pyes;
	  l1para->ms_zyes=zyes;
	  l1para->ms_yyes=yyes;
	
	}

      c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(5,l1para->niter*0.4), nsgtrace, MIN(l1para->highfreq*0.5f/deltaf,hfid),
		   lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight, l1para); 
      
      c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter), nsgtrace, MIN(l1para->highfreq/deltaf,hfid),    
		   lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight,l1para);  

    }
  

  datamode = 0;
  c_l1inv_TAUP_RV_MTOS(nsgtracep, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, 0);
  ARRAY_MATH(residual1_p, input_p, +, dataxt_p, nsgtracep, l1para->nsamp);
  
  datamode = 0;
  datum = 0;
  c_l1inv_TAUP_RV_CTOS(nsgtracep, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
  ARRAY_MATH(residual1_p, residual1_p, -, dataxt_p, nsgtracep, l1para->nsamp);

  /////////////////////////////////////////////////////////////////////////////////////////////  
  /////////////////////////////////////////////////////////////////////////////////////////////  

  datamode = 0;
  datum = DGRGCABLE90;
  c_l1inv_TAUP_RV_MTOS(nsgtracep, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
  ARRAY_MATH(output, output, -, dataxt_p, nsgtracep, l1para->nsamp);
  
  datamode = 0;
  datum = DGRGCABLE90;
  c_l1inv_TAUP_RV_CTOS(nsgtracep, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
  ARRAY_MATH(output, output, +, dataxt_p, nsgtracep, l1para->nsamp);
  
  /////////////////////////////////////////////////////////////////////////////////////////////  
  /////////////////////////////////////////////////////////////////////////////////////////////  

  if (l1para->ms_zyes==1)
    {
      datamode = 1;
      c_l1inv_TAUP_RV_MTOS(nsgtrace_az, hfid, lfid, npx, taup, dataxt_az, l1para, datamode, 0);
      ARRAY_MATH(residual1_az, input_az, +, dataxt_az, nsgtrace_az, l1para->nsamp);
      
      datamode = 1;
      datum = 0;
      c_l1inv_TAUP_RV_CTOS(nsgtrace_az, hfid, lfid, npx, taup, dataxt_az, l1para, datamode, datum);
      ARRAY_MATH(residual1_az, residual1_az, -, dataxt_az, nsgtrace_az, l1para->nsamp);

    }

  if (l1para->ms_yyes==1)
    {
      datamode = 2;
      c_l1inv_TAUP_RV_MTOS(nsgtrace_ay, hfid, lfid, npx, taup, dataxt_ay, l1para, datamode, 0);
      ARRAY_MATH(residual1_ay, input_ay, +, dataxt_ay, nsgtrace_ay, l1para->nsamp);
      
      datamode = 2;
      datum = 0;
      c_l1inv_TAUP_RV_CTOS(nsgtrace_ay, hfid, lfid, npx, taup, dataxt_ay, l1para, datamode, datum);
      ARRAY_MATH(residual1_ay, residual1_ay, -, dataxt_ay, nsgtrace_ay, l1para->nsamp);

    }

  npx = np_orig;

  myinit(pweight,npx,l1para->fftnr,1.0f);
  myinit(tweight,npx,l1para->fftnr,1.0f);

  myinit(taup,  npx, l1para->fftnr, 0.0f); 
  myinit(taup1, npx, l1para->fftnr, 0.0f); 

  for (ip = 0; ip < np_orig; ip++) indx[ip] = ip;

  calc_matrix_A_T_p(l1para, nsgtracep, offsetx, offsety, recz, lfid);

  if (l1para->lowguide==0)
    {

      l1para->ms_pyes=1;
      l1para->ms_zyes=0;
      l1para->ms_yyes=0;
      
      c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(3,l1para->niter*0.2), nsgtrace, MIN(l1para->lowstart/deltaf,hfid),
		   lfid0, npx, MIN(0.8f,sppow*0.5f), pweight, residual1_ponly, residual1_az, residual1_ay, taup, tweight, l1para); 
      
      c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(3,l1para->niter*0.2), nsgtrace, MIN((l1para->lowstart+5.0f)/deltaf,hfid),
		   lfid0, npx, MIN(0.8f,sppow*0.5f), pweight, residual1_ponly, residual1_az, residual1_ay, taup, tweight, l1para);
      
      l1para->ms_pyes=pyes;
      l1para->ms_zyes=zyes;
      l1para->ms_yyes=yyes;
      
    }
  else
    {
      c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(3,l1para->niter*0.2), nsgtrace, MIN(l1para->lowstart/deltaf,hfid),
		   lfid0, npx, MIN(0.8f,sppow*0.5f), pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para); 
      
      c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(3,l1para->niter*0.2), nsgtrace, MIN((l1para->lowstart+5.0f)/deltaf,hfid),
		   lfid0, npx, MIN(0.8f,sppow*0.5f), pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para); 
    }

  np_reduce = reorder(np_orig, l1para, indx, taup, 1.0f);
  npx = np_reduce;

  recalc_matrix_A_T_p(l1para, nsgtracep, offsetx, offsety, recz, indx, np_reduce, lfid);

  if ((l1para->outtype==3||l1para->outtype==5)&&l1para->zlowcut>=1.0f) //// <-- condition possibly different in production version v203.03
    {

      l1para->ms_pyes=1;
      l1para->ms_zyes=0;
      l1para->ms_yyes=0;
      
      c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(4,l1para->niter*0.2), nsgtrace, MIN(l1para->highfreq*0.5f/deltaf,hfid),
		   lfid, npx, sppow*0.5f, pweight, residual1_ponly, residual1_az, residual1_ay, taup1, tweight, l1para); 
      
      c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter*0.4), nsgtrace, MIN(l1para->highfreq/deltaf,hfid),    
		   lfid, npx, sppow*0.5f, pweight, residual1_ponly, residual1_az, residual1_ay, taup1, tweight, l1para);  

      l1para->ms_pyes=pyes;
      l1para->ms_zyes=zyes;
      l1para->ms_yyes=yyes;

      c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(4,l1para->niter*0.2), nsgtrace, MIN(l1para->highfreq*0.5f/deltaf,hfid),
		   lfid, npx, sppow*0.5f, pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para); 
      
      c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter*0.4), nsgtrace, MIN(l1para->highfreq/deltaf,hfid),    
		   lfid, npx, sppow*0.5f, pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para);  

      l1inv3d_taup_merge(l1para, hfid, npx, taup, taup1);    

    }

  else
    {

      c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(4,l1para->niter*0.2), nsgtrace, MIN(l1para->highfreq*0.5f/deltaf,hfid),
		   lfid, npx, sppow*0.5f, pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para); 
      
      c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter*0.4), nsgtrace, MIN(l1para->highfreq/deltaf,hfid),    
		   lfid, npx, sppow*0.5f, pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para);  

    }

  datamode = 0;
  datum = DGRGCABLE90;
  c_l1inv_TAUP_RV_MTOS(nsgtracep, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
  ARRAY_MATH(output, output, -, dataxt_p, nsgtracep, l1para->nsamp);
  
  datamode = 0;
  datum = DGRGCABLE90;
  c_l1inv_TAUP_RV_CTOS(nsgtracep, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
  ARRAY_MATH(output, output, +, dataxt_p, nsgtracep, l1para->nsamp);
  
  /////////////////////////////////////////////////////////////////////////////////////////////  
  /////////////////////////////////////////////////////////////////////////////////////////////  
  
  free(indx);

  flexible_free_array2d(dataxt_p);        //myfree (dataxt,nsgtrace);
  flexible_free_array2d(residual1_p);     //myfree (residual1,nsgtrace);
  flexible_free_array2d(residual1_ponly); //myfree (residual1,nsgtrace);

  flexible_free_array2d(dataxt_tmp); 

  if ( l1para->msdataz == YESYES ) 
    {
      flexible_free_array2d(dataxt_az);       //myfree (dataxt,nsgtrace);
      flexible_free_array2d(residual1_az);    //myfree (residual1,nsgtrace);
    }

  if ( l1para->msdatay == YESYES ) 
    {
      flexible_free_array2d(dataxt_ay);       //myfree (dataxt,nsgtrace);
      flexible_free_array2d(residual1_ay);    //myfree (residual1,nsgtrace);
    }

  flexible_free_array2d(taup);          //myfree (taup,   np_orig);
  flexible_free_array2d(taup1);        //myfree (taup,   np_orig);
  flexible_free_array2d(pweight);       //myfree (pweight,np_orig);
  flexible_free_array2d(tweight);       //myfree (tweight,np_orig);

  PFL_FREE(l1para->fftplan[0]); //kn
  PFL_FREE(l1para->fftplan[1]); //kn
}


//// Az/Vz to V : (outtype = 6)
//// Ay/Vy to V : (outtype = 7)
void ms3d_obliqinteg(l1inv_t *l1para, int nsgtrace, float *offsetx, float *offsety, float *recz,
		      float **input_p, float **input_az, float **input_ay, 
		      float **output, float sppow, int lfid)
{
  int npx, npy, nsgtracep, nsgtrace_az, nsgtrace_ay;
  int ip, itrc, isamp;
  float fnyq, deltaf;

  float maxamp=-1.0f;
  
  float **taup,**taup1;
  float **pweight,**tweight;

  float **dataxt_p, **residual1_p, **residual1_ponly;
  float **dataxt_az, **residual1_az;
  float **dataxt_ay, **residual1_ay;

  float **dataxt_tmp;

  int np_orig,np_reduce;
  int *indx;

  int hfid,lfid0,datamode,datum;

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
      for (isamp=0;isamp<l1para->nsamp;isamp++) output[itrc][isamp]=input_az[itrc][isamp];
    return;
   }

  if(l1para->zerotraceattr == 1)
    {
      nsgtracep   = l1para->nsgtracearr4d[0];
      nsgtrace_az = l1para->nsgtracearr4d[1];
      nsgtrace_ay = l1para->nsgtracearr4d[2];
    }
  else
    {
      nsgtracep   = nsgtrace;
      nsgtrace_az = nsgtrace;
      nsgtrace_ay = nsgtrace;
    }

  calc_matrix_A_T_p (l1para, nsgtracep, offsetx, offsety, recz, lfid);

  PFL_GET(&l1para->fftplan[0], PFL_PORTABLE_DEF); //kn
  PFL_GET(&l1para->fftplan[1], PFL_PORTABLE_DEF); //kn

  
  fnyq = 500000.0/l1para->srate;
  deltaf = fnyq/(float)(l1para->fftnc-1);

  lfid0 = 2.0f/deltaf; //// 2 Hz;

  ////index of the maximum/minimum frequency
  hfid = min(floor((float)(l1para->fftnc)*l1para->highfreq/fnyq),l1para->fftnc);

  npx = l1para->npxsave;
  npy = l1para->npysave;
  
  npx = npx*npy; 
  np_orig = npx;
 
  dataxt_p    = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));
  residual1_p = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));
  residual1_ponly = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));

  dataxt_tmp    = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));

  mycopy(input_p,dataxt_p,nsgtrace,l1para->nsamp);
  for (itrc=0;itrc<nsgtrace;itrc++) 
    for (isamp=l1para->nsamp;isamp<l1para->fftnr;isamp++) 
      dataxt_p[itrc][isamp] = 0.0;

  if ( l1para->msdataz == YESYES )
    {	
      dataxt_az    = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));
      residual1_az = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));

      mycopy(input_az,dataxt_az,nsgtrace,l1para->nsamp);
      for (itrc=0;itrc<nsgtrace;itrc++) 
	for (isamp=l1para->nsamp;isamp<l1para->fftnr;isamp++) 
	  dataxt_az[itrc][isamp] = 0.0;
      
      datamode = 1;
      ms_lowcutf(dataxt_az, nsgtrace, l1para, datamode);

      mycopy(dataxt_az,input_az,nsgtrace,l1para->nsamp);

    }

  if ( l1para->msdatay == YESYES )
    {	
      dataxt_ay    = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));
      residual1_ay = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));
      
      mycopy(input_ay,dataxt_ay,nsgtrace,l1para->nsamp);
      for (itrc=0;itrc<nsgtrace;itrc++) 
	for (isamp=l1para->nsamp;isamp<l1para->fftnr;isamp++) 
	  dataxt_ay[itrc][isamp] = 0.0;
      
      datamode = 2;
      ms_lowcutf(dataxt_ay, nsgtrace, l1para, datamode);
      
      mycopy(dataxt_ay,input_ay,nsgtrace,l1para->nsamp);
      
    }

  l1para->nsgtracepad   = ceil((float)(nsgtrace  /4.))*4;
  l1para->npxpad = ceil((float)(npx/4.))*4;
  
  taup   = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));
  taup1 = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));
  
  pweight = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));
  tweight = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));

  indx = (int *)calloc(np_orig,sizeof(int));

  myinit(pweight,npx,l1para->fftnr,1.0f);
  myinit(tweight,npx,l1para->fftnr,1.0f);

  myinit(output,nsgtrace,l1para->nsamp,0.0f);

  myinit(taup,  npx, l1para->fftnr, 0.0f); 
  myinit(taup1, npx, l1para->fftnr, 0.0f); 


  for (ip = 0; ip < np_orig; ip++) indx[ip] = ip;

  if (l1para->lowguide==0)
    {

      l1para->ms_pyes=1;  //// here recz = 1
      l1para->ms_zyes=0;
      l1para->ms_yyes=0;
      
      c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(5,l1para->niter*0.35), nsgtrace, MIN(l1para->lowstart/deltaf,hfid),
		   lfid0, npx, MIN(1.0f,sppow*0.5f), pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight, l1para);  
      
      c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter*0.35), nsgtrace, MIN((l1para->lowstart+5.0f)/deltaf,hfid),
		   lfid0, npx, MIN(1.0f,sppow*0.5f), pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight, l1para);
      
      l1para->ms_pyes=pyes;
      l1para->ms_zyes=zyes;
      l1para->ms_yyes=yyes;
      
    }
  else
    {
      
      c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(5,l1para->niter*0.35), nsgtrace, MIN(l1para->lowstart/deltaf,hfid),
		   lfid0, npx, MIN(1.0f,sppow*0.5f), pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight, l1para);  
      
      c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter*0.35), nsgtrace, MIN((l1para->lowstart+5.0f)/deltaf,hfid),
		   lfid0, npx, MIN(1.0f,sppow*0.5f), pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight, l1para);
    }

  np_reduce = reorder(np_orig, l1para, indx, taup, 0.5f);

  recalc_matrix_A_T_p(l1para, nsgtracep, offsetx, offsety, recz, indx, np_reduce, lfid);

  ////G2_logInfo("bsdg","nptotal=%d, np_reduce=%d",npx,np_reduce);fflush(0);
  npx = np_reduce;


  if (l1para->lowguide==0) //// calculate residual (p-only)
    {
      
      l1para->ms_pyes=1;
      l1para->ms_zyes=0;
      l1para->ms_yyes=0;
      
      c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(5,l1para->niter*0.4), nsgtrace, MIN(l1para->highfreq*0.5f/deltaf,hfid),
		       lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup1, tweight, l1para); 
      
      c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter), nsgtrace, MIN(l1para->highfreq/deltaf,hfid),    
		   lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup1, tweight,l1para);  
      
      
      datamode = 0;
      c_l1inv_TAUP_RV_MTOS(nsgtracep, hfid, lfid, npx, taup1, dataxt_tmp, l1para, datamode, 0);
      ARRAY_MATH(residual1_ponly, input_p, +, dataxt_tmp, nsgtracep, l1para->nsamp);
      
      
      datamode = 0;
      datum = 0;
      c_l1inv_TAUP_RV_CTOS(nsgtracep, hfid, lfid, npx, taup1, dataxt_tmp, l1para, datamode, datum);
      ARRAY_MATH(residual1_ponly, residual1_ponly, -, dataxt_tmp, nsgtracep, l1para->nsamp);
      
      l1para->ms_pyes=pyes;
      l1para->ms_zyes=zyes;
      l1para->ms_yyes=yyes;
      
    }
  
  //// only z-data (or y-data) and recz = 1 and zrfl = prfl:
  c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(5,l1para->niter*0.4), nsgtrace, MIN(l1para->highfreq*0.5f/deltaf,hfid),
	       lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight, l1para); 
  
  c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter), nsgtrace, MIN(l1para->highfreq/deltaf,hfid),    
	       lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight,l1para);  
  
  datamode = 0;
  c_l1inv_TAUP_RV_MTOS(nsgtracep, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, 0);
  ARRAY_MATH(residual1_p, input_p, +, dataxt_p, nsgtracep, l1para->nsamp);
  ARRAY_MATH(output, output, -, dataxt_p, nsgtracep, l1para->nsamp); //// <----
  
  datamode = 0;
  datum = 0;
  c_l1inv_TAUP_RV_CTOS(nsgtracep, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
  ARRAY_MATH(residual1_p, residual1_p, -, dataxt_p, nsgtracep, l1para->nsamp);
  ARRAY_MATH(output, output, +, dataxt_p, nsgtracep, l1para->nsamp); //// <----

  if (l1para->ms_zyes==1)
    {
      datamode = 1;
      c_l1inv_TAUP_RV_MTOS(nsgtrace_az, hfid, lfid, npx, taup, dataxt_az, l1para, datamode, 0);
      ARRAY_MATH(residual1_az, input_az, +, dataxt_az, nsgtrace_az, l1para->nsamp);
      
      datamode = 1;
      datum = 0;
      c_l1inv_TAUP_RV_CTOS(nsgtrace_az, hfid, lfid, npx, taup, dataxt_az, l1para, datamode, datum);
      ARRAY_MATH(residual1_az, residual1_az, -, dataxt_az, nsgtrace_az, l1para->nsamp);

    }

  if (l1para->ms_yyes==1)
    {
      datamode = 2;
      c_l1inv_TAUP_RV_MTOS(nsgtrace_ay, hfid, lfid, npx, taup, dataxt_ay, l1para, datamode, 0);
      ARRAY_MATH(residual1_ay, input_ay, +, dataxt_ay, nsgtrace_ay, l1para->nsamp);
      
      datamode = 2;
      datum = 0;
      c_l1inv_TAUP_RV_CTOS(nsgtrace_ay, hfid, lfid, npx, taup, dataxt_ay, l1para, datamode, datum);
      ARRAY_MATH(residual1_ay, residual1_ay, -, dataxt_ay, nsgtrace_ay, l1para->nsamp);

    }

  if( !(l1para->resscale<0.0f) ) // user choose to skip the fitting of weak energy if <0
    { 
      npx = np_orig;
  
      myinit(pweight,npx,l1para->fftnr,1.0f);
      myinit(tweight,npx,l1para->fftnr,1.0f);

      myinit(taup,  npx, l1para->fftnr, 0.0f); 
      myinit(taup1, npx, l1para->fftnr, 0.0f); 

      for (ip = 0; ip < np_orig; ip++) indx[ip] = ip;

      calc_matrix_A_T_p(l1para, nsgtracep, offsetx, offsety, recz, lfid);

      if (l1para->lowguide==0)
	{

	  l1para->ms_pyes=1;
	  l1para->ms_zyes=0;
	  l1para->ms_yyes=0;
      
	  c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(3,l1para->niter*0.2*l1para->resscale), nsgtrace, MIN(l1para->lowstart/deltaf,hfid),
		       lfid0, npx, MIN(0.8f,sppow*0.5f), pweight, residual1_ponly, residual1_az, residual1_ay, taup, tweight, l1para); 
      
	  c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(3,l1para->niter*0.2*l1para->resscale), nsgtrace, MIN((l1para->lowstart+5.0f)/deltaf,hfid),
		       lfid0, npx, MIN(0.8f,sppow*0.5f), pweight, residual1_ponly, residual1_az, residual1_ay, taup, tweight, l1para);
      
	  l1para->ms_pyes=pyes;
	  l1para->ms_zyes=zyes;
	  l1para->ms_yyes=yyes;
      
	}
      else
	{
	  c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(3,l1para->niter*0.2*l1para->resscale), nsgtrace, MIN(l1para->lowstart/deltaf,hfid),
		       lfid0, npx, MIN(0.8f,sppow*0.5f), pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para); 
      
	  c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(3,l1para->niter*0.2*l1para->resscale), nsgtrace, MIN((l1para->lowstart+5.0f)/deltaf,hfid),
		       lfid0, npx, MIN(0.8f,sppow*0.5f), pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para); 
	}

      np_reduce = reorder(np_orig, l1para, indx, taup, 1.0f);
      npx = np_reduce;

      recalc_matrix_A_T_p(l1para, nsgtracep, offsetx, offsety, recz, indx, np_reduce, lfid);


      c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(4,l1para->niter*0.2*l1para->resscale), nsgtrace, MIN(l1para->highfreq*0.5f/deltaf,hfid),
		   lfid, npx, sppow*0.5f, pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para); 
  
      c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter*0.4*l1para->resscale), nsgtrace, MIN(l1para->highfreq/deltaf,hfid),    
		   lfid, npx, sppow*0.5f, pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para);  

      datamode = 0;
      datum = 0;
      c_l1inv_TAUP_RV_CTOS(nsgtracep, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
      ARRAY_MATH(residual1_p, residual1_p, -, dataxt_p, nsgtracep, l1para->nsamp);  //// Not used at this moment
      ARRAY_MATH(output, output, +, dataxt_p, nsgtracep, l1para->nsamp); //// <----

      datamode = 0;
      c_l1inv_TAUP_RV_MTOS(nsgtracep, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, 0);
      ARRAY_MATH(residual1_p, residual1_p, +, dataxt_p, nsgtracep, l1para->nsamp); //// Not used at this moment
      ARRAY_MATH(output, output, -, dataxt_p, nsgtracep, l1para->nsamp); //// <---- 
    } // this ends the check on resscale

  if (l1para->outtype==6 && l1para->zmethod==2)
    {
      datamode = 1;
      c_l1inv_TAUP_RV_MTOS(nsgtrace_az, hfid, lfid, npx, taup, dataxt_az, l1para, datamode, 0);
      ARRAY_MATH(residual1_az, residual1_az, +, dataxt_az, nsgtrace_az, l1para->nsamp);
      
      datamode = 1;
      datum = 0;
      c_l1inv_TAUP_RV_CTOS(nsgtrace_az, hfid, lfid, npx, taup, dataxt_az, l1para, datamode, datum);
      ARRAY_MATH(residual1_az, residual1_az, -, dataxt_az, nsgtrace_az, l1para->nsamp);

      ARRAY_MATH(output, l1para->residual*residual1_az, +, output, nsgtrace_az, l1para->nsamp); 
    }
  if (l1para->outtype==7 && l1para->ymethod==2)
    {
      datamode = 2;
      c_l1inv_TAUP_RV_MTOS(nsgtrace_ay, hfid, lfid, npx, taup, dataxt_ay, l1para, datamode, 0);
      ARRAY_MATH(residual1_ay, residual1_ay, +, dataxt_ay, nsgtrace_ay, l1para->nsamp);
      
      datamode = 2;
      datum = 0;
      c_l1inv_TAUP_RV_CTOS(nsgtrace_ay, hfid, lfid, npx, taup, dataxt_ay, l1para, datamode, datum);
      ARRAY_MATH(residual1_ay, residual1_ay, -, dataxt_ay, nsgtrace_ay, l1para->nsamp);

      ARRAY_MATH(output, l1para->residual*residual1_ay, +, output, nsgtrace_ay, l1para->nsamp); 
    }
  
  free(indx);

  flexible_free_array2d(dataxt_p);        //myfree (dataxt,nsgtrace);
  flexible_free_array2d(residual1_p);     //myfree (residual1,nsgtrace);
  flexible_free_array2d(residual1_ponly); //myfree (residual1,nsgtrace);

  flexible_free_array2d(dataxt_tmp);  

  if ( l1para->msdataz == YESYES ) 
    {
      flexible_free_array2d(dataxt_az);       //myfree (dataxt,nsgtrace);
      flexible_free_array2d(residual1_az);    //myfree (residual1,nsgtrace);
    }

  if ( l1para->msdatay == YESYES ) 
    {
      flexible_free_array2d(dataxt_ay);       //myfree (dataxt,nsgtrace);
      flexible_free_array2d(residual1_ay);    //myfree (residual1,nsgtrace);
    }

  flexible_free_array2d(taup);          //myfree (taup,   np_orig);
  flexible_free_array2d(taup1);        //myfree (taup,   np_orig);
  flexible_free_array2d(pweight);       //myfree (pweight,np_orig);
  flexible_free_array2d(tweight);       //myfree (tweight,np_orig);

  PFL_FREE(l1para->fftplan[0]); //kn
  PFL_FREE(l1para->fftplan[1]); //kn
}

void ms3d_equivalfromp(l1inv_t *l1para, int nsgtrace, float *offsetx, float *offsety, float *recz,
		      float **input_p, float **input_az, float **input_ay, 
		      float **output, float sppow, int lfid)
{

  if (l1para->pmethod == 1 || l1para->pmethod == 3 || l1para->pmethod == 5)
    ms3d_equivalzfromp(l1para, nsgtrace, offsetx, offsety, recz, 
		       input_p, input_az, input_ay, output, 
		       sppow, lfid);
  else
    ms3d_equivalyfromp(l1para, nsgtrace, offsetx, offsety, recz, 
		      input_p, input_az, input_ay, output, 
		      sppow, lfid);

}

//// equivalent V,Az,Azup, etc.
void ms3d_equivalzfromp(l1inv_t *l1para, int nsgtrace, float *offsetx, float *offsety, float *recz,
		      float **input_p, float **input_az, float **input_ay, 
		      float **output, float sppow, int lfid)
{
  int npx, npy;
  int ip, itrc, isamp;
  float fnyq, deltaf;

  float maxamp=-1.0f;
  
  float **taup,**taup1;
  float **pweight,**tweight;

  float **dataxt_p, **residual1_p;
  float **dataxt_az, **residual1_az;
  float **dataxt_ay, **residual1_ay;

  int np_orig,np_reduce;
  int *indx;

  int hfid,lfid0,datamode,datum;

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
  
  npx = npx*npy; 
  np_orig = npx;
  
  dataxt_p    = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));
  residual1_p = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));
    
  mycopy(input_p,dataxt_p,nsgtrace,l1para->nsamp);
  for (itrc=0;itrc<nsgtrace;itrc++) 
    for (isamp=l1para->nsamp;isamp<l1para->fftnr;isamp++) 
      dataxt_p[itrc][isamp] = 0.0;

  l1para->nsgtracepad   = ceil((float)(nsgtrace  /4.))*4;
  l1para->npxpad = ceil((float)(npx/4.))*4;
  
  taup   = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));
  taup1 = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));
  
  pweight = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));
  tweight = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));

  indx = (int *)calloc(np_orig,sizeof(int));

  myinit(pweight,npx,l1para->fftnr,1.0f);
  myinit(tweight,npx,l1para->fftnr,1.0f);

  myinit(output,nsgtrace,l1para->nsamp,0.0f);

  myinit(taup,  npx, l1para->fftnr, 0.0f); 
  myinit(taup1, npx, l1para->fftnr, 0.0f); 

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
  
  datamode = 0;
  c_l1inv_TAUP_RV_MTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, 0);
  ARRAY_MATH(residual1_p, input_p, +, dataxt_p, nsgtrace, l1para->nsamp);

  if (l1para->pmethod == 5)
    {
      ARRAY_MATH(output, output, +, dataxt_p, nsgtrace, l1para->nsamp); //// <----
    }
  else if (l1para->pmethod == 1)
    {
      datamode = 1;
      l1para->zmethod = 0;
      c_l1inv_TAUP_RV_MTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, 0);
      ARRAY_MATH(output, output, -, dataxt_p, nsgtrace, l1para->nsamp); //// <----
    }

  datamode = 0;
  datum = 0;
  c_l1inv_TAUP_RV_CTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
  ARRAY_MATH(residual1_p, residual1_p, -, dataxt_p, nsgtrace, l1para->nsamp);

  if (l1para->pmethod == 5)
    {
      ARRAY_MATH(output, output, +, dataxt_p, nsgtrace, l1para->nsamp); //// <----
    }
  else if (l1para->pmethod == 1||l1para->pmethod == 3)
    {
      datamode = 1;
      l1para->zmethod = 0;
      datum = 0;
      c_l1inv_TAUP_RV_CTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
      ARRAY_MATH(output, output, +, dataxt_p, nsgtrace, l1para->nsamp); //// <----
    }

  if( !(l1para->resscale<0.0f) ) // user choose to skip the fitting of weak energy if <0
    {
      npx = np_orig;

      myinit(pweight,npx,l1para->fftnr,1.0f);
      myinit(tweight,npx,l1para->fftnr,1.0f);

      myinit(taup,  npx, l1para->fftnr, 0.0f); 
      myinit(taup1, npx, l1para->fftnr, 0.0f); 

      for (ip = 0; ip < np_orig; ip++) indx[ip] = ip;

      calc_matrix_A_T_p(l1para, nsgtrace, offsetx, offsety, recz, lfid);

      c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(3,l1para->niter*0.2*l1para->resscale), nsgtrace, MIN(l1para->lowstart/deltaf,hfid),
		   lfid0, npx, MIN(0.8f,sppow*0.5f), pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para); 
  
      c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(3,l1para->niter*0.2*l1para->resscale), nsgtrace, MIN((l1para->lowstart+5.0f)/deltaf,hfid),
		   lfid0, npx, MIN(0.8f,sppow*0.5f), pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para); 

      np_reduce = reorder(np_orig, l1para, indx, taup, 1.0f);
      npx = np_reduce;

      recalc_matrix_A_T_p(l1para, nsgtrace, offsetx, offsety, recz, indx, np_reduce, lfid);

      c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(4,l1para->niter*0.2*l1para->resscale), nsgtrace, MIN(l1para->highfreq*0.5f/deltaf,hfid),
		   lfid, npx, sppow*0.5f, pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para); 
  
      c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter*0.4*l1para->resscale), nsgtrace, MIN(l1para->highfreq/deltaf,hfid),    
		   lfid, npx, sppow*0.5f, pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para);  

      datamode = 0;
      datum = 0;
      c_l1inv_TAUP_RV_CTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
      ARRAY_MATH(residual1_p, residual1_p, -, dataxt_p, nsgtrace, l1para->nsamp);  //// Not used at this moment

      if (l1para->pmethod == 5)
	{
	  ARRAY_MATH(output, output, +, dataxt_p, nsgtrace, l1para->nsamp); //// <----
	}
      else if (l1para->pmethod == 1||l1para->pmethod == 3)
	{
	  datamode = 1;
	  l1para->zmethod = 0;
	  datum = 0;
	  c_l1inv_TAUP_RV_CTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
	  ARRAY_MATH(output, output, +, dataxt_p, nsgtrace, l1para->nsamp); //// <----
	}

      datamode = 0;
      c_l1inv_TAUP_RV_MTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, 0);
      ARRAY_MATH(residual1_p, residual1_p, +, dataxt_p, nsgtrace, l1para->nsamp); //// Not used at this moment

      if (l1para->pmethod == 5)
	{
	  ARRAY_MATH(output, output, +, dataxt_p, nsgtrace, l1para->nsamp); //// <----
	}
      else if (l1para->pmethod == 1)
	{
	  datamode = 1;
	  l1para->zmethod = 0;
	  c_l1inv_TAUP_RV_MTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, 0);
	  ARRAY_MATH(output, output, -, dataxt_p, nsgtrace, l1para->nsamp); //// <----
	}
    } // this ends the check on resscale

  free(indx);

  flexible_free_array2d(dataxt_p);        //myfree (dataxt,nsgtrace);
  flexible_free_array2d(residual1_p);     //myfree (residual1,nsgtrace);

  flexible_free_array2d(taup);         //myfree (taup,   np_orig);
  flexible_free_array2d(taup1);        //myfree (taup,   np_orig);
  flexible_free_array2d(pweight);      //myfree (pweight,np_orig);
  flexible_free_array2d(tweight);      //myfree (tweight,np_orig);

  PFL_FREE(l1para->fftplan[0]); //kn
  PFL_FREE(l1para->fftplan[1]); //kn
}

//// equivalent V,Ay,Ayup, etc.
void ms3d_equivalyfromp(l1inv_t *l1para, int nsgtrace, float *offsetx, float *offsety, float *recz,
		      float **input_p, float **input_az, float **input_ay, 
		      float **output, float sppow, int lfid)
{
  int npx, npy;
  int ip, itrc, isamp;
  float fnyq, deltaf;

  float maxamp=-1.0f;
  
  float **taup,**taup1;
  float **pweight,**tweight;

  float **dataxt_p, **residual1_p;
  float **dataxt_az, **residual1_az;
  float **dataxt_ay, **residual1_ay;

  int np_orig,np_reduce;
  int *indx;

  int hfid,lfid0,datamode,datum;

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
  
  npx = npx*npy; 
  np_orig = npx;
  
  dataxt_p    = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));
  residual1_p = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));
    
  mycopy(input_p,dataxt_p,nsgtrace,l1para->nsamp);
  for (itrc=0;itrc<nsgtrace;itrc++) 
    for (isamp=l1para->nsamp;isamp<l1para->fftnr;isamp++) 
      dataxt_p[itrc][isamp] = 0.0;

  l1para->nsgtracepad   = ceil((float)(nsgtrace  /4.))*4;
  l1para->npxpad = ceil((float)(npx/4.))*4;
  
  taup   = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));
  taup1 = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));
  
  pweight = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));
  tweight = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));

  indx = (int *)calloc(np_orig,sizeof(int));

  myinit(pweight,npx,l1para->fftnr,1.0f);
  myinit(tweight,npx,l1para->fftnr,1.0f);

  myinit(output,nsgtrace,l1para->nsamp,0.0f);

  myinit(taup,  npx, l1para->fftnr, 0.0f); 
  myinit(taup1, npx, l1para->fftnr, 0.0f); 

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
  
  datamode = 0;
  c_l1inv_TAUP_RV_MTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, 0);
  ARRAY_MATH(residual1_p, input_p, +, dataxt_p, nsgtrace, l1para->nsamp);

  if (l1para->pmethod == 6)
    {
      ARRAY_MATH(output, output, -, dataxt_p, nsgtrace, l1para->nsamp); //// <----
    }
  else if (l1para->pmethod == 2)
    {
      datamode = 2;
      l1para->ymethod = 0;
      c_l1inv_TAUP_RV_MTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, 0);
      ARRAY_MATH(output, output, -, dataxt_p, nsgtrace, l1para->nsamp); //// <----
    }

  datamode = 0;
  datum = 0;
  c_l1inv_TAUP_RV_CTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
  ARRAY_MATH(residual1_p, residual1_p, -, dataxt_p, nsgtrace, l1para->nsamp);

  if (l1para->pmethod == 6)
    {
      ARRAY_MATH(output, output, +, dataxt_p, nsgtrace, l1para->nsamp); //// <----
    }
  else if (l1para->pmethod == 2||l1para->pmethod == 4)
    {
      datamode = 2;
      l1para->ymethod = 0;
      datum = 0;
      c_l1inv_TAUP_RV_CTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
      ARRAY_MATH(output, output, +, dataxt_p, nsgtrace, l1para->nsamp); //// <----
    }

  if( !(l1para->resscale<0.0f) ) // user choose to skip the fitting of weak energy if <0
    {
      npx = np_orig;

      myinit(pweight,npx,l1para->fftnr,1.0f);
      myinit(tweight,npx,l1para->fftnr,1.0f);

      myinit(taup,  npx, l1para->fftnr, 0.0f); 
      myinit(taup1, npx, l1para->fftnr, 0.0f); 

      for (ip = 0; ip < np_orig; ip++) indx[ip] = ip;

      calc_matrix_A_T_p(l1para, nsgtrace, offsetx, offsety, recz, lfid);

      c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(3,l1para->niter*0.2*l1para->resscale), nsgtrace, MIN(l1para->lowstart/deltaf,hfid),
		   lfid0, npx, MIN(0.8f,sppow*0.5f), pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para); 
  
      c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(3,l1para->niter*0.2*l1para->resscale), nsgtrace, MIN((l1para->lowstart+5.0f)/deltaf,hfid),
		   lfid0, npx, MIN(0.8f,sppow*0.5f), pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para); 

      np_reduce = reorder(np_orig, l1para, indx, taup, 1.0f);
      npx = np_reduce;

      recalc_matrix_A_T_p(l1para, nsgtrace, offsetx, offsety, recz, indx, np_reduce, lfid);

      c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(4,l1para->niter*0.2*l1para->resscale), nsgtrace, MIN(l1para->highfreq*0.5f/deltaf,hfid),
		   lfid, npx, sppow*0.5f, pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para); 
  
      c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter*0.4*l1para->resscale), nsgtrace, MIN(l1para->highfreq/deltaf,hfid),    
		   lfid, npx, sppow*0.5f, pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para);  

      datamode = 0;
      datum = 0;
      c_l1inv_TAUP_RV_CTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
      ARRAY_MATH(residual1_p, residual1_p, -, dataxt_p, nsgtrace, l1para->nsamp);  //// Not used at this moment

      if (l1para->pmethod == 6)
	{
	  ARRAY_MATH(output, output, +, dataxt_p, nsgtrace, l1para->nsamp); //// <----
	}
      else if (l1para->pmethod == 2||l1para->pmethod == 4)
	{
	  datamode = 2;
	  l1para->ymethod = 0;
	  datum = 0;
	  c_l1inv_TAUP_RV_CTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
	  ARRAY_MATH(output, output, +, dataxt_p, nsgtrace, l1para->nsamp); //// <----
	}

      datamode = 0;
      c_l1inv_TAUP_RV_MTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, 0);
      ARRAY_MATH(residual1_p, residual1_p, +, dataxt_p, nsgtrace, l1para->nsamp); //// Not used at this moment

      if (l1para->pmethod == 6)
	{
	  ARRAY_MATH(output, output, -, dataxt_p, nsgtrace, l1para->nsamp); //// <----
	}
      else if (l1para->pmethod == 2)
	{
	  datamode = 2;
	  l1para->ymethod = 0;
	  c_l1inv_TAUP_RV_MTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, 0);
	  ARRAY_MATH(output, output, -, dataxt_p, nsgtrace, l1para->nsamp); //// <----
	}
    } // this ends the check on resscale

  free(indx);

  flexible_free_array2d(dataxt_p);        //myfree (dataxt,nsgtrace);
  flexible_free_array2d(residual1_p);     //myfree (residual1,nsgtrace);

  flexible_free_array2d(taup);         //myfree (taup,   np_orig);
  flexible_free_array2d(taup1);        //myfree (taup,   np_orig);
  flexible_free_array2d(pweight);      //myfree (pweight,np_orig);
  flexible_free_array2d(tweight);      //myfree (tweight,np_orig);

  PFL_FREE(l1para->fftplan[0]); //kn
  PFL_FREE(l1para->fftplan[1]); //kn
}

void calc_p(l1inv_t *l1para)
{

  int iq,ip,ipx,ipy,npx,npy,nq;
  float rq,px,py;
  float *pxsavetmp,*pysavetmp,*rqsavetmp;

  calc_px(l1para);
  calc_py(l1para);
  calc_rq(l1para);

  npx = l1para->npxsave;
  npy = l1para->npysave;
  nq  = l1para->nqsave;

  pxsavetmp = (float *) malloc(npx*npy*nq*sizeof(float));
  pysavetmp = (float *) malloc(npx*npy*nq*sizeof(float));
  rqsavetmp = (float *) malloc(npx*npy*nq*sizeof(float));

  ip = -1;

  for ( iq = 0;iq < nq; iq++ )
   {
  for ( ipy = 0;ipy < npy; ipy++ )
   {
     for ( ipx = 0;ipx < npx; ipx++ )
       {
	 ip = ip + 1;
	 
	 px = l1para->pxsave[ipx];
	 py = l1para->pysave[ipy];
	 rq = l1para->rqsave[iq];

	 pxsavetmp[ip] = px;
	 pysavetmp[ip] = py;
	 rqsavetmp[ip] = rq;
       }
   }

   }

  ip = -1;
  for ( iq = 0;iq < nq; iq++ )
   {
  for ( ipy = 0;ipy < npy; ipy++ )
    {
      for ( ipx = 0;ipx < npx; ipx++ )
	{
	  ip = ip + 1;
	  l1para->pxsave[ip] = pxsavetmp[ip];
	  l1para->pysave[ip] = pysavetmp[ip];
	  l1para->rqsave[ip] = rqsavetmp[ip];
	  
	}
    }
   }

  free(pxsavetmp);
  free(pysavetmp);
  free(rqsavetmp);
}

void calc_px(l1inv_t *l1para)
{
  int ip;
  float px,deltapxmin;

  int npx = l1para->npxorig;

  l1para->npxsave = npx;
  
  deltapxmin  = (l1para->pxmax-l1para->pxmin)/(float)(npx-1);
  
  for ( ip = 0;ip < npx; ip++ )
    {
      px = l1para->pxmin + (float)(ip)*deltapxmin;
      l1para->pxsave[ip] = px;
    }
  
  return;
  
    
}

void calc_py(l1inv_t *l1para)
{
  int ip;
  float py,deltapymin;

  int npy = l1para->npyorig;

  l1para->npysave = npy;
  
  deltapymin  = (l1para->pymax-l1para->pymin)/(float)(npy-1);
  
  for ( ip = 0;ip < npy; ip++ )
    {
      py = l1para->pymin + (float)(ip)*deltapymin;
      l1para->pysave[ip] = py;
    }

  return;

}

void calc_rq(l1inv_t *l1para)
{
  int iq;
  float rq,deltarq;

  int nq = l1para->nqorig;

  l1para->nqsave = nq;
  
  if (nq > 1) deltarq  = (l1para->rqmax-l1para->rqmin)/(float)(nq-1);
  else deltarq = 0.0f;

  for ( iq = 0;iq < nq; iq++ )
    {
      rq = l1para->rqmin + (float)(iq)*deltarq;
      l1para->rqsave[iq] = rq;
    }

  return;

}

void ms3d_deghosting_redatum(l1inv_t *l1para, int ntrcxy, float *offsetx, float *offsety, float *recz,
			     float **input_p, float **input_az, float **input_ay, 
			     float **outputi, float **output, float sppow, int lfid)
  
{
  
  int pyes=l1para->ms_pyes;
  int zyes=l1para->ms_zyes;
  int yyes=l1para->ms_yyes;
  int msdataz = l1para->msdataz;
  int msdatay = l1para->msdatay;
  int outtype = l1para->outtype;

  int datum = l1para->datum;
  int dgredatumflag = l1para->dgredatumflag;
  int fpredatum = l1para->fpredatum;
  int outpos = l1para->outpos;

  float zlowcut = l1para->zlowcut;
  float ylowcut = l1para->ylowcut;

  l1para->datum = CABLE90;
  l1para->dgredatumflag = 0;
  l1para->fpredatum = FPREDATUM_NO;
  l1para->outpos = 1;

  ms3d_deghosting  (l1para, ntrcxy, offsetx, offsety, recz, 
		    input_p, input_az, input_ay, outputi, 
		    sppow, lfid);
  
  l1para->ms_pyes = 1;
  l1para->ms_zyes = 0;
  l1para->ms_yyes = 0;
  l1para->msdataz = 0;
  l1para->msdatay = 0;
  l1para->outtype = 0;

  l1para->datum = datum;
  l1para->dgredatumflag = 1;
  l1para->fpredatum = FPREDATUM_CTOS;
  l1para->outpos = 4;
  
  l1para->zlowcut = 0.0f;
  l1para->ylowcut = 0.0f;

  ms3d_deghosting  (l1para, ntrcxy, offsetx, offsety, recz, 
		    outputi, input_az, input_ay, output, 
		    sppow, lfid);
  
  l1para->ms_pyes=pyes;
  l1para->ms_zyes=zyes;
  l1para->ms_yyes=yyes;
  l1para->msdataz = msdataz;
  l1para->msdatay = msdatay;
  l1para->outtype = outtype;

  l1para->datum = datum;
  l1para->dgredatumflag = dgredatumflag;
  l1para->fpredatum = fpredatum;
  l1para->outpos = outpos;

  l1para->zlowcut = zlowcut;
  l1para->ylowcut = ylowcut;
  
}

// CHECK : compare with jointsr3d_deghosting 

void jointmc_deghosting (l1inv_t *l1para, int nsgtrace, float *offsetx, float *offsety, float *recz,
		       float **input_p, float **input_az, float **input_ay, 
		       float **output, float sppow, int lfid)
{
  int npx, npy, nsgtracep, nsgtrace_az, nsgtrace_ay;
  int ip, itrc, isamp;
  float fnyq, deltaf;

  float maxamp=-1.0f;
  
  float **taup,**taup1;
  float **taup_filter;
  float **pweight,**tweight;

  float **dataxt_p, **residual1_p, **residual1_ponly;
  float **dataxt_az, **residual1_az;
  float **dataxt_ay, **residual1_ay;

  float **dataxt_tmp;

  int np_orig,np_reduce;
  int *indx;

  int hfid,lfid0,datamode,datum;

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

  if(l1para->zerotraceattr == 1)
    {
      nsgtracep   = l1para->nsgtracearr4d[0];
      nsgtrace_az = l1para->nsgtracearr4d[1];
      nsgtrace_ay = l1para->nsgtracearr4d[2];
    }
  else
    {
      nsgtracep   = nsgtrace;
      nsgtrace_az = nsgtrace;
      nsgtrace_ay = nsgtrace;
    }

  calc_matrix_A_T_p (l1para, nsgtracep, offsetx, offsety, recz, lfid);

  PFL_GET(&l1para->fftplan[0], PFL_PORTABLE_DEF); //kn
  PFL_GET(&l1para->fftplan[1], PFL_PORTABLE_DEF); //kn

  
  fnyq = 500000.0/l1para->srate;
  deltaf = fnyq/(float)(l1para->fftnc-1);

  lfid0 = 2.0f/deltaf; //// 2 Hz;

  ////index of the maximum/minimum frequency
  hfid = min(floor((float)(l1para->fftnc)*l1para->highfreq/fnyq),l1para->fftnc);

  npx = l1para->npxsave;
  npy = l1para->npysave;
  
  npx = npx*npy; 
  np_orig = npx;
  
  dataxt_p    = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));
  residual1_p = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));
  residual1_ponly = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));

  dataxt_tmp    = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));

 
  mycopy(input_p,dataxt_p,nsgtrace,l1para->nsamp);
  for (itrc=0;itrc<nsgtrace;itrc++) 
    for (isamp=l1para->nsamp;isamp<l1para->fftnr;isamp++) 
      dataxt_p[itrc][isamp] = 0.0;

  if ( l1para->msdataz == YESYES )
    {	
      dataxt_az    = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));
      residual1_az = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));

      mycopy(input_az,dataxt_az,nsgtrace,l1para->nsamp);
      for (itrc=0;itrc<nsgtrace;itrc++) 
	for (isamp=l1para->nsamp;isamp<l1para->fftnr;isamp++) 
	  dataxt_az[itrc][isamp] = 0.0;
      
      datamode = 1;
      if (l1para->iexternlc == 0) ms_lowcutf(dataxt_az, nsgtrace, l1para, datamode);

      mycopy(dataxt_az,input_az,nsgtrace,l1para->nsamp);

    }

  if ( l1para->msdatay == YESYES )
    {	
      dataxt_ay    = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));
      residual1_ay = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));
      
      mycopy(input_ay,dataxt_ay,nsgtrace,l1para->nsamp);
      for (itrc=0;itrc<nsgtrace;itrc++) 
	for (isamp=l1para->nsamp;isamp<l1para->fftnr;isamp++) 
	  dataxt_ay[itrc][isamp] = 0.0;
      
      datamode = 2;
      ms_lowcutf(dataxt_ay, nsgtrace, l1para, datamode);
      
      mycopy(dataxt_ay,input_ay,nsgtrace,l1para->nsamp);
      
    }

  l1para->nsgtracepad   = ceil((float)(nsgtrace  /4.))*4;
  l1para->npxpad = ceil((float)(npx/4.))*4;
  
  taup   = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));
  taup1 = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));
  taup_filter   = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));

  pweight = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));
  tweight = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));

  indx = (int *)calloc(np_orig,sizeof(int));

  myinit(pweight,npx,l1para->fftnr,1.0f);
  myinit(tweight,npx,l1para->fftnr,1.0f);

  myinit(output,nsgtrace,l1para->nsamp,0.0f);

  myinit(taup,  npx, l1para->fftnr, 0.0f); 
  myinit(taup1, npx, l1para->fftnr, 0.0f); 

  for (ip = 0; ip < np_orig; ip++) indx[ip] = ip;

  if (l1para->lowguide==0)
    {

      l1para->ms_pyes=1;
      l1para->ms_zyes=0;
      l1para->ms_yyes=0;
     
      c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(5,l1para->niterlow*0.35), nsgtrace, MIN(l1para->lowstart/deltaf,hfid),
		   lfid0, npx, MIN(1.0f,sppow*0.5f), pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight, l1para);  
      
      c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niterlow*0.35), nsgtrace, MIN((l1para->lowstart+5.0f)/deltaf,hfid),
		   lfid0, npx, MIN(1.0f,sppow*0.5f), pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight, l1para);
      
      l1para->ms_pyes=pyes;
      l1para->ms_zyes=zyes;
      l1para->ms_yyes=yyes;
      
    }
  else
    {
      c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(5,l1para->niter*0.35), nsgtrace, MIN(l1para->lowstart/deltaf,hfid),
		   lfid0, npx, MIN(1.0f,sppow*0.5f), pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight, l1para);  
      
      c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter*0.35), nsgtrace, MIN((l1para->lowstart+5.0f)/deltaf,hfid),
		   lfid0, npx, MIN(1.0f,sppow*0.5f), pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight, l1para);
    }

  np_reduce = reorder(np_orig, l1para, indx, taup, 0.5f);

  recalc_matrix_A_T_p(l1para, nsgtracep, offsetx, offsety, recz, indx, np_reduce, lfid);

  ////G2_logInfo("bsdg","nptotal=%d, np_reduce=%d",npx,np_reduce);fflush(0);
  npx = np_reduce;

  if ((l1para->outtype==3||l1para->outtype==5)&&(l1para->zlowcut>=1.0f||l1para->lzfilter>0))
    {
      l1para->ms_pyes=1;
      l1para->ms_zyes=0;
      l1para->ms_yyes=0;
      
      c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(5,l1para->niterlow*0.4), nsgtrace, MIN(l1para->highfreq_p*0.5f/deltaf,hfid),
		   lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup1, tweight, l1para); 
      
      c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niterlow), nsgtrace, MIN(l1para->highfreq_p/deltaf,hfid),    
		   lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup1, tweight,l1para);  
      
      if (l1para->lowguide==0)
	{
	  datamode = 0;
	  datum = 0;
          c_l1inv_TAUP_RV_STOS_MS(nsgtracep, hfid, lfid, npx, taup1, dataxt_tmp, l1para, datamode, datum, 0);
          ARRAY_MATH(residual1_ponly, input_p, -, dataxt_tmp, nsgtracep, l1para->nsamp);
	}

      l1para->ms_pyes=pyes;
      l1para->ms_zyes=zyes;
      l1para->ms_yyes=yyes;

      c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(5,l1para->niter*0.4), nsgtrace, MIN(l1para->highfreq*0.5f/deltaf,hfid),
		   lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight, l1para); 
      
      c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter), nsgtrace, MIN(l1para->highfreq/deltaf,hfid),    
		   lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight,l1para);  
      
      l1inv3d_taup_merge(l1para, hfid, npx, taup, taup1);    //// output is taup
      
    }
  else
    {

      if (l1para->lowguide==0)
	{

	  l1para->ms_pyes=1;
	  l1para->ms_zyes=0;
	  l1para->ms_yyes=0;
	  
	  c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(5,l1para->niter*0.4), nsgtrace, MIN(l1para->highfreq*0.5f/deltaf,hfid),
		       lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup1, tweight, l1para); 
	  
	  c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter), nsgtrace, MIN(l1para->highfreq/deltaf,hfid),    
		       lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup1, tweight,l1para);  
          datamode = 0;
	  datum = 0;
          c_l1inv_TAUP_RV_STOS_MS(nsgtracep, hfid, lfid, npx, taup1, dataxt_tmp, l1para, datamode, datum, 0);
          ARRAY_MATH(residual1_ponly, input_p, -, dataxt_tmp, nsgtracep, l1para->nsamp);

	  l1para->ms_pyes=pyes;
	  l1para->ms_zyes=zyes;
	  l1para->ms_yyes=yyes;
	
	}

      c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(5,l1para->niter*0.4), nsgtrace, MIN(l1para->highfreq*0.5f/deltaf,hfid),
		   lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight, l1para); 
      
      c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter), nsgtrace, MIN(l1para->highfreq/deltaf,hfid),    
		   lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight,l1para);  

    }
  
  mycopy (taup, taup_filter, npx, l1para->fftnr);
  if (l1para->ldipf > 0)
    apply_dipfilter3d (npx, indx, taup_filter, l1para);
  
  if (l1para->fpredatum==FPREDATUM_NO)
    {
      datamode = 0;
      datum = 0;
      c_l1inv_TAUP_RV_STOS_MS(nsgtracep, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum, 0);
      ARRAY_MATH(residual1_p, input_p, -, dataxt_p, nsgtracep, l1para->nsamp);

      datamode = 0;
      datum = l1para->datum;
      c_l1inv_TAUP_RV_STOS_MS(nsgtracep, hfid, lfid, npx, taup_filter, dataxt_p, l1para, datamode, datum, 1);

      if (l1para->jointmethod != JOINTDEGHOST)
	convolve_dataxt_with_target (nsgtracep, l1para->fftnr, dataxt_p, l1para); //// ** SRAY ** CHECK **

      ARRAY_MATH(output, output, +, dataxt_p, nsgtracep, l1para->nsamp);
    }
  else if (l1para->fpredatum==FPREDATUM_CTOS)
    {
      datamode = 0;
      datum = 0;
      c_l1inv_TAUP_RV_CTOS(nsgtracep, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
      ARRAY_MATH(residual1_p, input_p, -, dataxt_p, nsgtracep, l1para->nsamp);
      
      datamode = 0;
      datum = l1para->datum;
      c_l1inv_TAUP_RV_CTOS(nsgtracep, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
      ARRAY_MATH(output, output, +, dataxt_p, nsgtracep, l1para->nsamp);

    }

  if (l1para->ms_zyes==1)
    {
      datamode = 1;
      datum = 0;
      c_l1inv_TAUP_RV_STOS_MS(nsgtrace_az, hfid, lfid, npx, taup, dataxt_az, l1para, datamode, datum, 0);
      ARRAY_MATH(residual1_az, input_az, -, dataxt_az, nsgtrace_az, l1para->nsamp);
    }

  if (l1para->ms_yyes==1)
    {
      datamode = 2;
      datum = 0;
      c_l1inv_TAUP_RV_STOS_MS(nsgtrace_ay, hfid, lfid, npx, taup, dataxt_ay, l1para, datamode, datum, 0);
      ARRAY_MATH(residual1_ay, input_ay, -, dataxt_ay, nsgtrace_ay, l1para->nsamp);
    }

  if( !(l1para->resscale<0.0f) ) // user choose to skip the fitting of weak energy if <0
    {
      npx = np_orig;

      myinit(pweight,npx,l1para->fftnr,1.0f);
      myinit(tweight,npx,l1para->fftnr,1.0f);

      myinit(taup,  npx, l1para->fftnr, 0.0f); 
      myinit(taup1, npx, l1para->fftnr, 0.0f); 

      for (ip = 0; ip < np_orig; ip++) indx[ip] = ip;

      calc_matrix_A_T_p(l1para, nsgtracep, offsetx, offsety, recz, lfid);

      if (l1para->lowguide==0)
	{

	  l1para->ms_pyes=1;
	  l1para->ms_zyes=0;
	  l1para->ms_yyes=0;
      
	  c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(3,l1para->niter*0.2*l1para->resscale), nsgtrace, MIN(l1para->lowstart/deltaf,hfid),
		       lfid0, npx, MIN(0.8f,sppow*0.5f), pweight, residual1_ponly, residual1_az, residual1_ay, taup, tweight, l1para); 
      
	  c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(3,l1para->niter*0.2*l1para->resscale), nsgtrace, MIN((l1para->lowstart+5.0f)/deltaf,hfid),
		       lfid0, npx, MIN(0.8f,sppow*0.5f), pweight, residual1_ponly, residual1_az, residual1_ay, taup, tweight, l1para);
      
	  l1para->ms_pyes=pyes;
	  l1para->ms_zyes=zyes;
	  l1para->ms_yyes=yyes;
      
	}
      else
	{
	  c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(3,l1para->niter*0.2*l1para->resscale), nsgtrace, MIN(l1para->lowstart/deltaf,hfid),
		       lfid0, npx, MIN(0.8f,sppow*0.5f), pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para); 
      
	  c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(3,l1para->niter*0.2*l1para->resscale), nsgtrace, MIN((l1para->lowstart+5.0f)/deltaf,hfid),
		       lfid0, npx, MIN(0.8f,sppow*0.5f), pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para); 
	}

      np_reduce = reorder(np_orig, l1para, indx, taup, 1.0f);
      npx = np_reduce;

      recalc_matrix_A_T_p(l1para, nsgtracep, offsetx, offsety, recz, indx, np_reduce, lfid);

      if ((l1para->outtype==3||l1para->outtype==5)&&(l1para->zlowcut>=1.0f||l1para->lzfilter>0)) //// <-- condition possibly different in production version v203.03
	{

	  l1para->ms_pyes=1;
	  l1para->ms_zyes=0;
	  l1para->ms_yyes=0;
      
	  c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(4,l1para->niter*0.2*l1para->resscale), nsgtrace, MIN(l1para->highfreq_p*0.5f/deltaf,hfid),
		       lfid, npx, sppow*0.5f, pweight, residual1_ponly, residual1_az, residual1_ay, taup1, tweight, l1para); 
      
	  c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter*0.4*l1para->resscale), nsgtrace, MIN(l1para->highfreq_p/deltaf,hfid),    
		       lfid, npx, sppow*0.5f, pweight, residual1_ponly, residual1_az, residual1_ay, taup1, tweight, l1para);  

	  l1para->ms_pyes=pyes;
	  l1para->ms_zyes=zyes;
	  l1para->ms_yyes=yyes;

	  c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(4,l1para->niter*0.2*l1para->resscale), nsgtrace, MIN(l1para->highfreq*0.5f/deltaf,hfid),
		       lfid, npx, sppow*0.5f, pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para); 
      
	  c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter*0.4*l1para->resscale), nsgtrace, MIN(l1para->highfreq/deltaf,hfid),    
		       lfid, npx, sppow*0.5f, pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para);  

	  l1inv3d_taup_merge(l1para, hfid, npx, taup, taup1);    

	}

      else
	{

	  c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(4,l1para->niter*0.2*l1para->resscale), nsgtrace, MIN(l1para->highfreq*0.5f/deltaf,hfid),
		       lfid, npx, sppow*0.5f, pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para); 
      
	  c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter*0.4*l1para->resscale), nsgtrace, MIN(l1para->highfreq/deltaf,hfid),    
		       lfid, npx, sppow*0.5f, pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para);  

	}

      mycopy (taup, taup_filter, npx, l1para->fftnr);
      if (l1para->ldipf > 0)
	apply_dipfilter3d (npx, indx, taup_filter, l1para);

      if (l1para->fpredatum==FPREDATUM_NO)
	{
	  datamode = 0;
	  datum = 0;
	  c_l1inv_TAUP_RV_STOS_MS(nsgtracep, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum, 0);
	  ARRAY_MATH(residual1_p, residual1_p, -, dataxt_p, nsgtracep, l1para->nsamp);

	  datamode = 0;
	  datum = l1para->datum;
	  c_l1inv_TAUP_RV_STOS_MS(nsgtracep, hfid, lfid, npx, taup_filter, dataxt_p, l1para, datamode, datum, 1);

	  if (l1para->jointmethod != JOINTDEGHOST)
	    convolve_dataxt_with_target (nsgtracep, l1para->fftnr, dataxt_p, l1para); //// ** SRAY ** CHECK **
      
	  ARRAY_MATH(output, output, +, dataxt_p, nsgtracep, l1para->nsamp);
	}
      else if (l1para->fpredatum==FPREDATUM_CTOS)
	{
	  datamode = 0;
	  datum = 0;
	  c_l1inv_TAUP_RV_CTOS(nsgtracep, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
	  ARRAY_MATH(residual1_p, residual1_p, -, dataxt_p, nsgtracep, l1para->nsamp);
      
	  datamode = 0;
	  datum = l1para->datum;
	  c_l1inv_TAUP_RV_CTOS(nsgtracep, hfid, lfid, npx, taup, dataxt_p, l1para, datamode, datum);
	  ARRAY_MATH(output, output, +, dataxt_p, nsgtracep, l1para->nsamp);
	}
    } // this ends the check on resscale

  if (!((l1para->ms_pyes==0&&(l1para->zmethod==0))||(l1para->ms_pyes==0&&(l1para->ymethod==0))))
    {
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

      ARRAY_MATH(output, factor*l1para->residual*residual1_p, +, output, nsgtracep, l1para->nsamp); 
    }

  if (l1para->taupQC==2)
    {
      datamode = 1;
      c_l1inv_TAUP_RV_MTOS(nsgtrace_az, hfid, lfid, npx, taup, dataxt_az, l1para, datamode, 0);
      ARRAY_MATH(residual1_az, residual1_az, +, dataxt_az, nsgtrace_az, l1para->nsamp);
      
      datamode = 1;
      c_l1inv_TAUP_RV_CTOS(nsgtrace_az, hfid, lfid, npx, taup, dataxt_az, l1para, datamode, 0);
      ARRAY_MATH(residual1_az, residual1_az, -, dataxt_az, nsgtrace_az, l1para->nsamp);
    }

  if (l1para->taupQC==1) 
    {
      for (itrc=0;itrc<nsgtrace;itrc++){
	for (isamp=0;isamp<l1para->nsamp;isamp++){
	  output[itrc][isamp] = l1para->designorm_factor*residual1_p[itrc][isamp];
	}
      }
    }
  if (l1para->taupQC==2) 
    {
      for (itrc=0;itrc<nsgtrace;itrc++){
	for (isamp=0;isamp<l1para->nsamp;isamp++){
	  output[itrc][isamp] = residual1_az[itrc][isamp];
	}
      }
    }

  free(indx);

  flexible_free_array2d(dataxt_p);        //myfree (dataxt,nsgtrace);
  flexible_free_array2d(residual1_p);     //myfree (residual1,nsgtrace);
  flexible_free_array2d(residual1_ponly); //myfree (residual1,nsgtrace);

  flexible_free_array2d(dataxt_tmp); 

  if ( l1para->msdataz == YESYES ) 
    {
      flexible_free_array2d(dataxt_az);       //myfree (dataxt,nsgtrace);
      flexible_free_array2d(residual1_az);    //myfree (residual1,nsgtrace);
    }

  if ( l1para->msdatay == YESYES ) 
    {
      flexible_free_array2d(dataxt_ay);       //myfree (dataxt,nsgtrace);
      flexible_free_array2d(residual1_ay);    //myfree (residual1,nsgtrace);
    }

  flexible_free_array2d(taup);          //myfree (taup,   np_orig);
  flexible_free_array2d(taup1);        //myfree (taup,   np_orig);
  flexible_free_array2d(taup_filter);   //myfree (taup,   np_orig);
  flexible_free_array2d(pweight);       //myfree (pweight,np_orig);
  flexible_free_array2d(tweight);       //myfree (tweight,np_orig);

  PFL_FREE(l1para->fftplan[0]); //kn
  PFL_FREE(l1para->fftplan[1]); //kn
}


void ms3d_xghosting_wrapper(l1inv_t *l1para, int nsgtrace, float *offsetx, float *offsety, float *recz,
		       float **input_p, float **input_az, float **input_ay, 
		       float **output, float sppow, int lfid)
{

  int nsgtracep,itrc;

  if (l1para->outtype == VCROSSGHOST || l1para->outtype == VZTOVCR)
    {

      l1para->fpredatum = FPREDATUM_XG_Z;

      l1para->msdatap = 0;
      l1para->msdataz = 1;
      l1para->msdatay = 0;

      l1para->ms_pyes = 0;
      l1para->ms_zyes = 1;
      l1para->ms_yyes = 0;

	  // zrfl has become a pointer.
	  // In order to pass compilation, make the following change
      //l1para->zrfl = -l1para->zrfl; 
      l1para->zrfl[0] = -l1para->zrfl[0]; 

      ms3d_xghosting_z(l1para, nsgtrace, offsetx, offsety, recz, input_az, output, sppow, lfid);

      //l1para->zrfl = -l1para->zrfl; 
      l1para->zrfl[0] = -l1para->zrfl[0]; 

      l1para->msdatap = 1;
      l1para->msdataz = 1;
      l1para->msdatay = 0;

      l1para->ms_pyes = 1;
      l1para->ms_zyes = 1;
      l1para->ms_yyes = 0;

    }

  if (l1para->outtype == VZTOV)
    {
      
      l1para->fpredatum = FPREDATUM_NO;
      
      l1para->outtype = 6;

      l1para->ms_pyes = 0;
      l1para->ms_zyes = 1;
      l1para->ms_yyes = 0;
      
      //int zrfl = l1para->zrfl;
      int zrfl = l1para->zrfl[0];

      //l1para->zrfl = l1para->prfl;
      l1para->zrfl[0] = l1para->prfl;

      if(l1para->zerotraceattr == 1)
	nsgtracep   = l1para->nsgtracearr4d[0];
      else
	nsgtracep   = nsgtrace;
      
      float *pRecz     = (float *) calloc (nsgtracep, sizeof (float));

      for (itrc = 0; itrc < nsgtracep; itrc++)
	{
	  pRecz[itrc] = recz[itrc];
	  recz[itrc]  = 1.0f;
	}

      ms3d_obliqinteg(l1para, nsgtrace, offsetx, offsety, recz, input_p, input_az, input_ay, output, sppow, lfid);

      //l1para->zrfl = zrfl;
      l1para->zrfl[0] = zrfl;

      for (itrc = 0; itrc < nsgtracep; itrc++)
	recz[itrc] = pRecz[itrc];

      l1para->fpredatum = FPREDATUM_XG_Z;

      l1para->outtype = VZTOV;

      l1para->ms_pyes = 1;
      l1para->ms_zyes = 1;
      l1para->ms_yyes = 0;

      free(pRecz);

    }

  if (l1para->outtype == PCROSSGHOST)
    {

      l1para->fpredatum = FPREDATUM_CIN;

      l1para->msdatap = 1;
      l1para->msdataz = 0;
      l1para->msdatay = 0;

      l1para->ms_pyes = 1;
      l1para->ms_zyes = 0;
      l1para->ms_yyes = 0;

      l1para->prfl = -l1para->prfl; 

      bsdg3d_deghosting(l1para, nsgtrace, offsetx, offsety, recz, input_p, output, sppow, lfid);

      l1para->fpredatum = FPREDATUM_XG_Z;

      l1para->prfl = -l1para->prfl; 

      l1para->msdatap = 1;
      l1para->msdataz = 1;
      l1para->msdatay = 0;

      l1para->ms_pyes = 1;
      l1para->ms_zyes = 1;
      l1para->ms_yyes = 0;
 
    }

}


void ms3d_xghosting_z(l1inv_t *l1para, int nsgtrace, float *offsetx, float *offsety, float *recz,
		       float **input_az, float **output, float sppow, int lfid)
{
  int npx, npy, nsgtracep, nsgtrace_az;
  int ip, itrc, isamp;
  float fnyq, deltaf;

  float maxamp=-1.0f;
  
  float **taup,**taup1;
  float **pweight,**tweight;

  float **dataxt_az, **residual1_az;

  int np_orig,np_reduce;
  int *indx;

  int hfid,lfid0,datamode,datum;

  int zyes=l1para->ms_zyes;

  float **input_p      = NULL;
  float **dataxt_p     = NULL;
  float **residual1_p  = NULL;
  float **input_ay     = NULL;
  float **dataxt_ay    = NULL;
  float **residual1_ay = NULL;

  for (itrc=0;itrc<nsgtrace;itrc++)
     for (isamp=0;isamp<l1para->nsamp;isamp++)
       if (fabsf(input_az[itrc][isamp])>maxamp) maxamp=fabsf(input_az[itrc][isamp]);
  test_amp(maxamp);
  if (maxamp<1e-5)
   {
    for (itrc=0;itrc<nsgtrace;itrc++)
      for (isamp=0;isamp<l1para->nsamp;isamp++) output[itrc][isamp]=input_az[itrc][isamp];
    return;
   }

  if(l1para->zerotraceattr == 1)
    {
      nsgtracep   = l1para->nsgtracearr4d[0];
      nsgtrace_az = l1para->nsgtracearr4d[1];
    }
  else
    {
      nsgtracep   = nsgtrace;
      nsgtrace_az = nsgtrace;
    }

  calc_matrix_A_T_p (l1para, nsgtracep, offsetx, offsety, recz, lfid);

  PFL_GET(&l1para->fftplan[0], PFL_PORTABLE_DEF); //kn
  PFL_GET(&l1para->fftplan[1], PFL_PORTABLE_DEF); //kn

  
  fnyq = 500000.0/l1para->srate;
  deltaf = fnyq/(float)(l1para->fftnc-1);

  lfid0 = 2.0f/deltaf; //// 2 Hz;

  ////index of the maximum/minimum frequency
  hfid = min(floor((float)(l1para->fftnc)*l1para->highfreq/fnyq),l1para->fftnc);

  npx = l1para->npxsave;
  npy = l1para->npysave;
  
  npx = npx*npy; 
  np_orig = npx;
   
  dataxt_az    = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));
  residual1_az = (float**)flexible_array2d(nsgtrace, l1para->fftnr, sizeof(float));
  
  mycopy(input_az,dataxt_az,nsgtrace,l1para->nsamp);
  for (itrc=0;itrc<nsgtrace;itrc++) 
    for (isamp=l1para->nsamp;isamp<l1para->fftnr;isamp++) 
      dataxt_az[itrc][isamp] = 0.0;
  mycopy(dataxt_az,input_az,nsgtrace,l1para->nsamp);

  l1para->nsgtracepad   = ceil((float)(nsgtrace  /4.))*4;
  l1para->npxpad = ceil((float)(npx/4.))*4;
  
  taup   = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));
  taup1 = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));
  
  pweight = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));
  tweight = (float**)flexible_array2d(npx, l1para->fftnr, sizeof(float));

  indx = (int *)calloc(np_orig,sizeof(int));

  myinit(pweight,npx,l1para->fftnr,1.0f);
  myinit(tweight,npx,l1para->fftnr,1.0f);

  myinit(output,nsgtrace,l1para->nsamp,0.0f);

  myinit(taup,  npx, l1para->fftnr, 0.0f); 
  myinit(taup1, npx, l1para->fftnr, 0.0f); 

  for (ip = 0; ip < np_orig; ip++) indx[ip] = ip;

  c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(5,l1para->niter*0.35), nsgtrace, MIN(l1para->lowstart/deltaf,hfid),
	       lfid0, npx, MIN(1.0f,sppow*0.5f), pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight, l1para);  
  
  c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter*0.35), nsgtrace, MIN((l1para->lowstart+5.0f)/deltaf,hfid),
		   lfid0, npx, MIN(1.0f,sppow*0.5f), pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight, l1para);
  
  np_reduce = reorder(np_orig, l1para, indx, taup, 0.5f);

  recalc_matrix_A_T_p(l1para, nsgtracep, offsetx, offsety, recz, indx, np_reduce, lfid);

  ////G2_logInfo("bsdg","nptotal=%d, np_reduce=%d",npx,np_reduce);fflush(0);
  npx = np_reduce;

  c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(5,l1para->niter*0.4), nsgtrace, MIN(l1para->highfreq*0.5f/deltaf,hfid),
	       lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight, l1para); 
  
  c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter), nsgtrace, MIN(l1para->highfreq/deltaf,hfid),    
	       lfid, npx, sppow, pweight, dataxt_p, dataxt_az, dataxt_ay, taup, tweight,l1para);  

  if (l1para->outtype == VCROSSGHOST)
    {
      datamode = 1; 
      datum = 0;
      c_l1inv_TAUP_RV_CTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_az, l1para, datamode, datum);   
      ARRAY_MATH(residual1_az, input_az, -, dataxt_az, nsgtrace, l1para->nsamp);

      datamode = 0; 
      datum = 0;
      c_l1inv_TAUP_RV_CTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_az, l1para, datamode, datum); 
      ARRAY_MATH(output, output, +, dataxt_az, nsgtrace, l1para->nsamp);
      
      datamode = 0; 
      datum = 0;
      c_l1inv_TAUP_RV_MTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_az, l1para, datamode, datum);
      ARRAY_MATH(output, output, -, dataxt_az, nsgtrace, l1para->nsamp);
    }
  else if (l1para->outtype == VZTOVCR)
    {
      datamode = 1; 
      datum = 0;
      c_l1inv_TAUP_RV_CTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_az, l1para, datamode, datum);   
      ARRAY_MATH(residual1_az, input_az, -, dataxt_az, nsgtrace, l1para->nsamp);
      
      datamode = 0; 
      datum = 0;
      c_l1inv_TAUP_RV_CTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_az, l1para, datamode, datum);   
      ARRAY_MATH(output, output, +, dataxt_az, nsgtrace, l1para->nsamp);
    }

  if( !(l1para->resscale<0.0f) ) // user choose to skip the fitting of weak energy if <0
    { 
      
      npx = np_orig;

      myinit(pweight,npx,l1para->fftnr,1.0f);
      myinit(tweight,npx,l1para->fftnr,1.0f);

      myinit(taup,  npx, l1para->fftnr, 0.0f); 
      myinit(taup1, npx, l1para->fftnr, 0.0f); 

      for (ip = 0; ip < np_orig; ip++) indx[ip] = ip;

      calc_matrix_A_T_p(l1para, nsgtracep, offsetx, offsety, recz, lfid);


      c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(3,l1para->niter*0.2*l1para->resscale), nsgtrace, MIN(l1para->lowstart/deltaf,hfid),
		   lfid0, npx, MIN(0.8f,sppow*0.5f), pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para); 
  
      c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(3,l1para->niter*0.2*l1para->resscale), nsgtrace, MIN((l1para->lowstart+5.0f)/deltaf,hfid),
		   lfid0, npx, MIN(0.8f,sppow*0.5f), pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para); 
  

      np_reduce = reorder(np_orig, l1para, indx, taup, 1.0f);
      npx = np_reduce;

      recalc_matrix_A_T_p(l1para, nsgtracep, offsetx, offsety, recz, indx, np_reduce, lfid);

      c_l1inv_lsqr(MSDGINTERP_CG_GENERAL,l1para->tolerance, MAX(4,l1para->niter*0.2*l1para->resscale), nsgtrace, MIN(l1para->highfreq*0.5f/deltaf,hfid),
		   lfid, npx, sppow*0.5f, pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para); 
  
      c_l1inv_lsqr(MSDGINTERP_CG_WEIGHTED,l1para->tolerance, MAX(5,l1para->niter*0.4*l1para->resscale), nsgtrace, MIN(l1para->highfreq/deltaf,hfid),    
		   lfid, npx, sppow*0.5f, pweight, residual1_p, residual1_az, residual1_ay, taup, tweight, l1para);  

      if (l1para->outtype == VCROSSGHOST)
	{
	  datamode = 0; 
	  datum = 0;
	  c_l1inv_TAUP_RV_CTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_az, l1para, datamode, datum);
	  ARRAY_MATH(output, output, +, dataxt_az, nsgtrace, l1para->nsamp);
      
	  datamode = 0; 
	  datum = 0;
	  c_l1inv_TAUP_RV_MTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_az, l1para, datamode, datum);
	  ARRAY_MATH(output, output, -, dataxt_az, nsgtrace, l1para->nsamp);

	}
      else if (l1para->outtype == VZTOVCR)
	{
	  datamode = 1; 
	  datum = 0;
	  c_l1inv_TAUP_RV_CTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_az, l1para, datamode, datum);   
	  ARRAY_MATH(residual1_az, residual1_az, -, dataxt_az, nsgtrace, l1para->nsamp);
      
	  datamode = 0; 
	  datum = 0;
	  c_l1inv_TAUP_RV_CTOS(nsgtrace, hfid, lfid, npx, taup, dataxt_az, l1para, datamode, datum);   
	  ARRAY_MATH(output, output, +, dataxt_az, nsgtrace, l1para->nsamp);

	  ARRAY_MATH(output, l1para->residual*residual1_az, +, output, nsgtrace_az, l1para->nsamp);

	}
    } // this ends the check on resscale
  
  free(indx);

  flexible_free_array2d(dataxt_az);       //myfree (dataxt,nsgtrace);
  flexible_free_array2d(residual1_az);    //myfree (residual1,nsgtrace);
  
  flexible_free_array2d(taup);          //myfree (taup,   np_orig);
  flexible_free_array2d(taup1);        //myfree (taup,   np_orig);
  flexible_free_array2d(pweight);       //myfree (pweight,np_orig);
  flexible_free_array2d(tweight);       //myfree (tweight,np_orig);

  PFL_FREE(l1para->fftplan[0]); //kn
  PFL_FREE(l1para->fftplan[1]); //kn
}

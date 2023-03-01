#include <math.h>
#include <fftw3.h>
#include <time.h>
#include <pthread.h>
#include "string.h"
#include <stdlib.h>
#include "msdginterp.h"
#include "msdginterp.ci"
using namespace std;

void msdginterp_applydiverge(float * pRecv, master_t *msdginterp, int ntrace)
 {
   
   int isky,isamp,nsamp,trlen;
   float tscale;

   nsamp = msdginterp->nsamp;
   trlen = nsamp + NHEAD;

   for (isky=0;isky<ntrace;isky++)
     {
       pRecv[isky*trlen+NHEAD]=0.0f;
       for (isamp=1;isamp<nsamp;isamp++)
	 {
	   tscale = isamp*msdginterp->srate*1.0e-6;
	   tscale = powf(tscale,msdginterp->divalpha);
	   pRecv[isky*trlen+NHEAD+isamp] = pRecv[isky*trlen+NHEAD+isamp]*tscale; 
	 }
     }
 }

void msdginterp_removediverge(float * pRecv, master_t *msdginterp, int ntrace)
 {
   
   int isky,isamp,nsamp,trlen;
   float tscale;

   nsamp = msdginterp->nsamp;
   trlen = nsamp + NHEAD;

   for (isky=0;isky<ntrace;isky++)
     {
       pRecv[isky*trlen+NHEAD]=0.0f;
       for (isamp=1;isamp<nsamp;isamp++)
	 {
	   tscale = isamp*msdginterp->srate*1.0e-6;
	   tscale = powf(tscale,msdginterp->divalpha);
	   pRecv[isky*trlen+NHEAD+isamp] = pRecv[isky*trlen+NHEAD+isamp]/tscale; 
	 }
     }
 }

void
msdginterp_applyagc (float * pRecv, float *pScalar, master_t *msdginterp, int ntrace)
 {
   int isky,isamp,nsamp,trlen;
   float tscale;

   float *pInput;

   float *pScalar_agc;

   int gatelen;

   nsamp = msdginterp->nsamp;
   trlen = nsamp + NHEAD;

   pInput=(float *) calloc((size_t)(ntrace)*(size_t)nsamp,sizeof(float));
   pScalar_agc=(float *) calloc((size_t)(ntrace)*(size_t)nsamp,sizeof(float));

   for (isky=0;isky<ntrace;isky++)
     {
       for (isamp=0;isamp<nsamp;isamp++)
	 {
	   pInput[isky*nsamp+isamp] = pRecv[isky*trlen+NHEAD+isamp]; 
	 }
     }
   
   ////gatelen = MAX(3,(msdginterp->agcgatelen*1000/msdginterp->srate)+1);
   gatelen = MAX(3,(msdginterp->agcgatelen*1000/msdginterp->srate));

   agc_wrapper_(pInput,pScalar_agc,&ntrace,&nsamp,&gatelen);

   for (isky=0;isky<ntrace;isky++)
     {
       for (isamp=0;isamp<nsamp;isamp++)
	 {
	   pScalar[isky*trlen+NHEAD+isamp] = pScalar_agc[isky*nsamp+isamp]; 
	   
	   tscale = pScalar[isky*trlen+NHEAD+isamp];
	   pRecv[isky*trlen+NHEAD+isamp] = pRecv[isky*trlen+NHEAD+isamp]*tscale; 

	 }
     }

   free(pInput);   
   free(pScalar_agc);   

 }

void
msdginterp_removeagc (float * pRecv, float *pScalar, master_t *msdginterp, int ntrace)
 {
   int isky,isamp,nsamp,trlen;
   float tscale;

   nsamp = msdginterp->nsamp;
   trlen = nsamp + NHEAD;

   for (isky=0;isky<ntrace;isky++)
     {
       for (isamp=0;isamp<nsamp;isamp++)
	 {
	   tscale = pScalar[isky*trlen+NHEAD+isamp];
	   if (tscale!=0.0f) pRecv[isky*trlen+NHEAD+isamp] = pRecv[isky*trlen+NHEAD+isamp]/tscale; 
	 }
     }

 }

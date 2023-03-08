#include "msdginterp.h"
#include <complex>
#include <PFL_C.h>
#include "utility.h"
using namespace std;

float ** myallocatef(float **array, int dim1, int dim2)
{

	int i;

	array = (float **) malloc(dim1*sizeof(float *));

	for(i=0; i<dim1; i++) array[i] = (float *)calloc(dim2,sizeof(float));

	return array;

}

void myfree(float **array, int dim1)
{

	int i;

	for(i=0; i<dim1; i++) free(array[i]); 
	free(array); 

}

void myinit(float **array, int dim1, int dim2, int c) 
{

	int i,j;

	for(i=0; i<dim1; i++) 
	{
		for(j=0; j<dim2; j++) 
		{
			array[i][j] = c;
		}
	}

}

void mycopy(float **array1, float **array2, int dim1, int dim2)
{
	int i;
	for(i=0; i<dim1; i++) memcpy(array2[i],array1[i],dim2*sizeof(float));
}

void mymult(float **array1, float **array2, float **array3, int dim1, int dim2)
{
	int i,j;
	for(i=0; i<dim1; i++) 
	{
		for(j=0; j<dim2; j++) 
		{
			array3[i][j] = array1[i][j]*array2[i][j];
		}
	}
}

void allocate_deghost3d(l1inv_t *l1para, int ntrcxy) 
{

	int fftnr;
	int nsgtrace = ntrcxy;
	int nsgtracex = l1para->nsgtracexlastgate;
	int nsgtracey = l1para->nsgtraceyloc;
	int nsgtracepad   = ceil((float)(nsgtrace  /4.))*4;
	int npx, npy, nq, npxpad,npypad;

	if(l1para->clustertype == 0)
	{
		npx = MAX(ceil(((l1para->psample*l1para->chanitvx)*(float)(nsgtracex))*
					((l1para->pxmax-l1para->pxmin)*l1para->highfreq)), nsgtracex);
		l1para->npxorig = npx;

		npxpad = ceil((float)(npx/4.))*4;

		npy = MAX(ceil(((l1para->psample*l1para->chanitvy)*(float)(nsgtracey))*
					((l1para->pymax-l1para->pymin)*l1para->highfreq)), nsgtracey);
		l1para->npyorig = npy;

		npypad = ceil((float)(npy/4.))*4;

		//////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////
		npy = npy*l1para->nqorig;
		npypad = ceil((float)(npy/4.))*4;
		//////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////

	}

	if ((l1para->choosemethod == CHOOSEMETHOD_JOINTSR3D && 
				l1para->jointmethod != JOINTDEGHOST)
			||(l1para->choosemethod == CHOOSEMETHOD_MSDGI &&
				l1para->jointmc == YESYES &&
				l1para->jointmethod != JOINTDEGHOST))
		PFL_LENGTH(l1para->nsamporig+ l1para->nsampsrc,&fftnr);
	else
		PFL_LENGTH(l1para->nsamporig,&fftnr);

	fftnr+=2;

	//add by kyang try speed up
	l1para->orgnr =        fftnr-2;
	l1para->orgnc =        l1para->orgnr/2+1;
	l1para->hilbert_nt =   61;

	l1para->zwtilt = (float *) malloc(fftnr*sizeof(float));
	l1para->ywtilt = (float *) malloc(fftnr*sizeof(float));
	l1para->zlowf = (float *) malloc(fftnr*sizeof(float));
	l1para->ylowf = (float *) malloc(fftnr*sizeof(float));

	if(l1para->clustertype == 0)
	{
		l1para->ztop =  (float *) malloc(npx*npy*sizeof(float));
		l1para->ytop =  (float *) malloc(npx*npy*sizeof(float));
		l1para->freqfac = (float *) malloc(fftnr*sizeof(float));
		l1para->obliqcorr = (float *) malloc(MAX(nsgtrace,npx*npy)*sizeof(float));

		l1para->wtilt = (float *) malloc(fftnr*sizeof(float));  //// ** DO WE NEED THIS ?? **

		l1para->pxsave = (float *) malloc(npx*npy*sizeof(float));
		l1para->pysave = (float *) malloc(npx*npy*sizeof(float));
		l1para->rqsave = (float *) malloc(npx*npy*sizeof(float));

		l1para->AP_p =           (float**)flexible_array2d(npx*npy, nsgtracepad, 2*sizeof(float));
		l1para->AG_p =           (float**)flexible_array2d(npx*npy, nsgtracepad, 2*sizeof(float));
		l1para->AS_p =           (float**)flexible_array2d(npx*npy, nsgtracepad, 2*sizeof(float));
		l1para->TP_p =           (float**)flexible_array2d(nsgtrace, npxpad*npypad, 2*sizeof(float));
		l1para->TG_p =           (float**)flexible_array2d(nsgtrace, npxpad*npypad, 2*sizeof(float));
		l1para->TS_p =           (float**)flexible_array2d(nsgtrace, npxpad*npypad, 2*sizeof(float));

		l1para->fftNR =        ALIGN_FLOAT(l1para->orgnr+2);     //for compatable with fftnr
		l1para->fftNC =        ALIGN_FCPLX(l1para->orgnc+1);     //for compatable with fftnc

		l1para->pfp  =         (float**)flexible_array2d(npx*npy, l1para->fftNC, 2*sizeof(float));
		if(l1para->regparam > 0.0f)
			l1para->pfpinit  =     (float**)flexible_array2d(npx*npy, l1para->fftNC, 2*sizeof(float)); // LKL
		else
			l1para->pfpinit  = NULL;

		l1para->pfpin    =     (float**)flexible_array2d(npx*npy, l1para->fftNC, 2*sizeof(float));
		l1para->pfpadd   =     (float**)flexible_array2d(npx*npy, l1para->fftNC, 2*sizeof(float));
		l1para->pfpT =         (float**)flexible_array2d(l1para->fftNC, npx*npy, 2*sizeof(float));
		l1para->pfx  =         (float**)flexible_array2d(nsgtrace, l1para->fftNC, 2*sizeof(float));
		l1para->pfxT =         (float**)flexible_array2d(l1para->fftNC, nsgtrace, 2*sizeof(float));
                
		l1para->integration = (float *) malloc(l1para->fftNC*2*sizeof(float));
 
                posix_memalign((void **) &l1para->xGhost, 32,MAX(nsgtrace, npx*npy)*2*sizeof(float) );
                posix_memalign((void **) &l1para->xPrimary, 32,MAX(nsgtrace, npx*npy)*2*sizeof(float) );

		//l1para->xGhost =       (float*)malloc(MAX(nsgtrace, npx*npy)*2*sizeof(float));
		//l1para->xPrimary =     (float*)malloc(MAX(nsgtrace, npx*npy)*2*sizeof(float));

		if (l1para->flagl1tgate == 0)
		{
			l1para->hilbert_hw =   (float*)malloc((l1para->hilbert_nt*2+1)*sizeof(float));
			l1para->hilbert_conv = (float*)malloc((l1para->hilbert_nt*2+l1para->orgnr)*sizeof(float));
			l1para->hilbert_task = hilbert_init(l1para->hilbert_nt, l1para->orgnr, l1para->hilbert_hw);
		}
		//printf("l1para->fftNR=%d fftnr=%d, orgnr=%d\n", l1para->fftNR, fftnr, l1para->orgnr); fflush(0); 

		assert(l1para->fftNR>=fftnr);

		l1para->cg_d  = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));
		l1para->cg_r  = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));
		l1para->cg_wd = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));
		l1para->cg_Ad = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));

		if (l1para->cgmethod==1)
		{
			l1para->cg_s     = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));
			l1para->cg_As    = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));
			l1para->cg_rhat0 = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));
		}

		if(l1para->fpredatum==FPREDATUM_DGRG)
		{
			l1para->TPR_p =           (float**)flexible_array2d(nsgtrace, npxpad*npypad, 2*sizeof(float));
			l1para->TGR_p =           (float**)flexible_array2d(nsgtrace, npxpad*npypad, 2*sizeof(float));
		}
		if(l1para->choosemethod == CHOOSEMETHOD_JOINTSR3D || (l1para->choosemethod == CHOOSEMETHOD_MSDGI && 
					l1para->jointmc == YESYES))
		{
			l1para->xGhostS = (float*)malloc(MAX(nsgtrace, npx*npy)*2*sizeof(float));
			l1para->xGhostR = (float*)malloc(MAX(nsgtrace, npx*npy)*2*sizeof(float));
			if (l1para->computemode == 1)
				l1para->pSrcFilter = (float**)flexible_array2d(npx*npy, l1para->fftNC, 2*sizeof(float));

			if(l1para->jointmethod != JOINTDEGHOST)
			{
				l1para->pSrcFilterAz = (float**)flexible_array2d(npx*npy, l1para->fftNC, 2*sizeof(float));
				l1para->pTgtFilter = (float**)flexible_array2d(npx*npy, l1para->fftNC, 2*sizeof(float));
				l1para->p1dSrc = (float**)flexible_array2d(1, l1para->fftNC, 2*sizeof(float));
				l1para->p1dTgt = (float**)flexible_array2d(1, l1para->fftNC, 2*sizeof(float));
			}
		}
	}
	//end
}

void deallocate_deghost3d(l1inv_t *l1para) 
{

	free(l1para->zwtilt);
	free(l1para->ywtilt);
	free(l1para->zlowf);
	free(l1para->ylowf);

	if(l1para->clustertype == 0)
	{
		free(l1para->ztop);
		free(l1para->ytop);
		free(l1para->freqfac);
		free(l1para->obliqcorr);
		free(l1para->integration);

		free(l1para->wtilt);     //// ** DO WE NEED THIS ?? **

		free(l1para->pxsave);
		free(l1para->pysave);
		free(l1para->rqsave);

		flexible_free_array2d(l1para->AP_p);
		flexible_free_array2d(l1para->AG_p);
		flexible_free_array2d(l1para->AS_p);
		flexible_free_array2d(l1para->TP_p);
		flexible_free_array2d(l1para->TG_p);
		flexible_free_array2d(l1para->TS_p);

		flexible_free_array2d(l1para->pfp);
		if(l1para->regparam > 0.0f)flexible_free_array2d(l1para->pfpinit); // LKL
		flexible_free_array2d(l1para->pfpin);
		flexible_free_array2d(l1para->pfpadd);
		flexible_free_array2d(l1para->pfx);
		flexible_free_array2d(l1para->pfxT);
		flexible_free_array2d(l1para->pfpT);

		free(l1para->xGhost);
		free(l1para->xPrimary);

		if (l1para->flagl1tgate == 0)
		{
			free(l1para->hilbert_hw);
			free(l1para->hilbert_conv);
			hilbert_destroy(l1para->hilbert_task);
		}

		flexible_free_array2d(l1para->cg_d);
		flexible_free_array2d(l1para->cg_r);
		flexible_free_array2d(l1para->cg_wd);
		flexible_free_array2d(l1para->cg_Ad);

		if (l1para->cgmethod == 1)
		{
			flexible_free_array2d(l1para->cg_s); 
			flexible_free_array2d(l1para->cg_As);  
			flexible_free_array2d(l1para->cg_rhat0);  
		}

		if(l1para->fpredatum==FPREDATUM_DGRG)
		{
			flexible_free_array2d(l1para->TPR_p);
			flexible_free_array2d(l1para->TGR_p);
		}
		if(l1para->choosemethod == CHOOSEMETHOD_JOINTSR3D || (l1para->choosemethod == CHOOSEMETHOD_MSDGI && 
					l1para->jointmc == YESYES))
		{
			free(l1para->xGhostS);
			free(l1para->xGhostR);
			if (l1para->computemode == 1)
				flexible_free_array2d(l1para->pSrcFilter);

			if(l1para->jointmethod != JOINTDEGHOST)
			{
				flexible_free_array2d(l1para->pSrcFilterAz);
				flexible_free_array2d(l1para->pTgtFilter);
				flexible_free_array2d(l1para->p1dSrc);
				flexible_free_array2d(l1para->p1dTgt);
			}
		}
	}
}

void allocate_interp3d(l1inv_t *l1para, int ntrcxy, int ntrcxyout) 
{

	int fftnr;
	int nsgtrace = ntrcxy;
	int nsgtraceout = ntrcxyout;
	int nsgtracex = l1para->nsgtracexlastgate;
	int nsgtracey = l1para->nsgtraceyloc;
	int nsgtracepad   = ceil((float)(nsgtrace  /4.))*4;
	int npx,npxpad,npy,npypad;

	if(l1para->clustertype == 0)
	{
		npx = MAX(ceil(((l1para->psample*l1para->chanitvx)*(float)(nsgtracex))*
					((l1para->pxmax-l1para->pxmin)*l1para->highfreq)), nsgtracex);
		l1para->npxorig = npx;

		npxpad = ceil((float)(npx/4.))*4;

		npy = MAX(ceil(((l1para->psample*l1para->chanitvy)*(float)(nsgtracey))*
					((l1para->pymax-l1para->pymin)*l1para->highfreq)), nsgtracey);
		l1para->npyorig = npy;

		npypad = ceil((float)(npy/4.))*4;
	}

	PFL_LENGTH(l1para->nsamporig,&fftnr);
	fftnr+=2;

	//add by kyang try speed up
	l1para->orgnr =        fftnr-2;
	l1para->orgnc =        l1para->orgnr/2+1;
	l1para->hilbert_nt =   61;
	l1para->zwtilt = (float *) malloc(fftnr*sizeof(float));
	l1para->ywtilt = (float *) malloc(fftnr*sizeof(float));
	l1para->zlowf = (float *) malloc(fftnr*sizeof(float));
	l1para->ylowf = (float *) malloc(fftnr*sizeof(float));

	if(l1para->clustertype == 0)
	{
		l1para->ztop =  (float *) malloc(npx*npy*sizeof(float));
		l1para->ytop =  (float *) malloc(npx*npy*sizeof(float));

		l1para->freqfac = (float *) malloc(fftnr*sizeof(float));
		l1para->obliqcorr = (float *) malloc(MAX(nsgtrace,npx*npy)*sizeof(float));

		l1para->wtilt = (float *) malloc(fftnr*sizeof(float));  //// ** DO WE NEED THIS ?? **

		l1para->pxsave = (float *) malloc(npx*npy*sizeof(float));
		l1para->pysave = (float *) malloc(npx*npy*sizeof(float));
		l1para->rqsave = (float *) malloc(npx*npy*sizeof(float));      

		l1para->AP_p =           (float**)flexible_array2d(npx*npy, nsgtracepad, 2*sizeof(float));
		l1para->AG_p =           (float**)flexible_array2d(npx*npy, nsgtracepad, 2*sizeof(float));
		l1para->TP_p =           (float**)flexible_array2d(nsgtrace, npxpad*npypad, 2*sizeof(float));
		l1para->TG_p =           (float**)flexible_array2d(nsgtrace, npxpad*npypad, 2*sizeof(float));
		l1para->TS_p =           (float**)flexible_array2d(nsgtrace, npxpad*npypad, 2*sizeof(float));

		l1para->TP_interp3d =           (float**)flexible_array2d(nsgtraceout, npxpad*npypad, 2*sizeof(float));
		l1para->TS_interp3d =           (float**)flexible_array2d(nsgtraceout, npxpad*npypad, 2*sizeof(float));

		l1para->fftNR =        ALIGN_FLOAT(l1para->orgnr+2);     //for compatable with fftnr
		l1para->fftNC =        ALIGN_FCPLX(l1para->orgnc+1);     //for compatable with fftnc

		l1para->pfp  =         (float**)flexible_array2d(npx*npy, l1para->fftNC, 2*sizeof(float));
		l1para->pfpin  =         (float**)flexible_array2d(npx*npy, l1para->fftNC, 2*sizeof(float));
		l1para->pfpadd =         (float**)flexible_array2d(npx*npy, l1para->fftNC, 2*sizeof(float));
		l1para->pfpT =         (float**)flexible_array2d(l1para->fftNC, npx*npy, 2*sizeof(float));
		l1para->pfx  =         (float**)flexible_array2d(nsgtraceout, l1para->fftNC, 2*sizeof(float));
		l1para->pfxT =         (float**)flexible_array2d(l1para->fftNC, nsgtraceout, 2*sizeof(float));

		l1para->integration = (float *) malloc(l1para->fftNC*2*sizeof(float));

		l1para->xGhost =       (float*)malloc(MAX(nsgtraceout, npx*npy)*2*sizeof(float));
		l1para->xPrimary =     (float*)malloc(MAX(nsgtraceout, npx*npy)*2*sizeof(float));

		if (l1para->flagl1tgate == 0)
		{
			l1para->hilbert_nt =   61;
			l1para->hilbert_hw =   (float*)malloc((l1para->hilbert_nt*2+1)*sizeof(float));
			l1para->hilbert_conv = (float*)malloc((l1para->hilbert_nt*2+l1para->orgnr)*sizeof(float));
			l1para->hilbert_task = hilbert_init(l1para->hilbert_nt, l1para->orgnr, l1para->hilbert_hw);
		}

		//printf("l1para->fftNR=%d fftnr=%d, orgnr=%d\n", l1para->fftNR, fftnr, l1para->orgnr); fflush(0); 

		assert(l1para->fftNR>=fftnr);

		l1para->cg_d  = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));
		l1para->cg_r  = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));
		l1para->cg_wd = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));
		l1para->cg_Ad = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));

		if (l1para->cgmethod==1)
		{
			l1para->cg_s     = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));
			l1para->cg_As    = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));
			l1para->cg_rhat0 = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));
		}
	}

	//end
}

void deallocate_interp3d(l1inv_t *l1para) 
{

	free(l1para->zwtilt);
	free(l1para->ywtilt);
	free(l1para->zlowf);
	free(l1para->ylowf);

	if(l1para->clustertype == 0)
	{
		free(l1para->ztop);
		free(l1para->ytop);
		free(l1para->freqfac);
		free(l1para->obliqcorr);
		free(l1para->integration);

		free(l1para->wtilt);     //// ** DO WE NEED THIS ?? **

		free(l1para->pxsave);
		free(l1para->pysave);
		free(l1para->rqsave);

		flexible_free_array2d(l1para->AP_p);
		flexible_free_array2d(l1para->AG_p);
		flexible_free_array2d(l1para->TP_p);
		flexible_free_array2d(l1para->TG_p);
		flexible_free_array2d(l1para->TS_p);

		flexible_free_array2d(l1para->TP_interp3d);
		flexible_free_array2d(l1para->TS_interp3d);

		flexible_free_array2d(l1para->pfp);
		flexible_free_array2d(l1para->pfpin);
		flexible_free_array2d(l1para->pfpadd);
		flexible_free_array2d(l1para->pfx);
		flexible_free_array2d(l1para->pfxT);
		flexible_free_array2d(l1para->pfpT);

		free(l1para->xGhost);
		free(l1para->xPrimary);

		if (l1para->flagl1tgate == 0)
		{
			free(l1para->hilbert_hw);
			free(l1para->hilbert_conv);
			hilbert_destroy(l1para->hilbert_task);
		}

		flexible_free_array2d(l1para->cg_d);
		flexible_free_array2d(l1para->cg_r);
		flexible_free_array2d(l1para->cg_wd);
		flexible_free_array2d(l1para->cg_Ad);

		if (l1para->cgmethod == 1)
		{
			flexible_free_array2d(l1para->cg_s); 
			flexible_free_array2d(l1para->cg_As);  
			flexible_free_array2d(l1para->cg_rhat0);  
		}
	}
}

void allocate_fp4d(l1inv_t *l1para, int ntrcxy) 
{

	int fftnr, fftnc, idata;
	int nsgtrace = ntrcxy;
	int nsgtracex = l1para->nsgtracexlastgate;
	int nsgtracey = l1para->nsgtraceyloc;
	int nsgtracepad   = ceil((float)(nsgtrace  /4.))*4;
	int npx, npy,npxpad,npypad;

	if(l1para->clustertype == 0)
	{
		npx = MAX(ceil(((l1para->psample*l1para->chanitvx)*(float)(nsgtracex))*
					((l1para->pxmax-l1para->pxmin)*l1para->highfreq)), nsgtracex);
		l1para->npxorig = npx;

		npxpad = ceil((float)(npx/4.))*4;

		npy = MAX(ceil(((l1para->psample*l1para->chanitvy)*(float)(nsgtracey))*
					((l1para->pymax-l1para->pymin)*l1para->highfreq)), nsgtracey);
		l1para->npyorig = npy;

		npypad = ceil((float)(npy/4.))*4;
	}

	if (l1para->jointmc == YESYES)
		PFL_LENGTH(l1para->nsamporig+ l1para->nsampsrc,&fftnr);
	else
		PFL_LENGTH(l1para->nsamporig,&fftnr);

	fftnr+=2;
	fftnc = fftnr/2;
	l1para->fftnr = fftnr; 
	l1para->fftnc = fftnc; 

	//add by kyang try speed up
	l1para->orgnr =        fftnr-2;
	l1para->orgnc =        l1para->orgnr/2+1;
	l1para->hilbert_nt =   61;

	l1para->zwtilt = (float *) malloc(fftnr*sizeof(float));
	l1para->ywtilt = (float *) malloc(fftnr*sizeof(float));
	l1para->zlowf = (float *) malloc(fftnr*sizeof(float));
	l1para->ylowf = (float *) malloc(fftnr*sizeof(float));

	if(l1para->clustertype == 0)
	{
		l1para->ztop =  (float *) malloc(npx*npy*sizeof(float));
		l1para->ytop =  (float *) malloc(npx*npy*sizeof(float));
		l1para->freqfac = (float *) malloc(fftnr*sizeof(float));
		l1para->obliqcorr = (float *) malloc(MAX(nsgtrace,npx*npy)*sizeof(float));

		l1para->wtilt = (float *) malloc(fftnr*sizeof(float));  //// ** DO WE NEED THIS ?? **

		l1para->pxsave = (float *) malloc(npx*npy*sizeof(float));
		l1para->pysave = (float *) malloc(npx*npy*sizeof(float));
		l1para->rqsave = (float *) malloc(npx*npy*sizeof(float));

		l1para->AP_parray =           (float***)malloc(l1para->ndata*sizeof(float **));
		l1para->AG_parray =           (float***)malloc(l1para->ndata*sizeof(float **));
		l1para->AS_parray =           (float***)malloc(l1para->ndata*sizeof(float **));
		l1para->TP_parray =           (float***)malloc(l1para->ndata*sizeof(float **));
		l1para->TG_parray =           (float***)malloc(l1para->ndata*sizeof(float **));
		l1para->TS_parray =           (float***)malloc(l1para->ndata*sizeof(float **));

		for (idata=0;idata<l1para->ndata;idata++)
		{
			l1para->AP_parray[idata] = (float**)flexible_array2d(npx*npy, nsgtracepad, 2*sizeof(float));
			l1para->AG_parray[idata] = (float**)flexible_array2d(npx*npy, nsgtracepad, 2*sizeof(float));
			l1para->AS_parray[idata] = (float**)flexible_array2d(npx*npy, nsgtracepad, 2*sizeof(float));
			l1para->TP_parray[idata] = (float**)flexible_array2d(nsgtrace, npxpad*npypad, 2*sizeof(float));
			l1para->TG_parray[idata] = (float**)flexible_array2d(nsgtrace, npxpad*npypad, 2*sizeof(float));
			l1para->TS_parray[idata] = (float**)flexible_array2d(nsgtrace, npxpad*npypad, 2*sizeof(float));
		}

		l1para->fftNR =        ALIGN_FLOAT(l1para->orgnr+2);     //for compatable with fftnr
		l1para->fftNC =        ALIGN_FCPLX(l1para->orgnc+1);     //for compatable with fftnc

		l1para->pfp  =         (float**)flexible_array2d(npx*npy, l1para->fftNC, 2*sizeof(float));
		l1para->pfpin  =         (float**)flexible_array2d(npx*npy, l1para->fftNC, 2*sizeof(float));
		l1para->pfpinit=         (float**)flexible_array2d(npx*npy, l1para->fftNC, 2*sizeof(float));
		l1para->pfpadd =         (float**)flexible_array2d(npx*npy, l1para->fftNC, 2*sizeof(float));
		l1para->pfpT =         (float**)flexible_array2d(l1para->fftNC, npx*npy, 2*sizeof(float));
		l1para->pfx  =         (float**)flexible_array2d(nsgtrace, l1para->fftNC, 2*sizeof(float));
		l1para->pfxT =         (float**)flexible_array2d(l1para->fftNC, nsgtrace, 2*sizeof(float));

		l1para->integration = (float *) malloc(l1para->fftNC*2*sizeof(float));

		l1para->xGhost =       (float*)malloc(MAX(nsgtrace, npx*npy)*2*sizeof(float));
		l1para->xPrimary =     (float*)malloc(MAX(nsgtrace, npx*npy)*2*sizeof(float));


		l1para->hilbert_hw =   (float*)malloc((l1para->hilbert_nt*2+1)*sizeof(float));
		l1para->hilbert_conv = (float*)malloc((l1para->hilbert_nt*2+l1para->orgnr)*sizeof(float));
		l1para->hilbert_task = hilbert_init(l1para->hilbert_nt, l1para->orgnr, l1para->hilbert_hw);

		//printf("l1para->fftNR=%d fftnr=%d, orgnr=%d\n", l1para->fftNR, fftnr, l1para->orgnr); fflush(0); 

		assert(l1para->fftNR>=fftnr);

		l1para->cg_d  = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));
		l1para->cg_r  = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));
		l1para->cg_wd = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));
		l1para->cg_Ad = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));

		if (l1para->cgmethod==1)
		{
			l1para->cg_s     = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));
			l1para->cg_As    = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));
			l1para->cg_rhat0 = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));
		}

		if(l1para->fpredatum==FPREDATUM_DGRG)
		{
			l1para->TPR_p =           (float**)flexible_array2d(nsgtrace, npxpad*npypad, 2*sizeof(float));
			l1para->TGR_p =           (float**)flexible_array2d(nsgtrace, npxpad*npypad, 2*sizeof(float));
		}

		if(l1para->jointmc == YESYES)
		{

			l1para->AGS_parray =           (float***)malloc(l1para->ndata*sizeof(float **));
			l1para->AGR_parray =           (float***)malloc(l1para->ndata*sizeof(float **));
			l1para->TGS_parray =           (float***)malloc(l1para->ndata*sizeof(float **));
			l1para->TGR_parray =           (float***)malloc(l1para->ndata*sizeof(float **));

			for (idata=0;idata<l1para->ndata;idata++)
			{
				l1para->AGS_parray[idata] = (float**)flexible_array2d(npx*npy, nsgtracepad, 2*sizeof(float));
				l1para->AGR_parray[idata] = (float**)flexible_array2d(npx*npy, nsgtracepad, 2*sizeof(float));
				l1para->TGS_parray[idata] = (float**)flexible_array2d(nsgtrace, npxpad*npypad, 2*sizeof(float));
				l1para->TGR_parray[idata] = (float**)flexible_array2d(nsgtrace, npxpad*npypad, 2*sizeof(float));
			}

			l1para->xGhostS = (float*)malloc(MAX(nsgtrace, npx*npy)*2*sizeof(float));
			l1para->xGhostR = (float*)malloc(MAX(nsgtrace, npx*npy)*2*sizeof(float));
			if (l1para->computemode == 1)
			{
				l1para->pSrcFilter_array = (float***)malloc(l1para->ndata*sizeof(float **));
				for (idata=0;idata<l1para->ndata;idata++)
					l1para->pSrcFilter_array[idata] = (float**)flexible_array2d(npx*npy, l1para->fftNC, 2*sizeof(float));
			}

			if(l1para->jointmethod != JOINTDEGHOST)
			{
				l1para->pSrcFilterAz = (float**)flexible_array2d(npx*npy, l1para->fftNC, 2*sizeof(float));
				l1para->pTgtFilter = (float**)flexible_array2d(npx*npy, l1para->fftNC, 2*sizeof(float));
				l1para->p1dSrc = (float**)flexible_array2d(1, l1para->fftNC, 2*sizeof(float));
				l1para->p1dTgt = (float**)flexible_array2d(1, l1para->fftNC, 2*sizeof(float));
			}
		}
	}
	//end
}

void deallocate_fp4d(l1inv_t *l1para) 
{

	int idata;

	free(l1para->zwtilt);
	free(l1para->ywtilt);
	free(l1para->zlowf);
	free(l1para->ylowf);

	if(l1para->clustertype == 0)
	{
		free(l1para->ztop);
		free(l1para->ytop);
		free(l1para->freqfac);
		free(l1para->obliqcorr);
		free(l1para->integration);

		free(l1para->wtilt);     //// ** DO WE NEED THIS ?? **

		free(l1para->pxsave);
		free(l1para->pysave);
		free(l1para->rqsave);

		for (idata=0;idata<l1para->ndata;idata++)
		{
			flexible_free_array2d(l1para->AP_parray[idata]);
			flexible_free_array2d(l1para->AG_parray[idata]);
			flexible_free_array2d(l1para->AS_parray[idata]);
			flexible_free_array2d(l1para->TP_parray[idata]);
			flexible_free_array2d(l1para->TG_parray[idata]);
			flexible_free_array2d(l1para->TS_parray[idata]);
		}


		free(l1para->AP_parray);
		free(l1para->AG_parray);
		free(l1para->AS_parray);
		free(l1para->TP_parray);
		free(l1para->TG_parray);
		free(l1para->TS_parray);

		flexible_free_array2d(l1para->pfp);
		flexible_free_array2d(l1para->pfpin);
		flexible_free_array2d(l1para->pfpinit);
		flexible_free_array2d(l1para->pfpadd);
		flexible_free_array2d(l1para->pfx);
		flexible_free_array2d(l1para->pfxT);
		flexible_free_array2d(l1para->pfpT);

		free(l1para->xGhost);
		free(l1para->xPrimary);

		free(l1para->hilbert_hw);
		free(l1para->hilbert_conv);

		hilbert_destroy(l1para->hilbert_task);

		flexible_free_array2d(l1para->cg_d);
		flexible_free_array2d(l1para->cg_r);
		flexible_free_array2d(l1para->cg_wd);
		flexible_free_array2d(l1para->cg_Ad);

		if (l1para->cgmethod == 1)
		{
			flexible_free_array2d(l1para->cg_s); 
			flexible_free_array2d(l1para->cg_As);  
			flexible_free_array2d(l1para->cg_rhat0);  
		}

		if(l1para->fpredatum==FPREDATUM_DGRG)
		{
			flexible_free_array2d(l1para->TPR_p);
			flexible_free_array2d(l1para->TGR_p);
		}

		if(l1para->jointmc == YESYES)
		{
			free(l1para->xGhostS);
			free(l1para->xGhostR);

			if (l1para->computemode == 1)
			{
				for (idata=0;idata<l1para->ndata;idata++) 
					flexible_free_array2d(l1para->pSrcFilter_array[idata]);
				free(l1para->pSrcFilter_array);
			}

			if(l1para->jointmethod != JOINTDEGHOST)
			{
				flexible_free_array2d(l1para->pSrcFilterAz);
				flexible_free_array2d(l1para->pTgtFilter);
				flexible_free_array2d(l1para->p1dSrc);
				flexible_free_array2d(l1para->p1dTgt);
			}

			for (idata=0;idata<l1para->ndata;idata++) 
			{
				flexible_free_array2d(l1para->AGS_parray[idata]);
				flexible_free_array2d(l1para->AGR_parray[idata]);
				flexible_free_array2d(l1para->TGS_parray[idata]);
				flexible_free_array2d(l1para->TGR_parray[idata]);
			}

			free(l1para->AGS_parray);
			free(l1para->AGR_parray);
			free(l1para->TGS_parray);
			free(l1para->TGR_parray);

		}
	}
}


void allocate_dsdg3d(l1inv_t *l1para, int ntrcxy) 
{

	int fftnr;
	int nsgtrace = ntrcxy;
	int nsgtracex = l1para->nsgtracexlastgate;
	int nsgtracey = l1para->nsgtraceyloc;
	int nsgtracepad   = ceil((float)(nsgtrace  /4.))*4;

	int npx = MAX(ceil(((l1para->psample*l1para->chanitvx)*(float)(nsgtracex))*
				((l1para->pxmax-l1para->pxmin)*l1para->highfreq)), nsgtracex);
	l1para->npxorig = npx;

	int npxpad = ceil((float)(npx/4.))*4;

	int npy = MAX(ceil(((l1para->psample*l1para->chanitvy)*(float)(nsgtracey))*
				((l1para->pymax-l1para->pymin)*l1para->highfreq)), nsgtracey);
	l1para->npyorig = npy;

	int npypad = ceil((float)(npy/4.))*4;

	PFL_LENGTH(l1para->nsamporig,&fftnr);
	fftnr+=2;

	//add by kyang try speed up
	l1para->orgnr =        fftnr-2;
	l1para->orgnc =        l1para->orgnr/2+1;
	l1para->hilbert_nt =   61;

	l1para->zwtilt = (float *) malloc(fftnr*sizeof(float));
	l1para->ywtilt = (float *) malloc(fftnr*sizeof(float));
	l1para->zlowf = (float *) malloc(fftnr*sizeof(float));
	l1para->ylowf = (float *) malloc(fftnr*sizeof(float));

	if(l1para->clustertype == 0)
	{
		l1para->ztop =  (float *) malloc(npx*npy*sizeof(float));
		l1para->ytop =  (float *) malloc(npx*npy*sizeof(float));
		l1para->freqfac = (float *) malloc(fftnr*sizeof(float));
		l1para->obliqcorr  = (float *) malloc(MAX(nsgtrace,npx*npy)*sizeof(float));

		l1para->wtilt = (float *) malloc(fftnr*sizeof(float));  //// ** DO WE NEED THIS ?? **

		l1para->pxsave = (float *) malloc(npx*npy*sizeof(float));
		l1para->pysave = (float *) malloc(npx*npy*sizeof(float));
		l1para->rqsave = (float *) malloc(npx*npy*sizeof(float));    // zxue bug fix  

		l1para->AP_p =           (float**)flexible_array2d(npx*npy, nsgtracepad, 2*sizeof(float));
		l1para->AG_p =           (float**)flexible_array2d(npx*npy, nsgtracepad, 2*sizeof(float));
		l1para->AS_p =           (float**)flexible_array2d(npx*npy, nsgtracepad, 2*sizeof(float));
		l1para->TP_p =           (float**)flexible_array2d(nsgtrace, npxpad*npypad, 2*sizeof(float));
		l1para->TG_p =           (float**)flexible_array2d(nsgtrace, npxpad*npypad, 2*sizeof(float));
		l1para->TS_p =           (float**)flexible_array2d(nsgtrace, npxpad*npypad, 2*sizeof(float));

		l1para->fftNR =        ALIGN_FLOAT(l1para->orgnr+2);     //for compatable with fftnr
		l1para->fftNC =        ALIGN_FCPLX(l1para->orgnc+1);     //for compatable with fftnc

		l1para->pfp  =         (float**)flexible_array2d(npx*npy, l1para->fftNC, 2*sizeof(float));
		l1para->pfpin  =       (float**)flexible_array2d(npx*npy, l1para->fftNC, 2*sizeof(float));
		l1para->pfpadd =       (float**)flexible_array2d(npx*npy, l1para->fftNC, 2*sizeof(float));
		l1para->pfpT =         (float**)flexible_array2d(l1para->fftNC, npx*npy, 2*sizeof(float));
		l1para->pfpT2 =        (float**)flexible_array2d(l1para->fftNC, npx*npy, 2*sizeof(float));
		l1para->pfx  =         (float**)flexible_array2d(nsgtrace, l1para->fftNC, 2*sizeof(float));
		l1para->pfxT =         (float**)flexible_array2d(l1para->fftNC, nsgtrace, 2*sizeof(float));

		l1para->integration = (float *) malloc(l1para->fftNC*2*sizeof(float));

		l1para->xGhost =       (float*)malloc(MAX(nsgtrace, npx*npy)*2*sizeof(float));
		l1para->xPrimary =     (float*)malloc(MAX(nsgtrace, npx*npy)*2*sizeof(float));

		if (l1para->flagl1tgate == 0)
		{
			l1para->hilbert_hw =   (float*)malloc((l1para->hilbert_nt*2+1)*sizeof(float));
			l1para->hilbert_conv = (float*)malloc((l1para->hilbert_nt*2+l1para->orgnr)*sizeof(float));
			l1para->hilbert_task = hilbert_init(l1para->hilbert_nt, l1para->orgnr, l1para->hilbert_hw);
		}
		//printf("l1para->fftNR=%d fftnr=%d, orgnr=%d\n", l1para->fftNR, fftnr, l1para->orgnr); fflush(0); 

		assert(l1para->fftNR>=fftnr);

		l1para->cg_d  = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));
		l1para->cg_r  = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));
		l1para->cg_wd = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));
		l1para->cg_Ad = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));

		if (l1para->cgmethod==1)
		{
			l1para->cg_s     = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));
			l1para->cg_As    = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));
			l1para->cg_rhat0 = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));
		}

		if(l1para->fpredatum==FPREDATUM_DGRG)
		{
			l1para->TPR_p =           (float**)flexible_array2d(nsgtrace, npxpad*npypad, 2*sizeof(float));
			l1para->TGR_p =           (float**)flexible_array2d(nsgtrace, npxpad*npypad, 2*sizeof(float));
		}
	}
	//end
}

void deallocate_dsdg3d(l1inv_t *l1para) 
{

	free(l1para->zwtilt);
	free(l1para->ywtilt);
	free(l1para->zlowf);
	free(l1para->ylowf);

	if(l1para->clustertype == 0)
	{
		free(l1para->ztop);
		free(l1para->ytop);
		free(l1para->freqfac);
		free(l1para->obliqcorr);
		free(l1para->integration);

		free(l1para->wtilt);     //// ** DO WE NEED THIS ?? **

		free(l1para->pxsave);
		free(l1para->pysave);
		free(l1para->rqsave);

		flexible_free_array2d(l1para->AP_p);
		flexible_free_array2d(l1para->AG_p);
		flexible_free_array2d(l1para->AS_p);
		flexible_free_array2d(l1para->TP_p);
		flexible_free_array2d(l1para->TG_p);
		flexible_free_array2d(l1para->TS_p);

		flexible_free_array2d(l1para->pfp);
		flexible_free_array2d(l1para->pfpin);
		flexible_free_array2d(l1para->pfpadd);
		flexible_free_array2d(l1para->pfx);
		flexible_free_array2d(l1para->pfxT);
		flexible_free_array2d(l1para->pfpT);
		flexible_free_array2d(l1para->pfpT2);

		free(l1para->xGhost);
		free(l1para->xPrimary);

		if (l1para->flagl1tgate == 0)
		{
			free(l1para->hilbert_hw);
			free(l1para->hilbert_conv);
			hilbert_destroy(l1para->hilbert_task);
		}

		flexible_free_array2d(l1para->cg_d);
		flexible_free_array2d(l1para->cg_r);
		flexible_free_array2d(l1para->cg_wd);
		flexible_free_array2d(l1para->cg_Ad);

		if (l1para->cgmethod == 1)
		{
			flexible_free_array2d(l1para->cg_s); 
			flexible_free_array2d(l1para->cg_As);  
			flexible_free_array2d(l1para->cg_rhat0);  
		}

		if(l1para->fpredatum==FPREDATUM_DGRG)
		{
			flexible_free_array2d(l1para->TPR_p);
			flexible_free_array2d(l1para->TGR_p);
		}
	}
}


void allocate_l1inv(l1inv_t *l1para) 
{

	int fftnr;
	int nsgtrace = l1para->maxnsgtrace;
	int nsgtracepad   = ceil((float)(nsgtrace  /4.))*4;

	int npx = MAX(ceil(((l1para->psample*l1para->chanitv)*(float)(nsgtrace))*
				((l1para->pxmax-l1para->pxmin)*l1para->highfreq)), nsgtrace);
	l1para->npxorig = npx;

	int npxpad = ceil((float)(npx/4.))*4;

	int npy = 1;

	l1para->npyorig = npy;

	int npypad = ceil((float)(npy/4.))*4;

	//// ** SRAY ** CHECK ** NEW **
	if (l1para->choosemethod == CHOOSEMETHOD_JOINTSR2D &&
			l1para->jointmethod != JOINTDEGHOST)
		PFL_LENGTH(l1para->nsamporig+ l1para->nsampsrc,&fftnr);
	else
		PFL_LENGTH(l1para->nsamporig,&fftnr);

	fftnr+=2;

	//add by kyang try speed up
	l1para->orgnr =        fftnr-2;
	l1para->orgnc =        l1para->orgnr/2+1;
	l1para->hilbert_nt =   61;

	l1para->zwtilt = (float *) malloc(fftnr*sizeof(float));
	l1para->ywtilt = (float *) malloc(fftnr*sizeof(float));
	l1para->zlowf = (float *) malloc(fftnr*sizeof(float));
	l1para->ylowf = (float *) malloc(fftnr*sizeof(float));

	l1para->ztop =  (float *) malloc(npx*npy*sizeof(float));
	l1para->ytop =  (float *) malloc(npx*npy*sizeof(float));
	l1para->freqfac = (float *) malloc(fftnr*sizeof(float));
	l1para->obliqcorr = (float *) malloc(MAX(nsgtrace,npx*npy)*sizeof(float));

	l1para->wtilt = (float *) malloc(fftnr*sizeof(float));  //// ** DO WE NEED THIS ?? **

	l1para->pxsave = (float *) malloc(npx*npy*sizeof(float));
	l1para->pysave = (float *) malloc(npx*npy*sizeof(float));

	l1para->AP_p =           (float**)flexible_array2d(npx*npy, nsgtracepad, 2*sizeof(float));
	l1para->AG_p =           (float**)flexible_array2d(npx*npy, nsgtracepad, 2*sizeof(float));
	l1para->AS_p =           (float**)flexible_array2d(npx*npy, nsgtracepad, 2*sizeof(float));
	l1para->TP_p =           (float**)flexible_array2d(nsgtrace, npxpad*npypad, 2*sizeof(float));
	l1para->TG_p =           (float**)flexible_array2d(nsgtrace, npxpad*npypad, 2*sizeof(float));
	l1para->TS_p =           (float**)flexible_array2d(nsgtrace, npxpad*npypad, 2*sizeof(float));

	l1para->fftNR =        ALIGN_FLOAT(l1para->orgnr+2);     //for compatable with fftnr
	l1para->fftNC =        ALIGN_FCPLX(l1para->orgnc+1);     //for compatable with fftnc

	l1para->pfp  =         (float**)flexible_array2d(npx*npy, l1para->fftNC, 2*sizeof(float));
	if(l1para->regparam > 0.0f)
		l1para->pfpinit  =     (float**)flexible_array2d(npx*npy, l1para->fftNC, 2*sizeof(float)); // LKL
	else
		l1para->pfpinit  = NULL;
	l1para->pfpin    =     (float**)flexible_array2d(npx*npy, l1para->fftNC, 2*sizeof(float));
	l1para->pfpadd   =     (float**)flexible_array2d(npx*npy, l1para->fftNC, 2*sizeof(float));
	l1para->pfpT =         (float**)flexible_array2d(l1para->fftNC, npx*npy, 2*sizeof(float));
	l1para->pfx  =         (float**)flexible_array2d(nsgtrace, l1para->fftNC, 2*sizeof(float));
	l1para->pfxT =         (float**)flexible_array2d(l1para->fftNC, nsgtrace, 2*sizeof(float));

	l1para->integration = (float *) malloc(l1para->fftNC*2*sizeof(float));

	l1para->xGhost =       (float*)malloc(MAX(nsgtrace, npx*npy)*2*sizeof(float));
	l1para->xPrimary =     (float*)malloc(MAX(nsgtrace, npx*npy)*2*sizeof(float));

	if (l1para->flagl1tgate == 0)
	{
		l1para->hilbert_hw =   (float*)malloc((l1para->hilbert_nt*2+1)*sizeof(float));
		l1para->hilbert_conv = (float*)malloc((l1para->hilbert_nt*2+l1para->orgnr)*sizeof(float));
		l1para->hilbert_task = hilbert_init(l1para->hilbert_nt, l1para->orgnr, l1para->hilbert_hw);
	}
	//printf("l1para->fftNR=%d fftnr=%d, orgnr=%d\n", l1para->fftNR, fftnr, l1para->orgnr); fflush(0); 

	assert(l1para->fftNR>=fftnr);

	l1para->cg_d  = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));
	l1para->cg_r  = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));
	l1para->cg_wd = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));
	l1para->cg_Ad = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));

	if (l1para->cgmethod==1)
	{
		l1para->cg_s     = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));
		l1para->cg_As    = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));
		l1para->cg_rhat0 = (float**)flexible_array2d(npx*npy, fftnr, sizeof(float));
	}

	if(l1para->fpredatum==FPREDATUM_DGRG)
	{
		l1para->TPR_p =           (float**)flexible_array2d(nsgtrace, npxpad*npypad, 2*sizeof(float));
		l1para->TGR_p =           (float**)flexible_array2d(nsgtrace, npxpad*npypad, 2*sizeof(float));
	}
	if(l1para->choosemethod == CHOOSEMETHOD_JOINTSR2D)
	{
		l1para->AGS_p   = (float**)flexible_array2d(npx*npy, nsgtracepad, 2*sizeof(float));
		l1para->AGR_p   = (float**)flexible_array2d(npx*npy, nsgtracepad, 2*sizeof(float));
		l1para->TGS_p   = (float**)flexible_array2d(nsgtrace, npxpad*npypad, 2*sizeof(float));
		l1para->TGR_p   = (float**)flexible_array2d(nsgtrace, npxpad*npypad, 2*sizeof(float));
		l1para->xGhostS = (float*)malloc(MAX(nsgtrace, npx*npy)*2*sizeof(float));
		l1para->xGhostR = (float*)malloc(MAX(nsgtrace, npx*npy)*2*sizeof(float));
		if (l1para->computemode == 1)
			l1para->pSrcFilter = (float**)flexible_array2d(npx*npy, l1para->fftNC, 2*sizeof(float));
		if (l1para->jointmethod != JOINTDEGHOST)
		{
			l1para->pSrcFilterAz = (float**)flexible_array2d(npx*npy, l1para->fftNC, 2*sizeof(float));
			l1para->pTgtFilter   = (float**)flexible_array2d(npx*npy, l1para->fftNC, 2*sizeof(float));

			l1para->p1dSrc = (float**)flexible_array2d(1, l1para->fftNC, 2*sizeof(float));
			l1para->p1dTgt = (float**)flexible_array2d(1, l1para->fftNC, 2*sizeof(float));
		}
	}
}

void deallocate_l1inv(l1inv_t *l1para) 
{

	free(l1para->zwtilt);
	free(l1para->ywtilt);
	free(l1para->zlowf);
	free(l1para->ylowf);
	free(l1para->ztop);
	free(l1para->ytop);
	free(l1para->freqfac);
	free(l1para->obliqcorr);
	free(l1para->wtilt);
	free(l1para->integration);
	flexible_free_array2d(l1para->pfpadd);
	flexible_free_array2d(l1para->pfpin);

	free(l1para->pxsave);
	free(l1para->pysave);

	flexible_free_array2d(l1para->AP_p);
	flexible_free_array2d(l1para->AG_p);
	flexible_free_array2d(l1para->AS_p);
	flexible_free_array2d(l1para->TP_p);
	flexible_free_array2d(l1para->TG_p);
	flexible_free_array2d(l1para->TS_p);

	flexible_free_array2d(l1para->pfp);
	if(l1para->regparam > 0.0f) flexible_free_array2d(l1para->pfpinit); // LKL
	flexible_free_array2d(l1para->pfx);
	flexible_free_array2d(l1para->pfxT);
	flexible_free_array2d(l1para->pfpT);

	free(l1para->xGhost);
	free(l1para->xPrimary);

	if (l1para->flagl1tgate == 0)
	{
		free(l1para->hilbert_hw);
		free(l1para->hilbert_conv);
		hilbert_destroy(l1para->hilbert_task);
	}

	flexible_free_array2d(l1para->cg_d);
	flexible_free_array2d(l1para->cg_r);
	flexible_free_array2d(l1para->cg_wd);
	flexible_free_array2d(l1para->cg_Ad);

	if (l1para->cgmethod == 1)
	{
		flexible_free_array2d(l1para->cg_s); 
		flexible_free_array2d(l1para->cg_As);  
		flexible_free_array2d(l1para->cg_rhat0);  
	}

	if(l1para->fpredatum==FPREDATUM_DGRG)
	{
		flexible_free_array2d(l1para->TPR_p);
		flexible_free_array2d(l1para->TGR_p);
	}
	if(l1para->choosemethod == CHOOSEMETHOD_JOINTSR2D)
	{
		flexible_free_array2d(l1para->AGS_p);
		flexible_free_array2d(l1para->AGR_p);
		flexible_free_array2d(l1para->TGS_p);
		flexible_free_array2d(l1para->TGR_p);
		free(l1para->xGhostS);
		free(l1para->xGhostR);
		if (l1para->computemode == 1)
			flexible_free_array2d(l1para->pSrcFilter);
		if (l1para->jointmethod != JOINTDEGHOST)
		{
			flexible_free_array2d(l1para->pSrcFilterAz);
			flexible_free_array2d(l1para->pTgtFilter);

			flexible_free_array2d(l1para->p1dSrc);
			flexible_free_array2d(l1para->p1dTgt);
		}
	}

}

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

void calctgate_bsdg3d(int ntrc, int *wbid, int *tgate, l1inv_t *l1para);

void bsdg3d_deghostc (l1inv_t *l1para, float *offsetx, float *offsety, float *recz, int *wbid,
		    float *pOrigin_p, float *pDeghost, int ntrcxy)
{
        int i,k,isamp,itrc,fftnr,fftnc,orgnr,orgnc,ftrc,ltrc;
        float **winout, **input_p, sppow;
	int lfid, trlen;

	int *tgate,itgate,itbeg,itend,startid;
	float rtap,*lintapert;

	float **winout_tgate, **input_p_tgate;

	trlen = l1para->nsamp + NHEAD;
	l1para->nsamporig = l1para->nsamp;

	l1para->nsgtracewin = ntrcxy;

	input_p  = NULL;	
	winout = NULL;	
	input_p  = myallocatef (input_p,  ntrcxy, l1para->nsamp);
	winout = myallocatef (winout, ntrcxy, l1para->nsamp);

	if ( l1para->flagl1tgate == 1 )
	{
	    lintapert=(float *)malloc(l1para->nsamporig*sizeof(float));
	    for (i=0;i<l1para->nsamporig;i++) lintapert[i]=1.0f;
	    for (i=0;i<2*l1para->tgoverlap-1;i++) lintapert[i]=(float)(i+1)/(float)(2*l1para->tgoverlap);
	}

	for (itrc = 0; itrc < ntrcxy; itrc++)
		memcpy (input_p[itrc], pOrigin_p + itrc * trlen + NHEAD, l1para->nsamp * sizeof (float));

	if ( l1para->flagl1tgate == 0 )
	{

	    ////////////////////////////////////
	    allocate_deghost3d (l1para, ntrcxy); 
	    ////////////////////////////////////
    
	    sppow = l1para->swnear; 
	    
	    lfid = max (l1para->lowfar * l1para->srate / 500000.0f * l1para->fftnc, 4);
	    
	    calc_p (l1para);
	    
	    if (l1para->choosemethod == CHOOSEMETHOD_JOINTSR3D)
	      jointsr3d_deghosting  (l1para, ntrcxy, offsetx, offsety, recz, 
				     input_p, winout, 
				     sppow, lfid);
	    else
	      bsdg3d_deghosting  (l1para, ntrcxy, offsetx, offsety, recz, 
				  input_p, winout, 
				  sppow, lfid);
	    
	    for (itrc = 0; itrc < ntrcxy; itrc++)
	      memcpy (pDeghost + itrc * trlen + NHEAD, winout[itrc], l1para->nsamp * sizeof (float));
	    
	    ////////////////////////////////////
	    deallocate_deghost3d (l1para); 
	    ////////////////////////////////////
	    
	}
	else if ( l1para->flagl1tgate == 1 )
	{

	    tgate = (int *) calloc (l1para->ntgate, sizeof (int)); 

 	    calctgate_bsdg3d(ntrcxy, wbid, tgate, l1para);
	    
	    itgate = -1;
	    itend = -1111;
	    
	    while (itend < l1para->nsamporig-1)        
	    {
		
		itgate = itgate + 1;

		itbeg=MAX(0,tgate[itgate]-MAX(0,l1para->tgoverlap));
		itend=MIN(l1para->nsamporig-1,tgate[itgate+1]+(l1para->tgoverlap-1));
		if (tgate[itgate+1] == l1para->nsamporig-1) itend = l1para->nsamporig-1;

		l1para->nsamp = itend - itbeg + 1;              
		
		PFL_LENGTH(l1para->nsamp,&fftnr);
		fftnr+=2;
		fftnc = fftnr/2;

		l1para->fftnr = fftnr;
		l1para->fftnc = fftnc;

		orgnr =        fftnr-2;
		orgnc =        orgnr/2+1;
		
		l1para->hilbert_nt =   61;
		l1para->hilbert_hw =   (float*)malloc((l1para->hilbert_nt*2+1)*sizeof(float));
		l1para->hilbert_conv = (float*)malloc((l1para->hilbert_nt*2+orgnr)*sizeof(float));
		l1para->hilbert_task = hilbert_init(l1para->hilbert_nt, orgnr, l1para->hilbert_hw);
		
		////////////////////////////////////
		allocate_deghost3d (l1para, ntrcxy); 
		////////////////////////////////////
		
		sppow = l1para->swnear; 
		
		lfid = max (l1para->lowfar * l1para->srate / 500000.0f * l1para->fftnc, 4);
		
		calc_p (l1para);
		
		
 		//calc_matrix_A_T_p (l1para, ntrcxy, offsetx, offsety, recz, lfid);
		
		input_p_tgate  = NULL;
		winout_tgate   = NULL;
		
		input_p_tgate = myallocatef(input_p_tgate,ntrcxy,l1para->nsamp);
		winout_tgate  = myallocatef(winout_tgate, ntrcxy,l1para->nsamp);
		
		for (itrc=0;itrc<ntrcxy;itrc++)
		{
		    k = -1;
		    for (isamp=itbeg;isamp<=itend;isamp++)
		    {
			k = k + 1;
			input_p_tgate[itrc][k]= input_p[itrc][isamp];
		    }
		}

		if (l1para->choosemethod == CHOOSEMETHOD_JOINTSR3D)
		  jointsr3d_deghosting  (l1para, ntrcxy, offsetx, offsety, recz, 
				      input_p_tgate, winout_tgate, 
				      sppow, lfid);
		else
		  bsdg3d_deghosting  (l1para, ntrcxy, offsetx, offsety, recz, 
				      input_p_tgate, winout_tgate, 
				      sppow, lfid);
				
		free(l1para->hilbert_hw);
		free(l1para->hilbert_conv);
		hilbert_destroy(l1para->hilbert_task);

		ftrc = 0;
		ltrc = ntrcxy-1;


		if (itgate==0)
		{
		    for (itrc=ftrc;itrc<=ltrc;itrc++)
		      for (isamp=itbeg;isamp<=itend;isamp++)
			winout[itrc][isamp]=winout_tgate[itrc][isamp-itbeg];
		}
		else
		{
		    for (itrc=ftrc;itrc<=ltrc;itrc++)
		    {
			
			startid=itbeg;
			
			for (isamp=startid;isamp<=itend;isamp++)
			{
			    rtap=lintapert[isamp-startid];
			    winout[itrc][isamp]=winout_tgate[itrc][isamp-itbeg]*rtap+winout[itrc][isamp]*(1.0f-rtap);
			}
		    }
		}
		
		myfree(input_p_tgate, ntrcxy);
		myfree(winout_tgate,  ntrcxy);

		////////////////////////////////////
		deallocate_deghost3d (l1para); 
		////////////////////////////////////

	    }

	    l1para->nsamp = l1para->nsamporig;

	    for (itrc = 0; itrc < ntrcxy; itrc++)
	      memcpy (pDeghost + itrc * trlen + NHEAD, winout[itrc], l1para->nsamp * sizeof (float));

	    free(tgate);

	}

	myfree (input_p,  ntrcxy);
	myfree (winout,   ntrcxy);

	if ( l1para->flagl1tgate == 1 ) free(lintapert);

	return;
}


 
void calctgate_bsdg3d(int ntrc, int *wbid, int *tgate, l1inv_t *l1para)
 {


   //// we will call it inside the isgatex loop in str3d_deghost()

   int i,j;
   int nsamp;
   int wbidmax,sum;

   int ntgate,zntgate;

   int winsize_wb = 1000; //// in ms

   winsize_wb = winsize_wb*1000/(l1para->srate)+1;

   nsamp = l1para->nsamporig; 

   ntgate =  l1para->ntgate;
   zntgate = l1para->zntgate;

   int *zerot=(int *)malloc(ntgate*sizeof(int));

   sum=0;
   
   for (i=0;i<zntgate;i++)
     {
       zerot[i]=l1para->tgsampmin+(l1para->tgsampmax-l1para->tgsampmin)*i/zntgate;
       sum+=zerot[i];
     }
   zerot[zntgate-1]+=nsamp-sum;

   wbidmax = -1;
   for (i=0;i<ntrc;i++)
     {
       wbidmax = MAX(wbid[i],wbidmax);   
     }

   tgate[0]=0;
       
   tgate[1]=MIN(nsamp-1,wbidmax+winsize_wb);
       
   for (j=2;j<MIN(zntgate+1,ntgate);j++)
     {
       tgate[j]=MIN(nsamp-1, tgate[j-1]+zerot[j-2]);
       if (tgate[j]>nsamp-l1para->tgsampmin) tgate[j]=nsamp-1; 
     }
       
   for (j=MIN(zntgate+1,ntgate);j<ntgate;j++)
     {
       tgate[j]=nsamp-1;
     }

   tgate[ntgate-1]=nsamp-1;

   free(zerot);
   
 }

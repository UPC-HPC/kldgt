#include "msdginterp.h"
#include "mtmsdg.h"
#include "invmatrix.h"
#include <fftw3.h>
#include <PFL_C.h>
#include <complex>

using namespace std;

#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )

void calctgate_bsdg3d(int ntrc, int *wbid, int *tgate, l1inv_t *l1para);

void bsdg3d_deghostc_pre(l1inv_t * l1para, float *offsetx, float *offsety,
                     float *recz, int *wbid, float *pOrigin_p, int ntrcxy)


{
    l1para->nsamporig = l1para->nsamp;
    l1para->nsgtracewin = ntrcxy;

    // Allocate device memory
    mtbsdg_allocate_input(l1para);

    mtbsdg_allocate_offset(l1para, ntrcxy);

    mtmsdg_allocate_deghost3d(l1para, ntrcxy);

    // Copy data from host to device
    mtbsdg_h2d_deghost3d(l1para, offsetx, offsety, recz, wbid,
                         pOrigin_p, ntrcxy);

    return;

}

void bsdg3d_deghostcgpu(l1inv_t * l1para, int *wbid, int ntrcxy)
{
    float sppow;
    int lfid;
    float *d_input, *d_winout;
    float *d_winout_tgate, *d_input_tgate,*d_lintapert;
    int itgate, itbeg, itend, itrc, isamp, trlen, ntgate, zntgate,
        fftnr, fftnc, *tgate;
    head_t *src;

    d_input = l1para->d_input_p;
    d_winout = l1para->d_winout;

    if (l1para->flagl1tgate == 1)
    {
      d_lintapert = NULL;
      DEV_SAFE_CALL(devMalloc(l1para->gpuid, (void **) &d_lintapert,
			l1para->nsamporig * sizeof(float)));

      mtmsdg_myinitf(l1para->gpuid, d_lintapert, 2 * l1para->tgoverlap - 1,
		     1.0f / (float) (2 * l1para->tgoverlap),
		     1.0f / (float) (2 * l1para->tgoverlap));
      
      mtmsdg_myinitf(l1para->gpuid, d_lintapert + 2 * l1para->tgoverlap - 1,
		     l1para->nsamporig - 2 * l1para->tgoverlap + 1,
		     1.0f, 0.0f);
    }

    sppow = l1para->swnear;
    
    if (l1para->flagl1tgate == 0)
      {
	////////////////////////////////////
	allocate_deghost3d(l1para, ntrcxy);
	////////////////////////////////////

	lfid =
	  max(l1para->lowfar * l1para->srate / 500000.0f * l1para->fftnc, 4);
	
	mtmsdg_calc_p(l1para);

	l1para->finalflag = YESYES;

	if(l1para->choosemethod == CHOOSEMETHOD_JOINTSR3D)
	 mtbsdg_deghosting_jointsr3d(l1para, ntrcxy, d_input, d_winout, sppow, lfid);
	else
	  mtbsdg_deghosting(l1para, ntrcxy, d_input, d_winout, sppow, lfid);
	
	////////////////////////////////////
	deallocate_deghost3d(l1para);
	////////////////////////////////////

      }
    else if (l1para->flagl1tgate == 1)
      {
	tgate = (int *) calloc (l1para->ntgate, sizeof (int)); 
        
        calctgate_bsdg3d(ntrcxy, wbid, tgate, l1para);

	itgate = -1;
        itend = -1111;

        while (itend < l1para->nsamporig - 1)
        {
            float zero = 0;
            itgate++;

            itbeg = MAX(0, tgate[itgate] - MAX(0, l1para->tgoverlap));
            itend = MIN(l1para->nsamporig-1,tgate[itgate+1]+(l1para->tgoverlap-1));
            if (tgate[itgate + 1] == l1para->nsamporig - 1) itend = l1para->nsamporig-1;

            l1para->nsamp = itend - itbeg + 1;

            DEV_SAFE_CALL(devMalloc(l1para->gpuid, (void **) &d_input_tgate,
                              ntrcxy * l1para->nsamp * sizeof(float)));
            DEV_SAFE_CALL(devMalloc(l1para->gpuid, (void **) &d_winout_tgate,
                              ntrcxy * l1para->nsamp * sizeof(float)));
            DEV_SAFE_CALL(devMemset(l1para->gpuid, 0, d_input_tgate, &zero, sizeof(float),
                              ntrcxy * l1para->nsamp * sizeof(float), NULL));
            DEV_SAFE_CALL(devMemset(l1para->gpuid, 0, d_winout_tgate, &zero, sizeof(float),
                              ntrcxy * l1para->nsamp * sizeof(float), NULL));

           
            for (int ht = 0; ht < ntrcxy; ht++)
 	    {
	            DEV_SAFE_CALL(devMemcpyDtoDOffsetAsync(l1para->gpuid, 0, 
				d_input_tgate,
				d_input,
                                ht * l1para->nsamp * sizeof(float),
                                (itbeg + ht * l1para->nsamporig) * sizeof(float),
                                (itend - itbeg + 1) * sizeof(float),
                                NULL));

            }
	    if ((l1para->choosemethod == CHOOSEMETHOD_JOINTSR3D ||
		 (l1para->choosemethod == CHOOSEMETHOD_MSDGI && l1para->jointmc == YESYES))
		&& (l1para->jointmethod != JOINTDEGHOST))
		PFL_LENGTH(l1para->nsamp + l1para->nsampsrc, &fftnr);
		//fftnr = mtmsdg_fft_length(l1para->nsamp + l1para->nsampsrc);
	    else
		PFL_LENGTH(l1para->nsamp, &fftnr);
		//fftnr = mtmsdg_fft_length(l1para->nsamp);

	    fftnc = fftnr/ 2 + 1;
	    ////////////////////////////////////
	    allocate_deghost3d (l1para, ntrcxy); 
	    ////////////////////////////////////

            lfid =
                max(l1para->lowfar * l1para->srate / 500000.0f * fftnc, 4);

            mtmsdg_calc_p(l1para);

	    if(l1para->choosemethod == CHOOSEMETHOD_JOINTSR3D)
	      mtbsdg_deghosting_jointsr3d(l1para, ntrcxy, d_input_tgate, d_winout_tgate, sppow, lfid);
	    else
	      mtbsdg_deghosting(l1para, ntrcxy, d_input_tgate, d_winout_tgate, sppow, lfid);


            if (itgate == 0)
            {
		for(int ht=0; ht < ntrcxy; ht++)
                {
                    DEV_SAFE_CALL(devMemcpyDtoDOffsetAsync( l1para->gpuid, 0 ,
                                                            d_winout ,
                                                            d_winout_tgate,
                                                            (itbeg + ht* l1para->nsamporig)* sizeof(float),
                                                             ht *l1para->nsamp * sizeof(float),
                                                            (itend - itbeg + 1) * sizeof(float),
                                                            NULL));

                }

            }
            else
            {
                mtmsdg_mytaperadd(l1para->gpuid, d_winout + itbeg, d_winout_tgate,
                                  d_lintapert, itend - itbeg + 1, ntrcxy,
                                  l1para->nsamporig, l1para->nsamp);
            }

            DEV_SAFE_CALL(devFree(l1para->gpuid, d_input_tgate));
            DEV_SAFE_CALL(devFree(l1para->gpuid, d_winout_tgate));

        }       // // end of : while (itend < l1para->nsamporig-1) 
        l1para->nsamp = l1para->nsamporig;
	free(tgate);
        DEV_SAFE_CALL(devFree(l1para->gpuid, d_lintapert));

      }   // // end of : if ( l1para->flagl1tgate == 1 )

    return;
}

void bsdg3d_deghostc_post(l1inv_t * l1para,   float *deghost,
                     int ntrcxy, int tasksizemax)
{

  mtmsdg_copyback_deghost3d(l1para->gpuid, deghost, 
			    l1para->d_winout, 
			    ntrcxy,
			    l1para->nsamp, NHEAD);
  
  mtmsdg_deallocate_deghost3d(l1para);
  mtbsdg_deallocate_input(l1para);
  mtbsdg_deallocate_offset(l1para);
  
  return;
}



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
#define FFT_R2C_PP(plan, t, m, f, nf2)                                                    \
	for(int i=0; i<m; i++)                                                                  \
PFL_RCFFT(plan, (nf2-2),(float*)(t)[i],(PFLComplex*)(f)[i],1,-1,(nf2),(nf2),1.0f/((nf2)-2));

void calc_zfilt(l1inv_t *l1para);
void calc_matrix2d_A_T_p_joint(l1inv_t *l1para, int nsgtrace, float *offsetx, float *recz, int lfid) //// joint = jointdesignature option
{
	int fftnr, npx, npy, ip, ipx, ipy, itrc, isamp, ifreq, hpfl;
	float fnyq, deltaf, fsin, fcos, px=0.0f, py=0.0f, pxy=0.0f;
	float phase0, phase1, phase2, invwater=1.0/l1para->vwater, x,y;

	complex<float> ** __restrict AP,** __restrict AG, ** __restrict AS;
	complex<float> ** __restrict TP,** __restrict TG,** __restrict TS;
	complex<float> ** pSrcFilter, ** pSrcFilterAz, ** pTgtFilter;
	complex<float> ** p1dSrc = (complex<float> **) l1para->p1dSrc;
	complex<float> ** p1dTgt = (complex<float> **) l1para->p1dTgt;

	//PFL_LENGTH(l1para->nsamp,&fftnr);
	PFL_GET (&hpfl, PFL_PORTABLE_DEF);
	PFL_LENGTH(l1para->nsamp + l1para->nsampsrc,&fftnr);
	l1para->orgnr = fftnr;
	l1para->orgnc = fftnr/2+1;
	l1para->fftNR =        ALIGN_FLOAT(l1para->orgnr+2);     //for compatable with fftnr
	l1para->fftNC =        ALIGN_FCPLX(l1para->orgnc+1);     //for compatable with fftnc
	fftnr+=2;
	l1para->fftnr = fftnr;
	l1para->fftnc = fftnr/2;
	fnyq = 500000.0f/l1para->srate;
	deltaf = fnyq/(float)(l1para->fftnc-1);

	float twopif=deltaf*2.0f*3.1415927f;
	float *l_f = (float *) calloc (l1para->fftnc, sizeof (float));
	for (ifreq = 0; ifreq < l1para->fftnc; ifreq++) l_f[ifreq] = twopif * ifreq; 

	npx = l1para->npxorig;
	npy = 1;

	float deltap  = (l1para->pxmax-l1para->pxmin)/(float)(npx-1);  //// <-- BUG FIX

	AP = (complex<float> **)l1para->AP_p;
	AG = (complex<float> **)l1para->AG_p;
	AS = (complex<float> **)l1para->AS_p;

	TP = (complex<float> **)l1para->TP_p;
	TG = (complex<float> **)l1para->TG_p;
	TS = (complex<float> **)l1para->TS_p;

	if (l1para->choosemethod == CHOOSEMETHOD_JOINTSR2D) 
	{
		pSrcFilter   = (complex<float> **) l1para->pSrcFilter;
		pSrcFilterAz = (complex<float> **) l1para->pSrcFilterAz;
		pTgtFilter   = (complex<float> **) l1para->pTgtFilter;
	}

	head_t *src;
	float *gunvol = (float *) calloc (l1para->ntrcsrc, sizeof (float));
	float *psx    = (float *) calloc (l1para->ntrcsrc, sizeof (float));
	float *psy    = (float *) calloc (l1para->ntrcsrc, sizeof (float));
	float *psz    = (float *) calloc (l1para->ntrcsrc, sizeof (float));
	float shotx_avg, shoty_avg, shotz_avg, total_gunvol;
	float          **ptxsource = (float          **) flexible_array2d (l1para->ntrcsrc, l1para->fftnr, sizeof (float));
	float          **ptxtarget = (float          **) flexible_array2d (l1para->ntrcsrc, l1para->fftnr, sizeof (float));
	complex<float> **pfxsource = (complex<float> **) flexible_array2d (l1para->ntrcsrc, l1para->fftnc, sizeof (complex<float>));
	complex<float> **pfxtarget = (complex<float> **) flexible_array2d (l1para->ntrcsrc, l1para->fftnc, sizeof (complex<float>));

	float delaytot;
	int delaycount;

	int trlensig = NHEAD + l1para->nsampsrc;
	shotx_avg = shoty_avg = shotz_avg = total_gunvol = 0.0f;

	l1para->gunmin=10000;
	l1para->gunmax=0;

	float gunamptot = 0.0f;
	float shgun = 0.0f;
	float dpgun = 0.0f;

	for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
	{
		src = (head_t *) (l1para->psig + itrc * trlensig);
		psx[itrc] = src->shotx;
		psy[itrc] = src->shoty;
		psz[itrc] = src->shotz;
		shotx_avg += src->shotx;
		shoty_avg += src->shoty;
		shotz_avg += src->shotz;
		memset (ptxsource[itrc], 0, l1para->fftnr * sizeof (float));
		memset (ptxtarget[itrc], 0, l1para->fftnr * sizeof (float));
		memcpy (ptxsource[itrc], l1para->psig + itrc * trlensig + NHEAD,
				l1para->nsampsrc * sizeof (float));
		memcpy (ptxtarget[itrc], l1para->psig + itrc * trlensig + NHEAD,
				(l1para->tzerosrc + l1para->tzerotgt) * 1000.0f * 1000.0f / l1para->srate * sizeof (float));

		//for (isamp = 0; isamp < l1para->nsampsrc; isamp++)
		//  gunvol[itrc] = MAX (gunvol[itrc], fabsf (ptxsource[itrc][isamp]));
		if (l1para->designorm == DESIGNORM_MAX || l1para->designorm == DESIGNORM_NONE)
		{
			for (isamp = 0; isamp < l1para->nsampsrc; isamp++)
				gunvol[itrc] = MAX (gunvol[itrc], fabsf (ptxsource[itrc][isamp]));
		}
		else if (l1para->designorm == DESIGNORM_RMS)
		{
			for (isamp = 0; isamp < l1para->nsampsrc; isamp++)
				gunvol[itrc] += ptxsource[itrc][isamp] * ptxsource[itrc][isamp];
			gunvol[itrc] = sqrtf (gunvol[itrc] / l1para->nsampsrc);
		}
		total_gunvol += gunvol[itrc];

		l1para->gunmin=MIN(l1para->gunmin, src->shotz); 
		l1para->gunmax=MAX(l1para->gunmax, src->shotz); 

		gunamptot += src->gunamp;

		if (src->tdelay>0) 
		{
			delaycount++;
			delaytot+=src->tdelay;
			////printf("delaycount=%d,delaytot=%f\n",delaycount,delaytot);fflush(0);
		}
	}

	l1para->total_gunvol = total_gunvol;

	if (l1para->jointmethod==JOINT_SRC_DEGHOST)
	{
		for (itrc=0;itrc<nsgtrace;itrc++)
		{
			if (l1para->bsouttype==1||l1para->bsouttype==2)
				recz[itrc]  =  l1para->gunmin;
			else if (l1para->bsouttype==3||l1para->bsouttype==4)
				recz[itrc]  =  l1para->gunmax;
		}
	}

	for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
	{
		src = (head_t *) (l1para->psig + itrc * trlensig);
		if (src->shotz==l1para->gunmin) shgun += src->gunamp;
		else if (src->shotz==l1para->gunmax) dpgun += src->gunamp;
	}

	////printf("guntotal=%f,gundeep=%f,gunshallow=%f\n",l1para->totgun,l1para->shgun,l1para->dpgun);
	////fflush(0); 

	l1para->totgun = gunamptot;
	l1para->shgun = shgun;
	l1para->dpgun = dpgun;

	if (delaycount>0) l1para->gundelay=delaytot/(1000.0f*delaycount);
	else l1para->gundelay=0.0f;
	////printf("deep gun delay time = %f\n",l1para->gundelay);fflush(0);

	shotx_avg /= l1para->ntrcsrc;
	shoty_avg /= l1para->ntrcsrc;
	shotz_avg /= l1para->ntrcsrc;
	if (total_gunvol > 0.0f && l1para->designorm != DESIGNORM_NONE)
	{
		for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
		{
			for (isamp = 0; isamp < l1para->nsampsrc; isamp++)
			{
				ptxsource[itrc][isamp] /= total_gunvol;
				ptxtarget[itrc][isamp] /= total_gunvol;
			}
		}
	}
	if (total_gunvol > 0.0f)
		for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
			gunvol[itrc] /= total_gunvol;

	FFT_R2C_PP (hpfl, ptxsource, l1para->ntrcsrc, pfxsource, l1para->fftnr); 
	FFT_R2C_PP (hpfl, ptxtarget, l1para->ntrcsrc, pfxtarget, l1para->fftnr); 
	for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
	{
		for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
		{
			pfxsource[itrc][ifreq] *= (l1para->fftnr - 2);
			pfxtarget[itrc][ifreq] *= (l1para->fftnr - 2);
		}
	}

	ip = -1;
	//   for (ipy = 0; ipy < npy; ipy++)
	//   {
	for (ipx = 0; ipx < npx; ipx++)
	{
		ip++;

		px = l1para->pxmin + (float)(ip)*deltap; //// <-- BUG FIX
		py = 0.0f;

		pxy=sqrtf(px*px);
		fsin=MAX(-0.99999f,MIN(0.99999f,pxy*l1para->vwater));       
		fcos=sqrtf(1.0f-fsin*fsin);//cosine of 3D surface angle

		// For reverse tau-p operator and receiver ghost operator
		for (itrc = 0; itrc < nsgtrace; itrc++)
		{
			phase0 = twopif * px * offsetx[itrc];

			phase1 = twopif * MAX (5.0e-4 * l1para->tminscl, fcos * recz[itrc] * invwater);//half of ghost shift for reghosting
			phase2 = phase0 + phase1; //mirror cable
			phase1 = phase0 - phase1; //cable

			x = cosf(phase1); y = sinf(phase1); // for cable depth
			AP[ip][itrc] = complex<float> (x,  y);
			TP[itrc][ip] = complex<float> (x, -y);

			x = cosf(phase2); y = sinf(phase2); // for mirror cable depth
			AG[ip][itrc] = complex<float> (x,  y);
			TG[itrc][ip] = complex<float> (x, -y);

			x = cosf(phase0); y = sinf(phase0); // for reverse taup
			AS[ip][itrc] = complex<float> (x,  y);
			TS[itrc][ip] = complex<float> (x, -y);
		} // for ( itrc = 0;itrc < nsgtrace; itrc++ )

		// zxue
		// For source ghost operator
		for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
			pSrcFilter[ip][ifreq] = 0.0f;
		float tmp1, tmp2;
		if (l1para->datum == CABLE90)
			for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
			{
				src = (head_t *) (l1para->psig + itrc * (NHEAD + l1para->nsampsrc));
				float tsrcghost = MAX (5.0e-4 * l1para->tminscl, fcos * src->shotz * invwater);
				for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
				{
					tmp1 = l_f[ifreq] * 2.0f * tsrcghost;
					pSrcFilter[ip][ifreq] += gunvol[itrc] * (1.0f + l1para->srcrfl[ifreq] * complex<float> (cosf (tmp1), -sinf (tmp1))); 
				}
			}
		else
			for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
			{
				src = (head_t *) (l1para->psig + itrc * (NHEAD + l1para->nsampsrc));
				float tsrcghost = MAX (5.0e-4 * l1para->tminscl, fcos * src->shotz * invwater);
				for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
				{
					tmp1 = l_f[ifreq] * ( tsrcghost - src->tdelay);
					tmp2 = l_f[ifreq] * (-tsrcghost - src->tdelay);
					pSrcFilter[ip][ifreq] += gunvol[itrc] * (complex<float> (cosf (tmp1), sinf (tmp1)) + 
							l1para->srcrfl[ifreq] *  complex<float> (cosf (tmp2), sinf (tmp2)));
				}
			}

		// for source signature operator
		for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
			pSrcFilterAz[ip][ifreq] = pTgtFilter[ip][ifreq] = 0.0f;
		for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
		{
			head_t *src = (head_t *) (l1para->psig + itrc * (NHEAD + l1para->nsampsrc));
			float tdelay = -(psx[itrc] - shotx_avg) * px - (psy[itrc] - shoty_avg) * py - src->tdelay - l1para->tzerosrc;
			for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
			{
				tmp1 = l_f[ifreq] * tdelay;
				pSrcFilterAz[ip][ifreq] += pfxsource[itrc][ifreq] * complex<float> (cosf (tmp1), -sinf (tmp1));
				pTgtFilter[ip][ifreq]   += pfxtarget[itrc][ifreq] * complex<float> (cosf (tmp1), -sinf (tmp1));
			}
		}

		if (l1para->jointmethod == JOINTDESIG_SRC_DEGHOST || l1para->jointmethod == JOINTDESIG_SRCRCV_DEGHOST)
			for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
				pSrcFilter[ip][ifreq] *= pSrcFilterAz[ip][ifreq];
		else if (l1para->jointmethod == JOINTDESIGNATURE || l1para->jointmethod == JOINTDESIG_RCV_DEGHOST)
			for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
				pSrcFilter[ip][ifreq]  = pSrcFilterAz[ip][ifreq];
	} // for (ipx = 0; ipx < npx; ipx++)
	//   } // for (ipy = 0; ipy < npy; ipy++)

	if (l1para->lmerge1d == YESYES)
	{
		ip = (int) (-l1para->pxmin / deltap);
		assert (0 <= ip && ip < npx);
		memcpy (p1dSrc[0], pSrcFilter[ip], l1para->fftnc * 2 * sizeof (float));
		memcpy (p1dTgt[0], pTgtFilter[ip], l1para->fftnc * 2 * sizeof (float));
	}

	free(l_f);
	free (gunvol); free (psx); free (psy); free (psz);
	flexible_free_array2d (ptxsource);
	flexible_free_array2d (ptxtarget);
	flexible_free_array2d (pfxsource);
	flexible_free_array2d (pfxtarget);
	PFL_FREE (hpfl);
}

void calc_matrix2d_A_T_p(l1inv_t *l1para, int nsgtrace, float *offsetx, float *recz, float *shotz, int lfid)
{
	int fftnr, npx, npy;
	int ip, ipx, ipy, itrc;
	float fnyq, deltaf;
	float fsin, fcos, px=0.0f, py=0.0f, pxy=0.0f;
	float phase0, phase1, phase2, invwater=1.0/l1para->vwater;
	int ifreq;

	float _Complex ** __restrict AP,** __restrict AG, ** __restrict AS;
	float _Complex ** __restrict TP,** __restrict TG,** __restrict TS;
	float _Complex ** __restrict TPR,** __restrict TGR;
	float _Complex ** __restrict AGS,** __restrict AGR, ** __restrict TGS;
	complex<float> ** pSrcFilter;

	float *ptr,x,y;

	float *reczou;
	//  float *shotz;

	if (l1para->fpredatum == FPREDATUM_DGRG) reczou = &recz[nsgtrace];

	PFL_LENGTH(l1para->nsamp,&fftnr);
	l1para->orgnr = fftnr;
	l1para->orgnc = fftnr/2+1;
	l1para->fftNR =        ALIGN_FLOAT(l1para->orgnr+2);     //for compatable with fftnr
	l1para->fftNC =        ALIGN_FCPLX(l1para->orgnc+1);     //for compatable with fftnc
	fftnr+=2;
	l1para->fftnr = fftnr;
	l1para->fftnc = fftnr/2;

	fnyq = 500000.0f/l1para->srate;
	deltaf = fnyq/(float)(l1para->fftnc-1);

	float twopif=deltaf*2.0f*3.1415927f;
	float *l_f = (float *) calloc (l1para->fftnc, sizeof (float));
	for (ifreq = 0; ifreq < l1para->fftnc; ifreq++) l_f[ifreq] = twopif * ifreq; 

	npx = l1para->npxorig;
	npy = 1;

	float deltap  = (l1para->pxmax-l1para->pxmin)/(float)(npx-1);  //// <-- BUG FIX

	AP = (float _Complex **)l1para->AP_p;
	AG = (float _Complex **)l1para->AG_p;
	AS = (float _Complex **)l1para->AS_p;

	TP = (float _Complex **)l1para->TP_p;
	TG = (float _Complex **)l1para->TG_p;
	TS = (float _Complex **)l1para->TS_p;

	if (l1para->fpredatum == FPREDATUM_DGRG) 
	{
		TPR = (float _Complex **)l1para->TPR_p;
		TGR = (float _Complex **)l1para->TGR_p;
	}
	else if (l1para->choosemethod == CHOOSEMETHOD_JOINTSR2D) 
	{
		AGR = (float _Complex **)l1para->AGR_p;
		AGS = (float _Complex **)l1para->AGS_p;
		TGR = (float _Complex **)l1para->TGR_p;
		TGS = (float _Complex **)l1para->TGS_p;
		if (l1para->computemode == 1)
			pSrcFilter = (complex<float> **) l1para->pSrcFilter;
	}

	ip = -1;
	//  for ( ipy = 0;ipy < npy; ipy++ )
	//   {
	for ( ipx = 0;ipx < npx; ipx++ )
	{

		ip = ip + 1;

		px = l1para->pxmin + (float)(ip)*deltap; //// <-- BUG FIX

		////px = l1para->pxsave[ip];
		py = 0.0f;

		pxy=sqrtf(px*px+py*py);
		fsin=MAX(-0.99999f,MIN(0.99999f,pxy*l1para->vwater));       
		fcos=sqrtf(1.0f-fsin*fsin);//cosine of 3D surface angle

		for ( itrc = 0;itrc < nsgtrace; itrc++ )
		{
			phase0 = twopif*px*offsetx[itrc];

			phase1 = twopif*MAX(5.0e-4*l1para->tminscl,fcos*recz[itrc]*invwater);//half of ghost shift for reghosting
			phase2 = phase0+phase1; //mirror cable
			phase1 = phase0-phase1; //cable

			ptr = (float *) &AP[ip][itrc];
			x = cosf(phase1);
			y = sinf(phase1);
			ptr[0] = x;
			ptr[1] = y;

			ptr = (float *) &AG[ip][itrc];
			x = cosf(phase2);
			y = sinf(phase2);
			ptr[0] = x;
			ptr[1] = y;

			ptr = (float *) &TS[itrc][ip];
			x = cosf(phase0);
			y = - sinf(phase0);
			ptr[0] = x;
			ptr[1] = y;

			ptr = (float *) &AS[ip][itrc];
			x = cosf(phase0);
			y = sinf(phase0);
			ptr[0] = x;
			ptr[1] = y;

			ptr = (float *) &TP[itrc][ip];
			x = cosf(phase1);
			y = - sinf(phase1);
			ptr[0] = x;
			ptr[1] = y;

			ptr = (float *) &TG[itrc][ip];
			x = cosf(phase2);
			y = - sinf(phase2);
			ptr[0] = x;
			ptr[1] = y;

			if (l1para->fpredatum == FPREDATUM_DGRG)
			{
				phase1 = twopif*MAX(5.0e-4*l1para->tminscl,fcos*reczou[itrc]*invwater);//half of ghost shift for reghosting
				phase2 = phase0+phase1; //mirror cable
				phase1 = phase0-phase1; //cable

				ptr = (float *) &TPR[itrc][ip];
				x = cosf(phase1);
				y = - sinf(phase1);
				ptr[0] = x;
				ptr[1] = y;

				ptr = (float *) &TGR[itrc][ip];
				x = cosf(phase2);
				y = - sinf(phase2);
				ptr[0] = x;
				ptr[1] = y;
			}
		} // for ( itrc = 0;itrc < nsgtrace; itrc++ )

		float tmp1, tmp2;

		// zxue
		if (l1para->srdg == SHOTRCVRJNT2D && l1para->computemode == 1)
		{
			// for a source ghost filter
			float tsrcghost = MAX (5.0e-4 * l1para->tminscl, fcos * shotz[0] * invwater);
			if (l1para->datum == CABLE90)
			{
				for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
				{
					tmp1 = l_f[ifreq] * 2.0f * tsrcghost;
					pSrcFilter[ip][ifreq] = (1.0f + l1para->srcrfl[ifreq] * complex<float> (cosf (tmp1), -sinf (tmp1))); 
				}
			}
			else
			{
				for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
				{
					tmp1 = l_f[ifreq] * (tsrcghost);
					tmp2 = l_f[ifreq] * (-tsrcghost);
					pSrcFilter[ip][ifreq] = (complex<float> (cosf (tmp1), sinf (tmp1)) + 
							l1para->srcrfl[ifreq] *  complex<float> (cosf (tmp2), sinf (tmp2)));
				}
			}
		} // if (l1para->computemode == 1)
	}
	//   }
	free (l_f);
}

void calc_matrix_A_T_p(l1inv_t *l1para, int nsgtrace, float *offsetx, float *offsety, float *recz, int lfid)
{
	int fftnr, npx, npy, nq, hpfl;
	int ip, ipx, ipy, iq, itrc;
	float fnyq, deltaf;
	float fsin, fcos, rq, px=0.0f, py=0.0f, pxy=0.0f;
	float phase0, phase1, phase2, invwater=1.0/l1para->vwater;
	int ifreq,isamp;

	float _Complex ** __restrict AP,** __restrict AG, ** __restrict AS;
	float _Complex ** __restrict TP,** __restrict TG,** __restrict TS;
	float _Complex ** __restrict TPR,** __restrict TGR;
	float _Complex ** __restrict AGS,** __restrict AGR, ** __restrict TGS;
	complex<float> ** pSrcFilter, ** pSrcFilterAz;

	float *ptr,x,y;

	float *reczou;
	float *shotz;

	if (l1para->fpredatum == FPREDATUM_DGRG) reczou = &recz[nsgtrace];
	if (l1para->choosemethod == CHOOSEMETHOD_JOINTSR3D || (l1para->choosemethod == CHOOSEMETHOD_MSDGI && l1para->jointmc == YESYES)) 
		shotz = &recz[nsgtrace];
	if (l1para->choosemethod == CHOOSEMETHOD_FP4D  && l1para->jointmc == YESYES)
	{
		if (l1para->msdataz == YESYES) 
			shotz = &recz[(l1para->ndata+1)*nsgtrace]; 
		else
			shotz = &recz[l1para->ndata*nsgtrace]; 
	}
	if (l1para->choosemethod == CHOOSEMETHOD_JOINTSR3D && l1para->flag_recz_smooth == 1) reczou = &recz[2*nsgtrace];

	if ((l1para->choosemethod == CHOOSEMETHOD_MSDGI||l1para->choosemethod == CHOOSEMETHOD_FP4D)&& 
			l1para->jointmc == YESYES && l1para->jointmethod != JOINTDEGHOST) 
		PFL_LENGTH(l1para->nsamp + l1para->nsampsrc,&fftnr); 
	else
		PFL_LENGTH(l1para->nsamp,&fftnr);          

	l1para->orgnr = fftnr;
	l1para->orgnc = fftnr/2+1;
	l1para->fftNR =        ALIGN_FLOAT(l1para->orgnr+2);     //for compatable with fftnr
	l1para->fftNC =        ALIGN_FCPLX(l1para->orgnc+1);     //for compatable with fftnc
	fftnr+=2;
	l1para->fftnr = fftnr;
	l1para->fftnc = fftnr/2;

	fnyq = 500000.0f/l1para->srate;
	deltaf = fnyq/(float)(l1para->fftnc-1);

	float twopif=deltaf*2.0f*3.1415927f;
	float *l_f = (float *) calloc (l1para->fftnc, sizeof (float));
	for (ifreq = 0; ifreq < l1para->fftnc; ifreq++) l_f[ifreq] = twopif * ifreq; 

	npx = l1para->npxsave;
	npy = l1para->npysave;
	nq  = l1para->nqsave;

	AP = (float _Complex **)l1para->AP_p;
	AG = (float _Complex **)l1para->AG_p;
	AS = (float _Complex **)l1para->AS_p;

	TP = (float _Complex **)l1para->TP_p;
	TG = (float _Complex **)l1para->TG_p;
	TS = (float _Complex **)l1para->TS_p;

	if (l1para->fpredatum == FPREDATUM_DGRG) 
	{
		TPR = (float _Complex **)l1para->TPR_p;
		TGR = (float _Complex **)l1para->TGR_p;
	}
	if (l1para->choosemethod == CHOOSEMETHOD_JOINTSR3D || 
			((l1para->choosemethod == CHOOSEMETHOD_MSDGI||l1para->choosemethod == CHOOSEMETHOD_FP4D) && l1para->jointmc == YESYES)) 
	{
		AGS = (float _Complex **)l1para->AGS_p;
		AGR = (float _Complex **)l1para->AGR_p;
		TGS = (float _Complex **)l1para->TGS_p;
		TGR = (float _Complex **)l1para->TGR_p;
		if (l1para->computemode == 1) 
			pSrcFilter = (complex<float> **) l1para->pSrcFilter;
	}

	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////

	if ((l1para->choosemethod == CHOOSEMETHOD_MSDGI||l1para->choosemethod == CHOOSEMETHOD_FP4D) 
			&& l1para->jointmc == YESYES && l1para->jointmethod != JOINTDEGHOST) 
	{
		pSrcFilterAz = (complex<float> **) l1para->pSrcFilterAz;

		head_t *src;
		float *gunvol = (float *) calloc (l1para->ntrcsrc, sizeof (float));
		float *psx    = (float *) calloc (l1para->ntrcsrc, sizeof (float));
		float *psy    = (float *) calloc (l1para->ntrcsrc, sizeof (float));
		float *psz    = (float *) calloc (l1para->ntrcsrc, sizeof (float));
		float shotx_avg, shoty_avg, shotz_avg, total_gunvol;
		float          **ptxsource = (float          **) flexible_array2d (l1para->ntrcsrc, l1para->fftnr, sizeof (float));
		complex<float> **pfxsource = (complex<float> **) flexible_array2d (l1para->ntrcsrc, l1para->fftnc, sizeof (complex<float>));

		int trlensig = NHEAD + l1para->nsampsrc;

		shotx_avg = shoty_avg = shotz_avg = total_gunvol = 0.0f;
		for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
		{
			src = (head_t *) (l1para->psig + itrc * trlensig);
			psx[itrc] = src->shotx;
			psy[itrc] = src->shoty;
			psz[itrc] = src->shotz;
			shotx_avg += src->shotx;
			shoty_avg += src->shoty;
			shotz_avg += src->shotz;
			memset (ptxsource[itrc], 0, l1para->fftnr * sizeof (float));
			memcpy (ptxsource[itrc], l1para->psig + itrc * trlensig + NHEAD,
					l1para->nsampsrc * sizeof (float));

			int tmp1 = (l1para->tzerosrc + l1para->tzerotgt) * 1000.0f * 1000.0f / l1para->srate;
			int tmp2 = tmp1 + l1para->tzerotgttaper * 1000.0f * 1000.0f / l1para->srate;

			//for (isamp = 0; isamp < l1para->nsampsrc; isamp++)
			//  gunvol[itrc] = MAX (gunvol[itrc], fabsf (ptxsource[itrc][isamp]));
			if (l1para->designorm == DESIGNORM_MAX || l1para->designorm == DESIGNORM_NONE)
			{
				for (isamp = 0; isamp < l1para->nsampsrc; isamp++)
					gunvol[itrc] = MAX (gunvol[itrc], fabsf (ptxsource[itrc][isamp]));
			}
			else if (l1para->designorm == DESIGNORM_RMS)
			{
				for (isamp = 0; isamp < l1para->nsampsrc; isamp++)
					gunvol[itrc] += ptxsource[itrc][isamp] * ptxsource[itrc][isamp];
				gunvol[itrc] = sqrtf (gunvol[itrc] / l1para->nsampsrc);
			}
			total_gunvol += gunvol[itrc];
		}
		shotx_avg /= l1para->ntrcsrc;
		shoty_avg /= l1para->ntrcsrc;
		shotz_avg /= l1para->ntrcsrc;

		l1para->total_gunvol = total_gunvol;

		if (total_gunvol > 0.0f && l1para->designorm != DESIGNORM_NONE)
		{
			for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
			{
				for (isamp = 0; isamp < l1para->nsampsrc; isamp++)
				{
					ptxsource[itrc][isamp] /= total_gunvol;
				}
			}
		}
		if (total_gunvol > 0.0f)
			for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
				gunvol[itrc] /= total_gunvol;

		PFL_GET (&hpfl, PFL_PORTABLE_DEF);

		FFT_R2C_PP (hpfl, ptxsource, l1para->ntrcsrc, pfxsource, l1para->fftnr); 
		for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
		{
			for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
			{
				pfxsource[itrc][ifreq] *= (l1para->fftnr - 2);
			}
		}

		ip = -1;
		for (ipy = 0; ipy < npy; ipy++)
		{
			for (ipx = 0; ipx < npx; ipx++)
			{
				ip++;
				px = l1para->pxsave[ip];
				py = l1para->pysave[ip];

				fcos=MAX(-0.99999f,MIN(0.99999f,py*l1para->vwater));       

				pxy=sqrtf(px*px+py*py);
				fsin=MAX(-0.99999f,MIN(0.99999f,pxy*l1para->vwater));       
				fcos=sqrtf(1.0f-fsin*fsin);//cosine of 3D surface angle

				// For source ghost operator
				for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
					pSrcFilter[ip][ifreq] = 0.0f;

				// zxue
				complex<float> ctmp1, ctmp2, ctmp3, ctmp4;
				if (l1para->datum == CABLE90)
					for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
					{
						src = (head_t *) (l1para->psig + itrc * (NHEAD + l1para->nsampsrc));
						float tsrcghost = MAX (5.0e-4 * l1para->tminscl, fcos * src->shotz * invwater);
						ctmp1 = complex<float> (cosf (twopif * 2.0f * tsrcghost), -sinf (twopif * 2.0f * tsrcghost));
						ctmp2 = complex<float> (1.0f, 0.0f);
						for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
						{
							pSrcFilter[ip][ifreq] += gunvol[itrc] * (1.0f + l1para->srcrfl[ifreq] * ctmp2); 
							ctmp2 *= ctmp1;
						}
					}
				else
					for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
					{
						src = (head_t *) (l1para->psig + itrc * (NHEAD + l1para->nsampsrc));
						float tsrcghost = MAX (5.0e-4 * l1para->tminscl, fcos * src->shotz * invwater);
						ctmp1 = complex<float> (cosf (twopif * ( tsrcghost - src->tdelay)), sinf (twopif * ( tsrcghost - src->tdelay)));
						ctmp2 = complex<float> (1.0f, 0.0f);
						ctmp3 = complex<float> (cosf (twopif * (-tsrcghost - src->tdelay)), sinf (twopif * (-tsrcghost - src->tdelay)));
						ctmp4 = complex<float> (1.0f, 0.0f);
						for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
						{
							pSrcFilter[ip][ifreq] += gunvol[itrc] * (ctmp2 + l1para->srcrfl[ifreq] * ctmp4);
							ctmp2 *= ctmp1;
							ctmp4 *= ctmp3;
						}
					}

				// for source signature operator
				for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
					pSrcFilterAz[ip][ifreq] = 0.0f;
				for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
				{
					head_t *src = (head_t *) (l1para->psig + itrc * (NHEAD + l1para->nsampsrc));
					float tdelay = -(psx[itrc] - shotx_avg) * px - (psy[itrc] - shoty_avg) * py - src->tdelay - l1para->tzerosrc;
					ctmp1 = complex<float> (cosf (twopif * tdelay), -sinf (twopif * tdelay));
					ctmp2 = complex<float> (1.0f, 0.0f);
					for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
					{
						pSrcFilterAz[ip][ifreq] += pfxsource[itrc][ifreq] * ctmp2;
						ctmp2 *= ctmp1;
					}
				}

				if (l1para->jointmethod == JOINTDESIG_SRC_DEGHOST || l1para->jointmethod == JOINTDESIG_SRCRCV_DEGHOST)
					for (ifreq = 0; ifreq < l1para->fftnc; ifreq++) 
						pSrcFilter[ip][ifreq] *= pSrcFilterAz[ip][ifreq];  
				else if (l1para->jointmethod == JOINTDESIGNATURE || l1para->jointmethod == JOINTDESIG_RCV_DEGHOST)
					for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
						pSrcFilter[ip][ifreq]  = pSrcFilterAz[ip][ifreq];
			} // for (ipx = 0; ipx < npx; ipx++)
		} // for (ipy = 0; ipy < npy; ipy++)

		free (gunvol); free (psx); free (psy); free (psz);
		flexible_free_array2d (ptxsource);
		flexible_free_array2d (pfxsource);
		PFL_FREE (hpfl);
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	int thfid=MIN(l1para->fftnc*0.85,l1para->highfreq/deltaf);
	int yftap=thfid*0.1;
	int ftap=thfid*0.1;

	float *hftaper = (float *)calloc(l1para->fftnc,sizeof(float));

	for (ifreq=0;ifreq<l1para->fftnc;ifreq++)
	{
		if (ifreq<thfid-ftap)
			hftaper[ifreq]=1.0f;
		else if (ifreq<thfid+ftap)
			hftaper[ifreq]=(1.0+cosf(0.5f*PI*(ifreq-thfid+ftap+1)/ftap))*0.5f;
		else
			hftaper[ifreq]=0.0f;
	}

	yftap=l1para->ylowcut*(l1para->ylowtap*0.01)/deltaf;

	int ylfid=MAX(yftap+1,l1para->ylowcut/deltaf);

	calc_zfilt(l1para);

	for (ifreq=0;ifreq<l1para->fftnc;ifreq++)
	{
		if (l1para->pmethod>0)
			l1para->zwtilt[ifreq]=twopif*MAX(lfid,ifreq)*l1para->zlowf[ifreq]*hftaper[ifreq];
		else if (l1para->zmethod==0)
			l1para->zwtilt[ifreq]=twopif*MAX(lfid,ifreq)*l1para->zlowf[ifreq];
		else
			l1para->zwtilt[ifreq]=l1para->zlowf[ifreq];

		if (l1para->ylowcut<0.001f||l1para->pmethod>0)
			l1para->ylowf[ifreq]=1.0f;
		else
		{
			if (ifreq<ylfid-yftap)
				l1para->ylowf[ifreq]=0.0f;
			else if (ifreq<ylfid+yftap)
				l1para->ylowf[ifreq]=(1.0+cosf(0.5f*PI*(ylfid+yftap-ifreq)/yftap))*0.5f;
			else
				l1para->ylowf[ifreq]=1.0f;
		}
		if (l1para->pmethod>0)
			l1para->ywtilt[ifreq]=twopif*MAX(lfid,ifreq)*l1para->ylowf[ifreq]*hftaper[ifreq];
		else if (l1para->ymethod==0)
			l1para->ywtilt[ifreq]=twopif*MAX(lfid,ifreq)*l1para->ylowf[ifreq];//*hftaper[ifreq];
		else
			l1para->ywtilt[ifreq]=l1para->ylowf[ifreq];
	}

	float epx,epy; //// effective px,py (slope)

	ip = -1;
	for ( iq = 0;iq < nq; iq++ )
	{
		for ( ipy = 0;ipy < npy; ipy++ )
		{
			for ( ipx = 0;ipx < npx; ipx++ )
			{

				ip = ip + 1;

				px = l1para->pxsave[ip];
				py = l1para->pysave[ip];
				rq = l1para->rqsave[ip];

				//     fcos=MAX(-0.99999f,MIN(0.99999f,py*l1para->vwater));       
				//     l1para->ytop[ip]=-fcos;

				//     pxy=sqrtf(px*px+py*py);
				//     fsin=MAX(-0.99999f,MIN(0.99999f,pxy*l1para->vwater));       
				//     fcos=sqrtf(1.0f-fsin*fsin);//cosine of 3D surface angle

				//     l1para->ztop[ip]=MAX(l1para->zcosmin,fcos);

				if (l1para->choosemethod != CHOOSEMETHOD_JOINTSR3D ||
						(l1para->choosemethod == CHOOSEMETHOD_JOINTSR3D && l1para->computemode == 1))
				{

					for ( itrc = 0;itrc < nsgtrace; itrc++ )
					{
						//       phase0 = twopif*px*offsetx[itrc]+twopif*py*offsety[itrc];

						phase0 = twopif*px*offsetx[itrc]+twopif*py*offsety[itrc]+
							+twopif*rq*offsetx[itrc]*offsetx[itrc]+
							+twopif*rq*offsety[itrc]*offsety[itrc];

						epx = (2.0*offsetx[itrc]*rq+px);
						epy = (2.0*offsety[itrc]*rq+py);

						pxy=sqrtf(epx*epx+epy*epy);

						fcos=MAX(-0.99999f,MIN(0.99999f,epy*l1para->vwater));       
						l1para->ytop[ip]=-fcos;

						fsin=MAX(-0.99999f,MIN(0.99999f,pxy*l1para->vwater));       
						fcos=sqrtf(1.0f-fsin*fsin);//cosine of 3D surface angle
						l1para->ztop[ip]=MAX(l1para->zcosmin,fcos);

						phase1 = twopif*MAX(5.0e-4*l1para->tminscl,fcos*recz[itrc]*invwater);//half of ghost shift for reghosting
						phase2 = phase0+phase1; //mirror cable
						if (l1para->flag_recz_smooth==0)
							phase1 = phase0-phase1; //cable
						else
						{
							phase1 = twopif*MAX(5.0e-4*l1para->tminscl,fcos*reczou[itrc]*invwater);//half of primary shift for reghosting
							phase1 = phase0-phase1; //cable
						}

						ptr = (float *) &AP[ip][itrc];
						x = cosf(phase1);
						y = sinf(phase1);
						ptr[0] = x;
						ptr[1] = y;

						ptr = (float *) &AG[ip][itrc];
						x = cosf(phase2);
						y = sinf(phase2);
						ptr[0] = x;
						ptr[1] = y;

						ptr = (float *) &TS[itrc][ip];
						x = cosf(phase0);
						y = - sinf(phase0);
						ptr[0] = x;
						ptr[1] = y;

						ptr = (float *) &AS[ip][itrc];
						x = cosf(phase0);
						y = sinf(phase0);
						ptr[0] = x;
						ptr[1] = y;

						ptr = (float *) &TP[itrc][ip];
						x = cosf(phase1);
						y = - sinf(phase1);
						ptr[0] = x;
						ptr[1] = y;

						ptr = (float *) &TG[itrc][ip];
						x = cosf(phase2);
						y = - sinf(phase2);
						ptr[0] = x;
						ptr[1] = y;

						if (l1para->fpredatum == FPREDATUM_DGRG)
						{
							phase1 = twopif*MAX(5.0e-4*l1para->tminscl,fcos*reczou[itrc]*invwater);//half of ghost shift for reghosting
							phase2 = phase0+phase1; //mirror cable
							phase1 = phase0-phase1; //cable

							ptr = (float *) &TPR[itrc][ip];
							x = cosf(phase1);
							y = - sinf(phase1);
							ptr[0] = x;
							ptr[1] = y;

							ptr = (float *) &TGR[itrc][ip];
							x = cosf(phase2);
							y = - sinf(phase2);
							ptr[0] = x;
							ptr[1] = y;
						}

					} // for ( itrc = 0;itrc < nsgtrace; itrc++ )

					// zxue
					//single-sensor or ms joint deghosting + no designature
					// For jntsr, its deghosting + designature function was different from this one; this one only contains jonit deghosting.
					// For msjntsr and fp4d, its deghosting + designature was in the above; following was only joint deghosting.
					if (l1para->choosemethod == CHOOSEMETHOD_JOINTSR3D && l1para->computemode == 1 || 
							((l1para->choosemethod == CHOOSEMETHOD_MSDGI||l1para->choosemethod == CHOOSEMETHOD_FP4D) 
							 && l1para->jointmc == YESYES && l1para->jointmethod == JOINTDEGHOST))
					{

						// for a source ghost filter
						complex<float> ctmp1, ctmp2, ctmp3, ctmp4;
						float tsrcghost = MAX (5.0e-4 * l1para->tminscl, fcos * shotz[0] * invwater);

						if (l1para->datum == CABLE90)
						{
							ctmp1 = complex<float> (cosf (twopif * 2.0f * tsrcghost), -sinf (twopif * 2.0f * tsrcghost));
							ctmp2 = complex<float> (1.0f, 0.0f);

							for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
							{
								pSrcFilter[ip][ifreq] = (1.0f + l1para->srcrfl[ifreq] * ctmp2); 
								ctmp2 *= ctmp1;
							}
						}
						else
						{
							ctmp1 = complex<float> (cosf (twopif * tsrcghost), sinf (twopif * tsrcghost));
							ctmp2 = complex<float> (1.0f, 0.0f);
							ctmp3 = complex<float> (cosf (twopif * -tsrcghost), sinf (twopif * -tsrcghost));
							ctmp4 = complex<float> (1.0f, 0.0f);

							for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
							{
								pSrcFilter[ip][ifreq] = (ctmp2 + l1para->srcrfl[ifreq] * ctmp4);
								ctmp2 *= ctmp1;
								ctmp4 *= ctmp3;
							}
						}
					} // if (l1para->computemode == 1)
				} // if (l1para->choosemethod != CHOOSEMETHOD_JOINTSR3D)
				else //single sensor + ms sequential deghosting + ms seq/joint deghosting + designature
				{
					for ( itrc = 0;itrc < nsgtrace; itrc++ )
					{
						//       phase0 = twopif*px*offsetx[itrc]+twopif*py*offsety[itrc];

						phase0 = twopif*px*offsetx[itrc]+twopif*py*offsety[itrc]+
							+twopif*rq*offsetx[itrc]*offsetx[itrc]+
							+twopif*rq*offsety[itrc]*offsety[itrc];

						epx = (2.0*offsetx[itrc]*rq+px);
						epy = (2.0*offsety[itrc]*rq+py);

						pxy=sqrtf(epx*epx+epy*epy);

						fcos=MAX(-0.99999f,MIN(0.99999f,epy*l1para->vwater));       
						l1para->ytop[ip]=-fcos;

						fsin=MAX(-0.99999f,MIN(0.99999f,pxy*l1para->vwater));       
						fcos=sqrtf(1.0f-fsin*fsin);//cosine of 3D surface angle
						l1para->ztop[ip]=MAX(l1para->zcosmin,fcos);

						phase1 = twopif*MAX(5.0e-4*l1para->tminscl,fcos*recz[itrc]*invwater);//half of ghost shift for reghosting
						phase2 = twopif*MAX(5.0e-4*l1para->tminscl,fcos*shotz[itrc]*invwater);

						// shot depth to receiver depth
						x = cosf(phase0 - phase1 - phase2);
						y = sinf(phase0 - phase1 - phase2);
						ptr = (float *) &AP[ip][itrc];
						ptr[0] =  x; ptr[1] =  y;
						ptr = (float *) &TP[itrc][ip];
						ptr[0] =  x; ptr[1] = -y;

						// mirror shot to mirror receiver
						x = cosf(phase0 + phase1 + phase2);
						y = sinf(phase0 + phase1 + phase2);
						ptr = (float *) &AG[ip][itrc];
						ptr[0] =  x; ptr[1] =  y;
						ptr = (float *) &TG[itrc][ip];
						ptr[0] =  x; ptr[1] = -y;

						// shot to mirror receiver
						x = cosf(phase0 + phase1 - phase2);
						y = sinf(phase0 + phase1 - phase2);
						ptr = (float *) &AGR[ip][itrc];
						ptr[0] =  x; ptr[1] =  y;
						ptr = (float *) &TGR[itrc][ip];
						ptr[0] =  x; ptr[1] = -y;

						// mirror shot to receiver
						x = cosf(phase0 - phase1 + phase2);
						y = sinf(phase0 - phase1 + phase2);
						ptr = (float *) &AGS[ip][itrc];
						ptr[0] =  x; ptr[1] =  y;
						ptr = (float *) &TGS[itrc][ip];
						ptr[0] =  x; ptr[1] = -y;

						// surface to surface
						x = cosf(phase0);
						y = sinf(phase0);
						ptr = (float *) &AS[ip][itrc];
						ptr[0] =  x; ptr[1] =  y;
						ptr = (float *) &TS[itrc][ip];
						ptr[0] =  x; ptr[1] = -y;

					} // for ( itrc = 0;itrc < nsgtrace; itrc++ )
				} // if (l1para->choosemethod != CHOOSEMETHOD_JOINTSR3D)

			}
		}
	}

	free(hftaper);
	free(l_f);
}

void calc_matrix_A_T_p_joint(l1inv_t *l1para, int nsgtrace, float *offsetx, float *offsety, float *recz, int lfid) //// jointmethod!=JOINTDEGHOST
{
	int fftnr, npx, npy, ip, ipx, ipy, itrc, isamp, ifreq, hpfl;
	float fnyq, deltaf, fsin, fcos, px=0.0f, py=0.0f, pxy=0.0f;
	float phase0, phase1, phase2, invwater=1.0/l1para->vwater, x,y;

	complex<float> ** __restrict AP,** __restrict AG, ** __restrict AS;
	complex<float> ** __restrict TP,** __restrict TG,** __restrict TS;
	complex<float> ** pSrcFilter, ** pSrcFilterAz;
	complex<float> ** p1dSrc = (complex<float> **) l1para->p1dSrc;
	complex<float> ** p1dTgt = (complex<float> **) l1para->p1dTgt;

	//PFL_LENGTH(l1para->nsamp,&fftnr);
	PFL_GET (&hpfl, PFL_PORTABLE_DEF);
	PFL_LENGTH(l1para->nsamp + l1para->nsampsrc,&fftnr);
	l1para->orgnr = fftnr;
	l1para->orgnc = fftnr/2+1;
	l1para->fftNR =        ALIGN_FLOAT(l1para->orgnr+2);     //for compatable with fftnr
	l1para->fftNC =        ALIGN_FCPLX(l1para->orgnc+1);     //for compatable with fftnc
	fftnr+=2;
	l1para->fftnr = fftnr;
	l1para->fftnc = fftnr/2;
	fnyq = 500000.0f/l1para->srate;
	deltaf = fnyq/(float)(l1para->fftnc-1);

	float twopif=deltaf*2.0f*3.1415927f;
	float *l_f = (float *) calloc (l1para->fftnc, sizeof (float));
	for (ifreq = 0; ifreq < l1para->fftnc; ifreq++) l_f[ifreq] = twopif * ifreq; 

	npx = l1para->npxsave;
	npy = l1para->npysave;

	AP = (complex<float> **)l1para->AP_p;
	AG = (complex<float> **)l1para->AG_p;
	AS = (complex<float> **)l1para->AS_p;

	TP = (complex<float> **)l1para->TP_p;
	TG = (complex<float> **)l1para->TG_p;
	TS = (complex<float> **)l1para->TS_p;

	if (l1para->choosemethod == CHOOSEMETHOD_JOINTSR3D || l1para->jointmc == YESYES) //// ** SRAY ** CHECK ** : This should be redundant 
	{
		pSrcFilter   = (complex<float> **) l1para->pSrcFilter;
		pSrcFilterAz = (complex<float> **) l1para->pSrcFilterAz; 
	}

	head_t *src;
	float *gunvol = (float *) calloc (l1para->ntrcsrc, sizeof (float));
	float *psx    = (float *) calloc (l1para->ntrcsrc, sizeof (float));
	float *psy    = (float *) calloc (l1para->ntrcsrc, sizeof (float));
	float *psz    = (float *) calloc (l1para->ntrcsrc, sizeof (float));
	float shotx_avg, shoty_avg, shotz_avg, total_gunvol;
	float          **ptxsource = (float          **) flexible_array2d (l1para->ntrcsrc, l1para->fftnr, sizeof (float));
	float          **ptxtarget = (float          **) flexible_array2d (l1para->ntrcsrc, l1para->fftnr, sizeof (float));
	complex<float> **pfxsource = (complex<float> **) flexible_array2d (l1para->ntrcsrc, l1para->fftnc, sizeof (complex<float>));
	complex<float> **pfxtarget = (complex<float> **) flexible_array2d (l1para->ntrcsrc, l1para->fftnc, sizeof (complex<float>));
	int trlensig = NHEAD + l1para->nsampsrc;

	shotx_avg = shoty_avg = shotz_avg = total_gunvol = 0.0f;
	for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
	{
		src = (head_t *) (l1para->psig + itrc * trlensig);
		psx[itrc] = src->shotx;
		psy[itrc] = src->shoty;
		psz[itrc] = src->shotz;
		shotx_avg += src->shotx;
		shoty_avg += src->shoty;
		shotz_avg += src->shotz;
		memset (ptxsource[itrc], 0, l1para->fftnr * sizeof (float));
		memset (ptxtarget[itrc], 0, l1para->fftnr * sizeof (float));
		memcpy (ptxsource[itrc], l1para->psig + itrc * trlensig + NHEAD,
				l1para->nsampsrc * sizeof (float));
		//    memcpy (ptxtarget[itrc], l1para->psig + itrc * trlensig + NHEAD,
		//      (l1para->tzerosrc + l1para->tzerotgt) * 1000.0f * 1000.0f / l1para->srate * sizeof (float));
		int tmp1 = (l1para->tzerosrc + l1para->tzerotgt) * 1000.0f * 1000.0f / l1para->srate;
		int tmp2 = tmp1 + l1para->tzerotgttaper * 1000.0f * 1000.0f / l1para->srate;
		memcpy (ptxtarget[itrc], l1para->psig + itrc * trlensig + NHEAD, tmp2 * sizeof (float));
		for (isamp = tmp1; isamp < tmp2; isamp++)
			ptxtarget[itrc][isamp] *= powf (cosf (0.5f * 3.14159f * (isamp - tmp1 + 1) / (tmp2 - tmp1 + 1)), 2.0f);

		//for (isamp = 0; isamp < l1para->nsampsrc; isamp++)
		//  gunvol[itrc] = MAX (gunvol[itrc], fabsf (ptxsource[itrc][isamp]));
		if (l1para->designorm == DESIGNORM_MAX || l1para->designorm == DESIGNORM_NONE)
		{
			for (isamp = 0; isamp < l1para->nsampsrc; isamp++)
				gunvol[itrc] = MAX (gunvol[itrc], fabsf (ptxsource[itrc][isamp]));
		}
		else if (l1para->designorm == DESIGNORM_RMS)
		{
			for (isamp = 0; isamp < l1para->nsampsrc; isamp++)
				gunvol[itrc] += ptxsource[itrc][isamp] * ptxsource[itrc][isamp];
			gunvol[itrc] = sqrtf (gunvol[itrc] / l1para->nsampsrc);
		}
		total_gunvol += gunvol[itrc];
	}
	shotx_avg /= l1para->ntrcsrc;
	shoty_avg /= l1para->ntrcsrc;
	shotz_avg /= l1para->ntrcsrc;

	l1para->total_gunvol = total_gunvol;
	
	if (total_gunvol > 0.0f && l1para->designorm != DESIGNORM_NONE)
	{
		for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
		{
			for (isamp = 0; isamp < l1para->nsampsrc; isamp++)
			{
				ptxsource[itrc][isamp] /= total_gunvol;
				ptxtarget[itrc][isamp] /= total_gunvol;
			}
		}
	}
	if (total_gunvol > 0.0f)
		for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
			gunvol[itrc] /= total_gunvol;

	FFT_R2C_PP (hpfl, ptxsource, l1para->ntrcsrc, pfxsource, l1para->fftnr); 
	FFT_R2C_PP (hpfl, ptxtarget, l1para->ntrcsrc, pfxtarget, l1para->fftnr); 
	for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
	{
		for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
		{
			pfxsource[itrc][ifreq] *= (l1para->fftnr - 2);
			pfxtarget[itrc][ifreq] *= (l1para->fftnr - 2);
		}
	}

	ip = -1;
	for (ipy = 0; ipy < npy; ipy++)
	{
		for (ipx = 0; ipx < npx; ipx++)
		{
			ip++;
			px = l1para->pxsave[ip];
			py = l1para->pysave[ip];

			fcos=MAX(-0.99999f,MIN(0.99999f,py*l1para->vwater));       
			l1para->ytop[ip]=-fcos;

			pxy=sqrtf(px*px+py*py);
			fsin=MAX(-0.99999f,MIN(0.99999f,pxy*l1para->vwater));       
			fcos=sqrtf(1.0f-fsin*fsin);//cosine of 3D surface angle

			l1para->ztop[ip]=MAX(l1para->zcosmin,fcos);

			// For reverse tau-p operator and receiver ghost operator
			for (itrc = 0; itrc < nsgtrace; itrc++)
			{
				phase0 = twopif * px * offsetx[itrc] + twopif * py * offsety[itrc];

				phase1 = twopif * MAX (5.0e-4 * l1para->tminscl, fcos * recz[itrc] * invwater);//half of ghost shift for reghosting
				phase2 = phase0 + phase1; //mirror cable
				phase1 = phase0 - phase1; //cable

				x = cosf(phase1); y = sinf(phase1); // for cable depth
				AP[ip][itrc] = complex<float> (x,  y);
				TP[itrc][ip] = complex<float> (x, -y);

				x = cosf(phase2); y = sinf(phase2); // for mirror cable depth
				AG[ip][itrc] = complex<float> (x,  y);
				TG[itrc][ip] = complex<float> (x, -y);

				x = cosf(phase0); y = sinf(phase0); // for reverse taup
				AS[ip][itrc] = complex<float> (x,  y);
				TS[itrc][ip] = complex<float> (x, -y);
			} // for ( itrc = 0;itrc < nsgtrace; itrc++ )

			// zxue
			// For source ghost operator
			for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
				pSrcFilter[ip][ifreq] = 0.0f;

			complex<float> ctmp1, ctmp2, ctmp3, ctmp4;
			if (l1para->datum == CABLE90)
				for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
				{
					src = (head_t *) (l1para->psig + itrc * (NHEAD + l1para->nsampsrc));
					float tsrcghost = MAX (5.0e-4 * l1para->tminscl, fcos * src->shotz * invwater);
					ctmp1 = complex<float> (cosf (twopif * 2.0f * tsrcghost), -sinf (twopif * 2.0f * tsrcghost));
					ctmp2 = complex<float> (1.0f, 0.0f);
					for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
					{
						pSrcFilter[ip][ifreq] += gunvol[itrc] * (1.0f + l1para->srcrfl[ifreq] * ctmp2); 
						ctmp2 *= ctmp1;
					}
				}
			else
				for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
				{
					src = (head_t *) (l1para->psig + itrc * (NHEAD + l1para->nsampsrc));
					float tsrcghost = MAX (5.0e-4 * l1para->tminscl, fcos * src->shotz * invwater);
					ctmp1 = complex<float> (cosf (twopif * ( tsrcghost - src->tdelay)), sinf (twopif * ( tsrcghost - src->tdelay)));
					ctmp2 = complex<float> (1.0f, 0.0f);
					ctmp3 = complex<float> (cosf (twopif * (-tsrcghost - src->tdelay)), sinf (twopif * (-tsrcghost - src->tdelay)));
					ctmp4 = complex<float> (1.0f, 0.0f);
					for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
					{
						pSrcFilter[ip][ifreq] += gunvol[itrc] * (ctmp2 + l1para->srcrfl[ifreq] * ctmp4);
						ctmp2 *= ctmp1;
						ctmp4 *= ctmp3;
					}
				}

			// for source signature operator
			for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
				pSrcFilterAz[ip][ifreq] = 0.0f;
			for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
			{
				head_t *src = (head_t *) (l1para->psig + itrc * (NHEAD + l1para->nsampsrc));
				float tdelay = -(psx[itrc] - shotx_avg) * px - (psy[itrc] - shoty_avg) * py - src->tdelay - l1para->tzerosrc;
				ctmp1 = complex<float> (cosf (twopif * tdelay), -sinf (twopif * tdelay));
				ctmp2 = complex<float> (1.0f, 0.0f);
				for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
				{
					pSrcFilterAz[ip][ifreq] += pfxsource[itrc][ifreq] * ctmp2;
					ctmp2 *= ctmp1;
				}
			}

			if (l1para->jointmethod == JOINTDESIG_SRC_DEGHOST || l1para->jointmethod == JOINTDESIG_SRCRCV_DEGHOST)
				for (ifreq = 0; ifreq < l1para->fftnc; ifreq++) 
					pSrcFilter[ip][ifreq] *= pSrcFilterAz[ip][ifreq];  
			else if (l1para->jointmethod == JOINTDESIGNATURE || l1para->jointmethod == JOINTDESIG_RCV_DEGHOST)
				for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
					pSrcFilter[ip][ifreq]  = pSrcFilterAz[ip][ifreq];
		} // for (ipx = 0; ipx < npx; ipx++)
	} // for (ipy = 0; ipy < npy; ipy++)

	if (l1para->lmerge1d == YESYES)
	{
		float tmp1, tmp2;
		complex<float> ctmp1, ctmp2, ctmp3, ctmp4;
		for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
			p1dSrc[0][ifreq] = p1dTgt[0][ifreq] = 0.0f;

		for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
		{
			head_t *src = (head_t *) (l1para->psig + itrc * (NHEAD + l1para->nsampsrc));
			float tdelay = - src->tdelay - l1para->tzerosrc;
			ctmp1 = complex<float> (cosf (twopif * (-src->tdelay - l1para->tzerosrc)), -sinf (twopif * (-src->tdelay - l1para->tzerosrc)));
			ctmp2 = complex<float> (1.0f, 0.0f);
			for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
			{
				p1dSrc[0][ifreq] += pfxsource[itrc][ifreq] * ctmp2;
				ctmp2 *= ctmp1;
			}
		}
		// zxue
		if (l1para->jointmethod == JOINTDESIG_SRC_DEGHOST || l1para->jointmethod == JOINTDESIG_SRCRCV_DEGHOST)
		{
			complex<float> **p1dTmp = (complex<float> **) flexible_array2d (1, l1para->fftnc, sizeof (complex<float>));

			if (l1para->datum == CABLE90)
			{
				for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
				{
					head_t *src = (head_t *) (l1para->psig + itrc * (NHEAD + l1para->nsampsrc));
					float tsrcghost = src->shotz * invwater;
					ctmp1 = complex<float> (cosf (twopif * 2.0f * tsrcghost), -sinf (twopif * 2.0f * tsrcghost));
					ctmp2 = complex<float> (1.0f, 0.0f);
					for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
					{
						p1dTmp[0][ifreq] += gunvol[itrc] * (1.0f + l1para->srcrfl[ifreq] * ctmp2);
						ctmp2 *= ctmp1;
					}
				}
			}
			else
			{
				for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
				{
					head_t *src = (head_t *) (l1para->psig + itrc * (NHEAD + l1para->nsampsrc));
					float tsrcghost = src->shotz * invwater;
					ctmp1 = complex<float> (cosf (twopif * ( tsrcghost - src->tdelay)), sinf (twopif * ( tsrcghost - src->tdelay)));
					ctmp2 = complex<float> (1.0f, 0.0f);
					ctmp3 = complex<float> (cosf (twopif * (-tsrcghost - src->tdelay)), sinf (twopif * (-tsrcghost - src->tdelay)));
					ctmp2 = complex<float> (1.0f, 0.0f);
					for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
					{
						p1dTmp[0][ifreq] += gunvol[itrc] * (ctmp2 + l1para->srcrfl[ifreq] * ctmp4);
						ctmp2 *= ctmp1;
						ctmp4 *= ctmp3;
					}
				}
			}

			for (ifreq = 0; ifreq < l1para->fftnc; ifreq++) 
				p1dSrc[0][ifreq] *= p1dTmp[0][ifreq];  

			flexible_free_array2d (p1dTmp);
		}
		for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
		{
			head_t *src = (head_t *) (l1para->psig + itrc * (NHEAD + l1para->nsampsrc));
			float tdelay = - src->tdelay - l1para->tzerosrc;
			ctmp1 = complex<float> (cosf (twopif * (- src->tdelay - l1para->tzerosrc)), -sinf (twopif * (- src->tdelay - l1para->tzerosrc)));
			ctmp2 = complex<float> (1.0f, 0.0f);
			for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
			{
				p1dTgt[0][ifreq] += pfxtarget[itrc][ifreq] * ctmp2;
				ctmp2 *= ctmp1;
			}
		}
	}

	free(l_f);
	free (gunvol); free (psx); free (psy); free (psz);
	flexible_free_array2d (ptxsource);
	flexible_free_array2d (pfxsource);
	flexible_free_array2d (ptxtarget);
	flexible_free_array2d (pfxtarget);
	PFL_FREE (hpfl);
}

void recalc_matrix_A_T_p(l1inv_t *l1para, int nsgtrace, float *offsetx, float *offsety, float *recz, 
		int *ip_sorted, int np_reduce, int lfid)
{
	int ip, ipnew, itrc, hpfl;
	float fnyq, deltaf;
	float fsin, fcos, rq, px=0.0f, py=0.0f, pxy=0.0f;
	float phase0, phase1, phase2, invwater=1.0/l1para->vwater;
	int ifreq,isamp;

	float _Complex ** __restrict AP,** __restrict AG, ** __restrict AS;
	float _Complex ** __restrict TP,** __restrict TG,** __restrict TS;
	float _Complex ** __restrict TPR,** __restrict TGR;
	float _Complex ** __restrict AGS,** __restrict AGR, ** __restrict TGS;
	complex<float> ** pSrcFilter, ** pSrcFilterAz, ** pTgtFilter;

	float *ptr,x,y;

	float *reczou;
	float *shotz;

	if (l1para->fpredatum == FPREDATUM_DGRG) reczou = &recz[nsgtrace];
	if (l1para->choosemethod == CHOOSEMETHOD_JOINTSR3D || (l1para->choosemethod == CHOOSEMETHOD_MSDGI && l1para->jointmc == YESYES)) 
		shotz = &recz[nsgtrace];
	if (l1para->choosemethod == CHOOSEMETHOD_FP4D && l1para->jointmc == YESYES)
		if (l1para->msdataz == YESYES) 
			shotz = &recz[(l1para->ndata+1)*nsgtrace]; 
		else
			shotz = &recz[l1para->ndata*nsgtrace]; 
	if (l1para->choosemethod == CHOOSEMETHOD_JOINTSR3D && l1para->flag_recz_smooth == 1) reczou = &recz[2*nsgtrace];

	fnyq = 500000.0f/l1para->srate;
	deltaf = fnyq/(float)(l1para->fftnc-1);

	float twopif=deltaf*2.0f*3.1415927f;
	float *l_f = (float *) calloc (l1para->fftnc, sizeof (float));
	for (ifreq = 0; ifreq < l1para->fftnc; ifreq++) l_f[ifreq] = twopif * ifreq; 

	AP = (float _Complex **)l1para->AP_p;
	AG = (float _Complex **)l1para->AG_p;
	AS = (float _Complex **)l1para->AS_p;

	TP = (float _Complex **)l1para->TP_p;
	TG = (float _Complex **)l1para->TG_p;
	TS = (float _Complex **)l1para->TS_p;

	if (l1para->fpredatum == FPREDATUM_DGRG) 
	{
		TPR = (float _Complex **)l1para->TPR_p;
		TGR = (float _Complex **)l1para->TGR_p;
	}
	if (l1para->choosemethod == CHOOSEMETHOD_JOINTSR3D || 
			((l1para->choosemethod == CHOOSEMETHOD_MSDGI||l1para->choosemethod == CHOOSEMETHOD_FP4D) && l1para->jointmc == YESYES)) 
	{
		AGS = (float _Complex **)l1para->AGS_p;
		AGR = (float _Complex **)l1para->AGR_p;
		TGS = (float _Complex **)l1para->TGS_p;
		TGR = (float _Complex **)l1para->TGR_p;

		if (l1para->computemode == 1) 
			pSrcFilter = (complex<float> **) l1para->pSrcFilter;
	}

	////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////
	if ((l1para->choosemethod == CHOOSEMETHOD_MSDGI||l1para->choosemethod == CHOOSEMETHOD_FP4D)
			&& l1para->jointmc == YESYES && l1para->jointmethod != JOINTDEGHOST) 
	{
		pSrcFilterAz = (complex<float> **) l1para->pSrcFilterAz;
		pTgtFilter   = (complex<float> **) l1para->pTgtFilter;

		head_t *src;
		float *gunvol = (float *) calloc (l1para->ntrcsrc, sizeof (float));
		float *psx    = (float *) calloc (l1para->ntrcsrc, sizeof (float));
		float *psy    = (float *) calloc (l1para->ntrcsrc, sizeof (float));
		float *psz    = (float *) calloc (l1para->ntrcsrc, sizeof (float));
		float shotx_avg, shoty_avg, shotz_avg, total_gunvol;
		float          **ptxsource = (float          **) flexible_array2d (l1para->ntrcsrc, l1para->fftnr, sizeof (float));
		float          **ptxtarget = (float          **) flexible_array2d (l1para->ntrcsrc, l1para->fftnr, sizeof (float));
		complex<float> **pfxsource = (complex<float> **) flexible_array2d (l1para->ntrcsrc, l1para->fftnc, sizeof (complex<float>));
		complex<float> **pfxtarget = (complex<float> **) flexible_array2d (l1para->ntrcsrc, l1para->fftnc, sizeof (complex<float>));
		int trlensig = NHEAD + l1para->nsampsrc;

		shotx_avg = shoty_avg = shotz_avg = total_gunvol = 0.0f;
		for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
		{
			src = (head_t *) (l1para->psig + itrc * trlensig);
			psx[itrc] = src->shotx;
			psy[itrc] = src->shoty;
			psz[itrc] = src->shotz;
			shotx_avg += src->shotx;
			shoty_avg += src->shoty;
			shotz_avg += src->shotz;
			memset (ptxsource[itrc], 0, l1para->fftnr * sizeof (float));
			memset (ptxtarget[itrc], 0, l1para->fftnr * sizeof (float));
			memcpy (ptxsource[itrc], l1para->psig + itrc * trlensig + NHEAD,
					l1para->nsampsrc * sizeof (float));
			//    memcpy (ptxtarget[itrc], l1para->psig + itrc * trlensig + NHEAD,
			//      (l1para->tzerosrc + l1para->tzerotgt) * 1000.0f * 1000.0f / l1para->srate * sizeof (float));
			int tmp1 = (l1para->tzerosrc + l1para->tzerotgt) * 1000.0f * 1000.0f / l1para->srate;
			int tmp2 = tmp1 + l1para->tzerotgttaper * 1000.0f * 1000.0f / l1para->srate;
			memcpy (ptxtarget[itrc], l1para->psig + itrc * trlensig + NHEAD, tmp2 * sizeof (float));
			for (isamp = tmp1; isamp < tmp2; isamp++)
				ptxtarget[itrc][isamp] *= powf (cosf (0.5f * 3.14159f * (isamp - tmp1 + 1) / (tmp2 - tmp1 + 1)), 2.0f);

			//for (isamp = 0; isamp < l1para->nsampsrc; isamp++)
			//  gunvol[itrc] = MAX (gunvol[itrc], fabsf (ptxsource[itrc][isamp]));
			if (l1para->designorm == DESIGNORM_MAX || l1para->designorm == DESIGNORM_NONE)
			{
				for (isamp = 0; isamp < l1para->nsampsrc; isamp++)
					gunvol[itrc] = MAX (gunvol[itrc], fabsf (ptxsource[itrc][isamp]));
			}
			else if (l1para->designorm == DESIGNORM_RMS)
			{
				for (isamp = 0; isamp < l1para->nsampsrc; isamp++)
					gunvol[itrc] += ptxsource[itrc][isamp] * ptxsource[itrc][isamp];
				gunvol[itrc] = sqrtf (gunvol[itrc] / l1para->nsampsrc);
			}
			total_gunvol += gunvol[itrc];
		}
		shotx_avg /= l1para->ntrcsrc;
		shoty_avg /= l1para->ntrcsrc;
		shotz_avg /= l1para->ntrcsrc;
		if (total_gunvol > 0.0f && l1para->designorm != DESIGNORM_NONE)
		{
			for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
			{
				for (isamp = 0; isamp < l1para->nsampsrc; isamp++)
				{
					ptxsource[itrc][isamp] /= total_gunvol;
					ptxtarget[itrc][isamp] /= total_gunvol;
				}
			}
		}
		if (total_gunvol > 0.0f)
			for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
				gunvol[itrc] /= total_gunvol;

		PFL_GET (&hpfl, PFL_PORTABLE_DEF);

		FFT_R2C_PP (hpfl, ptxsource, l1para->ntrcsrc, pfxsource, l1para->fftnr); 
		FFT_R2C_PP (hpfl, ptxtarget, l1para->ntrcsrc, pfxtarget, l1para->fftnr); 
		for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
		{
			for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
			{
				pfxsource[itrc][ifreq] *= (l1para->fftnr - 2);
				pfxtarget[itrc][ifreq] *= (l1para->fftnr - 2);
			}
		}

		for (ipnew = 0; ipnew < np_reduce; ipnew++)
		{
			ip = ip_sorted[ipnew];

			px = l1para->pxsave[ip];
			py = l1para->pysave[ip];

			fcos=MAX(-0.99999f,MIN(0.99999f,py*l1para->vwater));       

			pxy=sqrtf(px*px+py*py);
			fsin=MAX(-0.99999f,MIN(0.99999f,pxy*l1para->vwater));       
			fcos=sqrtf(1.0f-fsin*fsin);//cosine of 3D surface angle

			// For source ghost operator
			for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
				pSrcFilter[ipnew][ifreq] = 0.0f;

			// zxue
			complex<float> ctmp1, ctmp2, ctmp3, ctmp4;
			if (l1para->datum == CABLE90)
				for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
				{
					src = (head_t *) (l1para->psig + itrc * (NHEAD + l1para->nsampsrc));
					float tsrcghost = MAX (5.0e-4 * l1para->tminscl, fcos * src->shotz * invwater);
					ctmp1 = complex<float> (cosf (twopif * 2.0f * tsrcghost), -sinf (twopif * 2.0f * tsrcghost));
					ctmp2 = complex<float> (1.0f, 0.0f);
					for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
					{
						pSrcFilter[ipnew][ifreq] += gunvol[itrc] * (1.0f + l1para->srcrfl[ifreq] * ctmp2); 
						ctmp2 *= ctmp1;
					}
				}
			else
				for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
				{
					src = (head_t *) (l1para->psig + itrc * (NHEAD + l1para->nsampsrc));
					float tsrcghost = MAX (5.0e-4 * l1para->tminscl, fcos * src->shotz * invwater);
					ctmp1 = complex<float> (cosf (twopif * ( tsrcghost - src->tdelay)), sinf (twopif * ( tsrcghost - src->tdelay)));
					ctmp2 = complex<float> (1.0f, 0.0f);
					ctmp3 = complex<float> (cosf (twopif * (-tsrcghost - src->tdelay)), sinf (twopif * (-tsrcghost - src->tdelay)));
					ctmp4 = complex<float> (1.0f, 0.0f);
					for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
					{
						pSrcFilter[ipnew][ifreq] += gunvol[itrc] * (ctmp2 + l1para->srcrfl[ifreq] * ctmp4);
						ctmp2 *= ctmp1;
						ctmp4 *= ctmp3;
					}
				}

			// for source signature operator
			for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
			{
				pSrcFilterAz[ipnew][ifreq] = 0.0f;
				pTgtFilter[ipnew][ifreq]   = 0.0f;
			}
			for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
			{
				head_t *src = (head_t *) (l1para->psig + itrc * (NHEAD + l1para->nsampsrc));
				float tdelay = -(psx[itrc] - shotx_avg) * px - (psy[itrc] - shoty_avg) * py - src->tdelay - l1para->tzerosrc;
				ctmp1 = complex<float> (cosf (twopif * tdelay), -sinf (twopif * tdelay));
				ctmp2 = complex<float> (1.0f, 0.0f);
				for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
				{
					pSrcFilterAz[ipnew][ifreq] += pfxsource[itrc][ifreq] * ctmp2;
					pTgtFilter[ipnew][ifreq]   += pfxtarget[itrc][ifreq] * ctmp2;
					ctmp2 *= ctmp1;
				}
			}

			if (l1para->jointmethod == JOINTDESIG_SRC_DEGHOST || l1para->jointmethod == JOINTDESIG_SRCRCV_DEGHOST)
				for (ifreq = 0; ifreq < l1para->fftnc; ifreq++) 
					pSrcFilter[ipnew][ifreq] *= pSrcFilterAz[ipnew][ifreq];  
			else if (l1para->jointmethod == JOINTDESIGNATURE || l1para->jointmethod == JOINTDESIG_RCV_DEGHOST)
				for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
					pSrcFilter[ipnew][ifreq]  = pSrcFilterAz[ipnew][ifreq];

		} 

		free (gunvol); free (psx); free (psy); free (psz);
		flexible_free_array2d (ptxsource);
		flexible_free_array2d (pfxsource);
		flexible_free_array2d (ptxtarget);
		flexible_free_array2d (pfxtarget);
		PFL_FREE (hpfl);
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////

	int thfid=MIN(l1para->fftnc*0.85,l1para->highfreq/deltaf);
	int yftap=thfid*0.1;
	int ftap=thfid*0.1;

	float *hftaper = (float *)calloc(l1para->fftnc,sizeof(float));

	for (ifreq=0;ifreq<l1para->fftnc;ifreq++)
	{
		if (ifreq<thfid-ftap)
			hftaper[ifreq]=1.0f;
		else if (ifreq<thfid+ftap)
			hftaper[ifreq]=(1.0+cosf(0.5f*PI*(ifreq-thfid+ftap+1)/ftap))*0.5f;
		else
			hftaper[ifreq]=0.0f;
	}

	calc_zfilt(l1para);

	yftap=l1para->ylowcut*(l1para->ylowtap*0.01)/deltaf;

	int ylfid=MAX(yftap+1,l1para->ylowcut/deltaf);

	for (ifreq=0;ifreq<l1para->fftnc;ifreq++)
	{
		if (l1para->pmethod>0)
			l1para->zwtilt[ifreq]=twopif*MAX(lfid,ifreq)*l1para->zlowf[ifreq]*hftaper[ifreq];
		else if (l1para->zmethod==0)
			l1para->zwtilt[ifreq]=twopif*MAX(lfid,ifreq)*l1para->zlowf[ifreq];
		else
			l1para->zwtilt[ifreq]=l1para->zlowf[ifreq];

		if (l1para->ylowcut<0.001f||l1para->pmethod>0)
			l1para->ylowf[ifreq]=1.0f;
		else
		{
			if (ifreq<ylfid-yftap)
				l1para->ylowf[ifreq]=0.0f;
			else if (ifreq<ylfid+yftap)
				l1para->ylowf[ifreq]=(1.0+cosf(0.5f*PI*(ylfid+yftap-ifreq)/yftap))*0.5f;
			else
				l1para->ylowf[ifreq]=1.0f;
		}
		if (l1para->pmethod>0)
			l1para->ywtilt[ifreq]=twopif*MAX(lfid,ifreq)*l1para->ylowf[ifreq]*hftaper[ifreq];
		else if (l1para->ymethod==0)
			l1para->ywtilt[ifreq]=twopif*MAX(lfid,ifreq)*l1para->ylowf[ifreq];//*hftaper[ifreq];
		else
			l1para->ywtilt[ifreq]=l1para->ylowf[ifreq];
	}

	float epx,epy; //// effective px (slope)

	for ( ipnew = 0;ipnew < np_reduce; ipnew++ )
	{

		ip = ip_sorted[ipnew];

		px = l1para->pxsave[ip];
		py = l1para->pysave[ip];
		rq = l1para->rqsave[ip];

		//     fcos=MAX(-0.99999f,MIN(0.99999f,py*l1para->vwater));       
		//     l1para->ytop[ipnew]=-fcos;

		//     pxy=sqrtf(px*px+py*py);
		//     fsin=MAX(-0.99999f,MIN(0.99999f,pxy*l1para->vwater));       
		//     fcos=sqrtf(1.0f-fsin*fsin);//cosine of 3D surface angle

		//     l1para->ztop[ipnew]=MAX(l1para->zcosmin,fcos);

		if (l1para->choosemethod != CHOOSEMETHOD_JOINTSR3D ||
				(l1para->choosemethod == CHOOSEMETHOD_JOINTSR3D && l1para->computemode == 1))
		{
			for ( itrc = 0;itrc < nsgtrace; itrc++ )
			{
				//       phase0 = twopif*px*offsetx[itrc]+twopif*py*offsety[itrc];

				phase0 = twopif*px*offsetx[itrc]+twopif*py*offsety[itrc]+
					+twopif*rq*offsetx[itrc]*offsetx[itrc]+
					+twopif*rq*offsety[itrc]*offsety[itrc];

				epx = (2.0*offsetx[itrc]*rq+px);
				epy = (2.0*offsety[itrc]*rq+py);

				pxy=sqrtf(epx*epx+epy*epy);

				fcos=MAX(-0.99999f,MIN(0.99999f,epy*l1para->vwater));       
				l1para->ytop[ipnew]=-fcos;

				fsin=MAX(-0.99999f,MIN(0.99999f,pxy*l1para->vwater));       
				fcos=sqrtf(1.0f-fsin*fsin);//cosine of 3D surface angle
				l1para->ztop[ipnew]=MAX(l1para->zcosmin,fcos);

				phase1 = twopif*MAX(5.0e-4*l1para->tminscl,fcos*recz[itrc]*invwater);//half of ghost shift for reghosting
				phase2 = phase0+phase1; //mirror cable
				if (l1para->flag_recz_smooth==0)
					phase1 = phase0-phase1; //cable
				else
				{
					phase1 = twopif*MAX(5.0e-4*l1para->tminscl,fcos*reczou[itrc]*invwater);//half of primary shift for reghosting
					phase1 = phase0-phase1; //cable
				}

				ptr = (float *) &AP[ipnew][itrc];
				x = cosf(phase1);
				y = sinf(phase1);
				ptr[0] = x;
				ptr[1] = y;

				ptr = (float *) &AG[ipnew][itrc];
				x = cosf(phase2);
				y = sinf(phase2);
				ptr[0] = x;
				ptr[1] = y;

				ptr = (float *) &TS[itrc][ipnew];
				x = cosf(phase0);
				y = - sinf(phase0);
				ptr[0] = x;
				ptr[1] = y;

				ptr = (float *) &AS[ipnew][itrc];
				x = cosf(phase0);
				y = sinf(phase0);
				ptr[0] = x;
				ptr[1] = y;

				ptr = (float *) &TP[itrc][ipnew];
				x = cosf(phase1);
				y = - sinf(phase1);
				ptr[0] = x;
				ptr[1] = y;

				ptr = (float *) &TG[itrc][ipnew];
				x = cosf(phase2);
				y = - sinf(phase2);
				ptr[0] = x;
				ptr[1] = y;

				if (l1para->fpredatum == FPREDATUM_DGRG)
				{
					phase1 = twopif*MAX(5.0e-4*l1para->tminscl,fcos*reczou[itrc]*invwater);//half of ghost shift for reghosting
					phase2 = phase0+phase1; //mirror cable
					phase1 = phase0-phase1; //cable

					ptr = (float *) &TPR[itrc][ipnew];
					x = cosf(phase1);
					y = - sinf(phase1);
					ptr[0] = x;
					ptr[1] = y;

					ptr = (float *) &TGR[itrc][ipnew];
					x = cosf(phase2);
					y = - sinf(phase2);
					ptr[0] = x;
					ptr[1] = y;
				}

			} // for ( itrc = 0;itrc < nsgtrace; itrc++ )

			// zxue
			if (l1para->choosemethod == CHOOSEMETHOD_JOINTSR3D && l1para->computemode == 1 || 
					((l1para->choosemethod == CHOOSEMETHOD_MSDGI||l1para->choosemethod == CHOOSEMETHOD_FP4D) 
					 && l1para->jointmc == YESYES && l1para->jointmethod == JOINTDEGHOST))
			{
				// for a source ghost filter
				complex<float> ctmp1, ctmp2, ctmp3, ctmp4;

				float tsrcghost = MAX (5.0e-4 * l1para->tminscl, fcos * shotz[0] * invwater);
				if (l1para->datum == CABLE90)
				{
					ctmp1 = complex<float> (cosf (twopif * 2.0f * tsrcghost), -sinf (twopif * 2.0f * tsrcghost));
					ctmp2 = complex<float> (1.0f, 0.0f);

					for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
					{
						pSrcFilter[ipnew][ifreq] = (1.0f + l1para->srcrfl[ifreq] * ctmp2); 
						ctmp2 *= ctmp1;
					}
				}
				else
				{
					ctmp1 = complex<float> (cosf (twopif * tsrcghost), sinf (twopif * tsrcghost));
					ctmp2 = complex<float> (1.0f, 0.0f);
					ctmp3 = complex<float> (cosf (twopif * -tsrcghost), sinf (twopif * -tsrcghost));
					ctmp4 = complex<float> (1.0f, 0.0f);

					for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
					{
						pSrcFilter[ipnew][ifreq] = (ctmp2 + l1para->srcrfl[ifreq] * ctmp4);
						ctmp2 *= ctmp1;
						ctmp4 *= ctmp3;
					}
				}
			} // if (l1para->computemode == 1)
		}
		else // if (l1para->choosemethod != CHOOSEMETHOD_JOINTSR3D)
		{
			for ( itrc = 0;itrc < nsgtrace; itrc++ )
			{
				//       phase0 = twopif*px*offsetx[itrc]+twopif*py*offsety[itrc];

				phase0 = twopif*px*offsetx[itrc]+twopif*py*offsety[itrc]+
					+twopif*rq*offsetx[itrc]*offsetx[itrc]+
					+twopif*rq*offsety[itrc]*offsety[itrc];

				epx = (2.0*offsetx[itrc]*rq+px);
				epy = (2.0*offsety[itrc]*rq+py);

				pxy=sqrtf(epx*epx+epy*epy);

				fcos=MAX(-0.99999f,MIN(0.99999f,epy*l1para->vwater));       
				l1para->ytop[ipnew]=-fcos;

				fsin=MAX(-0.99999f,MIN(0.99999f,pxy*l1para->vwater));       
				fcos=sqrtf(1.0f-fsin*fsin);//cosine of 3D surface angle
				l1para->ztop[ipnew]=MAX(l1para->zcosmin,fcos);

				phase1 = twopif*MAX(5.0e-4*l1para->tminscl,fcos*recz[itrc]*invwater);//half of ghost shift for reghosting
				phase2 = twopif*MAX(5.0e-4*l1para->tminscl,fcos*shotz[itrc]*invwater);

				// shot depth to receiver depth
				x = cosf(phase0 - phase1 - phase2);
				y = sinf(phase0 - phase1 - phase2);
				ptr = (float *) &AP[ipnew][itrc];
				ptr[0] =  x; ptr[1] =  y;
				ptr = (float *) &TP[itrc][ipnew];
				ptr[0] =  x; ptr[1] = -y;

				// mirror shot to mirror receiver
				x = cosf(phase0 + phase1 + phase2);
				y = sinf(phase0 + phase1 + phase2);
				ptr = (float *) &AG[ipnew][itrc];
				ptr[0] =  x; ptr[1] =  y;
				ptr = (float *) &TG[itrc][ipnew];
				ptr[0] =  x; ptr[1] = -y;

				// shot to mirror receiver
				x = cosf(phase0 + phase1 - phase2);
				y = sinf(phase0 + phase1 - phase2);
				ptr = (float *) &AGR[ipnew][itrc];
				ptr[0] =  x; ptr[1] =  y;
				ptr = (float *) &TGR[itrc][ipnew];
				ptr[0] =  x; ptr[1] = -y;

				// mirror shot to receiver
				x = cosf(phase0 - phase1 + phase2);
				y = sinf(phase0 - phase1 + phase2);
				ptr = (float *) &AGS[ipnew][itrc];
				ptr[0] =  x; ptr[1] =  y;
				ptr = (float *) &TGS[itrc][ipnew];
				ptr[0] =  x; ptr[1] = -y;

				// surface to surface
				x = cosf(phase0);
				y = sinf(phase0);
				ptr = (float *) &AS[ipnew][itrc];
				ptr[0] =  x; ptr[1] =  y;
				ptr = (float *) &TS[itrc][ipnew];
				ptr[0] =  x; ptr[1] = -y;
			} // for ( itrc = 0;itrc < nsgtrace; itrc++ )
		} // if (l1para->choosemethod != CHOOSEMETHOD_JOINTSR3D)

	}


	free(hftaper);
	free(l_f);
}

void recalc_matrix_A_T_p_joint(l1inv_t *l1para, int nsgtrace, float *offsetx, float *offsety, float *recz, 
		int *ip_sorted, int np_reduce, int lfid)
{
	int ip, ipnew, itrc, isamp, ifreq, hpfl;
	float fnyq, deltaf, fsin, fcos, px=0.0f, py=0.0f, pxy=0.0f;
	float phase0, phase1, phase2, invwater=1.0/l1para->vwater;

	complex<float> ** __restrict AP,** __restrict AG, ** __restrict AS;
	complex<float> ** __restrict TP,** __restrict TG,** __restrict TS;
	complex<float> ** pSrcFilter, ** pSrcFilterAz, ** pTgtFilter;

	float x,y;

	fnyq = 500000.0f/l1para->srate;
	deltaf = fnyq/(float)(l1para->fftnc-1);

	float twopif=deltaf*2.0f*3.1415927f;
	float *l_f = (float *) calloc (l1para->fftnc, sizeof (float));
	for (ifreq = 0; ifreq < l1para->fftnc; ifreq++) l_f[ifreq] = twopif * ifreq; 

	head_t *src;
	float *gunvol = (float *) calloc (l1para->ntrcsrc, sizeof (float));
	float *psx    = (float *) calloc (l1para->ntrcsrc, sizeof (float));
	float *psy    = (float *) calloc (l1para->ntrcsrc, sizeof (float));
	float *psz    = (float *) calloc (l1para->ntrcsrc, sizeof (float));
	float shotx_avg, shoty_avg, shotz_avg, total_gunvol;
	float          **ptxsource = (float          **) flexible_array2d (l1para->ntrcsrc, l1para->fftnr, sizeof (float));
	float          **ptxtarget = (float          **) flexible_array2d (l1para->ntrcsrc, l1para->fftnr, sizeof (float));
	complex<float> **pfxsource = (complex<float> **) flexible_array2d (l1para->ntrcsrc, l1para->fftnc, sizeof (complex<float>));
	complex<float> **pfxtarget = (complex<float> **) flexible_array2d (l1para->ntrcsrc, l1para->fftnc, sizeof (complex<float>));
	int trlensig = NHEAD + l1para->nsampsrc;

	shotx_avg = shoty_avg = shotz_avg = total_gunvol = 0.0f;
	for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
	{
		src = (head_t *) (l1para->psig + itrc * trlensig);
		psx[itrc] = src->shotx;
		psy[itrc] = src->shoty;
		psz[itrc] = src->shotz;
		shotx_avg += src->shotx;
		shoty_avg += src->shoty;
		shotz_avg += src->shotz;
		memset (ptxsource[itrc], 0, l1para->fftnr * sizeof (float));
		memset (ptxtarget[itrc], 0, l1para->fftnr * sizeof (float));
		memcpy (ptxsource[itrc], l1para->psig + itrc * trlensig + NHEAD,
				l1para->nsampsrc * sizeof (float));
		//    memcpy (ptxtarget[itrc], l1para->psig + itrc * trlensig + NHEAD,
		//      (l1para->tzerosrc + l1para->tzerotgt) * 1000.0f * 1000.0f / l1para->srate * sizeof (float));
		int tmp1 = (l1para->tzerosrc + l1para->tzerotgt) * 1000.0f * 1000.0f / l1para->srate;
		int tmp2 = tmp1 + l1para->tzerotgttaper * 1000.0f * 1000.0f / l1para->srate;
		memcpy (ptxtarget[itrc], l1para->psig + itrc * trlensig + NHEAD, tmp2 * sizeof (float));
		for (isamp = tmp1; isamp < tmp2; isamp++)
			ptxtarget[itrc][isamp] *= powf (cosf (0.5f * 3.14159f * (isamp - tmp1 + 1) / (tmp2 - tmp1 + 1)), 2.0f);

		//for (isamp = 0; isamp < l1para->nsampsrc; isamp++)
		//  gunvol[itrc] = MAX (gunvol[itrc], fabsf (ptxsource[itrc][isamp]));
		if (l1para->designorm == DESIGNORM_MAX || l1para->designorm == DESIGNORM_NONE)
		{
			for (isamp = 0; isamp < l1para->nsampsrc; isamp++)
				gunvol[itrc] = MAX (gunvol[itrc], fabsf (ptxsource[itrc][isamp]));
		}
		else if (l1para->designorm == DESIGNORM_RMS)
		{
			for (isamp = 0; isamp < l1para->nsampsrc; isamp++)
				gunvol[itrc] += ptxsource[itrc][isamp] * ptxsource[itrc][isamp];
			gunvol[itrc] = sqrtf (gunvol[itrc] / l1para->nsampsrc);
		}
		total_gunvol += gunvol[itrc];
	}
	shotx_avg /= l1para->ntrcsrc;
	shoty_avg /= l1para->ntrcsrc;
	shotz_avg /= l1para->ntrcsrc;

	if (total_gunvol > 0.0f && l1para->designorm != DESIGNORM_NONE)
	{
		for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
		{
			for (isamp = 0; isamp < l1para->nsampsrc; isamp++)
			{
				ptxsource[itrc][isamp] /= total_gunvol;
				ptxtarget[itrc][isamp] /= total_gunvol;
			}
		}
	}
	if (total_gunvol > 0.0f)
		for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
			gunvol[itrc] /= total_gunvol;

	PFL_GET (&hpfl, PFL_PORTABLE_DEF);

	FFT_R2C_PP (hpfl, ptxsource, l1para->ntrcsrc, pfxsource, l1para->fftnr); 
	FFT_R2C_PP (hpfl, ptxtarget, l1para->ntrcsrc, pfxtarget, l1para->fftnr); 
	for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
	{
		for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
		{
			pfxsource[itrc][ifreq] *= (l1para->fftnr - 2);
			pfxtarget[itrc][ifreq] *= (l1para->fftnr - 2);
		}
	}

	AP = (complex<float> **)l1para->AP_p;
	AG = (complex<float> **)l1para->AG_p;
	AS = (complex<float> **)l1para->AS_p;

	TP = (complex<float> **)l1para->TP_p;
	TG = (complex<float> **)l1para->TG_p;
	TS = (complex<float> **)l1para->TS_p;

	if (l1para->choosemethod == CHOOSEMETHOD_JOINTSR3D)
	{
		pSrcFilter   = (complex<float> **) l1para->pSrcFilter;
		pSrcFilterAz = (complex<float> **) l1para->pSrcFilterAz;
		pTgtFilter   = (complex<float> **) l1para->pTgtFilter;
	}

	for (ipnew = 0; ipnew < np_reduce; ipnew++)
	{
		ip = ip_sorted[ipnew];

		px = l1para->pxsave[ip];
		py = l1para->pysave[ip];

		fcos=MAX(-0.99999f,MIN(0.99999f,py*l1para->vwater));       
		l1para->ytop[ipnew]=-fcos;

		pxy=sqrtf(px*px+py*py);
		fsin=MAX(-0.99999f,MIN(0.99999f,pxy*l1para->vwater));       
		fcos=sqrtf(1.0f-fsin*fsin);//cosine of 3D surface angle

		l1para->ztop[ipnew]=MAX(l1para->zcosmin,fcos);

		// For reverse tau-p operator and receiver ghost operator
		for ( itrc = 0;itrc < nsgtrace; itrc++ )
		{
			phase0 = twopif * px * offsetx[itrc] + twopif * py * offsety[itrc];
			phase1 = twopif*MAX(5.0e-4*l1para->tminscl,fcos*recz[itrc]*invwater);//half of ghost shift for reghosting
			phase2 = phase0+phase1; //mirror cable
			phase1 = phase0-phase1; //cable

			x = cosf(phase1); y = sinf(phase1);
			AP[ipnew][itrc] = complex<float> (x,  y);
			TP[itrc][ipnew] = complex<float> (x, -y);

			x = cosf(phase2); y = sinf(phase2);
			AG[ipnew][itrc] = complex<float> (x,  y);
			TG[itrc][ipnew] = complex<float> (x, -y);

			x = cosf(phase0); y = sinf(phase0);
			AS[ipnew][itrc] = complex<float> (x,  y);
			TS[itrc][ipnew] = complex<float> (x, -y);
		} // for ( itrc = 0;itrc < nsgtrace; itrc++ )

		// For source ghost operator
		for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
			pSrcFilter[ipnew][ifreq] = 0.0f;

		// zxue
		complex<float> ctmp1, ctmp2, ctmp3, ctmp4;
		if (l1para->datum == CABLE90)
			for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
			{
				src = (head_t *) (l1para->psig + itrc * (NHEAD + l1para->nsampsrc));
				float tsrcghost = MAX (5.0e-4 * l1para->tminscl, fcos * src->shotz * invwater);

				ctmp1 = complex<float> (cosf (twopif * 2.0f * tsrcghost), -sinf (twopif * 2.0f * tsrcghost));
				ctmp2 = complex<float> (1.0f, 0.0f);
				for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
				{
					pSrcFilter[ipnew][ifreq] += gunvol[itrc] * (1.0f + l1para->srcrfl[ifreq] * ctmp2); 
					ctmp2 *= ctmp1;
				}
			}
		else
			for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
			{
				src = (head_t *) (l1para->psig + itrc * (NHEAD + l1para->nsampsrc));
				float tsrcghost = MAX (5.0e-4 * l1para->tminscl, fcos * src->shotz * invwater);
				ctmp1 = complex<float> (cosf (twopif * ( tsrcghost - src->tdelay)), sinf (twopif * ( tsrcghost - src->tdelay)));
				ctmp2 = complex<float> (1.0f, 0.0f);
				ctmp3 = complex<float> (cosf (twopif * (-tsrcghost - src->tdelay)), sinf (twopif * (-tsrcghost - src->tdelay)));
				ctmp4 = complex<float> (1.0f, 0.0f);
				for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
				{
					pSrcFilter[ipnew][ifreq] += gunvol[itrc] * (ctmp2 + l1para->srcrfl[ifreq] * ctmp4);
					ctmp2 *= ctmp1;
					ctmp4 *= ctmp3;
				}
			}


		// for source signature operator
		for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
		{
			pSrcFilterAz[ipnew][ifreq] = 0.0f;
			pTgtFilter[ipnew][ifreq]   = 0.0f;
		}
		for (itrc = 0; itrc < l1para->ntrcsrc; itrc++)
		{
			head_t *src = (head_t *) (l1para->psig + itrc * (NHEAD + l1para->nsampsrc));
			float tdelay = -(psx[itrc] - shotx_avg) * px - (psy[itrc] - shoty_avg) * py - src->tdelay - l1para->tzerosrc;

			ctmp1 = complex<float> (cosf (twopif * tdelay), -sinf (twopif * tdelay));
			ctmp2 = complex<float> (1.0f, 0.0f);
			for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
			{
				pSrcFilterAz[ipnew][ifreq] += pfxsource[itrc][ifreq] * ctmp2;
				pTgtFilter[ipnew][ifreq]   += pfxtarget[itrc][ifreq] * ctmp2;
				ctmp2 *= ctmp1;
			}
		}


		if (l1para->jointmethod == JOINTDESIG_SRC_DEGHOST || l1para->jointmethod == JOINTDESIG_SRCRCV_DEGHOST)
			for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
				pSrcFilter[ipnew][ifreq] *= pSrcFilterAz[ipnew][ifreq];
		else if (l1para->jointmethod == JOINTDESIGNATURE || l1para->jointmethod == JOINTDESIG_RCV_DEGHOST)
			for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
				pSrcFilter[ipnew][ifreq]  = pSrcFilterAz[ipnew][ifreq];
	} // for (ipnew = 0; ipnew < np_reduce; ipnew++)

	free(l_f);
	free (gunvol); free (psx); free (psy); free (psz);
	flexible_free_array2d (ptxsource);
	flexible_free_array2d (ptxtarget);
	flexible_free_array2d (pfxsource);
	flexible_free_array2d (pfxtarget);
	PFL_FREE (hpfl);
}


void calc_matrix_A_T_interp3d(l1inv_t *l1para, int nsgtrace, int nsgtraceout, 
		float *offsetx, float *offsety, float *recz, 
		float *offsetxou, float *offsetyou, float *reczou,
		int lfid)
{
	int fftnr, npx, npy;
	int ip, ipx, ipy, itrc;
	float fnyq, deltaf;
	float fsin, fcos, px=0.0f, py=0.0f, pxy=0.0f;
	float phase0, phase1, phase2, invwater=1.0/l1para->vwater;
	int ifreq;

	float _Complex ** __restrict AP,** __restrict AG;
	float _Complex ** __restrict TP,** __restrict TG,** __restrict TS;
	float _Complex ** __restrict TP_interp3d,** __restrict TS_interp3d;

	float *ptr,x,y;

	PFL_LENGTH(l1para->nsamp,&fftnr);
	l1para->orgnr = fftnr;
	l1para->orgnc = fftnr/2+1;
	l1para->fftNR =        ALIGN_FLOAT(l1para->orgnr+2);     //for compatable with fftnr
	l1para->fftNC =        ALIGN_FCPLX(l1para->orgnc+1);     //for compatable with fftnc
	fftnr+=2;
	l1para->fftnr = fftnr;
	l1para->fftnc = fftnr/2;

	fnyq = 500000.0f/l1para->srate;
	deltaf = fnyq/(float)(l1para->fftnc-1);

	float twopif=deltaf*2.0f*3.1415927f;

	npx = l1para->npxsave;
	npy = l1para->npysave;

	AP = (float _Complex **)l1para->AP_p;
	AG = (float _Complex **)l1para->AG_p;

	TP = (float _Complex **)l1para->TP_p;
	TG = (float _Complex **)l1para->TG_p;
	TS = (float _Complex **)l1para->TS_p;

	TP_interp3d = (float _Complex **)l1para->TP_interp3d;
	TS_interp3d = (float _Complex **)l1para->TS_interp3d;

	int thfid=MIN(l1para->fftnc*0.85,l1para->highfreq/deltaf);
	int zftap=thfid*0.1;
	int yftap=thfid*0.1;
	int ftap=thfid*0.1;

	float *hftaper = (float *)calloc(l1para->fftnc,sizeof(float));

	for (ifreq=0;ifreq<l1para->fftnc;ifreq++)
	{
		if (ifreq<thfid-ftap)
			hftaper[ifreq]=1.0f;
		else if (ifreq<thfid+ftap)
			hftaper[ifreq]=(1.0+cosf(0.5f*PI*(ifreq-thfid+ftap+1)/ftap))*0.5f;
		else
			hftaper[ifreq]=0.0f;
	}

	zftap=l1para->zlowcut*(l1para->zlowtap*0.01)/deltaf;
	yftap=l1para->ylowcut*(l1para->ylowtap*0.01)/deltaf;

	int zlfid=MAX(zftap+1,l1para->zlowcut/deltaf);
	int ylfid=MAX(yftap+1,l1para->ylowcut/deltaf);

	for (ifreq=0;ifreq<l1para->fftnc;ifreq++)
	{
		if (l1para->zlowcut<0.001f||l1para->pmethod>0)
			l1para->zlowf[ifreq]=1.0f;
		else
		{
			if (ifreq<zlfid-zftap)
				l1para->zlowf[ifreq]=0.0f;
			else if (ifreq<zlfid+zftap)
				l1para->zlowf[ifreq]=(1.0+cosf(0.5f*PI*(zlfid+zftap-ifreq)/zftap))*0.5f;
			else
				l1para->zlowf[ifreq]=1.0f;
		}
		if (l1para->pmethod>0)
			l1para->zwtilt[ifreq]=twopif*MAX(lfid,ifreq)*l1para->zlowf[ifreq]*hftaper[ifreq];
		else if (l1para->zmethod==0)
			l1para->zwtilt[ifreq]=twopif*MAX(lfid,ifreq)*l1para->zlowf[ifreq];
		else
			l1para->zwtilt[ifreq]=l1para->zlowf[ifreq];

		if (l1para->ylowcut<0.001f||l1para->pmethod>0)
			l1para->ylowf[ifreq]=1.0f;
		else
		{
			if (ifreq<ylfid-yftap)
				l1para->ylowf[ifreq]=0.0f;
			else if (ifreq<ylfid+yftap)
				l1para->ylowf[ifreq]=(1.0+cosf(0.5f*PI*(ylfid+yftap-ifreq)/yftap))*0.5f;
			else
				l1para->ylowf[ifreq]=1.0f;
		}
		if (l1para->pmethod>0)
			l1para->ywtilt[ifreq]=twopif*MAX(lfid,ifreq)*l1para->ylowf[ifreq]*hftaper[ifreq];
		else if (l1para->ymethod==0)
			l1para->ywtilt[ifreq]=twopif*MAX(lfid,ifreq)*l1para->ylowf[ifreq];//*hftaper[ifreq];
		else
			l1para->ywtilt[ifreq]=l1para->ylowf[ifreq];
	}


	ip = -1;
	for ( ipy = 0;ipy < npy; ipy++ )
	{
		for ( ipx = 0;ipx < npx; ipx++ )
		{

			ip = ip + 1;

			px = l1para->pxsave[ip];
			py = l1para->pysave[ip];

			fcos=MAX(-0.99999f,MIN(0.99999f,py*l1para->vwater));       
			l1para->ytop[ip]=-fcos;

			pxy=sqrtf(px*px+py*py);
			fsin=MAX(-0.99999f,MIN(0.99999f,pxy*l1para->vwater));       
			fcos=sqrtf(1.0f-fsin*fsin);//cosine of 3D surface angle

			l1para->ztop[ip]=MAX(l1para->zcosmin,fcos);

			for ( itrc = 0;itrc < nsgtrace; itrc++ )
			{
				phase0 = twopif*px*offsetx[itrc]+twopif*py*offsety[itrc];

				phase1 = twopif*MAX(5.0e-4*l1para->tminscl,fcos*recz[itrc]*invwater);//half of ghost shift for reghosting
				phase2 = phase0+phase1; //mirror cable
				phase1 = phase0-phase1; //cable

				ptr = (float *) &AP[ip][itrc];
				x = cosf(phase1);
				y = sinf(phase1);
				ptr[0] = x;
				ptr[1] = y;

				ptr = (float *) &AG[ip][itrc];
				x = cosf(phase2);
				y = sinf(phase2);
				ptr[0] = x;
				ptr[1] = y;

				ptr = (float *) &TS[itrc][ip];
				x = cosf(phase0);
				y = - sinf(phase0);
				ptr[0] = x;
				ptr[1] = y;

				ptr = (float *) &TP[itrc][ip];
				x = cosf(phase1);
				y = - sinf(phase1);
				ptr[0] = x;
				ptr[1] = y;

				ptr = (float *) &TG[itrc][ip];
				x = cosf(phase2);
				y = - sinf(phase2);
				ptr[0] = x;
				ptr[1] = y;


			}

			for ( itrc = 0;itrc < nsgtraceout; itrc++ )
			{
				phase0 = twopif*px*offsetxou[itrc]+twopif*py*offsetyou[itrc];

				phase1 = twopif*MAX(5.0e-4*l1para->tminscl,fcos*reczou[itrc]*invwater);//half of ghost shift for reghosting
				phase2 = phase0+phase1; //mirror cable
				phase1 = phase0-phase1; //cable

				ptr = (float *) &TS_interp3d[itrc][ip];
				x = cosf(phase0);
				y = - sinf(phase0);
				ptr[0] = x;
				ptr[1] = y;

				ptr = (float *) &TP_interp3d[itrc][ip];
				x = cosf(phase1);
				y = - sinf(phase1);
				ptr[0] = x;
				ptr[1] = y;

			}

		}
	}

	free(hftaper);
}


void recalc_matrix_A_T_interp3d(l1inv_t *l1para, int nsgtrace, int nsgtraceout,
		float *offsetx, float *offsety, float *recz, 
		float *offsetxou, float *offsetyou, float *reczou,
		int *ip_sorted, int np_reduce, int lfid)
{
	int ip, ipnew, itrc;
	float fnyq, deltaf;
	float fsin, fcos, px=0.0f, py=0.0f, pxy=0.0f;
	float phase0, phase1, phase2, invwater=1.0/l1para->vwater;
	int ifreq;

	float _Complex ** __restrict AP,** __restrict AG;
	float _Complex ** __restrict TP,** __restrict TG,** __restrict TS;
	float _Complex ** __restrict TP_interp3d,** __restrict TS_interp3d;

	float *ptr,x,y;

	fnyq = 500000.0f/l1para->srate;
	deltaf = fnyq/(float)(l1para->fftnc-1);

	float twopif=deltaf*2.0f*3.1415927f;

	AP = (float _Complex **)l1para->AP_p;
	AG = (float _Complex **)l1para->AG_p;

	TP = (float _Complex **)l1para->TP_p;
	TG = (float _Complex **)l1para->TG_p;
	TS = (float _Complex **)l1para->TS_p;

	TP_interp3d = (float _Complex **)l1para->TP_interp3d;
	TS_interp3d = (float _Complex **)l1para->TS_interp3d;

	int thfid=MIN(l1para->fftnc*0.85,l1para->highfreq/deltaf);
	int zftap=thfid*0.1;
	int yftap=thfid*0.1;
	int ftap=thfid*0.1;

	float *hftaper = (float *)calloc(l1para->fftnc,sizeof(float));

	for (ifreq=0;ifreq<l1para->fftnc;ifreq++)
	{
		if (ifreq<thfid-ftap)
			hftaper[ifreq]=1.0f;
		else if (ifreq<thfid+ftap)
			hftaper[ifreq]=(1.0+cosf(0.5f*PI*(ifreq-thfid+ftap+1)/ftap))*0.5f;
		else
			hftaper[ifreq]=0.0f;
	}

	zftap=l1para->zlowcut*(l1para->zlowtap*0.01)/deltaf;
	yftap=l1para->ylowcut*(l1para->ylowtap*0.01)/deltaf;

	int zlfid=MAX(zftap+1,l1para->zlowcut/deltaf);
	int ylfid=MAX(yftap+1,l1para->ylowcut/deltaf);

	for (ifreq=0;ifreq<l1para->fftnc;ifreq++)
	{
		if (l1para->zlowcut<0.001f||l1para->pmethod>0)
			l1para->zlowf[ifreq]=1.0f;
		else
		{
			if (ifreq<zlfid-zftap)
				l1para->zlowf[ifreq]=0.0f;
			else if (ifreq<zlfid+zftap)
				l1para->zlowf[ifreq]=(1.0+cosf(0.5f*PI*(zlfid+zftap-ifreq)/zftap))*0.5f;
			else
				l1para->zlowf[ifreq]=1.0f;
		}
		if (l1para->pmethod>0)
			l1para->zwtilt[ifreq]=twopif*MAX(lfid,ifreq)*l1para->zlowf[ifreq]*hftaper[ifreq];
		else if (l1para->zmethod==0)
			l1para->zwtilt[ifreq]=twopif*MAX(lfid,ifreq)*l1para->zlowf[ifreq];
		else
			l1para->zwtilt[ifreq]=l1para->zlowf[ifreq];

		if (l1para->ylowcut<0.001f||l1para->pmethod>0)
			l1para->ylowf[ifreq]=1.0f;
		else
		{
			if (ifreq<ylfid-yftap)
				l1para->ylowf[ifreq]=0.0f;
			else if (ifreq<ylfid+yftap)
				l1para->ylowf[ifreq]=(1.0+cosf(0.5f*PI*(ylfid+yftap-ifreq)/yftap))*0.5f;
			else
				l1para->ylowf[ifreq]=1.0f;
		}
		if (l1para->pmethod>0)
			l1para->ywtilt[ifreq]=twopif*MAX(lfid,ifreq)*l1para->ylowf[ifreq]*hftaper[ifreq];
		else if (l1para->ymethod==0)
			l1para->ywtilt[ifreq]=twopif*MAX(lfid,ifreq)*l1para->ylowf[ifreq];//*hftaper[ifreq];
		else
			l1para->ywtilt[ifreq]=l1para->ylowf[ifreq];
	}


	for ( ipnew = 0;ipnew < np_reduce; ipnew++ )
	{

		ip = ip_sorted[ipnew];

		px = l1para->pxsave[ip];
		py = l1para->pysave[ip];

		fcos=MAX(-0.99999f,MIN(0.99999f,py*l1para->vwater));       
		l1para->ytop[ipnew]=-fcos;

		pxy=sqrtf(px*px+py*py);
		fsin=MAX(-0.99999f,MIN(0.99999f,pxy*l1para->vwater));       
		fcos=sqrtf(1.0f-fsin*fsin);//cosine of 3D surface angle

		l1para->ztop[ipnew]=MAX(l1para->zcosmin,fcos);

		for ( itrc = 0;itrc < nsgtrace; itrc++ )
		{
			phase0 = twopif*px*offsetx[itrc]+twopif*py*offsety[itrc];

			phase1 = twopif*MAX(5.0e-4*l1para->tminscl,fcos*recz[itrc]*invwater);//half of ghost shift for reghosting
			phase2 = phase0+phase1; //mirror cable
			phase1 = phase0-phase1; //cable

			ptr = (float *) &AP[ipnew][itrc];
			x = cosf(phase1);
			y = sinf(phase1);
			ptr[0] = x;
			ptr[1] = y;

			ptr = (float *) &AG[ipnew][itrc];
			x = cosf(phase2);
			y = sinf(phase2);
			ptr[0] = x;
			ptr[1] = y;

			ptr = (float *) &TS[itrc][ipnew];
			x = cosf(phase0);
			y = - sinf(phase0);
			ptr[0] = x;
			ptr[1] = y;

			ptr = (float *) &TP[itrc][ipnew];
			x = cosf(phase1);
			y = - sinf(phase1);
			ptr[0] = x;
			ptr[1] = y;

			ptr = (float *) &TG[itrc][ipnew];
			x = cosf(phase2);
			y = - sinf(phase2);
			ptr[0] = x;
			ptr[1] = y;

		}

		for ( itrc = 0;itrc < nsgtraceout; itrc++ )
		{
			phase0 = twopif*px*offsetxou[itrc]+twopif*py*offsetyou[itrc];

			phase1 = twopif*MAX(5.0e-4*l1para->tminscl,fcos*reczou[itrc]*invwater);//half of ghost shift for reghosting
			phase2 = phase0+phase1; //mirror cable
			phase1 = phase0-phase1; //cable

			ptr = (float *) &TS_interp3d[itrc][ipnew];
			x = cosf(phase0);
			y = - sinf(phase0);
			ptr[0] = x;
			ptr[1] = y;

			ptr = (float *) &TP_interp3d[itrc][ipnew];
			x = cosf(phase1);
			y = - sinf(phase1);
			ptr[0] = x;
			ptr[1] = y;

		}
	}

	free(hftaper);
}

void calc_zfilt(l1inv_t *l1para)
{
	int fftnr;
	if(l1para->clustertype == 0)
		fftnr = l1para->fftnr;
	else
		fftnr = l1para->fftnr + 2;
	float fnyq, deltaf,maxamp = 0.0f;
	int ifreq,hpfl;
	float realpart, imaginary, rms;

	fnyq = 500000.0f/l1para->srate;
	deltaf = fnyq/(float)(l1para->fftnc-1);

	float twopif=deltaf*2.0f*3.1415927f;
	int thfid=MIN(l1para->fftnc*0.85,l1para->highfreq/deltaf);
	int zftap=thfid*0.1;
	zftap=l1para->zlowcut*(l1para->zlowtap*0.01)/deltaf;
	int zlfid=MAX(zftap+1,l1para->zlowcut/deltaf);

	if ((l1para->zlowcut<0.001f && l1para->lzfilter == 0)||l1para->pmethod>0)
	  {
	    for (ifreq=0;ifreq<l1para->fftnc;ifreq++)
	      l1para->zlowf[ifreq]=1.0f;
	  }
	else
	  {
	    if(l1para->lzfilter > 0)
	      {
		float tdelayzfilt  = l1para->tzerozfilt;
		float **zfilter = (float **) flexible_array2d (1, fftnr, sizeof (float)); 
		complex<float> **pfxzfilter = (complex<float> **) flexible_array2d (1, l1para->fftnc, sizeof (complex<float>));
		
		PFL_GET (&hpfl, PFL_PORTABLE_DEF);
		memcpy (zfilter[0], l1para->zfilter, l1para->lzfilter * sizeof (float));

		FFT_R2C_PP (hpfl, zfilter, 1, pfxzfilter, fftnr); 
		for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
		  {
		    pfxzfilter[0][ifreq] *= (fftnr - 2);
					pfxzfilter[0][ifreq] *=  exp (complex<float> (0.0f, twopif * ifreq * tdelayzfilt));
		  }
		
		if(l1para->qczfilt == 1 && l1para->myrank ==1) printf(" Z filter: \n");
		
		for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
		  {
		    realpart = real(pfxzfilter[0][ifreq]) * real(pfxzfilter[0][ifreq]);
		    imaginary = imag(pfxzfilter[0][ifreq])* imag(pfxzfilter[0][ifreq]);
		    rms = sqrtf(realpart+imaginary);
		    l1para->zlowf[ifreq] = rms;
		  }
		for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
		  maxamp = MAX (maxamp, fabsf (l1para->zlowf[ifreq]));
		
		for (ifreq = 0; ifreq < l1para->fftnc; ifreq++)
		  {
		    l1para->zlowf[ifreq] /= maxamp;
		    if(l1para->qczfilt == 1 && l1para->myrank ==1) printf("%f  %f  \n",deltaf*ifreq,l1para->zlowf[ifreq]);
		  }
		
		l1para->qczfilt = 0;
		
		flexible_free_array2d (zfilter);
		flexible_free_array2d (pfxzfilter);
		PFL_FREE (hpfl);
	      }
	    else
	      {
		for (ifreq=0;ifreq<l1para->fftnc;ifreq++){
		  if (ifreq<zlfid-zftap)
		    l1para->zlowf[ifreq]=0.0f;
		  else if (ifreq<zlfid+zftap)
		    l1para->zlowf[ifreq]=(1.0+cosf(0.5f*PI*(zlfid+zftap-ifreq)/zftap))*0.5f;
		  else
		    l1para->zlowf[ifreq]=1.0f;
		}
	      }
	  }

}

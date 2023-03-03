#include "msdginterp.h"
#include <sys/types.h>
#include <unistd.h>
#include <vector>
using namespace std;

void data_distribute_prep (master_t *msdginterp);
void vsp_distribute_prep  (master_t *msdginterp, float *fshotx, float *fshoty);
void copyin_head(head_t *head, master_t *msdginterp, float *fshotx, float *fshoty, float *frecx, float *frecy);

void func_defsource_read(master_t *msdginterp, float *pdefSource);
void func_insert_defsource(master_t *msdginterp, int fsignature);

void
data_readin_deghost3d_pdata(int forigin, master_t *msdginterp)
{
  int i, j, scope, status;
  int vclen, pklen;
  int *src, *pvect, *pbuff;
  
  float *fobnxpos,*fobnypos;

  int sratems;

  head_t *head;
  size_t nwrite;

  int init,final;
  int NWRITE_SIZE_MAX,nbcast,nrem;

  CapiID oinst;

  ////////////////////////////////
  oinst = msdginterp->dbio->oinst;
  ////////////////////////////////

  NWRITE_SIZE_MAX = MIN(4 * 1024, msdginterp->ninst);
  nbcast = msdginterp->ninst / NWRITE_SIZE_MAX;
  nrem = msdginterp->ninst - nbcast * NWRITE_SIZE_MAX;
  if ( nrem != 0 ) nbcast = nbcast + 1;

  msdginterp->sl.v = static_cast<int*>(malloc(msdginterp->ninst*sizeof(int)));
  msdginterp->xl.v = static_cast<int*>(malloc(msdginterp->ninst*sizeof(int)));  

  float *fshotx = static_cast<float*>(malloc(msdginterp->ninst*sizeof(float)));
  float *fshoty = static_cast<float*>(malloc(msdginterp->ninst*sizeof(float)));
  float *frecx  = static_cast<float*>(malloc(msdginterp->ninst*sizeof(float)));
  float *frecy  = static_cast<float*>(malloc(msdginterp->ninst*sizeof(float)));
  msdginterp->shotz.v = static_cast<int*>(malloc(msdginterp->ninst*sizeof(int)));
  msdginterp->recz.v  = static_cast<int*>(malloc(msdginterp->ninst*sizeof(int)));

  if (msdginterp->l1para.flag_recz_smooth==1)
    msdginterp->recz_smooth.v  = static_cast<int*>(malloc(msdginterp->ninst*sizeof(int)));

  fobnxpos=(float*)msdginterp->obnxpos.v;
  fobnypos=(float*)msdginterp->obnypos.v;
  
  sratems = msdginterp->srate / 1000;

  int first, tlast, nakey, nbkey;
  if (msdginterp->highfreq==0) msdginterp->highfreq=500000/msdginterp->srate;
  nbkey = 1000;
  nakey = (int)ceilf(msdginterp->ninst*1.0/nbkey);
  if (msdginterp->ninst > nakey*nbkey) G2_logCritical(msdginterp->logName, "Number of input instances=%d. Must be <= %d x %d", msdginterp->ninst, nakey, nbkey);

  vclen = msdginterp->vlen;
  pklen = NHEAD + vclen;
 
  pvect = static_cast<int*>(malloc((vclen) * nbkey*sizeof(int)));
  pbuff = static_cast<int*>(malloc((pklen) * nbkey*sizeof(int)));

  if(msdginterp->l1para.zerotraceattr == 1)
    {
      msdginterp->zerotrace.v = static_cast<int*>(malloc(msdginterp->ninst*sizeof(int)));
      memset(msdginterp->zerotrace.v, 0, msdginterp->ninst*sizeof(int));
    }

  lseek (forigin, 0, SEEK_SET);

//////////////////////////////////////////
//       read in normal data          //
//////////////////////////////////////////
   scope = 0;
   msdginterp->wbd.v = static_cast<int*>(calloc(msdginterp->ninst,sizeof(int)));
   // to me cpp casting doesn't make sense with C malloc/calloc -- haghakha
   msdginterp->regattr.rv = (float*) calloc(msdginterp->ninst,sizeof(float)); 
      
   for(j=1; j<=nbcast; j++)
     {
       
       init = (j - 1) * NWRITE_SIZE_MAX + 1;
       final = j * NWRITE_SIZE_MAX;
       if ( j == nbcast )
	 {
	   if ( nrem != 0 )
	     {
	       final = init + nrem - 1;
	    }
	 }

       if ( cCapi_readByName( oinst, msdginterp->shotz.name, msdginterp->shotz.v+(size_t)(init-1), init, final ) != FCAPI_SUCCESS)
         G2_logCritical(msdginterp->logName, "Fail to read shotz attribute %s from input", msdginterp->shotz.name);
      
       if ( cCapi_readByName( oinst, msdginterp->recz.name,  msdginterp->recz.v +(size_t)(init-1), init, final ) != FCAPI_SUCCESS)
         G2_logCritical(msdginterp->logName, "Fail to read recz attribute %s from input", msdginterp->recz.name);

       if (msdginterp->l1para.flag_recz_smooth==1)
	 {
	   if ( cCapi_readByName( oinst, msdginterp->recz_smooth.name,  msdginterp->recz_smooth.v +(size_t)(init-1), init, final ) != FCAPI_SUCCESS)
              G2_logCritical(msdginterp->logName, "Fail to read recz_smooth attribute %s from input", msdginterp->recz_smooth.name);
	 }

       if ( cCapi_readByName( oinst, msdginterp->sl.name,    msdginterp->sl.v   +(size_t)(init-1), init, final ) != FCAPI_SUCCESS)
         G2_logCritical(msdginterp->logName, "Fail to read subline attribute %s from input", msdginterp->sl.name);
       
       if ( cCapi_readByName( oinst, msdginterp->xl.name,    msdginterp->xl.v   +(size_t)(init-1), init, final ) != FCAPI_SUCCESS)
         G2_logCritical(msdginterp->logName, "Fail to read crossline attribute %s from input", msdginterp->xl.name);
       
       // SHOTX
       if ( cCapi_readConvertByName( oinst, msdginterp->shotx.name, fshotx+(size_t)(init-1), FCAPI_FLOAT32, init, final ) != FCAPI_SUCCESS)
         G2_logCritical(msdginterp->logName, "Fail to read shotx attribute %s from input", msdginterp->shotx.name);
       
       // SHOTY
       if ( cCapi_readConvertByName( oinst, msdginterp->shoty.name, fshoty+(size_t)(init-1), FCAPI_FLOAT32, init, final ) != FCAPI_SUCCESS)
         G2_logCritical(msdginterp->logName, "Fail to read shoty attribute %s from input", msdginterp->shoty.name);
       
       // RECX
       if ( cCapi_readConvertByName( oinst, msdginterp->recx.name,  frecx+(size_t)(init-1), FCAPI_FLOAT32, init, final ) != FCAPI_SUCCESS)
         G2_logCritical(msdginterp->logName, "Fail to read recx attribute %s from input", msdginterp->recx.name);
       
       // RECY
       if ( cCapi_readConvertByName( oinst, msdginterp->recy.name,  frecy+(size_t)(init-1),  FCAPI_FLOAT32, init, final ) != FCAPI_SUCCESS)
         G2_logCritical(msdginterp->logName, "Fail to read recy attribute %s from input", msdginterp->recy.name);


       if (j==1) G2_logInfo(msdginterp->logName, "msdginterp->wbduse=%d",msdginterp->wbduse);
       if (msdginterp->wbduse==YESYES)
	 {
	   if ( cCapi_readByName( oinst, msdginterp->wbd.name,   msdginterp->wbd.v+(size_t)(init-1),   init, final ) != FCAPI_SUCCESS)
             G2_logCritical(msdginterp->logName, "Fail to read water bottom attribute %s from input", msdginterp->wbd.name);
	   if (j==1) G2_logInfo(msdginterp->logName, "Water bottom attribute read in!");
	 }
      
       // reading regularization parameter
       if (j==1) G2_logInfo(msdginterp->logName, "msdginterp->regattruse=%d",msdginterp->regattruse);
       if (msdginterp->regattruse==YESYES){
         if ( cCapi_readConvertByName( oinst, msdginterp->regattr.name,   msdginterp->regattr.rv+(size_t)(init-1),FCAPI_FLOAT32,   init, final ) != FCAPI_SUCCESS)
           G2_logCritical(msdginterp->logName, "Fail to read  Regularization attribute %s from input", msdginterp->regattr.name);
         if (j==1) G2_logInfo(msdginterp->logName, "Regularization attribute read in!");
       }

       //zerotraceattr
       if(msdginterp->l1para.zerotraceattr == 1)
	 {
	   if ( cCapi_readByName( oinst, msdginterp->zerotrace.name, msdginterp->zerotrace.v+(size_t)(init-1), init, final ) != FCAPI_SUCCESS)
             G2_logCritical(msdginterp->logName, "Fail to read zerotrace attribute %s from input", msdginterp->zerotrace.name);
	   if (j==1) G2_logInfo(msdginterp->logName, "Zerotraceattr attribute read in!");
	 }

     }//// end of : for(j=1; j<=nbcast; j++)

   for(i=0; i<nakey; i++)
    {
      first = i*nbkey + 1;
      tlast = MIN((i+1)*nbkey, msdginterp->ninst);

      if ( cCapi_readByName( oinst, "TRACE", pvect, first, tlast ) != FCAPI_SUCCESS)
        G2_logCritical(msdginterp->logName, "Fail to read TRACE vector from input");
      
      for(j=0; j<tlast-first+1; j++)
	{
	  head = (head_t*)(pbuff+j*pklen);
	  memset(head, 0, pklen*sizeof(int));
	  head->tid = first+j;
	  head->pval = msdginterp->ipky.v[first-1+j];
	  head->sval = msdginterp->isky.v[first-1+j];
	  head->ipky = head->pval;
	  head->isky = head->sval;
	  
	  copyin_head(head, msdginterp, fshotx, fshoty, frecx, frecy);
	  src = (int*)(pvect + j*vclen);
	  memcpy((int*)head+NHEAD, src, vclen*sizeof(int));
	}
      nwrite = write(forigin, pbuff, (tlast-first+1)*pklen*sizeof(int));
      if (nwrite != (tlast-first+1)*pklen*sizeof(int)) 
        G2_logCritical(msdginterp->logName, "Fail to write data from input. Should write %zu, but wrote only %zu", (tlast-first+1)*pklen*sizeof(int), nwrite );
    }

   if (msdginterp->dgmethod == DGMETHOD_VSP) vsp_distribute_prep (msdginterp,fshotx,fshoty); 

   G2_logInfo(msdginterp->logName, "finished reading the traces");

   free(pvect);                //// <---
   free(pbuff);                //// <--- 
   free(msdginterp->wbd.v);    //// <---
   free(msdginterp->regattr.rv);

   free(msdginterp->sl.v);       //// <---
   free(msdginterp->xl.v);       //// <---
   if(msdginterp->l1para.zerotraceattr == 1) free(msdginterp->zerotrace.v);

   free(fshotx);      //// <---
   free(fshoty);      //// <---
   free(frecx);       //// <---
   free(frecy);       //// <---

   free(msdginterp->shotz.v);  //// <---
   free(msdginterp->recz.v);   //// <---

   if (msdginterp->l1para.flag_recz_smooth==1)
     free(msdginterp->recz_smooth.v);   //// <---

}

void
copyin_head(head_t *head, master_t *msdginterp, float *fshotx, float *fshoty, float *frecx, float *frecy)
{
  
  int tid; 
  float *frecz_smooth;

  if (msdginterp->l1para.interp3d==YESYES) tid = head->tidinp-1;
  else tid = head->tid-1;

  assert(tid>=0 && tid<msdginterp->ninst);  //// <---

  float *fshotz=(float*)msdginterp->shotz.v;
  float *frecz=(float*)msdginterp->recz.v;

  if (msdginterp->l1para.flag_recz_smooth==1) frecz_smooth=(float*)msdginterp->recz_smooth.v;

  head->shotx = fshotx[tid];
  head->shoty = fshoty[tid];
  head->shotz = fshotz[tid];
  head->recx =  frecx[tid];
  head->recy =  frecy[tid];
  head->recz =  frecz[tid];
  if (msdginterp->l1para.flag_recz_smooth==1) head->recz_smooth =  frecz_smooth[tid];
  head->offs =  msdginterp->offs.v[tid];
  head->offset = head->offs;
  // if msdginterp->regattruse==NONO this value is 0 
  head->regparam = msdginterp->regattr.rv[tid];

  if ((msdginterp->l1para.deghost3d == YESYES)&&(head->recz<=0||head->shotz<=0))
   {
    G2_logCritical(msdginterp->logName, "Zero or negative receiver/shot depth detected (shotz,recz)=(%f,%f)", head->shotz,head->recz);
   }

  if (msdginterp->l1para.flag_recz_smooth==1)
    {
      if ((msdginterp->l1para.deghost3d == YESYES)&&(head->recz_smooth<=0))
	{
	  G2_logCritical(msdginterp->logName, "Zero or negative receiver depth2 detected (reczout)=(%f)", head->recz_smooth);
	}
    }

  float *fwbd=(float *)msdginterp->wbd.v;
  
  if (msdginterp->wbd.type==FCAPI_FLOAT32)
   {
    head->wbid = MAX(2,MIN(msdginterp->nsamp-1,round((fwbd[tid]+msdginterp->wbshift)*1000.0f/msdginterp->srate)));
    if (head->wbid>msdginterp->wdepthmax) msdginterp->wdepthmax=head->wbid;
    if (head->wbid<msdginterp->wdepthmin) msdginterp->wdepthmin=head->wbid;
   }
  else if (msdginterp->wbd.type==FCAPI_INT32)
   {
    head->wbid = MAX(2,MIN(msdginterp->nsamp-1,round((msdginterp->wbd.v[tid]+msdginterp->wbshift)*1000.0f/msdginterp->srate)));
    if (head->wbid>msdginterp->wdepthmax) msdginterp->wdepthmax=head->wbid;
    if (head->wbid<msdginterp->wdepthmin) msdginterp->wdepthmin=head->wbid;
   }

  if ( msdginterp->wbduse==YESYES&&head->wbid<0)
   {
    G2_logCritical(msdginterp->logName, "Zero or negative water depth detected!");
   }

  head->xl = msdginterp->xl.v[tid];
  head->sl = msdginterp->sl.v[tid];

  float *fobnxpos=(float*)msdginterp->obnxpos.v;
  float *fobnypos=(float*)msdginterp->obnypos.v;
        
  head->obnxpos = fobnxpos[tid];
  head->obnypos = fobnypos[tid];
  
  head->cableid = msdginterp->cableid.v[tid];
  head->newcableid = msdginterp->newcableid.v[tid];
  head->cableidseq = msdginterp->cableidseq.v[tid];

  if(msdginterp->l1para.zerotraceattr == 1)
    head->zeroattr = msdginterp->zerotrace.v[tid];
}
  
void data_distribute_prep (master_t *msdginterp)
{
  int ipky, itrc, itrcx, itrcy, xpos, ypos;
  int isbegx, isendx, isbegy, isendy, jtrc, ifac;
  int obnixmax, obniymax, nsgtracex, nsgtracey, nsgoverlapx, nsgoverlapy,iwinincy,iwinincx;
  int ntaskmax = 200000, ntaskinc = 200000, ntaskprv = 0;
  
  vector<vector<vector<int> > > trcloc;
  float *fobnxpos = (float *) msdginterp->obnxpos.v, *fobnypos = (float *) msdginterp->obnypos.v, rtap;
  
  msdginterp->l1para.nsgtracexlastgate = 0;
  
  obnixmax          = msdginterp->l1para.obnixmax; 
  obniymax          = msdginterp->l1para.obniymax; 
  nsgtracex         = msdginterp->l1para.nsgtracex; 
  nsgtracey         = msdginterp->l1para.nsgtracey; 
  nsgoverlapx       = msdginterp->l1para.nsgoverlapx; 
  nsgoverlapy       = msdginterp->l1para.nsgoverlapy; 
  msdginterp->tasksizemax = 4 * (nsgtracex  * nsgtracey); 

  if (msdginterp->dgmethod == DGMETHOD_OBN)
    msdginterp->tasksizemax = 10 * (nsgtracex  * nsgtracey); 

  if ( msdginterp->l1para.choosemethod == CHOOSEMETHOD_FP4D )
    ifac = msdginterp->l1para.ndata*2;
  else
    ifac = 1;

  msdginterp->taskin          = (int *)   malloc (ntaskmax * msdginterp->tasksizemax * sizeof (int));
  msdginterp->taskou          = (int *)   malloc (ntaskmax * msdginterp->tasksizemax * sizeof (int));
  msdginterp->taskweight      = (float *) malloc (ntaskmax * msdginterp->tasksizemax * sizeof (float));
  msdginterp->ntaskpershot    = (int *)   malloc (msdginterp->nipky * sizeof (int));
  msdginterp->ntrcinpertask   = (int *)   malloc (ntaskmax * sizeof (int));
  msdginterp->ntrcoupertask   = (int *)   malloc (ntaskmax * sizeof (int));
  msdginterp->shotidpertask   = (int *)   malloc (ntaskmax * sizeof (int));
  msdginterp->ampscalar       = (float *) calloc (msdginterp->nipky, sizeof (float));
  msdginterp->fileidpershot   = (int *)   malloc (msdginterp->nipky * ifac * sizeof (int));
  msdginterp->filenamepershot = (char **) malloc (msdginterp->nipky * ifac * sizeof (char *));
  
  memset (msdginterp->taskin,        -1, ntaskmax * msdginterp->tasksizemax * sizeof (int));
  memset (msdginterp->taskou,        -1, ntaskmax * msdginterp->tasksizemax * sizeof (int));
  memset (msdginterp->taskweight,    -1, ntaskmax * msdginterp->tasksizemax * sizeof (float));
  memset (msdginterp->fileidpershot, -1, msdginterp->nipky * ifac * sizeof (int));
  
  msdginterp->ntask = 0;

  trcloc.resize (obnixmax);
  for (itrc = 0; itrc < obnixmax; itrc++)
    trcloc[itrc].resize (obniymax);
  
  for (ipky = 0; ipky < msdginterp->nipky; ipky++)
    {

      for (itrcx = 0; itrcx < obnixmax; itrcx++)
	for (itrcy = 0; itrcy < obniymax; itrcy++)
	  trcloc[itrcx][itrcy].resize (0);
      for (itrc = msdginterp->ftrc[ipky]; itrc < msdginterp->ftrc[ipky] + msdginterp->nsky[ipky]; itrc++)
	{
	  xpos = (int) fobnxpos[itrc];
	  ypos = (int) fobnypos[itrc];
	  if (!(0 <= xpos && xpos < obnixmax)) G2_logCritical(msdginterp->logName, "xpos=%d is not in the range 0 - %d", xpos, obnixmax);
	  if (!(0 <= ypos && ypos < obniymax)) G2_logCritical(msdginterp->logName, "ypos=%d is not in the range 0 - %d", ypos, obniymax);
	  trcloc[xpos][ypos].push_back (itrc);
	} // for (itrc = 0; itrc < msdginterp->nsky[ipky]; itrc++)
      
      isbegy = 0;
      isendy = -1;
      iwinincy = nsgtracey-nsgoverlapy;
      iwinincx = nsgtracex-nsgoverlapx;
      
      while (isendy < obniymax-1)
	{
	  isendy = MIN (obniymax-1,isbegy +nsgtracey-1);
	  isbegy = MAX (0,isendy-nsgtracey+1);
	  isbegx = 0;
	  isendx= -1;
	  
	  while (isendx < obnixmax-1)
	    {
	      isendx = MIN (obnixmax-1,isbegx +nsgtracex-1);
	      isbegx = MAX (0,isendx-nsgtracex+1);

	      itrc = 0;
	      for (itrcy = isbegy; itrcy <= isendy; itrcy++)
		{
		  for (itrcx = isbegx; itrcx <= isendx; itrcx++)
		    {
		      for (jtrc = 0; jtrc < trcloc[itrcx][itrcy].size (); jtrc++)
			{
			  rtap  = MIN (1.0f, (float) (itrcx - isbegx + 1.0f) / 
				       (float) (nsgoverlapx ));
			  rtap *= MIN (1.0f, (float) (isendx+1 - itrcx) / 
				       (float) (nsgoverlapx));
			  rtap *= MIN (1.0f, (float) (itrcy - isbegy + 1.0f) / 
				       (float) (nsgoverlapy));
			  rtap *= MIN (1.0f, (float) (isendy+1 - itrcy) / 
				       (float) (nsgoverlapy));

			  msdginterp->taskin[msdginterp->ntask * msdginterp->tasksizemax + itrc] = trcloc[itrcx][itrcy][jtrc];
			  msdginterp->taskou[msdginterp->ntask * msdginterp->tasksizemax + itrc] = trcloc[itrcx][itrcy][jtrc];
			  msdginterp->taskweight[msdginterp->ntask * msdginterp->tasksizemax + itrc++] = rtap;
			  assert (0.0f < rtap  && rtap <= 1.0f);
			  assert (itrc <= msdginterp->tasksizemax);
			}
		    }
		}
	      msdginterp->ntrcinpertask[msdginterp->ntask] = itrc;
	      msdginterp->ntrcoupertask[msdginterp->ntask] = itrc;
	      msdginterp->shotidpertask[msdginterp->ntask] = ipky;
	      
	      msdginterp->l1para.nsgtracexlastgate = MAX (msdginterp->l1para.nsgtracexlastgate, isendx - isbegx+1);
	      
	      if (++(msdginterp->ntask) >= ntaskmax) 
		{
		  ntaskmax += ntaskinc;
		  msdginterp->taskin = (int *) realloc 
		    (msdginterp->taskin, ntaskmax * msdginterp->tasksizemax * sizeof (int));
		  msdginterp->taskou = (int *) realloc 
		    (msdginterp->taskou, ntaskmax * msdginterp->tasksizemax * sizeof (int));
		  msdginterp->taskweight = (float *) realloc 
		    (msdginterp->taskweight, ntaskmax * msdginterp->tasksizemax * sizeof (float));
		  memset (&msdginterp->taskin[(ntaskmax - ntaskinc) * msdginterp->tasksizemax], -1, 
			  ntaskinc * msdginterp->tasksizemax * sizeof (int));
		  memset (&msdginterp->taskou[(ntaskmax - ntaskinc) * msdginterp->tasksizemax], -1, 
			  ntaskinc * msdginterp->tasksizemax * sizeof (int));
		  memset (&msdginterp->taskweight[(ntaskmax - ntaskinc) * msdginterp->tasksizemax], -1, 
			  ntaskinc * msdginterp->tasksizemax * sizeof (float));
		  
		  msdginterp->ntrcinpertask = (int *) realloc 
		    (msdginterp->ntrcinpertask, ntaskmax * sizeof (int));
		  msdginterp->ntrcoupertask = (int *) realloc 
		    (msdginterp->ntrcoupertask, ntaskmax * sizeof (int));
		  msdginterp->shotidpertask = (int *) realloc 
		    (msdginterp->shotidpertask, ntaskmax * sizeof (int));
		}
	      isbegx = isbegx+iwinincx;
	    } // while (isendx < obnixmax)
	  isbegy = isbegy+iwinincy;
	} // while (isendy < obniymax)
      
      if (ipky == 0)
	{
	  ntaskprv = 0;
	  msdginterp->ntaskpershot[ipky] = msdginterp->ntask;
	}
      else
	{
	  ntaskprv += msdginterp->ntaskpershot[ipky - 1];
	  msdginterp->ntaskpershot[ipky] = msdginterp->ntask - ntaskprv;
	}
      

      
    } // for (ipky = 0; ipky < msdginterp->nipky; ipky++)
  
  msdginterp->taskin     = (int *)   realloc (msdginterp->taskin,     
					      msdginterp->ntask * msdginterp->tasksizemax * sizeof (int));
  msdginterp->taskou     = (int *)   realloc (msdginterp->taskou,     
					      msdginterp->ntask * msdginterp->tasksizemax * sizeof (int));
  msdginterp->taskweight = (float *) realloc (msdginterp->taskweight, 
					      msdginterp->ntask * msdginterp->tasksizemax * sizeof (float));
  
  msdginterp->ntrcinpertask = (int *) realloc (msdginterp->ntrcinpertask, msdginterp->ntask * sizeof (int));
  msdginterp->ntrcoupertask = (int *) realloc (msdginterp->ntrcoupertask, msdginterp->ntask * sizeof (int));
  msdginterp->shotidpertask = (int *) realloc (msdginterp->shotidpertask, msdginterp->ntask * sizeof (int));
  
  return;
}


void vsp_distribute_prep (master_t *bsdg, float *fshotx, float *fshoty)
{
	int ipky, isky, itrc, itrcx, itrcy, xpos, ypos, nsgatex, nsgatey, isgatex, isgatey, jsgatex, jsgatey;
	int isbegx, isendx, isbegy, isendy;
	int obnixmax, obniymax, nsgtracex, nsgtracey, nsgoverlapx, nsgoverlapy;
	int ntaskmax = 200000, ntaskinc = 200000, ntaskprv = 0;
	int ***trcloc;
    float xmin = fshotx[0], xmax = fshotx[0], ymin = fshoty[0], ymax = fshoty[0], rtap;

	bsdg->l1para.nsgtracexlastgate = 0;

	nsgtracex         = bsdg->l1para.nsgtracex; 
	nsgtracey         = bsdg->l1para.nsgtracey; 
	nsgoverlapx       = bsdg->l1para.nsgoverlapx; 
	nsgoverlapy       = bsdg->l1para.nsgoverlapy; 
	bsdg->tasksizemax = 15 * nsgtracex * nsgtracey;

	bsdg->taskin          = (int *)   malloc (ntaskmax * bsdg->tasksizemax * sizeof (int));
	bsdg->taskou          = (int *)   malloc (ntaskmax * bsdg->tasksizemax * sizeof (int));
	bsdg->taskweight      = (float *) malloc (ntaskmax * bsdg->tasksizemax * sizeof (float));
	//bsdg->ntaskpershot    = (int *)   malloc (bsdg->nipky * sizeof (int));
	bsdg->ntaskpershot    = (int *)   malloc (1 * sizeof (int));
	bsdg->ntrcinpertask   = (int *)   malloc (ntaskmax * sizeof (int));
	bsdg->ntrcoupertask   = (int *)   malloc (ntaskmax * sizeof (int));
	bsdg->shotidpertask   = (int *)   malloc (ntaskmax * sizeof (int));
	//bsdg->ampscalar       = (float *) calloc (bsdg->nipky, sizeof (float));
	//bsdg->fileidpershot   = (int *)   malloc (bsdg->nipky * sizeof (int));
	//bsdg->filenamepershot = (char **) malloc (bsdg->nipky * sizeof (char *));
	bsdg->ampscalar       = (float *) calloc (1, sizeof (float));
	bsdg->fileidpershot   = (int *)   malloc (1 * sizeof (int));
	bsdg->filenamepershot = (char **) malloc (1 * sizeof (char *));

	memset (bsdg->taskin,        -1, ntaskmax * bsdg->tasksizemax * sizeof (int));
	memset (bsdg->taskou,        -1, ntaskmax * bsdg->tasksizemax * sizeof (int));
	memset (bsdg->taskweight,    -1, ntaskmax * bsdg->tasksizemax * sizeof (float));
	//memset (bsdg->fileidpershot, -1, bsdg->nipky * sizeof (int));
	memset (bsdg->fileidpershot, -1, 1 * sizeof (int));

    for (itrc = 0; itrc < bsdg->ninst; itrc++)
    {
        xmin = MIN (xmin, fshotx[itrc]);
        xmax = MAX (xmax, fshotx[itrc]);
        ymin = MIN (ymin, fshoty[itrc]);
        ymax = MAX (ymax, fshoty[itrc]);
    }
    obnixmax = (int) ((xmax - xmin) / bsdg->l1para.chanitvx + 1);
    obniymax = (int) ((ymax - ymin) / bsdg->l1para.chanitvy + 1);
    nsgatex = MAX (1, (int) (obnixmax / nsgtracex)); 
    nsgatey = MAX (1, (int) (obniymax / nsgtracey)); 

    G2_logInfo (bsdg->logName, "%s: shotx range (%f, %f)", __func__, xmin, xmax);
    G2_logInfo (bsdg->logName, "%s: shoty range (%f, %f)", __func__, ymin, ymax);
    G2_logInfo (bsdg->logName, "%s: obnixmax %d obniymax %d nsgatex %d nsgatey %d", __func__, 
                obnixmax, obniymax, nsgatex, nsgatey);

	bsdg->ntask = 0;

	trcloc = (int ***) malloc (obnixmax * sizeof (int **));
	for (itrcx = 0; itrcx < obnixmax; itrcx++)
    {
	    trcloc[itrcx] = (int **) malloc (obniymax * sizeof (int *));
	    for (itrcy = 0; itrcy < obniymax; itrcy++)
            trcloc[itrcx][itrcy] = (int *) malloc (60 * sizeof (int));
    }

	//for (ipky = 0; ipky < bsdg->nipky; ipky++)
	for (ipky = 0; ipky < 1; ipky++)
	{
		for (itrcx = 0; itrcx < obnixmax; itrcx++)
		    for (itrcy = 0; itrcy < obniymax; itrcy++)
			    memset (trcloc[itrcx][itrcy], -1, 60 * sizeof (int));

		//for (itrc = bsdg->ftrc[ipky]; itrc < bsdg->ftrc[ipky] + bsdg->nsky[ipky]; itrc++)
		for (itrc = 0; itrc < bsdg->ninst; itrc++)
		{
			xpos = (int) ((fshotx[itrc] - xmin) / bsdg->l1para.chanitvx);
			ypos = (int) ((fshoty[itrc] - ymin) / bsdg->l1para.chanitvy);
	                if (!(0 <= xpos && xpos < obnixmax)) G2_logCritical(bsdg->logName, "xpos=%d is not in the range 0 - %d", xpos, obnixmax);
	                if (!(0 <= ypos && ypos < obniymax)) G2_logCritical(bsdg->logName, "ypos=%d is not in the range 0 - %d", ypos, obniymax);

                        for (isky = 0; isky < 60; isky++)
                        {
                            if (trcloc[xpos][ypos][isky] == -1)
                            {
            			        trcloc[xpos][ypos][isky] = itrc;
                                break;
                            }
                        }
                        if (isky == 60) G2_logWarn(bsdg->logName, "Number of input traces in 1 grid is more than 60");
		} // for (itrc = 0; itrc < bsdg->nsky[ipky]; itrc++)

	    isbegy = 0;
	    isendy = -1;
	    int iwinincy = nsgtracey-nsgoverlapy;
	    int iwinincx = nsgtracex-nsgoverlapx;
		
	    while (isendy < obniymax-1)
	      {
		      isendy = MIN (obniymax-1,isbegy +nsgtracey-1);
		      isbegy = MAX (0,isendy-nsgtracey+1);
		      isbegx = 0;
		      isendx= -1;
		
		      while (isendx < obnixmax-1)
		           {
		                isendx = MIN (obnixmax-1,isbegx +nsgtracex-1);
		                isbegx = MAX (0,isendx-nsgtracex+1);

				itrc = 0;
				for (itrcy = isbegy; itrcy <= isendy; itrcy++)
				{
					for (itrcx = isbegx; itrcx <= isendx; itrcx++)
					{
					    for (isky = 0; isky < 60; isky++)
                        {
						    if (trcloc[itrcx][itrcy][isky] >= 0)
						    {
							    rtap  = MIN (1.0f, (float) (itrcx - isbegx + 1.0f) / 
								    (float) (nsgoverlapx));
							    rtap *= MIN (1.0f, (float) (isendx + 1 - itrcx) / 
								    (float) (nsgoverlapx));
							    rtap *= MIN (1.0f, (float) (itrcy - isbegy + 1.0f) / 
								    (float) (nsgoverlapy));
							    rtap *= MIN (1.0f, (float) (isendy + 1 - itrcy) / 
								    (float) (nsgoverlapy));

							    bsdg->taskin[bsdg->ntask * bsdg->tasksizemax + itrc] = 
							      trcloc[itrcx][itrcy][isky];
							    bsdg->taskou[bsdg->ntask * bsdg->tasksizemax + itrc] = 
							      trcloc[itrcx][itrcy][isky];
							    bsdg->taskweight[bsdg->ntask * bsdg->tasksizemax + itrc++] = rtap;
							    assert (0.0f < rtap  && rtap <= 1.0f);
							    assert (itrc <= bsdg->tasksizemax);
						    }
                        }
					}
				}
				bsdg->ntrcinpertask[bsdg->ntask] = itrc;
				bsdg->ntrcoupertask[bsdg->ntask] = itrc;
				bsdg->shotidpertask[bsdg->ntask] = ipky;

				bsdg->l1para.nsgtracexlastgate = MAX (bsdg->l1para.nsgtracexlastgate, isendx - isbegx);

				if (++(bsdg->ntask) >= ntaskmax) 
				{
					ntaskmax += ntaskinc;
					bsdg->taskin = (int *) realloc 
						(bsdg->taskin, ntaskmax * bsdg->tasksizemax * sizeof (int));
					bsdg->taskou = (int *) realloc 
						(bsdg->taskou, ntaskmax * bsdg->tasksizemax * sizeof (int));
					bsdg->taskweight = (float *) realloc 
						(bsdg->taskweight, ntaskmax * bsdg->tasksizemax * sizeof (float));
					memset (&bsdg->taskin[(ntaskmax - ntaskinc) * bsdg->tasksizemax], -1, 
						ntaskinc * bsdg->tasksizemax * sizeof (int));
					memset (&bsdg->taskou[(ntaskmax - ntaskinc) * bsdg->tasksizemax], -1, 
						ntaskinc * bsdg->tasksizemax * sizeof (int));
					memset (&bsdg->taskweight[(ntaskmax - ntaskinc) * bsdg->tasksizemax], -1, 
						ntaskinc * bsdg->tasksizemax * sizeof (float));

					bsdg->ntrcinpertask = (int *) realloc 
						(bsdg->ntrcinpertask, ntaskmax * sizeof (int));
					bsdg->ntrcoupertask = (int *) realloc 
						(bsdg->ntrcoupertask, ntaskmax * sizeof (int));
					bsdg->shotidpertask = (int *) realloc 
						(bsdg->shotidpertask, ntaskmax * sizeof (int));
				}

				isbegx = isbegx+iwinincx;
			} // while (isendx < obnixmax)
		      isbegy = isbegy+iwinincy;
	        } // while (isendy < obniymax)


		if (ipky == 0)
		{
			ntaskprv = 0;
			bsdg->ntaskpershot[ipky] = bsdg->ntask;
		}
		else
		{
			ntaskprv += bsdg->ntaskpershot[ipky - 1];
			bsdg->ntaskpershot[ipky] = bsdg->ntask - ntaskprv;
		}


	} // for (ipky = 0; ipky < bsdg->nipky; ipky++)

	bsdg->taskin     = (int *)   realloc (bsdg->taskin,     bsdg->ntask * bsdg->tasksizemax * sizeof (int));
	bsdg->taskou     = (int *)   realloc (bsdg->taskou,     bsdg->ntask * bsdg->tasksizemax * sizeof (int));
	bsdg->taskweight = (float *) realloc (bsdg->taskweight, bsdg->ntask * bsdg->tasksizemax * sizeof (float));

	bsdg->ntrcinpertask = (int *) realloc (bsdg->ntrcinpertask, bsdg->ntask * sizeof (int));
	bsdg->ntrcoupertask = (int *) realloc (bsdg->ntrcoupertask, bsdg->ntask * sizeof (int));
	bsdg->shotidpertask = (int *) realloc (bsdg->shotidpertask, bsdg->ntask * sizeof (int));

	for (itrcx = 0; itrcx < obnixmax; itrcx++)
    {
	    for (itrcy = 0; itrcy < obniymax; itrcy++)
		    free (trcloc[itrcx][itrcy]);
        free (trcloc[itrcx]);
    }
	free (trcloc);

	return;
}

void check_consistency_input_and_notional (master_t *msdginterp)
{
    CapiID srcinst = msdginterp->dbio->srcinst;
    Int32 status, datatype, length, samplerate, keynumber;
    Int64 *groupheader, *uniqheader, itrc, ntrcsrc, ipky = 1, nipkysrc = 1;

    status = cCapi_hasAttribute (srcinst, msdginterp->ipky.name);
    if (status != FCAPI_SUCCESS)
        G2_logCritical (msdginterp->logName, "Notional source does not have [%s] attribute", msdginterp->ipky.name);

    status = cCapi_getAttributeByName (srcinst, msdginterp->ipky.name, &datatype, &length, &samplerate, &keynumber);
    if (status != FCAPI_SUCCESS || keynumber != 1) // fail if ipky.name is not a group name for notional data
        G2_logCritical (msdginterp->logName, "The [%s] attribute is not a group name for notional", msdginterp->ipky.name);

    if ( cCapi_getNInstances (srcinst, &ntrcsrc) != FCAPI_SUCCESS)
      G2_logCritical(msdginterp->logName, "Fail to get number of instances from source wavelet");

    groupheader = (Int64 *) calloc (ntrcsrc, sizeof (Int64));
    if ( cCapi_readConvertByName (srcinst, msdginterp->ipky.name, groupheader, FCAPI_INT64, 1, ntrcsrc) != FCAPI_SUCCESS)
      G2_logCritical(msdginterp->logName, "Fail to read primary key %s from source wavelet", msdginterp->ipky.name);

    for (itrc = 1; itrc < ntrcsrc; itrc++)
    {
        if (groupheader[itrc] > groupheader[itrc - 1])
            nipkysrc++;
        else if (groupheader[itrc] < groupheader[itrc - 1])
            G2_logCritical (msdginterp->logName, "Data order of notional source does not have an increasing order.");
    }
    G2_logInfo (msdginterp->logName, "Number of input groups = %d Number of notional groups = %lld", msdginterp->nipky, nipkysrc);
    if (msdginterp->lshotbyshot == YESYES && msdginterp->ldefaultsrc == 0 && msdginterp->nipky != nipkysrc)
        G2_logCritical (msdginterp->logName, "Number of input groups (%d) does not match with number of notional groups (%lld)", msdginterp->nipky, nipkysrc);
    if (msdginterp->lshotbyshot == NONO && nipkysrc != 1)
        G2_logCritical (msdginterp->logName, "Number of notional groups (%lld) must be 1 for shotbyshot = NO", nipkysrc);

    uniqheader = (Int64 *) calloc (nipkysrc, sizeof (Int64));
    uniqheader[0] = groupheader[0];
    for (itrc = 1; itrc < ntrcsrc; itrc++)
        if (groupheader[itrc] > groupheader[itrc - 1])
            uniqheader[ipky++] = groupheader[itrc];

    if (msdginterp->lshotbyshot == YESYES && msdginterp->ldefaultsrc == 0)
        for (ipky = 0; ipky < msdginterp->nipky; ipky++)
            if (msdginterp->pkyno[ipky] != uniqheader[ipky])
                G2_logCritical (msdginterp->logName, "Group headers do not match between input at %s = %d and notional source at %s = %lld.", 
                    msdginterp->ipky.name, msdginterp->pkyno[ipky], msdginterp->ipky.name, uniqheader[ipky]);

    free (uniqheader);
    free (groupheader);
}

void data_readin_deghost3d_signature (int fsignature, master_t *msdginterp)
{
    check_consistency_input_and_notional (msdginterp);

    // read notional source and target wavelet
    Int64 ntrcsrc;
    Int32 status, nipky_notion, dtype, nsampsrc, nsamptgt, sratesrc, sratetgt, keynum, trlen,nipky_real, isamp;
    Float32 *pShotx, *pShoty, *pShotz, *pTdelay, *pGunamp, *pTarget;
    CapiID srcinst = msdginterp->dbio->srcinst, targetinst = msdginterp->dbio->targetinst; 
    CapiID defsrcinst = msdginterp->dbio->defsrcinst;

    // get source wavelet //
    if ( cCapi_getNInstances (srcinst, &ntrcsrc) != FCAPI_SUCCESS)
      G2_logCritical(msdginterp->logName, "Fail to get the number of instances from source wavelet");

    if(msdginterp->dataorderfp == 0) nipky_real = msdginterp->nipky;
    else nipky_real = msdginterp->pkyindex[msdginterp->nipky-1];
    
    if (msdginterp->lshotbyshot == YESYES && msdginterp->ldefaultsrc == 0 
	&& ntrcsrc % nipky_real != 0)
    {
        G2_logCritical (msdginterp->logName, "Number of traces for source wavelet must be a multiple of number of shots!!!!\n#src = %lld #shots = %d", ntrcsrc, nipky_real);
    }
    if (msdginterp->ldefaultsrc == 0)
      msdginterp->l1para.ntrcsrc = (msdginterp->lshotbyshot) ? (int) ntrcsrc / nipky_real : (int) ntrcsrc;
    else
      {
	Int64 ntrcdefsrc;
	if ( cCapi_getNInstances (defsrcinst, &ntrcdefsrc) != FCAPI_SUCCESS)
          G2_logCritical(msdginterp->logName, "Fail to get the number of instances from default source wavelet");
	msdginterp->l1para.ntrcsrc = (int)ntrcdefsrc;
      }

    if (msdginterp->ldefaultsrc == 0)
      nipky_notion = (msdginterp->lshotbyshot) ? nipky_real : 1;
    else
      {
	Int64 irg = 0;
	Int64 irgs = 0;

	if ( cCapi_getGroupsAndSizesByName(srcinst,msdginterp->ipky.name,&irg,&irgs,FCAPI_INCREASING) != FCAPI_SUCCESS)
          G2_logCritical(msdginterp->logName, "Fail to get the group and sizes of source wavelet");

	nipky_notion=I8_n(irg);
      }
    
    if ( cCapi_getAttributeByName (srcinst, "TRACE", &dtype, &nsampsrc, &sratesrc, &keynum) != FCAPI_SUCCESS)
      G2_logCritical(msdginterp->logName, "Fail to read TRACE vector from source wavelet");
    G2_logInfo (msdginterp->logName, "SOURCE: #traces %d #sample %d sample rate %d #notional %d NHEAD %d",
                                     msdginterp->l1para.ntrcsrc, nsampsrc, sratesrc, nipky_notion, NHEAD);
    if ( sratesrc != msdginterp->srate )
      G2_logCritical(msdginterp->logName, "Must have source srate %d = input srate %d", 
                     sratesrc, msdginterp->srate); 

    trlen = NHEAD + nsampsrc;
    msdginterp->l1para.nsampsrc = nsampsrc;

    // get target wavelet //
    if (msdginterp->l1para.target_type == TARGETTYPE_WAVELET)
    {
        if ( cCapi_getNInstances (targetinst, &ntrcsrc) != FCAPI_SUCCESS)
          G2_logCritical(msdginterp->logName, "Fail to get the number of instances from target wavelet");
        if (ntrcsrc != 1)
        {
          G2_logCritical (msdginterp->logName, "Target wavelet must have 1 trace!!!!");
        }
        if ( cCapi_getAttributeByName (targetinst, "TRACE", &dtype, &nsamptgt, &sratetgt, &keynum) != FCAPI_SUCCESS)
          G2_logCritical(msdginterp->logName, "Fail to read TRACE vector from target wavelet");
        G2_logInfo(msdginterp->logName, "TARGET: #traces  1 #sample %d sample rate %d", nsamptgt, sratetgt);

        if ( !(nsampsrc == nsamptgt && sratesrc <= msdginterp->srate && 
            sratesrc == sratetgt && sratetgt == msdginterp->srate) )
          G2_logCritical(msdginterp->logName, 
                         "Must have source wavelet nsamp (=%d) equal to target nsamp (=%d). Must have source srate %d = target srate %d = input srate %d", 
                         nsampsrc, nsamptgt, sratesrc, sratetgt, msdginterp->srate); 
    }
    else
    {
        nsamptgt = nsampsrc;
    }

    pShotx  = (Float32 *) calloc (msdginterp->l1para.ntrcsrc * nipky_notion, sizeof (Float32));
    pShoty  = (Float32 *) calloc (msdginterp->l1para.ntrcsrc * nipky_notion, sizeof (Float32));
    pShotz  = (Float32 *) calloc (msdginterp->l1para.ntrcsrc * nipky_notion, sizeof (Float32));
    pTdelay = (Float32 *) calloc (msdginterp->l1para.ntrcsrc * nipky_notion, sizeof (Float32));
    pGunamp = (Float32 *) calloc (msdginterp->l1para.ntrcsrc * nipky_notion, sizeof (Float32));
    pTarget = (Float32 *) calloc (nsamptgt, sizeof (Float32));

    msdginterp->psignature = (Float32 *) calloc ((msdginterp->l1para.ntrcsrc + 1) * trlen, sizeof (Float32));

    status = cCapi_hasAttribute (srcinst, msdginterp->srcshotx.name);
    if (status != FCAPI_SUCCESS)
        G2_logCritical (msdginterp->logName, "[%s] from notional source does not exist", msdginterp->srcshotx.name);
    status = cCapi_hasAttribute (srcinst, msdginterp->srcshoty.name);
    if (status != FCAPI_SUCCESS)
        G2_logCritical (msdginterp->logName, "[%s] from notional source does not exist", msdginterp->srcshoty.name);
    status = cCapi_hasAttribute (srcinst, msdginterp->srcshotz.name);
    if (status != FCAPI_SUCCESS)
        G2_logCritical (msdginterp->logName, "[%s] from notional source does not exist", msdginterp->srcshotz.name);

    if ( cCapi_readConvertByName (srcinst, msdginterp->srcshotx.name, pShotx,  FCAPI_FLOAT32, 1, msdginterp->l1para.ntrcsrc * nipky_notion) != FCAPI_SUCCESS)
      G2_logCritical(msdginterp->logName, "Fail to read shotx attribute %s from source wavelet", msdginterp->srcshotx.name);
    if ( cCapi_readConvertByName (srcinst, msdginterp->srcshoty.name, pShoty,  FCAPI_FLOAT32, 1, msdginterp->l1para.ntrcsrc * nipky_notion) != FCAPI_SUCCESS)
      G2_logCritical(msdginterp->logName, "Fail to read shoty attribute %s from source wavelet", msdginterp->srcshoty.name);
    if ( cCapi_readConvertByName (srcinst, msdginterp->srcshotz.name, pShotz,  FCAPI_FLOAT32, 1, msdginterp->l1para.ntrcsrc * nipky_notion) != FCAPI_SUCCESS)
      G2_logCritical(msdginterp->logName, "Fail to read shotz attribute %s from source wavelet", msdginterp->srcshotz.name);
    if (!EQSTR (msdginterp->gunamp.name, "NONE"))
    {
        if ( cCapi_readConvertByName (srcinst, msdginterp->gunamp.name, pGunamp, FCAPI_FLOAT32, 1, msdginterp->l1para.ntrcsrc * nipky_notion) != FCAPI_SUCCESS)
          G2_logCritical(msdginterp->logName, "Fail to read gun amp attribute %s from source wavelet");
    }
    if (!EQSTR (msdginterp->tdelay.name, "NONE"))
    {
        if ( cCapi_readConvertByName (srcinst, msdginterp->tdelay.name, pTdelay, FCAPI_FLOAT32, 1, msdginterp->l1para.ntrcsrc * nipky_notion) != FCAPI_SUCCESS)
          G2_logCritical(msdginterp->logName, "Fail to read delay time attribute %s from source wavelet", msdginterp->tdelay.name);
    }

    if (msdginterp->l1para.target_type == TARGETTYPE_WAVELET)
    {
        if ( cCapi_readConvertByName (targetinst, "TRACE", pTarget, FCAPI_FLOAT32, 1, 1) != FCAPI_SUCCESS)
          G2_logCritical(msdginterp->logName, "Fail to read TRACE vector from target wavelet");
    }

    Int32 nbkey = 1000, nakey = (Int32) ceilf (msdginterp->l1para.ntrcsrc * nipky_notion / (float) nbkey);
    if (msdginterp->l1para.ntrcsrc * nipky_notion > nakey * nbkey)
      G2_logCritical(msdginterp->logName, "Number of instances in source wavelet must be <= %s x %s", msdginterp->l1para.ntrcsrc*nipky_notion, nakey, nbkey);
    Float32 *pSource = (Float32 *) calloc (nbkey * nsampsrc, sizeof (Float32));
    Float32 *pBuffer = (Float32 *) calloc (nbkey * trlen,    sizeof (Float32));
    Float32 *des;
    head_t *src;

    int tapend1 = MAX ((int) (100 * 1000 / msdginterp->srate), 10);
    int tapend2 = MAX ((int) (20  * 1000 / msdginterp->srate),  2);
    for (int ikey = 0; ikey < nakey; ikey++)
    {
        int first = ikey * nbkey + 1;
        int tlast = MIN ((ikey + 1) * nbkey, msdginterp->l1para.ntrcsrc * nipky_notion);
        if ( cCapi_readConvertByName (srcinst, "TRACE", pSource, FCAPI_FLOAT32, first, tlast) != FCAPI_SUCCESS)
          G2_logCritical(msdginterp->logName, "Fail to read TRACE vector from source wavelet");

        //// linear taper 100 ms at the end of trace length to avoid artifacts in autocorrelation QC ////
        for (int itrc = 0; itrc < tlast - first + 1; itrc++)
        {
            for (int isamp = msdginterp->l1para.nsampsrc - tapend1; isamp < msdginterp->l1para.nsampsrc - tapend2; isamp++)
                pSource[itrc * msdginterp->l1para.nsampsrc + isamp] *= 
                    (float) (msdginterp->l1para.nsampsrc - isamp - tapend2 - 1) / (tapend1 - tapend2);
                //float rtap = (float) (msdginterp->l1para.nsampsrc - isamp - tapend2 - 1) / (tapend1 - tapend2);
            for (int isamp = msdginterp->l1para.nsampsrc - tapend2; isamp < msdginterp->l1para.nsampsrc; isamp++)
                pSource[itrc * msdginterp->l1para.nsampsrc + isamp] = 0.0f;
        }
        //// linear taper 100 ms at the end of trace length to avoid artifacts in autocorrelation QC ////

        for (int j = 0; j < tlast - first + 1; j++)
        {
            src = (head_t *) (pBuffer + j * trlen);
            src->shotx  = pShotx [first - 1 + j];
            src->shoty  = pShoty [first - 1 + j];
            src->shotz  = pShotz [first - 1 + j];
            src->tdelay = pTdelay[first - 1 + j] * 0.001f; //unit is second
            src->gunamp = pGunamp[first - 1 + j];

            des = (float *) (pSource + j * nsampsrc);
            memcpy ((float *) src + NHEAD, des, nsampsrc * sizeof (float));
        }

	//If shotbyshot, write source to disk file. Otherwise copy it to buffer psignature
	if (msdginterp->lshotbyshot == YESYES)
        {
	    ssize_t nwrite = write (fsignature, pBuffer, (size_t) (tlast - first + 1) * trlen * sizeof (float));
	    if (nwrite != (tlast - first + 1) * trlen * sizeof (float))
              G2_logCritical(msdginterp->logName, "Fail write source wavelet to file. Should write %zu, but wrote only %zu", 
                                                  (tlast - first + 1) * trlen * sizeof (float), nwrite);
        }
	else
        {
	    memcpy (msdginterp->psignature+(first-1)*trlen, pBuffer, (size_t) (tlast - first + 1) * trlen * sizeof (float));
        }
    }

    if (msdginterp->ldefaultsrc == 1) func_insert_defsource(msdginterp, fsignature);

    //write target trace to buffer psignature
    memcpy (msdginterp->psignature + msdginterp->l1para.ntrcsrc * trlen + NHEAD,
        pTarget, (size_t) nsamptgt * sizeof (float));

    free (pSource); free (pBuffer); free (pTarget); free (pShotx); free (pShoty); free (pShotz);
    free (pGunamp); free (pTdelay);
}

void func_insert_defsource(master_t *msdginterp, int fsignature)
{

  int i,j,nipkysrc,status;

  CapiID srcinst = msdginterp->dbio->srcinst;
  CapiID oinst = msdginterp->dbio->oinst;

  //Get unique list of skey
  Int64 irg = 0;
  Int64 irgs = 0;

  int     trlensig  = NHEAD + msdginterp->l1para.nsampsrc;
  Float32 *pdefSource = (Float32 *) calloc (msdginterp->l1para.ntrcsrc * 1 * trlensig,    sizeof (Float32));

  func_defsource_read(msdginterp,pdefSource);

  int *pkeyval    = (int *) malloc((size_t)msdginterp->nipky*(size_t)sizeof(int));
  int *pkeyvalsrc = (int *) malloc((size_t)msdginterp->nipky*(size_t)sizeof(int));
  int *map        = (int *) malloc((size_t)msdginterp->nipky*(size_t)sizeof(int));
  int *missing    = (int *) malloc((size_t)msdginterp->nipky*(size_t)sizeof(int));

  if ( cCapi_getGroupsAndSizesByName(oinst,msdginterp->ipky.name,&irg,&irgs,FCAPI_INCREASING) != FCAPI_SUCCESS)
    G2_logCritical(msdginterp->logName, "Fail to get groups and sizes from input");

  for(i=0; i<msdginterp->nipky; i++) {
    pkeyval[i] = I8_geti(irg);
  }
  S8_remove(irg);
  S8_remove(irgs);

  irg = 0;
  irgs = 0;

  if ( cCapi_getGroupsAndSizesByName(srcinst,msdginterp->ipky.name,&irg,&irgs,FCAPI_INCREASING) != FCAPI_SUCCESS)
    G2_logCritical(msdginterp->logName, "Fail to get groups and sizes from source wavelet");

  nipkysrc=I8_n(irg);

  for(i=0; i<nipkysrc; i++) {
    pkeyvalsrc[i] = I8_geti(irg);
  }
  S8_remove(irg);
  S8_remove(irgs);

  for(i=0; i<msdginterp->nipky; i++) 
    {
      missing[i] = 1;
      map[i] = -1;
      for(j=0; j<nipkysrc; j++) 
	{
	  
	  if (pkeyvalsrc[j] == pkeyval[i])
	    {
	      map[i] = j;
	      missing[i] = 0;
	    }

	}
    }

  ssize_t buff_size = msdginterp->l1para.ntrcsrc * trlensig * sizeof (float);
  off_t   buff_loc;
  Float32 *pbuff    = (Float32 *) malloc (buff_size);
  ssize_t nread,nwrite;     

  char *oname_s =  (char *) malloc (LMSDGINTERP_LNAME*sizeof(char));
  int fsignature1 = opentmpfile("msdginterpsignature1", oname_s,  LMSDGINTERP_LNAME);

  for(i=0; i<msdginterp->nipky; i++) 
    {
      if (missing[i] == 0)
	{
	  j = map[i];
	  buff_loc  = msdginterp->l1para.ntrcsrc * trlensig * sizeof (float) * j;

	  nread     = pread (fsignature, pbuff, buff_size, buff_loc);
	  if (nread != buff_size) 
            G2_logCritical(msdginterp->logName, "Fail to read from source signature file. Should read %zu, but read only %zu", buff_size, nread);
	}
      else
	memcpy (pbuff, pdefSource, buff_size);

      buff_loc  = msdginterp->l1para.ntrcsrc * trlensig * sizeof (float) * i;
      
      nwrite     = pwrite (fsignature1, pbuff, buff_size, buff_loc);
      if (nwrite != buff_size)
        G2_logCritical(msdginterp->logName, "Fail to write source signature to file. Should write %zu, but wrote only %zu", buff_size, nwrite);
  
    }

  for(i=0; i<msdginterp->nipky; i++) 
    {
      buff_loc  = msdginterp->l1para.ntrcsrc * trlensig * sizeof (float) * i;
     
      nread     = pread  (fsignature1, pbuff, buff_size, buff_loc);

      nwrite    = pwrite (fsignature,  pbuff, buff_size, buff_loc);
    }

  removetmpfile(fsignature1,oname_s);
  free(oname_s);

  free (pdefSource);
  free (pbuff);
  free(pkeyval);
  free(pkeyvalsrc);
  free(map);
  free(missing);

}


void func_defsource_read(master_t *msdginterp, float *pdefSource)
{

  int i,j,itrc,status;

  Float32 *des;
  head_t *src;

  Float32 *pShotx, *pShoty, *pShotz, *pTdelay, *pGunamp;
  
  int nipky_notion = 1; //// for the default source signature

  CapiID defsrcinst = msdginterp->dbio->defsrcinst;

  pShotx  = (Float32 *) calloc (msdginterp->l1para.ntrcsrc * nipky_notion, sizeof (Float32));
  pShoty  = (Float32 *) calloc (msdginterp->l1para.ntrcsrc * nipky_notion, sizeof (Float32));
  pShotz  = (Float32 *) calloc (msdginterp->l1para.ntrcsrc * nipky_notion, sizeof (Float32));
  pTdelay = (Float32 *) calloc (msdginterp->l1para.ntrcsrc * nipky_notion, sizeof (Float32));
  pGunamp = (Float32 *) calloc (msdginterp->l1para.ntrcsrc * nipky_notion, sizeof (Float32));
  
  if ( cCapi_readConvertByName (defsrcinst, msdginterp->srcshotx.name, pShotx,  FCAPI_FLOAT32, 1, msdginterp->l1para.ntrcsrc * 1) != FCAPI_SUCCESS)
    G2_logCritical(msdginterp->logName, "Fail to read shotx attribute %s from default source wavelet", msdginterp->srcshotx.name);
  if ( cCapi_readConvertByName (defsrcinst, msdginterp->srcshoty.name, pShoty,  FCAPI_FLOAT32, 1, msdginterp->l1para.ntrcsrc * 1) != FCAPI_SUCCESS)
    G2_logCritical(msdginterp->logName, "Fail to read shoty attribute %s from default source wavelet", msdginterp->srcshoty.name);
  if ( cCapi_readConvertByName (defsrcinst, msdginterp->srcshotz.name, pShotz,  FCAPI_FLOAT32, 1, msdginterp->l1para.ntrcsrc * 1) != FCAPI_SUCCESS)
    G2_logCritical(msdginterp->logName, "Fail to read shotz attribute %s from default source wavelet", msdginterp->srcshotz.name);
  if (!EQSTR (msdginterp->gunamp.name, "NONE"))
    {
      if ( cCapi_readConvertByName (defsrcinst, msdginterp->gunamp.name, pGunamp, FCAPI_FLOAT32, 1, msdginterp->l1para.ntrcsrc * 1) != FCAPI_SUCCESS)
        G2_logCritical(msdginterp->logName, "Fail to read gun amp attribute %s from default source wavelet", msdginterp->gunamp.name);
    }
  if (!EQSTR (msdginterp->tdelay.name, "NONE"))
    {
      if ( cCapi_readConvertByName (defsrcinst, msdginterp->tdelay.name, pTdelay, FCAPI_FLOAT32, 1, msdginterp->l1para.ntrcsrc * 1) != FCAPI_SUCCESS)
        G2_logCritical(msdginterp->logName, "Fail to read delay time attribute %s from default source wavelet", msdginterp->tdelay.name);
    }

  Float32 *pSource      = (Float32 *) calloc (msdginterp->l1para.ntrcsrc * 1 * msdginterp->l1para.nsampsrc, sizeof (Float32));
  
  if ( cCapi_readConvertByName (defsrcinst,    "TRACE", pSource, FCAPI_FLOAT32, 1, msdginterp->l1para.ntrcsrc) != FCAPI_SUCCESS)
    G2_logCritical(msdginterp->logName, "Fail to read TRACE vector from default source wavelet");

  int tapend1 = MAX ((int) (100 * 1000 / msdginterp->srate), 10);
  int tapend2 = MAX ((int) (20  * 1000 / msdginterp->srate),  2);
  
  //// linear taper 100 ms at the end of trace length to avoid artifacts in autocorrelation QC ////
  for (int itrc = 0; itrc < msdginterp->l1para.ntrcsrc; itrc++)
    {
      for (int isamp = msdginterp->l1para.nsampsrc - tapend1; isamp < msdginterp->l1para.nsampsrc - tapend2; isamp++)
	pSource[itrc * msdginterp->l1para.nsampsrc + isamp] *= 
	  (float) (msdginterp->l1para.nsampsrc - isamp - tapend2 - 1) / (tapend1 - tapend2);
     
      for (int isamp = msdginterp->l1para.nsampsrc - tapend2; isamp < msdginterp->l1para.nsampsrc; isamp++)
	pSource[itrc * msdginterp->l1para.nsampsrc + isamp] = 0.0f;
    }
  //// linear taper 100 ms at the end of trace length to avoid artifacts in autocorrelation QC ////

  int trlen = NHEAD + msdginterp->l1para.nsampsrc;

  for (int j = 0; j < msdginterp->l1para.ntrcsrc; j++)
    {
      src = (head_t *) (pdefSource + j * trlen);
      src->shotx  = pShotx [j];
      src->shoty  = pShoty [j];
      src->shotz  = pShotz [j];
      src->tdelay = pTdelay[j] * 0.001f; //unit is second
      src->gunamp = pGunamp[j];
      
      des = (float *) (pSource + j * msdginterp->l1para.nsampsrc);
      memcpy ((float *) src + NHEAD, des, msdginterp->l1para.nsampsrc * sizeof (float));
    }

  free (pSource);

  free (pShotx); free (pShoty); free (pShotz);
  free (pGunamp); free (pTdelay);


}

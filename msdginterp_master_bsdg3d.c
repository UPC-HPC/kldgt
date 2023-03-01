#include "msdginterp.h"
#include <time.h>
#include <string.h>
#include <errno.h>
#include <sys/mman.h>
#include "bsdg_utils.h"
#include "TagDefs.hpp"
#include "G2WorkProgress.ci"

void next_task (int itask, int *ppky, int fagcscl, int forigin, int forigin_az, int forigin_ay, int fagc2d,
		int fagc2d_az,int fagc2d_ay,
		int fsignature,float *pOrigin, float *pTemp, float *pTemp_az, float *pTemp_ay, master_t *msdginterp);

void display_progress (int count, int ntask, float *progress, float increment, 
		       int task, int rank, int send_or_receive, char *logname);
void save_task_to_disk (int tag, float *pDeghost, master_t *msdginterp);
void *write_complete_gather (void *arguments);
void *large_alloc (size_t size, int *fd);
void large_free (void *ptr, size_t size, int fd);
void myread (int fd, void *buffer, size_t size, off_t offset);
void mywrite (int fd, void *buffer, size_t size, off_t offset);

struct thrd_args {
  master_t *msdginterp;
  int fagcscl;
  int fagc2d;
  int fagc2d_az;
  int fagc2d_ay;
};

void msdginterp_master_bsdg3d (int forigin, int fsignature,
			       master_t *msdginterp, MPI_Comm MPI_LEAD_WORLD, int leadrank, int leadsize)
{
  int slaverank, itask = 0, recv_count = 0, ipky = -1, fdTemp, fagc2d, fagc2d_az, fagc2d_ay;
  float *pOrigin, *pDeghost, *pTemp, *pTemp_az, *pTemp_ay;
  float fagcscl, send_complete = 0.0, recv_complete = 0.0, inc_complete = 5.0;
  char *onameagc = (char *) malloc (LMSDGINTERP_LNAME * sizeof (char));
  char *onameagc2d = (char *) malloc (LMSDGINTERP_LNAME * sizeof (char));
  char *onameagc2d_az = (char *) malloc (LMSDGINTERP_LNAME * sizeof (char));
  char *onameagc2d_ay = (char *) malloc (LMSDGINTERP_LNAME * sizeof (char));
  	
  MPI_Status mpistatus;

  int forigin_az = -1;
  int forigin_ay = -1;

  msdginterp->nsamporig = msdginterp->nsamp;

  if (msdginterp->nshiftdata > 0)
    msdginterp->nsamp = msdginterp->nsamp + msdginterp->nshiftdata;
  
  size_t trlen = NHEAD + msdginterp->nsamp;
  size_t trlenorig = NHEAD + msdginterp->nsamporig;

  fagcscl = opentmpfile ("fagcscl", onameagc, LMSDGINTERP_LNAME);
  fagc2d = opentmpfile ("fagc2d", onameagc2d, LMSDGINTERP_LNAME);
  fagc2d_az = opentmpfile ("fagc2d_az", onameagc2d_az, LMSDGINTERP_LNAME);
  fagc2d_ay = opentmpfile ("fagc2d_ay", onameagc2d_ay, LMSDGINTERP_LNAME);

  lseek (forigin,  0, SEEK_SET);

  pthread_t output_thrd;
  struct thrd_args args;
  args.msdginterp    = msdginterp;   
  args.fagcscl = fagcscl;
  args.fagc2d = fagc2d;
  args.fagc2d_az = fagc2d_az;
  args.fagc2d_ay = fagc2d_ay;
  pthread_create (&output_thrd, NULL, &write_complete_gather, &args);

  if (msdginterp->dgmethod == DGMETHOD_VSP)
    {
      msdginterp->nipky = 1;
      msdginterp->nsky[0] = msdginterp->ninst;
      msdginterp->nisky   = msdginterp->ninst;
    }

  pTemp    = (float *) large_alloc (msdginterp->nisky * trlen * sizeof (float), &fdTemp);

  size_t buffsize = msdginterp->l1para.ndata * msdginterp->tasksizemax * trlen;
  //Read the source data too
  if (msdginterp->l1para.choosemethod == CHOOSEMETHOD_JOINTSR3D &&
      msdginterp->l1para.jointmethod != JOINTDEGHOST && msdginterp->lshotbyshot == YESYES)
    buffsize += msdginterp->l1para.ntrcsrc * (NHEAD + msdginterp->l1para.nsampsrc); 
  
  pOrigin  = (float *) malloc (buffsize * sizeof (float));
  pDeghost = (float *) malloc (msdginterp->tasksizemax * trlen * sizeof (float));

  G2startphase(&msdginterp->wpid, msdginterp->ntask,"Data send/recv and deghosting",GXPHASE,1);   
  if (  msdginterp->g2progress)
    G2_logState("bsdg", G2LOG_STATE_PROCESSING, "started bsdg processing");

  //Free nodes that will not be used
  for (slaverank = msdginterp->ntask+1;slaverank<leadsize;slaverank++)
    {
      //Send end of data message and free the node 
      MPI_Send (pOrigin, buffsize, MPI_FLOAT, slaverank, msdginterp->ntask, MPI_LEAD_WORLD);
      cl_free_node(slaverank);
    }

  G2updatenumerictag(&msdginterp->wpid, TAG_TASKS, msdginterp->ntask);

  G2_logInfo(msdginterp->logName, "Starting to read in and send %d %s (%6d tasks) to slaves!", 
                        	  msdginterp->nipky, msdginterp->ipky.name, msdginterp->ntask);

  if (msdginterp->l1para.agc2dflag == YESYES)
    {
      msdginterp->l1para.scalaragc2d = (float *) calloc ((size_t) msdginterp->nsamp, sizeof (float));
      agc2d_sgatexy (msdginterp, forigin, fagcscl); 
      free (msdginterp->l1para.scalaragc2d);
    }

  lseek (fagcscl, 0, SEEK_SET);
  lseek (fagc2d, 0, SEEK_SET);
  lseek (fagc2d_az, 0, SEEK_SET);
  lseek (fagc2d_ay, 0, SEEK_SET);

  for (slaverank = 1; slaverank < MIN (leadsize, msdginterp->ntask + 1); slaverank++)
    {
      display_progress (itask, msdginterp->ntask, &send_complete, inc_complete, itask, slaverank, 1, msdginterp->logName);

      next_task (itask, &ipky, fagcscl, forigin, forigin_az, forigin_ay, fagc2d,fagc2d_az,fagc2d_ay, fsignature, pOrigin, pTemp, pTemp_az, pTemp_ay, msdginterp);

      MPI_Send (pOrigin, buffsize, MPI_FLOAT, slaverank, itask++, MPI_LEAD_WORLD);
    }

  while (itask < msdginterp->ntask)
    {
      // receive data from slave //
      MPI_Recv (pDeghost, msdginterp->tasksizemax * trlen, MPI_FLOAT, MPI_ANY_SOURCE,
		MPI_ANY_TAG, MPI_LEAD_WORLD, &mpistatus);

      if ((recv_count+1) * 100.0f / msdginterp->ntask >=  recv_complete){
	G2updateprogress(&msdginterp->wpid, recv_count+1, GXPHASE, 1);
        if (  msdginterp->g2progress)
          G2_logProgress("bsdg", G2LOG_STATE_PROCESSING, "bsdg processing", msdginterp->ntask, recv_count+1, "processing windows", 1, 1, "bsdg");
      }
      
      display_progress (++recv_count, msdginterp->ntask, &recv_complete, inc_complete, 
			mpistatus.MPI_TAG, mpistatus.MPI_SOURCE, 0, msdginterp->logName);
      

      // send data to next available slave //
      slaverank = mpistatus.MPI_SOURCE;
      display_progress (itask, msdginterp->ntask, &send_complete, inc_complete, itask, slaverank, 1, msdginterp->logName);

      next_task (itask, &ipky, fagcscl, forigin, forigin_az, forigin_ay, fagc2d,fagc2d_az,fagc2d_ay, fsignature, pOrigin, pTemp, pTemp_az, pTemp_ay, msdginterp);

      MPI_Send (pOrigin, buffsize, MPI_FLOAT, slaverank, itask++, MPI_LEAD_WORLD);

      // save received data to disk //
      save_task_to_disk (mpistatus.MPI_TAG, pDeghost, msdginterp);
    }

  for (slaverank = 1; slaverank < MIN (leadsize, msdginterp->ntask + 1); slaverank++)
    {
      // receive data from slave //
      MPI_Recv (pDeghost, msdginterp->tasksizemax * trlen, MPI_FLOAT, MPI_ANY_SOURCE,
		MPI_ANY_TAG, MPI_LEAD_WORLD, &mpistatus);

      if ((recv_count+1) * 100.0f / msdginterp->ntask >=  recv_complete){
	G2updateprogress(&msdginterp->wpid, recv_count+1, GXPHASE, 1);
        if (  msdginterp->g2progress)
           G2_logProgress("bsdg", G2LOG_STATE_PROCESSING, "bsdg processing", msdginterp->ntask, recv_count+1, "processing windows", 1, 1, "bsdg");
      }
      
      display_progress (++recv_count, msdginterp->ntask, &recv_complete, inc_complete, 
			mpistatus.MPI_TAG, mpistatus.MPI_SOURCE, 0, msdginterp->logName);
      
      //Send end of message 
      MPI_Send (pOrigin, buffsize, MPI_FLOAT, mpistatus.MPI_SOURCE, itask, MPI_LEAD_WORLD);
      cl_free_node(mpistatus.MPI_SOURCE);

      // save received data to disk //
      save_task_to_disk (mpistatus.MPI_TAG, pDeghost, msdginterp);
    }

  G2endphase(&msdginterp->wpid,GXPHASE,1);  
  if (  msdginterp->g2progress)
     G2_logState("bsdg", G2LOG_STATE_PROCESSING, "finished bsdg processing");

  pthread_join (output_thrd, NULL);

  removetmpfile (fagcscl, onameagc);
  removetmpfile (fagc2d, onameagc2d);
  removetmpfile (fagc2d_az, onameagc2d_az);
  removetmpfile (fagc2d_ay, onameagc2d_ay);
  free (onameagc);
  free (onameagc2d);
  free (onameagc2d_az);
  free (onameagc2d_ay);

  free (pOrigin);
  free (pDeghost);
  large_free (pTemp, msdginterp->nisky * trlen * sizeof (float), fdTemp);

  return;
}

void agc2d_sgatexy(master_t *msdginterp, int fOrigin, int fScalar) 
{
  
  int i, j, l, isky, jsky, ntracewin, itrc, isamp, trlen, pklen, nsamp;

  int isgatex,isbegx,isendx,nsgatex;

  int startidx,endidx,outovlapx;

  int obnkeyx,obnkeyy;

  int wbmin,ntocopy;
  size_t nread,nwrite,n2read, nx2read, n2write,nx2write, nRemain, nDoneSize, batchsize;

  float *pOrigin,*pPack,*outputscl,*trace, *temp;
  float **input,**winoutscl;
  int *mapkey,*mapkeyx, *wbval, *ishift;;

  float rstartidx,rendidx,rtap,rval;

  head_t *src, *des;

  size_t offset,offsetou;

  l1inv_t *l1para;
  
  l1para = &msdginterp->l1para;
  l1para->nsamp = msdginterp->nsamp;
  l1para->srate = msdginterp->srate;

  pklen = NHEAD + msdginterp->vlen;
  trlen = NHEAD + msdginterp->nsamp;

  float obnxposmin = l1para->obnxposmin;
  float obnxposmax = l1para->obnxposmax;
  float obnyposmin = l1para->obnyposmin;
  float obnyposmax = l1para->obnyposmax;
  int   obnixmax   = l1para->obnixmax;
  
  nsgatex=MAX(1,(int)((obnixmax)/l1para->winxagc2d));    
  outovlapx = l1para->ovlapxagc2d;

  trace =   (float*)malloc(trlen*sizeof(float));

  offset   = 0;
  offsetou = 0;
  
  if(l1para->agc2dmethod ==2)
    {
      wbval  = (int *) malloc (msdginterp->nisky * sizeof (float));
      temp   = (float *) malloc (msdginterp->nsamp * sizeof (float));
      ishift = (int *) malloc (msdginterp->nisky * sizeof (float));
    }

  for(i=0; i<msdginterp->nipky; i++)
   {

     pPack =     (float*)calloc(msdginterp->nisky*pklen, sizeof(float));
     pOrigin =   (float*)calloc(msdginterp->nisky*trlen, sizeof(float));
     outputscl = (float*)calloc(msdginterp->nisky*trlen, sizeof(float));

     batchsize = 536870911.0 / pklen; // maximum number of traces that can be read in 1 go
     nRemain   = msdginterp->nsky[i];
     nDoneSize = 0;
     while (nRemain > 0) { // read in batches
       nx2read = (nRemain > batchsize)? batchsize : nRemain;
       n2read = nx2read*(size_t)pklen*sizeof(float);
       nread = pread(fOrigin, pPack + nDoneSize, n2read, (offset+nDoneSize)*sizeof(float));
       if (nread != n2read)
         G2_logCritical(msdginterp->logName, "Fail to read input data. Should read %zu, but read only %zu. Num. pkey = %d, current num. skey = %d", 
                                             n2read,nread, msdginterp->nipky, msdginterp->nsky[i]);
       nDoneSize += nx2read*(size_t)pklen;
       nRemain   -= nx2read;
     } 

     wbmin = 1.0e+8;
     for(j=0; j<msdginterp->nsky[i]; j++)
       {
	 src = (head_t*)((float*)pPack + j*pklen);
	 
	 des = (head_t*)((float*)pOrigin + (j)*trlen);
	 dbunpack(&msdginterp->pack, src, des, NHEAD);
	 
	 //Calculate minimum WB value for this gather
	 if(l1para->agc2dmethod ==2)
	   {
	     wbval[j] = src->wbid;
	     wbmin = MIN(wbval[j],wbmin);
	   }
       }
     if(l1para->agc2dmethod ==2)
       {
	 //Shift input data upward to wbmin
	 for(isky=0; isky<msdginterp->nsky[i]; isky++)
	   {
	     ishift[isky] = wbval[isky]-wbmin;
	     ntocopy = msdginterp->nsamp-ishift[isky];
	     memcpy(temp,  pOrigin+isky*trlen+NHEAD+ishift[isky], ntocopy*sizeof(float));

	     //For last samples copy over data in mirror fashion
	     j = 1;
	     for (isamp=ntocopy;isamp<msdginterp->nsamp;isamp++)
	       {
		 if( ntocopy-j >= 0)
		   temp[isamp] = temp[ ntocopy-j ];
		 else
		   temp[isamp] = temp[2*ntocopy];
		 j++;
	       }
	     memcpy(pOrigin+isky*trlen+NHEAD,temp,msdginterp->nsamp*sizeof(float));
	   }
       }
     else wbmin = 0;

     for (isgatex=0;isgatex<nsgatex;isgatex++)
       {
	 
	 isbegx=MAX(0,isgatex*l1para->winxagc2d-outovlapx);
	 isendx=MIN(obnixmax-1,(isgatex+1)*l1para->winxagc2d+outovlapx-1);
	 
	 if (isgatex==nsgatex-1) isendx=obnixmax-1;
	 if (isgatex==0) isbegx=0;
	 
	 jsky = -1;
	 
	 for (isky=0;isky<msdginterp->nsky[i];isky++)
	   {
	     
	     src=(head_t*)(pOrigin+isky*trlen);
	     obnkeyx=src->obnxpos;
	     
	     if (obnkeyx >=isbegx && obnkeyx<=isendx) jsky = jsky + 1;
    
	   } 

	 ntracewin = jsky+1;

	 winoutscl = NULL;
	 winoutscl = myallocatef(winoutscl,ntracewin,msdginterp->nsamp);
     
	 input = NULL;	
	 input = myallocatef(input,ntracewin,msdginterp->nsamp);

	 mapkey  = (int *) malloc(ntracewin*sizeof(int));
	 mapkeyx = (int *) malloc(ntracewin*sizeof(int));

	 jsky = -1;
	 ntocopy = msdginterp->nsamp-wbmin;
	 for (isky=0;isky<msdginterp->nsky[i];isky++)
	   {
	     
	     memcpy (trace,  pOrigin+isky*trlen, trlen*sizeof(float));

	     src=(head_t*)(trace);
	     obnkeyx=src->obnxpos;
	     
	     if (obnkeyx >=isbegx && obnkeyx<=isendx)
	       {
		 jsky = jsky + 1;
		 
		 mapkey[jsky] =  isky;
		 mapkeyx[jsky] = obnkeyx;
		 
		 for (l=0;l<ntocopy;l++)
		   { 
		     input[jsky][l]=trace[NHEAD+l+wbmin];
		   }
	       } 
	   } 

	 l1para->nsamp = ntocopy;
	 agc2d(l1para,ntracewin,&input[0],&winoutscl[0]);
	 l1para->nsamp =msdginterp->nsamp;
	 	 
	 if (isgatex==0)
	   {
	     for (jsky=0;jsky<ntracewin;jsky++)
	       {

		 isky = mapkey[jsky];
		 for (isamp=0;isamp<msdginterp->nsamp;isamp++)
		   {
		     if(isamp<wbmin)
		       outputscl[isky*trlen+NHEAD+isamp]=winoutscl[jsky][0];
		     else
		       outputscl[isky*trlen+NHEAD+isamp]=winoutscl[jsky][isamp-wbmin];
		   }
	       }
	   }
	 else
	   {
	     
	     startidx=isbegx;
	     endidx  =isbegx+2*outovlapx-2;
	     rstartidx = (float)startidx*l1para->obngridx; 
	     rendidx =   (float)endidx  *l1para->obngridx; 
	     
	     for (jsky=0;jsky<ntracewin;jsky++)
	       {
		 isky = mapkey[jsky];
		 obnkeyx = mapkeyx[jsky];
		 if (obnkeyx >=rstartidx &&obnkeyx<=rendidx)
		   {
		     ////rtap=lintaperx[jsky-startidx];
		     rtap=(float)(obnkeyx-rstartidx+l1para->obngridx)/(float)((2*outovlapx)*l1para->obngridx);
		 
		   }
		 else rtap = 1.0f;
		 for (isamp=0;isamp<msdginterp->nsamp;isamp++)
		   {
		      if(isamp<wbmin)
		       outputscl[isky*trlen+NHEAD+isamp]=winoutscl[jsky][0]*rtap+
			 outputscl[isky*trlen+NHEAD+isamp]*(1.0f-rtap);
		     else
		       outputscl[isky*trlen+NHEAD+isamp]=winoutscl[jsky][isamp-wbmin]*rtap+
			 outputscl[isky*trlen+NHEAD+isamp]*(1.0f-rtap);
		   }
		   
	       }
	   }
	 
	 free(mapkey);	 
	 free(mapkeyx);	 

	 for (j=0;j<ntracewin;j++)
	   {
	     free(input[j]);
	     free(winoutscl[j]);
	   }
	 
	 free(input);
	 free(winoutscl);

       }

     if(l1para->agc2dmethod ==2)
       {
	 //Shift the scalars down to original position. 
	 for(isky=0; isky<msdginterp->nsky[i]; isky++)
	   {
	     memset (temp,0 , msdginterp->nsamp * sizeof (float));
	     ntocopy = msdginterp->nsamp-ishift[isky];
	     memcpy(temp+ishift[isky],  outputscl+isky*trlen+NHEAD, ntocopy*sizeof(float));
	     
	     //For all samples above WB, use scalar at WB
	     rval =temp[wbval[isky]];
	     for (j=0;j<wbval[isky];j++)
	       temp[j] = rval;
	     
	     memcpy(outputscl+isky*trlen+NHEAD,temp,msdginterp->nsamp*sizeof(float));
	   }
       }

     batchsize = 536870911.0 / trlen; // maximum number of traces that can be read in 1 go
     nRemain   = msdginterp->nsky[i];
     nDoneSize = 0;
     while (nRemain > 0) { // write in batches
       nx2write = (nRemain > batchsize)? batchsize : nRemain;
       n2write = nx2write * (size_t) trlen*sizeof(float);
       nwrite = pwrite(fScalar,outputscl+nDoneSize,n2write,(offsetou+nDoneSize)*sizeof(float));
       if (nwrite != n2write)
         G2_logCritical(msdginterp->logName, "Fail to write to output. Should write %zu, but wrote only %zu. Num. pkey = %d, current num. skey = %d", 
                                             n2write, nwrite, msdginterp->nipky, msdginterp->nsky[i]);
       nDoneSize += nx2write*(size_t) trlen;
       nRemain   -= nx2write;
     }

     offset =   (size_t)(msdginterp->nsky[i])*(size_t)pklen + offset;
     offsetou = (size_t)(msdginterp->nsky[i])*(size_t)trlen + offsetou;
     
     free(outputscl);
     free(pPack);
     free(pOrigin);

   }//// end of : for(i=0; i<msdginterp->nipky; i++)

  free(trace);
  if(l1para->agc2dmethod ==2)
    {
       free(wbval);
       free(temp);
       free(ishift);
    }
}

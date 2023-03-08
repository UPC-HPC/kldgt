#include <stddef.h>
#include <stdbool.h>
#include "/usr/include/sched.h"
#include "DB.h"
#include "size.h"

#include "demult.h"
#include "G2WorkProgress.ci"
#include "dev_wrapper.h"

#define CHOOSEMETHOD_MSDGI      1 
#define CHOOSEMETHOD_FP4D       2 
#define CHOOSEMETHOD_BSDG3D     3 
#define CHOOSEMETHOD_JOINTSR3D  4 
#define CHOOSEMETHOD_BSDG2D     5
#define CHOOSEMETHOD_JOINTSR2D  6 
#define CHOOSEMETHOD_DAS        7

#define SMTHTYPE_NO         0
#define SMTHTYPE_MEDIAN     1
#define SMTHTYPE_GAUSSIAN   2

#define MASKQCNO            0
#define MASKQCRATIO         1
#define MASKQCMSK           2

#define WNOISE              0
#define DBCTRL              1
#define NOTCHCTRL           2

#define FPMETHOD_INTERPOLATE    1 
#define FPMETHOD_FILLIN         2
#define FPMETHOD_REGULARIZE     3
#define FPMETHOD_ANTIALIAS      4

#define IOTYPE_SDS        1
#define IOTYPE_TANGO      2
#define IOTYPE_CAPI       3
#define IOTYPE_CONST      4 //for velocities
#define IOTYPE_FDM        5

#define CABLE               0
#define SURFACE             1

#define DEGHOST             0

#define FIXED               0
#define VARIABLE            1

#define NONO                0
#define YESYES              1
#define YESNO               2

#define SHOT                1
#define RECEIVER            2
#define SHOTRCVR            3
#define SHOTRCVRJNT2D       4

#define CABLE90             0
#define SURFACE0            1
#define SURFACE90           2
#define DGRGCABLE90         3

#define DSDGOUT_PZ          0
#define DSDGOUT_P           1
#define DSDGOUT_Z           2
#define DSDGOUT_Pnotch      3
#define DSDGOUT_Pupnotch    4
#define DSDGOUT_Znotch      5
#define DSDGOUT_Zupnotch    6
#define DSDGOUT_Ppeak       7 
#define DSDGOUT_Puppeak     8
#define DSDGOUT_Zpeak       9
#define DSDGOUT_Zuppeak    10
#define DSDGOUT_ZtoV       11

#define FPREDATUM_NO      0
#define FPREDATUM_CTOS    1
#define FPREDATUM_MTOS    2
#define FPREDATUM_DGRG    3
#define FPREDATUM_SIN     4
#define FPREDATUM_CIN     5
#define FPREDATUM_STOC    6
#define FPREDATUM_STOM    7
#define FPREDATUM_XG_Z    8

#define PCROSSGHOST       111
#define VCROSSGHOST       222
#define PVSUM             333
#define VZTOV             444
#define VZTOVCR           555

#define DGMETHOD_3DSTR     1
#define DGMETHOD_OBN       2
#define DGMETHOD_FP        3
#define DGMETHOD_BOOTSTRAP 4
#define DGMETHOD_RESDEG    5
#define DGMETHOD_VSP       6

#define LMSDGINTERP_SNAME          128     //short name
#define LMSDGINTERP_LNAME          512     //long name
#define LBSDG_SNAME          128     //short name
#define LBSDG_LNAME          512     //long name

#define RDCTYPE_UPDATE      0

#define VDIM_T3VL2D         2
#define VDIM_T3VL3S         30
#define VDIM_T3VL3X         31

#define OFFMODE_OFFSET     0
#define OFFMODE_OFFSETX    1
#define OFFMODE_OFFSETY    2
#define OFFMODE_OFFSETC   3
#define MAXSECDATA 10

#define FKTYPE_QC        0
#define FKTYPE_MULT      1
#define FKTYPE_CABLE     2
#define FKTYPE_SURFACE   3
#define FKTYPE_CTOS      4
#define FKTYPE_MTOS      5
#define FKTYPE_STOC      6 
#define FKTYPE_STOM      7
#define FK_REGTYPE_RCV   0
#define FK_REGTYPE_SHOT  1
#define FK_REGTYPE_BOTH  2

#define MSIGN_POSITIVE      1
#define MSIGN_NEGATIVE     -1

#define JOINTDEGHOST               10
#define JOINTDESIGNATURE           11
#define JOINTDESIG_SRC_DEGHOST     12
#define JOINTDESIG_RCV_DEGHOST     13
#define JOINTDESIG_SRCRCV_DEGHOST  14
#define JOINT_SRC_DEGHOST          15

#define TARGETTYPE_WAVELET     0
#define TARGETTYPE_SPIKE       1
#define TARGETTYPE_GAPDECON    2

#define DESIGNORM_NONE     0
#define DESIGNORM_MAX      1
#define DESIGNORM_RMS      2

#define MAX_FFT_PLAN	   1

typedef struct ATTR_T attr_t;
typedef struct FATTR_T fattr_t;
typedef struct DBIO_T dbio_t;
typedef struct HEAD_T head_t;

typedef struct ST_MASTER master_t;

typedef struct ST_L1INV l1inv_t;

typedef struct TAUP_STRUCT taup_struct;

typedef struct ST_L1INV_FPINTERP fpinterp_t;

typedef struct ST_PACKSS packss_t;


struct ST_PACKSS { int nbyte, nscal, nplot, nb4pk, pkwin; };
//typedef struct FREQBD_T { float lowcut, lowpass, highpass, highcut, nyquist; } freqbd_t;

#define LEN_FILE_NAME   128
#define MAX_HOSTS       100
#define TANGODS         1
#define WZDATA          2

#define LFPINTERP_SNAME                 128     //short name
#define LFPINTERP_LNAME                 512     //long name

///////////////////////////////////////////////////////////////////////

enum {ANGLE_FILTER_OUTPUT_TOTAL, ANGLE_FILTER_OUTPUT_RESIDUAL, ANGLE_FILTER_OUTPUT_SIGNAL};

struct GPU_FFT_T
{
    int lfft, nbatch;
    int bufFlag; //if need to allocate buffer
    void *buf;
    void *fftplan;
    dev_fft_type type;
    int usage;
};

struct DBIO_T
{
  int flag;

  dbinfo_t odbinfo, mdbinfo, ddbinfo, tdbinfo, adbinfo, sdbinfo, srcdbinfo, targetdbinfo;
  dbinfo_t desigdbinfo1, desigdbinfo2;
  dbinfo_t ddbinfo2[MAXSECDATA],  odbinfo2[MAXSECDATA];
  
  char  oProj [LMSDGINTERP_LNAME], oDset [LMSDGINTERP_LNAME];
  char  mProj [LMSDGINTERP_LNAME], mDset [LMSDGINTERP_LNAME];
  char  srcProj [LMSDGINTERP_LNAME], srcDset [LMSDGINTERP_LNAME];
  char  targetProj [LMSDGINTERP_LNAME], targetDset [LMSDGINTERP_LNAME];
  char  dProj [LMSDGINTERP_LNAME], dDset [LMSDGINTERP_LNAME];

  char  desigProj[LMSDGINTERP_LNAME];
  char  desigDset1[LMSDGINTERP_LNAME], desigDset2[LMSDGINTERP_LNAME];
  char  prProj [LMSDGINTERP_LNAME], prDset [LMSDGINTERP_LNAME];
  char  sgProj [LMSDGINTERP_LNAME], sgDset [LMSDGINTERP_LNAME];
  char  rgProj [LMSDGINTERP_LNAME], rgDset [LMSDGINTERP_LNAME];
  char  srgProj [LMSDGINTERP_LNAME], srgDset [LMSDGINTERP_LNAME];

  char  tProj [LMSDGINTERP_LNAME], tDset [LMSDGINTERP_LNAME];
  char  aProj [LMSDGINTERP_LNAME], aDset [LMSDGINTERP_LNAME];
  char  sProj [LMSDGINTERP_LNAME], sDset [LMSDGINTERP_LNAME];
  char  vProj [LMSDGINTERP_LNAME], vDset [LMSDGINTERP_LNAME];

  char  dsdgou1Proj [LMSDGINTERP_LNAME], dsdgou1Dset [LMSDGINTERP_LNAME];
  char  dsdgou2Proj [LMSDGINTERP_LNAME], dsdgou2Dset [LMSDGINTERP_LNAME];
  char  dsdgou3Proj [LMSDGINTERP_LNAME], dsdgou3Dset [LMSDGINTERP_LNAME];
  char  dsdgou4Proj [LMSDGINTERP_LNAME], dsdgou4Dset [LMSDGINTERP_LNAME];

  char recvdfile[LMSDGINTERP_LNAME],gfilter[LMSDGINTERP_LNAME];
  char dirheader[LMSDGINTERP_LNAME], desigmethod[LMSDGINTERP_LNAME], desigoutput[LMSDGINTERP_LNAME];
  char  nhost[2];

  float vextrap;
  float nottzero, tartzero, wnoise, freqhi, limithi,maxalpha;
  int smooth;



  int next;         //for tango velocity interpolation.
  //=====================================================
  int datatype; // 1 tango dataset; 2 wz data 
  int icid, ocid;

  //jhong added for WZ data
  char iWZfile [FS_SZFILE];
  char oWZfile [FS_SZFILE];
//  int  nohost;

  int iotype, vltype, otype, mtype, srctype, targettype, dtype, ttype, atype, stype;
  int dsdgou1type, dsdgou2type, dsdgou3type, dsdgou4type;
  CapiID oinst, minst, srcinst, targetinst, dinst, tinst, ainst, sinst,prinst, sginst, rginst, srginst; // FCAPI instances
  CapiID desig1inst, desig2inst;
  CapiID zfiltinst;
  CapiID srcinst2[MAXSECDATA];
  CapiID defsrcinst;
  CapiID oinst2[MAXSECDATA],dinst2[MAXSECDATA];
  CapiID dsdgou1inst, dsdgou2inst, dsdgou3inst, dsdgou4inst;
  long int vinst, nmoinst;
  int dim;
  
  char  oIdent [LMSDGINTERP_LNAME], oVers [LMSDGINTERP_LNAME], oSite [LMSDGINTERP_LNAME];
  char  mIdent [LMSDGINTERP_LNAME], mVers [LMSDGINTERP_LNAME], mSite [LMSDGINTERP_LNAME];
  char  dIdent [LMSDGINTERP_LNAME], dVers [LMSDGINTERP_LNAME], dSite [LMSDGINTERP_LNAME];
  char  srcIdent [LMSDGINTERP_LNAME], srcVers [LMSDGINTERP_LNAME], srcSite [LMSDGINTERP_LNAME];
  char  targetIdent [LMSDGINTERP_LNAME], targetVers [LMSDGINTERP_LNAME], targetSite [LMSDGINTERP_LNAME];

  char  dsdgou1Ident [LMSDGINTERP_LNAME], dsdgou1Vers [LMSDGINTERP_LNAME], dsdgou1Site [LMSDGINTERP_LNAME];
  char  dsdgou2Ident [LMSDGINTERP_LNAME], dsdgou2Vers [LMSDGINTERP_LNAME], dsdgou2Site [LMSDGINTERP_LNAME];
  char  dsdgou3Ident [LMSDGINTERP_LNAME], dsdgou3Vers [LMSDGINTERP_LNAME], dsdgou3Site [LMSDGINTERP_LNAME];
  char  dsdgou4Ident [LMSDGINTERP_LNAME], dsdgou4Vers [LMSDGINTERP_LNAME], dsdgou4Site [LMSDGINTERP_LNAME];

  char  prIdent [LMSDGINTERP_LNAME], prVers [LMSDGINTERP_LNAME], prSite [LMSDGINTERP_LNAME];
  char  sgIdent [LMSDGINTERP_LNAME], sgVers [LMSDGINTERP_LNAME], sgSite [LMSDGINTERP_LNAME];
  char  rgIdent [LMSDGINTERP_LNAME], rgVers [LMSDGINTERP_LNAME], rgSite [LMSDGINTERP_LNAME];
  char  srgIdent [LMSDGINTERP_LNAME], srgVers [LMSDGINTERP_LNAME], srgSite [LMSDGINTERP_LNAME];

  char  tIdent [LMSDGINTERP_LNAME], tVers [LMSDGINTERP_LNAME], tSite [LMSDGINTERP_LNAME];
  char  aIdent [LMSDGINTERP_LNAME], aVers [LMSDGINTERP_LNAME], aSite [LMSDGINTERP_LNAME];
  char  sIdent [LMSDGINTERP_LNAME], sVers [LMSDGINTERP_LNAME], sSite [LMSDGINTERP_LNAME];
  char  vIdent [LMSDGINTERP_LNAME], vVers [LMSDGINTERP_LNAME], vSite [LMSDGINTERP_LNAME];
  
  char  oPath [6*LMSDGINTERP_LNAME];
  char  oPath2 [MAXSECDATA] [6*LMSDGINTERP_LNAME];
  char  mPath [6*LMSDGINTERP_LNAME];
  char  dPath [6*LMSDGINTERP_LNAME];
  char  dPath2 [MAXSECDATA][6*LMSDGINTERP_LNAME]; 
  char  srcPath [6*LMSDGINTERP_LNAME];
  char  srcPath2 [MAXSECDATA][6*LMSDGINTERP_LNAME];
  char  defsrcPath [6*LMSDGINTERP_LNAME];
  char  targetPath [6*LMSDGINTERP_LNAME];
  char  zfiltPath [6*LMSDGINTERP_LNAME];
  char  desig1Path [6*LMSDGINTERP_LNAME];
  char  desig2Path [6*LMSDGINTERP_LNAME];

  char  prPath [6*LMSDGINTERP_LNAME];
  char  sgPath [6*LMSDGINTERP_LNAME];
  char  rgPath [6*LMSDGINTERP_LNAME];
  char  srgPath [6*LMSDGINTERP_LNAME];

  // LKL: add this to handle conversion between vectors and local files during dosetup and dogroup phase
  char  prVecname [LMSDGINTERP_LNAME]; /* for fpbootstrap multi-component */
  char  sgVecname [LMSDGINTERP_LNAME];
  char  rgVecname [LMSDGINTERP_LNAME];
  char  srgVecname[LMSDGINTERP_LNAME];
 
  /* for msdgi 3 types of input vectors / datasets */
  /* pressure data always have trace as the input vector */
  char  mVecname [LMSDGINTERP_LNAME]; /* dataz    */
  char  tVecname [LMSDGINTERP_LNAME]; /* datay    */

  char  tPath [6*LMSDGINTERP_LNAME];
  char  aPath [6*LMSDGINTERP_LNAME];
  char  sPath [6*LMSDGINTERP_LNAME];
  char  vPath [6*LMSDGINTERP_LNAME];

  char  dsdgou1Path [6*LMSDGINTERP_LNAME];
  char  dsdgou2Path [6*LMSDGINTERP_LNAME];
  char  dsdgou3Path [6*LMSDGINTERP_LNAME];
  char  dsdgou4Path [6*LMSDGINTERP_LNAME];

};

struct FATTR_T
{
  int         mark;
  int         id;
  int         type;
  char        name[LMSDGINTERP_LNAME];
  int         *v;           //attribute value
  int         *key, nkey;
  float       min, max,*rv;
};

struct ATTR_T
{
  int         mark;
  int         id;
  int         type;
  char        name[LMSDGINTERP_LNAME];
  int         *v;           //attribute value
  int         *key, nkey;
  int         min, max;
};


struct ST_L1INV_FPINTERP
{
    int    nsampfpintp,itbeg_fpintp,itend_fpintp;
    int    flagtgate, zerotime, method;
    int    srate, nsamp, nsamporig, nxwin;
    int    nxoverlap, nsky, nxpad;
    int    ntrc, nthread, niter;
    int    fftnr, fftnc, ntwin, ntlap, nloop;
    int    orgnr, orgnc, fftNR, fftNC;
    int    npxpad, nxwinpad, ninterp, naltp, keeporig;
    int    maxhole, nbeg, nend, sgatelen, sgatelen_reg;
    int    npall, np_reduce, npx, npy, memMode;
    float  twopif;
    float  lownear, lowfar, highfreq, pxmin, pxmax, highalias;
    float  chanitv, psample, vwater;
    float  swnear, swfar, tolerance;
    float  wbshift, ampscale, freqlen;
    int    nskynew, rmslen, spcntwin, spcntlap, imode;
    int    nywin, nyoverlap, npypad, nywinpad, nswinpad, nppad, mode3d;
    int    highcut_flt;
    float  pyminout,pymaxout, pymin, pymax, xchanitv, ychanitv, percent, percent2, pyminorg, pymaxorg;
    int    nxwinmax, nywinmax, maxout, nxnymax, ncab_intp;
    float  cablespacing, chanspacing, pxhi, pyhi, highcut;
    int    laaf, ltvp, langlemute, tvgatenum[10], flt_order, lresidual, lwzbin;
    float  tvmaxp[10], tvmaxf[10], tvratio[10], tvmigdx, tvmigdy, migdx, migdy, abratio; 
    float  anglemin, anglemax, obliqpower, threshold, deltapx, deltapy;
    fftwf_plan fftplan[2];
    float *xGhost_r, *xPrimary_r,*xGhost_i, *xPrimary_i;
    float **pfp_r, **pfpT_r, **pfx_r, **pfxT_r;
    float **pfp_i, **pfpT_i, **pfx_i, **pfxT_i;
    float **A11_r, **A12_r, **A21_r, **A22_r, **T11_r, **T12_r, **T21_r,**T22_r;
    float **A11_i, **A12_i, **A21_i, **A22_i, **T11_i, **T12_i, **T21_i,**T22_i;
    float *xPrimary, **aaa, **bbb, **ccc, **ddd;
    float **pfp, **pfpT, **pfx, **pfxT;
    float **A11, **A12, **T11, **T12;
    float **taup, **taup1, **taupall, **sparse, **weight;
    float *d_offset, *d_offsetout, *d_offx, *d_offy, *d_offxout,
        *d_offyout, *d_ztop;
    float *d_pfp, *d_pfpT, *d_pfx, *d_pfxT;
    float *d_xPrimary, *d_bbb, *d_ccc, *d_ddd;
    float *d_Ain, *d_Tin, *d_Asort, *d_Tsort, *d_Areg, *d_Treg;
    float *d_taup, *d_taupall, *d_sparse, *d_weight;
    float *d_hilbert_hw;
    float *d_lpfilter;
    int *d_wbid;
    size_t pfxTmem;
    int fftnr_orig, fftnc_orig;
    int fftnr_next, fftnc_next;
    int downratio, downratio_next;
    int sort_flag, filter_flag;
    int hilbert_radius;
    int fft_plan_num;
};

struct ST_L1INV
{
  int    srate,nsamp,nsgtrace,nsgtaper, computemode, srdg; 
  int    nsgoverlap,nsky,jtdgitpflag,fp4didata; 
  int    ntrc,ntrcall,datum,nthread,niter,niterlow,niter2,outpos,lowguide;
  int    fftnr,fftnc,fpredatum,offmode,dsdgout,outtype;
  int    npxpad,nsgtracepad;
  float  lownear,lowfar,highfreq,highfreq_p,lowstart,pxmin, pxmax,tminscl,residual,resscale; 
  float  chanitv,psample,vwater,reflectivity;  
  float  swnear,swfar,tolerance,regparam,pvamp,pvphase;
  int    srcnum,rcvnum,znum,ynum,nfreq2;
  float  srcfreq[20],srcrefl[20],rcvfreq[20],rcvrefl[20],zfreq[20],zrefl[20],yfreq[20],yrefl[20];
  float  *srcrfl,*d_srcrfl,prfl,*prfl2,*d_prfl,*zrfl,*d_zrfl,*yrfl,*d_yrfl,*d_reflect; 
  float  plowcut,zlowcut,zlowmerge,zphase,zampscl,zcosmin,zlowtap,ztminscl;
  float  ylowcut,yphase,yampscl,ycosmin,ylowtap,ppower,zpower;
  int    zmethod,zobliq,zinteg,pmethod;
  int    ymethod,yobliq,yinteg;
  int    ms_pyes,ms_zyes,ms_yyes;
  int    clustertype; //=0 for CPU, =1 for GPU
  float *ztop,*ytop,*wtilt,*zwtilt,*ywtilt,*zlowf,*ylowf;

  float *freqfac,*obliqcorr,*integration;

  //add by kyang testing speed up
  //fftwf_plan fftplan[2];
  int fftplan[2];
  int zerotraceattr;
  int lzfilter, qczfilt, myrank;
  float *zfilter;
  int orgnr, orgnc, fftNR, fftNC;         //for vectorization
  float **AP_p, **AG_p, **AS_p, **TP_p, **TG_p, **TS_p;
  float ***AP_parray, ***AG_parray, ***AS_parray, ***TP_parray, ***TG_parray, ***TS_parray, ***pSrcFilter_array;
  float **TPR_p, **TGR_p;
  float **AGS_p, **AGR_p, **TGS_p, *xGhostS, *xGhostR;
  float ***AGS_parray, ***AGR_parray, ***TGS_parray, ***TGR_parray;
  float **TP_interp3d, **TS_interp3d;
  float **pfx, **pfp, **pfpin, **pfpinit, **pfpadd, **pfpT, **pfpT2, **pfxT, **pSrcFilter, **pSrcFilterAz, **pTgtFilter;
  float **p1dSrc, **p1dTgt;
  float *xGhost, *xPrimary;
  int hilbert_nt; void *hilbert_task; float *hilbert_hw, *hilbert_conv;
  float **cg_d, **cg_wd, **cg_r, **cg_Ad;
  float **cg_s,**cg_As,**cg_rhat0; // new BCG method
  //end

  float *pxsave, *pysave;

  int npx, npy, npxsave, npysave,npxorig, npyorig;        //should calculate once only

  int nsamporig;

  float dpx;

  ////////////////////////////////////////////
  ////////////////////////////////////////////
  int obnixmax,obniymax;
  float obnxposmin,obnxposmax;
  float obnyposmin,obnyposmax;
  float obngridx,obngridy;
  int nsgtracex,nsgtracey;
  int nsgoverlapx,nsgoverlapy;

  float pymin,pymax,pyminout,pymaxout;
  int npypad;
  float chanitvx,chanitvy;
  int maipky,mzipky;
  int nipkyou;

  int nsgtracewin,nsgtracewinout;
  float npperc;
  int nsgtraceyloc;
  int nsgtracexlastgate;
  ////////////////////////////////////////////
  ////////////////////////////////////////////

  //////////////////////////////////////////
  //////////////////////////////////////////
  int ndata;

  float *recz_p, *offsetx_p, *offsety_p;
  float *recz_az,*offsetx_az,*offsety_az;
  float *recz_ay,*offsetx_ay,*offsety_ay;

  //////////////////////////////////////////////////////////////////////////
  ////////////////////////// 3D INTERPOLATION //////////////////////////////
  //////////////////////////////////////////////////////////////////////////

  int interp3d,deghost3d;
  float *offsetx_p_intp, *offsety_p_intp, *recz_p_intp;
  int ncab_intp,noffsetou; 

  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  int jointmc;

  int msdatap,msdataz,msdatay;

  int filtorder;

  float highfreqtap;

  int msdgi_option;

  ///GPU releted
  int memMode;
  float twopif;
  float sratescl;
  float deltapx, deltapy;
  size_t pfxTmem;
  float *d_ztop, *d_ytop;
  float *d_zlowf, *d_ylowf;
  float *d_zwtilt, *d_ywtilt;
  float *d_pfp, *d_pfpin, *d_pfpadd, *d_pfpT,*d_pfx,*d_pfxT,*d_pfpinit;
  float *d_recz,*d_reczou;
  float *d_offset,*d_offsetou;
  float *d_obliqcorr;
  float *d_freqfac;
  float *d_hilbert_hw;
  float *d_cg_d, *d_cg_wd, *d_cg_r, *d_cg_Ad;
  float *d_cg_s, *d_cg_As, *d_cg_rhat0; // new BCG method

  int *d_wbid;  
  float *d_input_p, *d_input_az, *d_input_ay, *d_winout, *d_reswinout;
  
  //desig option
  float *d_shotsrc, *d_gunvol, *d_pfxsource, *d_pfxtarget;


  int taupQC;

  int dgredatumflag; 
  int finalflag;

  int ntgate,zntgate,tgmin,tgmax,tgoverlap,tglength,flagl1tgate;
  int tgsamp,tgsampmin,tgsampmax;

  int choosemethod, reusemem;

  int *nsgtracearr4d;
  int addrestores;
  float sweight2fac,fraccomm;

  float pxlo[10], pxhi[10], flocut[10], fhicut[10];
  int ldipf;

  int agc2dflag;
  float *scalaragc2d;
  float winagc2d, taptagc2d, sclmaxagc2d;
  int winxagc2d,ovlapxagc2d,agc2dmethod,agcqc;

  /////////////////////////////////////////////////////////////////
  ////////////////////////// BSDG 2D //////////////////////////////
  int nxpad,nsampl1orig,nsampl1,nsgatetot;
  int switchreczshotz;
  int *tgate;
  int dataghost,acquis3d;
  int flagapplc;
  int zerotime,win0,dwbid,l1ismin,flagfastl1inv;
  int imaxholesize,sparsenormflag,flagghostout;
  int xgsampmax,xgsampmin;
  int xgoverlapmin,xgoverlapmax;
  int    *zoverlap;
  int pchan,pvobliq,pvmethod,pvinteg;
  int *pkyno;
  int tidalheightsflag,tidalheightsiter,nlag,nlead,nlags,ipickwinend,ipickwinstt,seastresampfac,primaryredatout,thtsmoothlen;
  int nsampcorr,sratecorr,maxnsgtrace;


  float deltat_3dl1,wnnear_3dl1,wnfar_3dl1,thetacut3dl1,dcut3dl1;
  float lowcutf,swscale;
  float maxholesize;
  float highstart;
  int *isbeg,*isend,*isgatesubgath;

  float tcorrlag,tcorrlead,thtmax,pickwinstt,pickwinend,corrtimewin,offstop;
  int npxfw;
  float *taper2;
  fftwf_plan *fftplanorig;
  fftwf_plan *fftplanrs;

  int fft_plan_num;

  int cgmethod;
  /////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////
  float tzerosrc, tzerotgt, tzerozfilt, tzerotgttaper,tzerosrc2[MAXSECDATA]; 
  float *psig;   
  int jointmethod, ntrcsrc, nsampsrc, sratesrc, ntrcsrc2[MAXSECDATA];
  float gunmin,gunmax,gundelay;
  int bsouttype;
  float totgun,shgun,dpgun;


  int ipn,ipn_type;
  float ftap;
  int flag_recz_smooth, target_type, lmerge1d; 
  float mrg1dfreq, mrg1dmaxg;
  int iprocessing_nthres, designorm, choosemethod_msdgi_pseudo;
  int iexternlc;
 
  float rqmin,rqmax,rqinc,rqminsave,rqmaxsave,rqincsave;
  float  pxmin_linear, pxmax_linear, pymin_linear, pymax_linear;	//ezhang: added to differentiate from parabolic skew p terms
  float pxmin_parab, pxmax_parab, pymin_parab, pymax_parab;     //ezhang: added in px and py specifically for parabolic transform
  int nqorig,nqsave;
  float idoparabtrans;						//ezhang: changed from integer to float to allow for taper for parabolic transform
  float *rqsave;
  float offxmax_parab;
  float taper_parab;
 
  // related to RMS qc
  bool rmsqc;
  char rmsqcfile[LMSDGINTERP_LNAME];
  
  int targetconvopt;

  /* ************** DAS-related ************** */
  bool applyobliq;
  float nobliq2d;
  bool isVariableRecVel;
  float recVelAverage;
  float cutOffAngle;
  float CosineCutOffAngle;
  bool applyTaper;
  float taperBegAngle;
  float taperEndAngle;
  float taperBegValue;
  float taperEndValue;
  float compensationTaper[101];  // Fixed tapering function interval

  /* Angle filtering */
  bool applyAngleFilter;
  int nAngleFilters;
  float lowCutTaperBegAngle[99];
  float lowCutTaperBegValue[99];
  float lowCutTaperEndAngle[99];
  float lowCutTaperEndValue[99];
  float highCutTaperBegAngle[99];
  float highCutTaperBegValue[99];
  float highCutTaperEndAngle[99];
  float highCutTaperEndValue[99];
  float lowCutTapers[99][101];
  float highCutTapers[99][101];
  int angleFilterOutputOption;
  float lowCutTaperBegAngleIterator, lowCutTaperBegValueIterator;
  float lowCutTaperEndAngleIterator, lowCutTaperEndValueIterator;
  float highCutTaperBegAngleIterator, highCutTaperBegValueIterator;
  float highCutTaperEndAngleIterator, highCutTaperEndValueIterator;
  float *lowCutTapersIterator, *highCutTapersIterator;

  /* Positive/negative angle separation */
  bool applyPosNegSeparation;
  float posNegSeparationMode;
  float posNegSeparationTaperLimit;
  float posNegSeparationTapers[2][101];
  float *posNegSeparationTaperIterator;
  /* ************** DAS-related ************** */
  int keeporig;  
  float total_gunvol,designorm_factor;

};

struct TAUP_STRUCT
{
  int         shiftmethod;
  int         tpresampfac;
  int         nsampnew;
};

struct ST_MASTER
{
  bool  g2progress;                  
  int datum;
  char    logName[LEN_LONG]; // logging name for module

  bool rmsqc; //whether to print out convergence curve of rms error in CG iterations
  char rmsqcdir[LMSDGINTERP_LNAME]; // where to write the convergence curve 

  int dgtype,dgmethod,wbduse, regattruse; //deghost type,stacktype, deghost method
  int ninput;
  float reflectivity,tolerance,regparam,resscale,swnear,swfar;
  int niter,niterlow,niter2;
  float lowfreq,highfreq,highfreq_p,deltat,vconst;
  float shotdepth;
  float vwater,wbshift,ampscale,ampscale_p,ampscale_az,ampscale_ay;
  float chanitv,psample,vscale;
  int fftnr,fftnc,fftln,wdepthmax,wdepthmin;
  struct ST_PACKSS      pack,mpack,tpack;
  struct ST_PACKSS      desigpack1,desigpack2;

  int nbas;
  int vlen,mvlen,tvlen;
  int desigvlen1, desigvlen2;
  int mem_MB, dat_MB;               //number of core and memory

  attr_t shootdir;
  attr_t opky, osky,ipky,isky, iskyintp, intpflag, cableidintp, cabtrintp, offs,oreg;
  attr_t vpky, vsky,sl,xl;
  fattr_t shotx,shoty,shotz,pz,recx,recy,recz,recz_smooth,offsetx,offsety,offsetz,ampth,recz_bsdg, gunamp, tdelay;
  fattr_t srcshotx, srcshoty, srcshotz, regattr;
  attr_t wbd;
  fattr_t recvel;
  
  int aipky,zipky,iipky,nipky,maipky,mzipky;
  int aisky,zisky,iisky,nisky,maisky,mzisky;
  int avnipky,exipky, zerotraceattr;
  int *nsky,*pkyno,*ftrc,*pkyindex;

  dbio_t  *dbio; 


  double          rms;
  char            pjnm[SX_PROJ];

  int             type;           //project type  (db|wz)
  int             vtype;          //velocity type (vl|fdm)

  char jsname[JSFILE_LEN];
  char pjname[SIZE_PROJ_NAME];

  int nsamp;
  int srate;
  int ninst; 

  int designsamp1, designsamp2;
  int sratedesig1, sratedesig2;
  int ninstdesig1, ninstdesig2;
 
  G2WorkProgress_ID wpid; 

  l1inv_t l1para; 

  taup_struct *taupparams; 

  fpinterp_t fpinterppara; 

  long int mastermem,slavemem;

  float chint,dpx,pxmin,pxmax;
  int npx;
  
  
  int nsamporig,sratenew,srateorig;

  float df;

  /////////////////////////////////////
  /////////////////////////////////////
  attr_t obnpkey,obnskey;
  fattr_t obnxpos,obnypos;
  fattr_t  zerotrace;
  attr_t cabtr;
  /////////////////////////////////////
  /////////////////////////////////////


  //////////////////////////////////////////////////
  //////////////////////////////////////////////////
  
  int numgatey,maipky3d,mzipky3d;
  int *cableno,nipky3d,nisky3d,*ftrc3d,*nsky3d,**nsky3d_tmp;

  int *obnixmax3d,*obniymax3d;
  float *obnxposmin3d,*obnxposmax3d;
  float *obngridx3d;

  attr_t cableid,newcableid,cableidseq,pkeyindx;

  attr_t ztwbd;

  //////////////////////////////////////////////////
  //////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  ////////////////////////// 3D INTERPOLATION //////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  int ninstintp;
  int *ncabpershot,*ncabpershotorig;
  int *pkyno3d;
  int ncab_intp; //// == number of interpolated cables in between 2 cables 

  //// tidinp, intpflag, iskyintp, cableidintp, cabtrintp added in HEAD_T and ST_MASTER
  //// robnkeyy added in HEAD_T
  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  ////////////////////////////////// FPINTERP //////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  int nattrout,lshotongrid,lshotbyshot,luseprereg,lwbduse,ldefaultsrc;
  float xorig,yorig;
  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////

  int msdatap,msdataz,msdatay;
  
  dmt_t *dmt;

  int flt_order;

  int ntask, tasksizemax, tasksizemax_interp3d;
  int *taskin, *taskou, *ntaskpershot, *ntrcinpertask, *ntrcoupertask, *shotidpertask, *fileidpershot;
  float *ampscalar, *taskweight;
  char **filenamepershot;

  int taupQC;

  int *isbegxwinpertask,*isendxwinpertask;
  int *isbegywinpertask,*isendywinpertask;
  int *isendypertask;
  int *numtracexpertask;

  int debug;

  int *nskyou,*ftrcou;

  int partialnmoflag;
  int nmooverlapt,nmowinlen,nmowinlenmin;
  
  float geomdx,geomdy,geomx0,geomy0;

  int *vtgate0,*vtgate1;
  int maxvtgate1;
  int taskid,dwbid,numtracex,nsky3dloc;

  float divalpha;
  int invdivflag;

  int invagcflag;
  int agcgatelen;

  int dataorderdg3d; 

  /////////////////////////////////////////////////////////////////
  ////////////////////////// BSDG 2D //////////////////////////////
  int deringf,mthread,onkeep,wtype,ndip,constflag,dbzero,dbfull,tosduse;
  int skeylen,pkeylen,skeysg,pkeysg,offflag,arrszcorrflag,msign,misstrcflag,horzattrsflag,horzfmt;
  int smthratio,sigp,sigs,smthtype,isSmthLocal, shortlen,longlen;
  int dataorderfp,fpghostou,flagapplc,bosduse,dpduse;
  int nsgtap,nxpad,nsgtrace,nsgoverlap;
  int addresidue,taupmethod,resampflag,maskqc,splitmergeflag,nmor1dflag,predcntype,recuse;
  int fxouttapt,fxghostou,nwbtap,nwbtapx,rwbring,l1mergefreq,l1all,zerotime,flagmcredatum,flagl1inv,flagl1invmc,wmarker;
  int sparsenormflag,flagghostout,deadtrflagresdeg;
  int ntgate,zntgate,tgoverlap,tgsamp,tgsampmin,tgsampmax,tgmin,tgmax,tglength,flagl1tgate,flagl2inv,win0,wbfx;
  int nfxwin,tapt,tapx,l1l2tap,l1tgmin,l1ovlap,srdg,l2onepass,flagobnl1inv;
  int resampfac,pchan,nsampnew;
  int npgrmin,npgrmax,znpxwin,npxwin,sgoverlap,naltp,imaxholesize,nadapt;
  int evnpickingflag,nipky3dlocout;
  int nnotches,winstart,rphase;
  int horiwidth,ltaper,nwindow,nnotchattrflag,tgateattrflag,itapzonemax;
  int isgateyi,isgateyf,npxfw;
  int *isbeg,*isend,*isgatesubgath;
  int nsgatetot,maxnsgtrace;
  int nhorzdef[NHORIZMAX];
  float lowdip,highdip,diptaper;
  int ipx100p,ipx100n;
  int ipxp0,ipxn0;
  attr_t tosd,bosd,dpd;

  float sgdist,lowlimit,wbzone,drngbw,extrat,extratmin,extratmax;
  float reczconst,nfreq,pkeylap,skeylap,nearalign,faralign;
  float intperr,lfaddbck,lowdbdp,lowdbsh,highdbdp,highdbsh,highrefdp,highrefsh,l1noise;
  int halfwidth,taperlen,npts,mxiter,fktype;
  int shaping, ntpad,switchreczshotz;

  float offmax,offmin,l1time,theta_wbc,maxboost,lowfreqboost,offyl2,hpcutf,lowcutf;
  float lfshallow,lfdeep;
  float lowscalesh,lowscaledp,notchscalesh,notchscaledp,dnotchscale;
  float minnoise,maxnoise,lowp,highp,recerror,shterror,wnnear,wnfar,wndeep;
  float minnoise1,maxnoise1;
  float filtertype;
  float maxholesize,freql1l2;
  float ds,narray,strratio;

  int *l1isminarr,*l1ismaxarr,*listofcablocal;
  int l1ismin,l1ismax,nsampl1;
  attr_t corrs,corre;
  fattr_t horz1,horz2,horz3;        // for horizon attributes
  attr_t nnotchattr,minlenattr,maxlenattr,datadeadtr;
  float *psignature,*psignature_sbys; 
  /////////////////////////////////////////////////////////////////
  ////////////////////////// DESIGNATURE //////////////////////////////

  int desigflag;
  fattr_t desigshotx,desigshoty,desigshotz;
  int shootdirpar;
  int nshiftdata;
  /////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////

  int nounderboost;
  int flowprepflag; 
};

struct gthdcn_t
{
  master_t *bsdg;
  float *hConv1,*hConv2,*pFFT1,*pFFT2,*Gmin,*Gmax,*fptmp1,*fptmp2,*Ddata1,*Ddata2,*pWorkin,*pWorkout,*pWorkoutgth,*pWorkoutgthm;
  float *tgood,*agood,*tmpphase,*pOrigin,*pMirror,*l_f,*pPack1,*pPack2,*pPack3,*taper,*taper2,*gdelay;
  float *tmp,*recxz,resrecz;
  int *vtgate,*voverlap,zntgate,*freqref;
  fftwf_plan plan[2];
  int ithread,ipky,nisky,ntgate,nsgate,nsgate2,fsgate,lsgate,fsgate2,lsgate2,rank,lfid,hfid,ihfcut,*toff;
  float *osFFT1,*osPack1,*oshConv1,*freqboost;
  float *pdesig1,*pdesig2;
 } ;
typedef struct gthdcn_t gthdcn_thread_t;
struct gthdcn_fw
{
  master_t *bsdg;

  float *pOrigintp,*pOrigintpall,*pOrigin;
  int ithread,nthread,isgate,itp,lsgate;

  float *costable2,*xkey,*ykey,*px,*py;
  int *pxgate,*ipxbeg,*ipxend;
  fftwf_plan plan[2];  
  fftwf_plan plan8[2];  
  fftwf_plan plan_rs[2];  
  fftwf_plan plan_rs8[2];  
  int rank;
 } ;
typedef struct gthdcn_fw gthdcn_thread_fw;
//#pragma pack(push, 1)		/* push current alignment to stack */
struct HEAD_T
{
  int flag;             //flag, by bit
  int sl, xl, live,datadeadtr,ztwbid,mwbid,tosid,bosid,dpid;   	//subline & crossline
  int tid, tidinp, intpflag, cableidintp, cabtrintp, hid, wid,wbid; //trace id, used to ident in TANGO DB
  float shotx,shoty,shotz,px,pz,recx,recy,recz,recz_smooth,offsetx,offsety,offsetz,recz_bsdg, regparam;
  float recvel; // For spatially variable velocities at receiver locations 
  int shootdir;
  float horz1,horz2,horz3;  // for horizon attributes
  int ipky, isky, iskyintp, offs, oreg;       //primary and secondary key index   (0:n-1)
  float offsx, offsy;  
  int pval, sval;       //shot and chan value   (attrs)
//  int wina[2];          //start of data window
//  int winz[2];          //end of data window
  float offset;          //shot/receiver depth and offset
  int l1wbid;
  int nnotchattr;
  /////////////////////////////////////
  /////////////////////////////////////
  int obnkeyx,obnkeyy,obnpkey,obnskey;
  float obnxpos,obnypos;
  /////////////////////////////////////
  /////////////////////////////////////

  /////////////////////////////////////
  /////////////////////////////////////
  int cableid,newcableid,cableidseq,pkeyindx,ntrcoupertask,zeroattr;
  /////////////////////////////////////
  /////////////////////////////////////
  int minlenattr,maxlenattr;
  float weight;
  float nmoweight2,nmoweight1;
  int numtracex,nsky3dloc;
  int corrstart,corrend;
  float theta_wb;
  float ampth;
  float gunamp, tdelay;
};
//#pragma pack(pop)		/* restore normal alignment from stack */

#define NHEAD       ((int)(sizeof(head_t)/sizeof(int)))

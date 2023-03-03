#include "msdginterp.h"
#include <complex>
using namespace std;

extern "C" {
  void agc2df90_(float *winagc2d, float *taptagc2d, float *sclmaxagc2d, 
		    int *nsamp, int *isi, int *nsgtrace, float *scalaragc2d, float *input);
}

void agc2d (l1inv_t *l1para, int nsgtrace, float **input, float **winoutscl)
{
  if(nsgtrace == 0) return;

  int itrc,isamp;
  float winagc2d, taptagc2d, sclmaxagc2d;
  int isi,nsamp;

  float *scalaragc2d;
  float *inputf90;

  scalaragc2d = l1para->scalaragc2d;

  memset(scalaragc2d,0,l1para->nsamp*sizeof(float));

  isi = l1para->srate;
  nsamp = l1para->nsamp;
  winagc2d = l1para->winagc2d;
  taptagc2d = l1para->taptagc2d;
  sclmaxagc2d = l1para->sclmaxagc2d;

  inputf90 = (float *) malloc((size_t)nsamp*(size_t)nsgtrace*sizeof(float));
  for (itrc=0;itrc<nsgtrace;itrc++)
     for (isamp=0;isamp<nsamp;isamp++)
       inputf90[nsamp*itrc+isamp] = input[itrc][isamp];

  agc2df90_(&winagc2d,&taptagc2d,&sclmaxagc2d,&nsamp,&isi,&nsgtrace,scalaragc2d,inputf90);

  for (itrc=0;itrc<nsgtrace;itrc++)
    for (isamp=0;isamp<nsamp;isamp++)
      {
	winoutscl[itrc][isamp] = scalaragc2d[isamp];
      }

  free(inputf90);

}

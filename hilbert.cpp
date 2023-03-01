#include "hilbert.h"

hilbert::hilbert()
{
    nt = 61;
    hw = (float*)malloc((nt*2+1)*sizeof(float));
    hw[nt]=0.0;
    float inv=3.14159/nt;
    float twodpi=2.0/3.14159;
    for(int i = 1; i<=nt; i++){
        float taper = .54+0.46*cosf(i*inv);
	hw[nt+i]    = -taper*(i%2)*twodpi/i;
	hw[nt-i]    = -hw[nt+i];
    }
}

void hilbert::hilbertit(float* a, float* b, int n)
{
    int lx  = nt*2+1;
    int ifx =-nt;
    int ly  = n;
    int ify = 0;
    int lz  = n;
    int ifz = 0;

    int ilx = ifx+lx-1;
    int ily = ify+ly-1;
    int ilz = ifz+lz-1;
    for(int i= ifz; i<=ilz; i++){
        int jlow = i-ily;
	if(jlow<ifx) jlow = ifx;
	int jhigh = i-ify;
	if(jhigh>ilx) jhigh = ilx;
	float sum = 0.0;
	for(int j=jlow; j<=jhigh; j++){
	    sum += hw[j-ifx]*a[i-j-ify];
	}
	b[i-ifz]=sum;
    }
}

void hilbert::amplitude(float* a, float* b, int n)
{
    hilbertit( a,  b,  n);
    for(int i=0;i<n;i++) b[i]=sqrt(b[i]*b[i]+a[i]*a[i]);
}

hilbert::~hilbert()
{
  if(hw) free(hw); hw=NULL;
}

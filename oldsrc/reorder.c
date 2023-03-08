#include "msdginterp.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "invmatrix.h"

int calc_reorder_indx(int dim1, int dim2, float **taup, float threshold, int *indx); //// taup[npx][dim2]
void reorder_taup(int dim1, int dim2, float **array, int *indx); //// array[npx][dim2]

void bsdg_indexx_(int *n,float *arr,int *indx);

int reorder(int npx, l1inv_t *l1para, int *indx, float **taup_p, float pscale) 
 {

   int npxnew = calc_reorder_indx(npx, l1para->fftnr, taup_p, l1para->npperc*pscale, indx);

   reorder_taup(npx, l1para->fftnr, taup_p,  indx);

   return npxnew;   

 }

int calc_reorder_indx(int dim1, int dim2, float **taup, float threshold, int *indx) //// taup[npx][dim2]
 {
   
   int i,j;

   double ampl2;

   double energymax = 0.0;

   float *energy;

   int dim1new;

   energy = calloc(dim1,sizeof(double));

   energymax = -100000.0f;
   for (i=0;i<dim1;i++)
     { 
       energy[i] = 0.0f; 
       for (j=0;j<dim2;j++)
	 {  
	   ampl2 = taup[i][j];
	   ampl2 = ampl2 * ampl2;
	   energy[i] = energy[i] + ampl2;
	 }
       energy[i] = -sqrt(energy[i]);  //// to get the output from bsdg_indexx in descending order of magnitude 
       energymax = MAX(energymax,-energy[i]);
     }
   
   bsdg_indexx_(&dim1,energy,indx);

   for(i=0; i<dim1; i++)indx[i] = indx[i]-1;  //// FORTRAN TO C

   dim1new = MAX(1, MIN((int)(dim1*threshold), dim1));
   
   free(energy);

   return dim1new;   

 }


void reorder_taup(int dim1, int dim2, float **array, int *indx) //// array[npx][fftnr]
 {

   int i,k;
   float **arraysave;
   
   arraysave = (float **) malloc(dim1*sizeof(float *));
   for(i=0; i<dim1; i++) arraysave[i] = (float *)calloc(dim2,sizeof(float));
   
   for(i=0; i<dim1; i++) memcpy(arraysave[i],array[i],dim2*sizeof(float));
   
   for(i=0; i<dim1; i++)
     {
       k = indx[i];
       memcpy(array[i],arraysave[k],dim2*sizeof(float));
     }
   
   for(i=0; i<dim1; i++) free(arraysave[i]); 
   free(arraysave); 

   return;

 }



void merge_indices_dsdg3d(l1inv_t *l1para, float **taup, float **taup_res, float **taup_merge, 
			  int npx1, int npx2, int np_orig, 
			  int *indx1, int *indx2)
 {

   int isamp,ip,iporig;

   myinit(taup_merge, np_orig, l1para->fftnr, 0); 

   for (ip=0;ip<npx1;ip++)
     {
       iporig = indx1[ip];
       for (isamp=0;isamp<l1para->fftnr;isamp++)
	 {
	   taup_merge[iporig][isamp] = taup [ip][isamp];
	 }
     }

   for (ip=0;ip<npx2;ip++)
     {
       iporig = indx2[ip];
       for (isamp=0;isamp<l1para->fftnr;isamp++)
	 {
	   taup_merge[iporig][isamp] = taup_merge[iporig][isamp] + taup_res[ip][isamp];
	 }
     }
   
   for (ip=0;ip<np_orig;ip++)
     {
       for (isamp=0;isamp<l1para->fftnr;isamp++)
	 {
	   taup[ip][isamp] = taup_merge[ip][isamp];
	 }
     }

   return;
   
 }



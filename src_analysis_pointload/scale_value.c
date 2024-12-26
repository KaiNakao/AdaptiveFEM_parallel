#include "omp.h"

void pack_scalevalue_(nn, mpibuf_f, mpibuf_i)
     int * restrict nn;
     short * restrict mpibuf_f;
     short * restrict mpibuf_i;
{
   int i;
#pragma omp for
   for(i=0; i<npc*12*(*nn); i++){
      mpibuf_i[i] = mpibuf_f[i];
   }
}

void unpack_scalevalue_(nn, mpibuf_i, mpibuf_f)
     int * restrict nn;
     short * restrict mpibuf_i;
     short * restrict mpibuf_f;
{
   int i;
#pragma omp for
   for(i=0; i<npc*12*(*nn); i++){
      mpibuf_f[i] = mpibuf_i[i];
   }
}


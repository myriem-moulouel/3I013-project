/* -------------------- */
/* --- sigmadelta.c --- */
/* -------------------- */

/*
* Copyright (c) 2004 - 2021, Lionel Lacassagne, All rights reserved
* University of Paris Sud, Laboratoire de Recherche en Informatique
* Sorbonne University - Laboratoire d'Informatique de Paris 6
*/

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#include "nrtype.h"
#include "nrdef.h"
#include "nrutil.h"
#include "macro.h"

#include "sigmadelta.h"

// -----------------------------------------------------------------------------------------
void SigmaDelta_Step0_line(uint8 *I, uint8 *M, uint8 *O, uint8 *V, uint8 *E, int j0, int j1)
// -----------------------------------------------------------------------------------------
{
    uint8 m, i, o, v, e;
    for(int j=j0; j<=j1; j++){
        m=load1(M,j);
        i=load1(I,j);
        o=load1(O,j);
        v=load1(V,j);
        e=load1(E,j);

        if(m<i){
            m = m + 1;
            o = i - m;
        }else{
            if(m>i){
                m = m - 1;
                o = m - i;
            }else{
                o = 0;
            }
        }
        if(v<o){
            v = v + 1;
        }else{
            if(v>o){
                v = v - 1;
            }
        }
        uint8 v = min2(v,SD_VMAX);
        v = max2(v,SD_VMIN);
        v = v;
        if(o<v){
            e = 0;
        }else{
            e = 1;
        }
        store1(M,j,m);
        store1(I,j,i);
        store1(O,j,i);
        store1(V,j,v);
        store1(E,j,e);
    }
}
// ------------------------------------------------------------------------------------------------
void SigmaDelta_1Step_line(uint8 *I, uint8 *M, uint8 *O, uint8 *V, uint8 *E, int k, int j0, int j1)
// ------------------------------------------------------------------------------------------------
{
    uint8 m, i, o, v, e;
    for(int j=j0; j<=j1; j++){
        m=load1(M,j);
        i=load1(I,j);
        o=load1(O,j);
        v=load1(V,j);
        e=load1(E,j);

        if(m<i){
            m = m + 1;
            o = i - m;
        }else{
            if(m>i){
                m = m - 1;
                o = m - i;
            }else{
                o = 0;
            }
        }
        if(v<k*o){
            v = v + 1;
        }else{
            if(v>k*o){
                v = v - 1;
            }
        }
        uint8 v = min2(v,SD_VMAX);
        v = max2(v,SD_VMIN);
        v = v;
        if(o<v){
            e = 0;
        }else{
            e = 1;
        }
        store1(M,j,m);
        store1(I,j,i);
        store1(O,j,i);
        store1(V,j,v);
        store1(E,j,e);
    }
}
// ---------------------------------------------------------------------------------------------------------
void SigmaDelta_Step0(uint8 **I, uint8 **M, uint8 **O, uint8 **V, uint8 **E, int i0, int i1, int j0, int j1)
// ---------------------------------------------------------------------------------------------------------
{
    for(int i=i0; i<=i1;i++){
        SigmaDelta_Step0_line(I[i],M[i],O[i],V[i],E[i],j0,j1);
    }
}
// ----------------------------------------------------------------------------------------------------------------
void SigmaDelta_1Step(uint8 **I, uint8 **M, uint8 **O, uint8 **V, uint8 **E, int k, int i0, int i1, int j0, int j1)
// ----------------------------------------------------------------------------------------------------------------
{
    for(int i=i0; i<=i1;i++){
        SigmaDelta_1Step_line(I[i],M[i],O[i],V[i],E[i],k,j0,j1);
    }
}

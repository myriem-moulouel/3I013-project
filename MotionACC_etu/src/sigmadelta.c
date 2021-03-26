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

#include "sigmadelta.h"

// -----------------------------------------------------------------------------------------
void SigmaDelta_Step0_line(uint8 *I, uint8 *M, uint8 *O, uint8 *V, uint8 *E, int j0, int j1)
// -----------------------------------------------------------------------------------------
{
    uint8 mj, ij, oj, vj, ej;
    for(int j=j0; j<=j1; j++){
        mj=load1(M,j);
        ij=load1(I,j);
        oj=load1(O,j);
        vj=load1(V,j);
        ej=load1(E,j);

        if(mj<ij){
            mj = mj + 1;
            oj = ij - mj;
        }else{
            if(mj>ij){
                mj = mj - 1;
                oj = mj - ij;
            }else{
                oj = 0;
            }
        }
        if(vj<oj){
            vj = vj + 1;
        }else{
            if(vj>oj){
                vj = vj - 1;
            }
        }
        uint8 v = min2(vj,SD_VMAX);
        v = max2(v,SD_VMIN);
        vj = v;
        if(oj<vj){
            ej = 0;
        }else{
            ej = 1;
        }
        store1(M,j,mj);
        store1(I,j,ij);
        store1(O,j,ij);
        store1(V,j,vj);
        store1(E,j,ej);
    }
}
// ------------------------------------------------------------------------------------------------
void SigmaDelta_1Step_line(uint8 *I, uint8 *M, uint8 *O, uint8 *V, uint8 *E, int k, int j0, int j1)
// ------------------------------------------------------------------------------------------------
{
    uint8 mj, ij, oj, vj, ej;
    for(int j=j0; j<=j1; j++){
        mj=load1(M,j);
        ij=load1(I,j);
        oj=load1(O,j);
        vj=load1(V,j);
        ej=load1(E,j);

        if(mj<ij){
            mj = mj + 1;
            oj = ij - mj;
        }else{
            if(mj>ij){
                mj = mj - 1;
                oj = mj - ij;
            }else{
                oj = 0;
            }
        }
        if(vj<k*oj){
            vj = vj + 1;
        }else{
            if(vj>k*oj){
                vj = vj - 1;
            }
        }
        uint8 v = min2(vj,SD_VMAX);
        v = max2(v,SD_VMIN);
        vj = v;
        if(oj<vj){
            ej = 0;
        }else{
            ej = 1;
        }
        store1(M,j,mj);
        store1(I,j,ij);
        store1(O,j,ij);
        store1(V,j,vj);
        store1(E,j,ej);
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

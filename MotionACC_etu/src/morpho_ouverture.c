/* -------------------------- */
/* --- morpho_ouverture.c --- */
/* -------------------------- */

/*
 * Copyright (c) 2020 - 2021, Lionel Lacassagne, All rights reserved
  
 
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#include "nrtype.h"
#include "nrdef.h"
#include "nrutil.h"
//#include "sequence.h"

#include "swp.h"
#include "morpho_min.h"
#include "morpho_max.h"

// -------------------------------------------------------------------------------
void line_ouverture3_ui8matrix_fusion(uint8 **X, int i, int j0, int j1, uint8 **Y)
// -------------------------------------------------------------------------------
{
    uint8 x00, x01, x02;
    uint8 x10, x11, x12;
    uint8 x20, x21, x22;

    uint8 min00, min01, min02;
    uint8 min10, min11, min12;
    uint8 min20, min21, min22;

    uint8 maxi;

    for(int j=j0; j<=j1; j++){
        x00 = load2(X,i-2,j-2);     x01 = load2(X,i-2,j-1);     x02 = load2(X,i-2, j );
        x10 = load2(X,i-1,j-2);     x11 = load2(X,i-1,j-1);     x12 = load2(X,i-1, j );
        x20 = load2(X, i ,j-2);     x21 = load2(X, i ,j-1);     x22 = load2(X, i , j );
        min9(x00,x01,x02,x10,x11,x12,x20,x21,x22,min00);
        x00 = x01;                  x01 = x02;                  x02 = load2(X,i-2,j+1);
        x10 = x11;                  x11 = x12;                  x12 = load2(X,i-1,j+1);
        x20 = x21;                  x21 = x22;                  x22 = load2(X, i ,j+1);
        min9(x00,x01,x02,x10,x11,x12,x20,x21,x22,min01);
        x00 = x01;                  x01 = x02;                  x02 = load2(X,i-2,j+2);
        x10 = x11;                  x11 = x12;                  x12 = load2(X,i-1,j+2);
        x20 = x21;                  x21 = x22;                  x22 = load2(X, i ,j+2);
        min9(x00,x01,x02,x10,x11,x12,x20,x21,x22,min02);


        x00 = load2(X,i-1,j-2);     x01 = load2(X,i-1,j-1);     x02 = load2(X,i-1, j );
        x10 = load2(X, i ,j-2);     x11 = load2(X, i ,j-1);     x12 = load2(X, i , j );
        x20 = load2(X,i+1,j-2);     x21 = load2(X,i+1,j-1);     x22 = load2(X,i+1, j );
        min9(x00,x01,x02,x10,x11,x12,x20,x21,x22,min10);
        x00 = x01;                  x01 = x02;                  x02 = load2(X,i-1,j+1);
        x10 = x11;                  x11 = x12;                  x12 = load2(X, i ,j+1);
        x20 = x21;                  x21 = x22;                  x22 = load2(X,i+1,j+1);
        min9(x00,x01,x02,x10,x11,x12,x20,x21,x22,min11);
        x00 = x01;                  x01 = x02;                  x02 = load2(X,i-1,j+2);
        x10 = x11;                  x11 = x12;                  x12 = load2(X, i ,j+2);
        x20 = x21;                  x21 = x22;                  x22 = load2(X,i+1,j+2);
        min9(x00,x01,x02,x10,x11,x12,x20,x21,x22,min12);

        
        x00 = load2(X, i ,j-2);     x01 = load2(X, i ,j-1);     x02 = load2(X, i , j );
        x10 = load2(X,i+1,j-2);     x11 = load2(X,i+1,j-1);     x12 = load2(X,i+1, j );
        x20 = load2(X,i+2,j-2);     x21 = load2(X,i+2,j-1);     x22 = load2(X,i+2, j );
        min9(x00,x01,x02,x10,x11,x12,x20,x21,x22,min20);
        x00 = x01;                  x01 = x02;                  x02 = load2(X, i ,j+1);
        x10 = x11;                  x11 = x12;                  x12 = load2(X,i+1,j+1);
        x20 = x21;                  x21 = x22;                  x22 = load2(X,i+2,j+1);
        min9(x00,x01,x02,x10,x11,x12,x20,x21,x22,min21);
        x00 = x01;                  x01 = x02;                  x02 = load2(X, i ,j+2);
        x10 = x11;                  x11 = x12;                  x12 = load2(X,i+1,j+2);
        x20 = x21;                  x21 = x22;                  x22 = load2(X,i+2,j+2);
        min9(x00,x01,x02,x10,x11,x12,x20,x21,x22,min22);
        


        max9(min00,min01,min02,min10,min11,min12,min20,min21,min22,maxi);
        store2(Y,i,j,maxi);
    }
}
// ----------------------------------------------------------------------------------------
void line_ouverture3_ui8matrix_fusion_ilu5_red(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ----------------------------------------------------------------------------------------
{
    /*
    int r = (j1-j0+3)%5;
    for(int j=j0-1; j<j1+1-r; j+=5){
        line_min3_ui8matrix_ilu3_red(X,i-1, j ,j+4,Z);
        line_min3_ui8matrix_ilu3_red(X, i , j ,j+4,Z);
        line_min3_ui8matrix_ilu3_red(X,i+1, j ,j+4,Z);

        line_max3_ui8matrix_ilu3_red(Z,i  , j ,j+4,Y);
    }
    switch (r){
    case 1:
        line_min3_ui8matrix_ilu3_red(X,i-1,j1-1,j1+1,Z);
        line_min3_ui8matrix_ilu3_red(X, i ,j1-1,j1+1,Z);
        line_min3_ui8matrix_ilu3_red(X,i+1,j1-1,j1+1,Z);

        line_max3_ui8matrix_ilu3_red(Z,i  ,j1-1,j1+1,Y);
    case 2:
        line_min3_ui8matrix_ilu3_red(X,i-1,j1-2,j1+1,Z);
        line_min3_ui8matrix_ilu3_red(X, i ,j1-2,j1+1,Z);
        line_min3_ui8matrix_ilu3_red(X,i+1,j1-2,j1+1,Z);

        line_max3_ui8matrix_ilu3_red(Z,i  ,j1-2,j1+1,Y);
    case 3:
        line_min3_ui8matrix_ilu3_red(X,i-1,j1-3,j1+1,Z);
        line_min3_ui8matrix_ilu3_red(X, i ,j1-3,j1+1,Z);
        line_min3_ui8matrix_ilu3_red(X,i+1,j1-3,j1+1,Z);

        line_max3_ui8matrix_ilu3_red(Z,i  ,j1-3,j1+1,Y);
    case 4:
        line_min3_ui8matrix_ilu3_red(X,i-1,j1-4,j1+1,Z);
        line_min3_ui8matrix_ilu3_red(X, i ,j1-4,j1+1,Z);
        line_min3_ui8matrix_ilu3_red(X,i+1,j1-4,j1+1,Z);

        line_max3_ui8matrix_ilu3_red(Z,i  ,j1-4,j1+1,Y);
    }
    */
}
// ---------------------------------------------------------------------------------------------
void line_ouverture3_ui8matrix_fusion_ilu5_elu2_red(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ---------------------------------------------------------------------------------------------
{
    /*
    int r = (j1-j0+3)%5;
    for(int j=j0-1; j<j1+1-r; j+=5){
        line_min3_ui8matrix_ilu3_red(X,i-1, j ,j+4,Z);
        line_min3_ui8matrix_ilu3_red(X, i , j ,j+4,Z);
        line_min3_ui8matrix_ilu3_red(X,i+1, j ,j+4,Z);
        line_min3_ui8matrix_ilu3_red(X,i+2, j ,j+4,Z);

        line_max3_ui8matrix_ilu3_red(Z,i  , j ,j+4,Y);
        line_max3_ui8matrix_ilu3_red(Z,i+1, j ,j+4,Y);
    }
    switch (r){
    case 1:
        line_min3_ui8matrix_ilu3_red(X,i-1,j1,j1+1,Z);
        line_min3_ui8matrix_ilu3_red(X, i ,j1,j1+1,Z);
        line_min3_ui8matrix_ilu3_red(X,i+1,j1,j1+1,Z);
        line_min3_ui8matrix_ilu3_red(X,i+2,j1,j1+1,Z);

        line_max3_ui8matrix_ilu3_red(Z,i  ,j1,j1+1,Y);

        line_max3_ui8matrix_ilu3_red(Z,i+1,j1,j1+1,Y);
    case 2:
        line_min3_ui8matrix_ilu3_red(X,i-1,j1-1,j1+1,Z);
        line_min3_ui8matrix_ilu3_red(X, i ,j1-1,j1+1,Z);
        line_min3_ui8matrix_ilu3_red(X,i+1,j1-1,j1+1,Z);
        line_min3_ui8matrix_ilu3_red(X,i+2,j1-1,j1+1,Z);

        line_max3_ui8matrix_ilu3_red(Z,i  ,j1-1,j1+1,Y);

        line_max3_ui8matrix_ilu3_red(Z,i+1,j1-1,j1+1,Y);
    case 3:
        line_min3_ui8matrix_ilu3_red(X,i-1,j1-2,j1+1,Z);
        line_min3_ui8matrix_ilu3_red(X, i ,j1-2,j1+1,Z);
        line_min3_ui8matrix_ilu3_red(X,i+1,j1-2,j1+1,Z);
        line_min3_ui8matrix_ilu3_red(X,i+2,j1-2,j1+1,Z);

        line_max3_ui8matrix_ilu3_red(Z,i  ,j1-2,j1+1,Y);

        line_max3_ui8matrix_ilu3_red(Z,i+1,j1-2,j1+1,Y);
    case 4:
        line_min3_ui8matrix_ilu3_red(X,i-1,j1-3,j1+1,Z);
        line_min3_ui8matrix_ilu3_red(X, i ,j1-3,j1+1,Z);
        line_min3_ui8matrix_ilu3_red(X,i+1,j1-3,j1+1,Z);
        line_min3_ui8matrix_ilu3_red(X,i+2,j1-3,j1+1,Z);

        line_max3_ui8matrix_ilu3_red(Z,i  ,j1-3,j1+1,Y);

        line_max3_ui8matrix_ilu3_red(Z,i+1,j1-3,j1+1,Y);
    }
    */
}
// ----------------------------------------------------------------------------------------------------
void line_ouverture3_ui8matrix_fusion_ilu5_elu2_red_factor(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ----------------------------------------------------------------------------------------------------
{
}
// -----------------------------------------------------------------------------------------
void line_ouverture3_ui8matrix_fusion_ilu15_red(uint8 **X, int i, int j0, int j1, uint8 **Y)
// -----------------------------------------------------------------------------------------
{
}
// ---------------------------------------------------------------------------------------------
void ouverture3_ui8matrix_basic(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y, uint8 **Z)
// ---------------------------------------------------------------------------------------------
{
    /*
    printf("X\n");
    for(int i=i0; i<=i1; i++) {
        for(int j=j0; j<=j1; j++){
            printf("%d  ", X[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    */
    min3_ui8matrix_basic(X, i0-1, i1+1, j0-1, j1+1, Y);
    max3_ui8matrix_basic(Y, i0,   i1,   j0,   j1,   Z);
    
    /*printf("Z\n");
    for(int i=i0; i<=i1; i++) {
        for(int j=j0; j<=j1; j++){
            printf("%d  ", Z[i][j]);
        }
        printf("\n");
    }
    printf("\n");*/
    
}
// -----------------------------------------------------------------------------------
void ouverture3_ui8matrix_fusion(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// -----------------------------------------------------------------------------------
{
    for(int i=i0; i<=i1; i++){
        line_ouverture3_ui8matrix_fusion(X,i,j0,j1,Y);
    }
}
// --------------------------------------------------------------------------------------------
void ouverture3_ui8matrix_fusion_ilu5_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// --------------------------------------------------------------------------------------------
{
    /*
    for(int i=i0; i<=i1; i++){
        line_ouverture3_ui8matrix_fusion_ilu5_red(X,i,j0,j1,Z,Y);
    }
    printf("Z\n");
    for(int i=i0; i<=i1; i++) {
        for(int j=j0; j<=j1; j++){
            printf("%d  ", Y[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    */
}
// -------------------------------------------------------------------------------------------------
void ouverture3_ui8matrix_fusion_ilu5_elu2_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// -------------------------------------------------------------------------------------------------
{
    /*
    int r = (i1 - i0)%2;
    for(int i=i0; i<=i1-r; i+=2){
        line_ouverture3_ui8matrix_fusion_ilu5_elu2_red(X,i,j0,j1,Z,Y);
    }
    if(r==1){
        line_ouverture3_ui8matrix_fusion_ilu5_red(X,i1,j0,j1,Z,Y);
    }
    */
}
// --------------------------------------------------------------------------------------------------------
void ouverture3_ui8matrix_fusion_ilu5_elu2_red_factor(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// --------------------------------------------------------------------------------------------------------
{
}
// ---------------------------------------------------------------------------------------------
void ouverture3_ui8matrix_fusion_ilu15_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// ---------------------------------------------------------------------------------------------
{
}
// ------------------------------------------------------------------------------------------------------
void ouverture3_ui8matrix_pipeline_basic(uint8 **X, int i0, int i1, int j0, int j1, uint8 **T, uint8 **Y)
// ------------------------------------------------------------------------------------------------------
{
    for(int i=i0-1; i<=i0; i++){
        line_min3_ui8matrix_basic(X,i,j0,j1,T);
    }

    for(int i=i0; i<=i1; i++){
        line_min3_ui8matrix_basic(X,i+1,j0,j1,T);

        line_max3_ui8matrix_basic(T,i,j0,j1,Y);
    }
}
// ----------------------------------------------------------------------------------------------------
void ouverture3_ui8matrix_pipeline_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **T, uint8 **Y)
// ----------------------------------------------------------------------------------------------------
{
    for(int i=i0-1; i<=i0; i++){
        line_min3_ui8matrix_red(X,i,j0,j1,T);
    }

    for(int i=i0; i<=i1; i++){
        line_min3_ui8matrix_red(X,i+1,j0,j1,T);

        line_max3_ui8matrix_red(T,i,j0,j1,Y);
    }
}
// ---------------------------------------------------------------------------------------------------------
void ouverture3_ui8matrix_pipeline_ilu3_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **T, uint8 **Y)
// ---------------------------------------------------------------------------------------------------------
{
    for(int i=i0-1; i<=i0; i++){
        line_min3_ui8matrix_ilu3_red(X,i,j0,j1,T);
    }

    for(int i=i0; i<=i1; i++){
        line_min3_ui8matrix_ilu3_red(X,i+1,j0,j1,T);

        line_max3_ui8matrix_ilu3_red(T,i,j0,j1,Y);
    }
}
// ---------------------------------------------------------------------------------------------------------
void ouverture3_ui8matrix_pipeline_elu2_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **T, uint8 **Y)
// ---------------------------------------------------------------------------------------------------------
{
    for(int i=i0-1; i<=i0; i++){
        line_min3_ui8matrix_elu2_red(X,i,j0,j1,T);
    }    

    int r = (i1-i0+1)%2;
    for(int i=i0; i<i1-r; i+=2){
        line_min3_ui8matrix_elu2_red(X,i+1,j0,j1,T);
        line_max3_ui8matrix_elu2_red(T, i ,j0,j1,Y); 
    }

    if(r==1){
        line_min3_ui8matrix_ilu3_red(X,i1+1,j0,j1,T);

        line_max3_ui8matrix_ilu3_red(T,i1,j0,j1,Y);
    }
}
// ----------------------------------------------------------------------------------------------------------------
void ouverture3_ui8matrix_pipeline_elu2_red_factor(uint8 **X, int i0, int i1, int j0, int j1, uint8 **T, uint8 **Y)
// ----------------------------------------------------------------------------------------------------------------
{
    for(int i=i0-1; i<=i0; i++){
        line_min3_ui8matrix_elu2_red_factor(X,i,j0,j1,T);
    }    

    int r = (i1-i0+1)%2;
    for(int i=i0; i<i1-r; i+=2){
        line_min3_ui8matrix_elu2_red_factor(X,i+1,j0,j1,T);
        line_max3_ui8matrix_elu2_red_factor(T, i ,j0,j1,Y); 
    }

    if(r==1){
        line_min3_ui8matrix_ilu3_red(X,i1+1,j0,j1,T);

        line_max3_ui8matrix_ilu3_red(T,i1,j0,j1,Y);
    }
}
// --------------------------------------------------------------------------------------------------------------
void ouverture3_ui8matrix_pipeline_ilu3_elu2_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **T, uint8 **Y)
// --------------------------------------------------------------------------------------------------------------
{
    for(int i=i0-1; i<=i0; i++){
        line_min3_ui8matrix_ilu3_elu2_red(X,i,j0,j1,T);
    }    

    int r = (i1-i0+1)%2;
    for(int i=i0; i<i1-r; i+=2){
        line_min3_ui8matrix_ilu3_elu2_red(X,i+1,j0,j1,T);
        line_max3_ui8matrix_ilu3_elu2_red(T, i ,j0,j1,Y); 
    }

    if(r==1){
        line_min3_ui8matrix_ilu3_red(X,i1+1,j0,j1,T);

        line_max3_ui8matrix_ilu3_red(T,i1,j0,j1,Y);
    }
}
// ---------------------------------------------------------------------------------------------------------------------
void ouverture3_ui8matrix_pipeline_ilu3_elu2_red_factor(uint8 **X, int i0, int i1, int j0, int j1, uint8 **T, uint8 **Y)
// ---------------------------------------------------------------------------------------------------------------------
{
    for(int i=i0-1; i<=i0; i++){
        line_min3_ui8matrix_ilu3_elu2_red_factor(X,i,j0,j1,T);
    }    

    int r = (i1-i0+1)%2;
    for(int i=i0; i<i1-r; i+=2){
        line_min3_ui8matrix_ilu3_elu2_red_factor(X,i+1,j0,j1,T);
        line_max3_ui8matrix_ilu3_elu2_red_factor(T, i ,j0,j1,Y); 
    }

    if(r==1){
        line_min3_ui8matrix_ilu3_red(X,i1+1,j0,j1,T);

        line_max3_ui8matrix_ilu3_red(T,i1,j0,j1,Y);
    }
}
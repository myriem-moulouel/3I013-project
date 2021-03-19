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
    /*
    line_min3_ui8matrix_rot(X,i-1,j0-1,j1+1,Z);
    line_min3_ui8matrix_rot(X,i  ,j0-1,j1+1,Z);
    line_min3_ui8matrix_rot(X,i+1,j0-1,j1+1,Z);

    line_max3_ui8matrix_rot(Z,i,j0,j1,Y);
    */
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
{/*
    for(int i=i0; i<=i1; i++){
        line_ouverture3_ui8matrix_fusion(X,i,j0,j1,Z,Y);
    }
    */
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
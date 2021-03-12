/* -------------------- */
/* --- morpho_min.c --- */
/* -------------------- */

/*
 * Copyright (c) 2004 - 2013, Lionel Lacassagne, All rights reserved
 * University of Paris Sud, Laboratoire de Recherche en Informatique 
 * Creation: 2004-05-18 : creation, from an idea of Antoine Manzanera
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#include "nrtype.h"
#include "nrdef.h"
#include "nrutil.h"

#define load1(X,i) X[i]
#define load2(X,i,j) X[i][j]
#define store2(x,i,j,Y) do{ Y[i][j]=x; }while(0)
#define min2(a,b) (a<b)?a:b
#define min3(a,b,c,d) do{\
    if(a<b){        \
        if(a<c){    \
            d=a;    \
        }else{      \
            d=c;    \
        }           \
    }else{          \
        if(b<c){    \
            d=b;    \
        }else{      \
            d=c;    \
        }           \
    }               \
}while(0)           \

#define min9(a,b,c,d,e,f,g,h,i,j) do{               \
    uint8 m1, m2, m3;                               \
    min3(a,b,c,m1); min3(d,e,f,m2); min3(g,h,i,m3); \
    min3(m1,m2,m3,j);                               \
}while(0);                                          \

// ------------------------------------------------------------------------
void line_min3_ui8matrix_basic(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ------------------------------------------------------------------------
{
    for(int j=j0; j<j1; j++){
        min9(X[i-1][j-1], X[i-1][j], X[i-1][j+1],
             X[i  ][j-1], X[i  ][j], X[i  ][j+1],
             X[i+1][j-1], X[i+1][j], X[i+1][j+1],Y[i][j]);
    }
}
// ----------------------------------------------------------------------
void line_min3_ui8matrix_reg(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ----------------------------------------------------------------------
{
    //avec scalarisation
    uint8 x11, x12, x13;
    uint8 x21, x22, x23;
    uint8 x31, x32, x33;
    uint8 y;

    for(int j = j0; j<=j1; j++){
        x11 = load2(X,i-1,j-1); x12 = load2(X,i-1,j); x13 = load2(X,i-1,j+1);
        x21 = load2(X, i ,j-1); x22 = load2(X, i ,j); x23 = load2(X, i ,j+1);
        x31 = load2(X,i+1,j-1); x32 = load2(X,i+1,j); x33 = load2(X,i+1,j+1);
        min9(x11, x12, x13,
             x21, x22, x23,
             x31, x32, x33, y);
        store2(y, i, j, Y);
    }
}
// ----------------------------------------------------------------------
void line_min3_ui8matrix_rot(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ----------------------------------------------------------------------
{
    //avec scalarisation et rotation
    uint8 x11, x12, x13;
    uint8 x21, x22, x23;
    uint8 x31, x32, x33;
    uint8 y;

    x12 = load2(X,i-1,j0-1); x13 = load2(X,i-1,j0);
    x22 = load2(X, i ,j0-1); x23 = load2(X, i ,j0);
    x32 = load2(X,i+1,j0-1); x33 = load2(X,i+1,j0);

    for(int j=j0; j<=j1; j++){
        x11 = x12; x12 = x13; x13 = load2(X,i-1,j+1);
        x21 = x22; x22 = x23; x23 = load2(X, i ,j+1);
        x31 = x32; x32 = x33; x33 = load2(X,i+1,j+1);
        min9(x11, x12, x13,
             x21, x22, x23,
             x31, x32, x33, y);
        store2(y,i,j,Y);
    }
}
// ----------------------------------------------------------------------
void line_min3_ui8matrix_red(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ----------------------------------------------------------------------
{
    //par colonne !
    uint8 a0, b0, c0;
    uint8 a1, b1, c1;
    uint8 a2, b2, c2;
    uint8 ra, rb, rc;
    uint8 y;
    //reduction de la colonne a
    a0 = load2(X,i-1,j0-1);
    a1 = load2(X, i ,j0-1);
    a2 = load2(X,i+1,j0-1);
    min3(a0,a1,a2,ra);
    //reduction de la colonne b
    b0 = load2(X,i-1, j0 );
    b1 = load2(X, i , j0 );
    b2 = load2(X,i+1, j0 );
    min3(b0,b1,b2,rb);
    
    for (int j=j0; j<=j1 ; j++) {
        
        //reduction de la colonne c
        c0 = load2(X,i-1,j+1);
        c1 = load2(X, i ,j+1);
        c2 = load2(X,i+1,j+1);
        min3(c0,c1,c2,rc);
        min3(ra,rb,rc,y);
        //operateur ligne
        store2(y,i,j,Y) ;
        ra = rb ;
        rb = rc ;
    }
}
// -----------------------------------------------------------------------
void line_min3_ui8matrix_ilu3(uint8 **X, int i, int j0, int j1, uint8 **Y)
// -----------------------------------------------------------------------
{
    //avec scalarisation et rotation
    uint8 *x1 = load1(X,i-1); uint8 x11, x12, x13;
    uint8 *x2 = load1(X,  i); uint8 x21, x22, x23;
    uint8 *x3 = load1(X,i+1); uint8 x31, x32, x33;
    uint8 y;
    x12 = load1(x1,j0-1); x13 = load1(x1,j0);
    x22 = load1(x2,j0-1); x23 = load1(x2,j0);
    x32 = load1(x3,j0-1); x33 = load1(x3,j0);

    int r = (j1-j0+1)%3;
    for(int j=j0; j<j1-r; j+=3){
        x11 = x12; x12 = x13; x13 = load1(x1,j+1);
        x21 = x22; x22 = x23; x23 = load1(x2,j+1);
        x31 = x32; x32 = x33; x33 = load1(x3,j+1);
        min9(x11, x12, x13,
             x21, x22, x23,
             x31, x32, x33, y);
        store2(y,i,j,Y);

        x11 = x12; x12 = x13; x13 = load1(x1,j+2);
        x21 = x22; x22 = x23; x23 = load1(x2,j+2);
        x31 = x32; x32 = x33; x33 = load1(x3,j+2);
        min9(x11, x12, x13,
             x21, x22, x23,
             x31, x32, x33, y);
        store2(y,i,j+1,Y);

        x11 = x12; x12 = x13; x13 = load1(x1,j+3);
        x21 = x22; x22 = x23; x23 = load1(x2,j+3);
        x31 = x32; x32 = x33; x33 = load1(x3,j+3);
        min9(x11, x12, x13,
             x21, x22, x23,
             x31, x32, x33, y);
        store2(y,i,j+2,Y);
    }
    
    if(r==1){
        x11 = x12; x12 = x13; x13 = load1(x1,j1+1);
        x21 = x22; x22 = x23; x23 = load1(x2,j1+1);
        x31 = x32; x32 = x33; x33 = load1(x3,j1+1);
        min9(x11, x12, x13,
             x21, x22, x23,
             x31, x32, x33, y);
        store2(y,i,j1,Y);
    }else{
        if(r==2){
            x11 = x12; x12 = x13; x13 = x1[ j1 ];
            x21 = x22; x22 = x23; x23 = x2[ j1 ];
            x31 = x32; x32 = x33; x33 = x3[ j1 ];
            min9(x11, x12, x13,
                 x21, x22, x23,
                 x31, x32, x33, y);
            store2(y,i,j1-1,Y);

            x11 = x12; x12 = x13; x13 = x1[j1+1];
            x21 = x22; x22 = x23; x23 = x2[j1+1];
            x31 = x32; x32 = x33; x33 = x3[j1+1];
            min9(x11, x12, x13,
                 x21, x22, x23,
                 x31, x32, x33, y);
            store2(y,i,j1,Y);
        }
    }
}
// ---------------------------------------------------------------------------
void line_min3_ui8matrix_ilu3_red(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ---------------------------------------------------------------------------
{
    //par colonne !
    uint8 a0, b0, c0;
    uint8 a1, b1, c1;
    uint8 a2, b2, c2;
    uint8 ra, rb, rc;
    uint8 y;
    
    int r = (j1-j0+1)%3;

    for (int j=j0; j<j1-r ; j+=3) {

        //reduction colonne a //reduction colonne b     //reduction colonne c
        a0 = load2(X,i-1,j-1);       b0 = load2(X,i-1, j );       c0 = load2(X,i-1,j+1);
        a1 = load2(X, i ,j-1);       b1 = load2(X, i , j );       c1 = load2(X, i ,j+1);
        a2 = load2(X,i+1,j-1);       b2 = load2(X,i+1, j );       c2 = load2(X,i+1,j+1);
        min3( a0, a1, a2, ra);       min3( b0, b1, b2, rb);       min3( c0, c1, c2, rc);

        min3( ra, rb, rc, y);
        //operateur ligne
        store2(y,i,j,Y) ;
        
        //reduction colonne a //reduction colonne b     //reduction colonne c
        a0 = load2(X,i-1, j );       b0 = load2(X,i-1,j+1);       c0 = load2(X,i-1,j+2);
        a1 = load2(X, i , j );       b1 = load2(X, i ,j+1);       c1 = load2(X, i ,j+2);
        a2 = load2(X,i+1, j );       b2 = load2(X,i+1,j+1);       c2 = load2(X,i+1,j+2);
        min3( a0, a1, a2, ra);       min3( b0, b1, b2, rb);       min3( c0, c1, c2, rc);

        min3( ra, rb, rc, y);
        //operateur ligne
        store2(y,i,j+1,Y) ;
        
        //reduction colonne a //reduction colonne b     //reduction colonne c
        a0 = load2(X,i-1,j+1);       b0 = load2(X,i-1,j+2);       c0 = load2(X,i-1,j+3);
        a1 = load2(X, i ,j+1);       b1 = load2(X, i ,j+2);       c1 = load2(X, i ,j+3);
        a2 = load2(X,i+1,j+1);       b2 = load2(X,i+1,j+2);       c2 = load2(X,i+1,j+3);
        min3( a0, a1, a2, ra);       min3( b0, b1, b2, rb);       min3( c0, c1, c2, rc);

        min3( ra, rb, rc, y);
        //operateur ligne
        store2(y,i,j+2,Y) ;
    }
    if(r==1){
        //reduction colonne a //reduction colonne b     //reduction colonne c
        a0 = load2(X,i-1,j1-1);       b0 = load2(X,i-1,j1);       c0 = load2(X,i-1,j1+1);
        a1 = load2(X, i ,j1-1);       b1 = load2(X, i ,j1);       c1 = load2(X, i ,j1+1);
        a2 = load2(X,i+1,j1-1);       b2 = load2(X,i+1,j1);       c2 = load2(X,i+1,j1+1);
        min3( a0, a1, a2, ra);       min3( b0, b1, b2, rb);       min3( c0, c1, c2, rc);

        min3( ra, rb, rc, y);
        //operateur ligne
        store2(y,i,j1,Y) ;
    }else{
        if(r==2){
            //reduction colonne a //reduction colonne b     //reduction colonne c
            a0 = load2(X,i-1,j1-2);       b0 = load2(X,i-1,j1-1);       c0 = load2(X,i-1, j1 );
            a1 = load2(X, i ,j1-2);       b1 = load2(X, i ,j1-1);       c1 = load2(X, i , j1 );
            a2 = load2(X,i+1,j1-2);       b2 = load2(X,i+1,j1-1);       c2 = load2(X,i+1, j1 );
            min3( a0, a1, a2, ra);       min3( b0, b1, b2, rb);       min3( c0, c1, c2, rc);

            min3( ra, rb, rc, y);
            //operateur ligne
            store2(y,i,j1-1,Y) ;

           //reduction colonne a //reduction colonne b     //reduction colonne c
            a0 = load2(X,i-1,j1-1);       b0 = load2(X,i-1,j1);       c0 = load2(X,i-1,j1+1);
            a1 = load2(X, i ,j1-1);       b1 = load2(X, i ,j1);       c1 = load2(X, i ,j1+1);
            a2 = load2(X,i+1,j1-1);       b2 = load2(X,i+1,j1);       c2 = load2(X,i+1,j1+1);
            min3( a0, a1, a2, ra);       min3( b0, b1, b2, rb);       min3( c0, c1, c2, rc);

            min3( ra, rb, rc, y);
            //operateur ligne
            store2(y,i,j1,Y) ;
        }
    }
}
// ---------------------------------------------------------------------------
void line_min3_ui8matrix_elu2_red(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ---------------------------------------------------------------------------
{
    //par colonne !
    uint8 a0, b0, c0;
    uint8 a1, b1, c1;
    uint8 a2, b2, c2;
    uint8 ra, rb, rc;
    uint8 y;
    
    for (int j=j0; j<=j1 ; j++) {
        //reduction colonne a //reduction colonne b     //reduction colonne c
        a0 = load2(X,i-1,j-1);       b0 = load2(X,i-1, j );       c0 = load2(X,i-1,j+1);
        a1 = load2(X, i ,j-1);       b1 = load2(X, i , j );       c1 = load2(X, i ,j+1);
        a2 = load2(X,i+1,j-1);       b2 = load2(X,i+1, j );       c2 = load2(X,i+1,j+1);
        min3( a0, a1, a2, ra);       min3( b0, b1, b2, rb);       min3( c0, c1, c2, rc);
        min3( ra, rb, rc, y);
        //operateur ligne
        store2(y,i,j,Y) ;
        //reduction colonne a //reduction colonne b     //reduction colonne c
        a0 = load2(X, i ,j-1);       b0 = load2(X, i , j );       c0 = load2(X, i ,j+1);
        a1 = load2(X,i+1,j-1);       b1 = load2(X,i+1, j );       c1 = load2(X,i+1,j+1);
        a2 = load2(X,i+2,j-1);       b2 = load2(X,i+2, j );       c2 = load2(X,i+2,j+1);
        min3( a0, a1, a2, ra);       min3( b0, b1, b2, rb);       min3( c0, c1, c2, rc);
        min3( ra, rb, rc, y);
        //operateur ligne
        store2(y,i+1,j,Y) ;

    }
}
// ----------------------------------------------------------------------------------
void line_min3_ui8matrix_elu2_red_factor(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ----------------------------------------------------------------------------------
{
    //version elu2_red avec rotation de registres
    //par colonne !
    uint8 a0, b0, c0;
    uint8 a1, b1, c1;
    uint8 a2, b2, c2;
    uint8 ra1, rb1, rc1;
    uint8 ra2, rb2, rc2;
    uint8 y;
    // pour i
    a0 = load2(X,i-1,j0-1);      b0 = load2(X,i-1, j0 );
    a1 = load2(X, i ,j0-1);      b1 = load2(X, i , j0 );
    a2 = load2(X,i+1,j0-1);      b2 = load2(X,i+1, j0 );
    min3( a0, a1, a2, ra1);      min3( b0, b1, b2, rb1);
    // pour i+1
    a0 = a1;                     b0 = b1;
    a1 = a2;                     b1 = b2;
    a2 = load2(X,i+2,j0-1);      b2 = load2(X,i+2, j0 );
    min3( a0, a1, a2, ra2);      min3( b0, b1, b2, rb2);

    for (int j=j0; j<=j1 ; j++) {
        
        //reduction de la colonne c
        c0 = load2(X,i-1,j+1);
        c1 = load2(X, i ,j+1);
        c2 = load2(X,i+1,j+1);
        min3( c0, c1, c2, rc1);
        min3(ra1, rb1, rc1, y);
        //operateur ligne
        store2(y,i,j,Y) ;
        ra1 = rb1 ;
        rb1 = rc1 ;
        
        //reduction de la colonne c
        c0 = c1;
        c1 = c2;
        c2 = load2(X,i+2,j+1);
        min3( c0, c1, c2, rc2);
        min3(ra2, rb2, rc2, y);
        //operateur ligne
        store2(y,i+1,j,Y) ;
        ra2 = rb2 ;
        rb2 = rc2 ;
    }
}
// --------------------------------------------------------------------------------
void line_min3_ui8matrix_ilu3_elu2_red(uint8 **X, int i, int j0, int j1, uint8 **Y)
// --------------------------------------------------------------------------------
{
}
// ---------------------------------------------------------------------------------------
void line_min3_ui8matrix_ilu3_elu2_red_factor(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ---------------------------------------------------------------------------------------
{
    //par colonne !
    uint8 a0, b0, c0, d0, e0;
    uint8 a1, b1, c1, d1, e1;
    uint8 a2, b2, c2, d2, e2;
    uint8 ra1, rb1, rc1, rd1, re1;
    uint8 ra2, rb2, rc2, rd2, re2;
    uint8 y0, y1, y2;
    // pour i
    a0 = load2(X,i-1,j0-1);      b0 = load2(X,i-1, j0 );
    a1 = load2(X, i ,j0-1);      b1 = load2(X, i , j0 );
    a2 = load2(X,i+1,j0-1);      b2 = load2(X,i+1, j0 );
    min3( a0, a1, a2, ra1);      min3( b0, b1, b2, rb1);
    // pour i+1
    a0 = a1;                    b0 = b1;
    a1 = a2;                    b1 = b2;
    a2 = load2(X,i+2,j0-1);     b2 = load2(X,i+2, j0 );
    min3( a0, a1, a2, ra2);     min3( b0, b1, b2, rb2);
    
    int r = (j1-j0+1)%3;

    for (int j=j0; j<j1-r ; j+=3) {
        //reduction de la colonne c
        c0 = load2(X,i-1,j+1);          d0 = load2(X,i-1,j+2);          e0 = load2(X,i-1,j+3);
        c1 = load2(X, i ,j+1);          d1 = load2(X, i ,j+2);          e1 = load2(X, i ,j+3);
        c2 = load2(X,i+1,j+1);          d2 = load2(X,i+1,j+2);          e2 = load2(X,i+1,j+3);
        min3( c0, c1, c2,rc1);          min3( d0, d1, d2,rd1);          min3( e0, e1, e2,re1);

        min3(ra1,rb1,rc1, y0);
        store2(y0,i,j,Y) ;                      
        ra1 = rb1;
        rb1 = rc1;
        rc1 = rd1;

        min3(ra1,rb1,rc1, y1);
        store2(y1,i,j+1,Y) ;
        ra1 = rb1;
        rb1 = rc1;
        rc1 = re1;

        min3(ra1,rb1,rc1, y2);
        store2(y2,i,j+2,Y) ;
        ra1 = rb1;
        rb1 = rc1;



        //reduction de la colonne c
        c0 = c1;                  d0 = d1;                e0 = e1;
        c1 = c2;                  d1 = d2;                e1 = e2;
        c2 = load2(X,i+2,j+1);    d2 = load2(X,i+2,j+2);  e2 = load2(X,i+2,j+3);
        min3( c0, c1, c2,rc2);    min3( d0, d1, d2,rd2);  min3( e0, e1, e2,re2);

        min3(ra2,rb2,rc2, y0);
        store2(y0,i+1,j,Y) ;                      
        ra2 = rb2;
        rb2 = rc2;
        rc2 = rd2;

        min3(ra2,rb2,rc2, y1);
        store2(y1,i+1,j+1,Y) ;
        ra2 = rb2;
        rb2 = rc2;
        rc2 = re2;
        
        min3(ra2,rb2,rc2, y2);
        store2(y2,i+1,j+2,Y) ;
        ra2 = rb2;
        rb2 = rc2;
        
    }
    if(r==1){
        //reduction de la colonne c
        c0 = load2(X,i-1,j1+1);
        c1 = load2(X, i ,j1+1);
        c2 = load2(X,i+1,j1+1);
        min3( c0, c1, c2,rc1);
        min3(ra1,rb1,rc1, y0);
        //operateur ligne
        store2(y0,i,j1,Y) ;

        //reduction de la colonne c
        c0 = c1;
        c1 = c2;
        c2 = load2(X,i+2,j1+1);
        min3( c0, c1, c2,rc2);
        min3(ra2,rb2,rc2, y0);
        //operateur ligne
        store2(y0,i+1,j1,Y) ;

    }else{
        if(r==2){
            //reduction [de] la colonne c
            c0 = load2(X,i-1, j1 );         d0 = load2(X,i-1,j1+1);
            c1 = load2(X, i , j1 );         d1 = load2(X, i ,j1+1);
            c2 = load2(X,i+1, j1 );         d2 = load2(X,i+1,j1+1);
            min3( c0, c1, c2, rc1);         min3( d0, d1, d2, rd1);
                     
            min3(ra1,rb1,rc1, y0);
            store2(y0,i,j1-1,Y) ;

            ra1 = rb1;
            rb1 = rc1;
            rc1 = rd1;

            min3(ra1,rb1,rc1, y1);
            store2(y1,i,j1,Y) ;

            c0 = c1;                        d0 = d1;
            c1 = c2;                        d1 = d2;
            c2 = load2(X,i+2, j1 );         d2 = load2(X,i+2,j1+1);
            min3( c0, c1, c2, rc2);         min3( d0, d1, d2, rd2);

            min3(ra2,rb2,rc2, y0);
            store2(y0,i+1,j1-1,Y) ; 

            ra2 = rb2;
            rb2 = rc2;
            rc2 = rd2;
            
            min3(ra2,rb2,rc2, y1);
            store2(y1,i+1,j1,Y);
        }
    }
}
// ----------------------------------------------------------------------------
void min3_ui8matrix_basic(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// ----------------------------------------------------------------------------
{
    for(int i=i0; i<=i1; i++) {
        line_min3_ui8matrix_basic(X, i, j0, j1, Y);
    }

}
// --------------------------------------------------------------------------
void min3_ui8matrix_reg(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// --------------------------------------------------------------------------
{
    for(int i=i0; i<=i1; i++) {
        line_min3_ui8matrix_reg(X, i, j0, j1, Y);
    }
}
// --------------------------------------------------------------------------
void min3_ui8matrix_rot(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// --------------------------------------------------------------------------
{
    for(int i=i0; i<=i1; i++) {
        line_min3_ui8matrix_rot(X, i, j0, j1, Y);
    }
}
// --------------------------------------------------------------------------
void min3_ui8matrix_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// --------------------------------------------------------------------------
{
    for(int i=i0; i<=i1; i++) {
        line_min3_ui8matrix_red(X, i, j0, j1, Y);
    }
}
// ---------------------------------------------------------------------------
void min3_ui8matrix_ilu3(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// ---------------------------------------------------------------------------
{
    for(int i=i0; i<=i1; i++) {
        line_min3_ui8matrix_ilu3(X, i, j0, j1, Y);
    }
}
// -------------------------------------------------------------------------------
void min3_ui8matrix_ilu3_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// -------------------------------------------------------------------------------
{
    for(int i=i0; i<=i1; i++) {
        line_min3_ui8matrix_ilu3_red(X, i, j0, j1, Y);
    }
}
// -------------------------------------------------------------------------------
void min3_ui8matrix_elu2_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// -------------------------------------------------------------------------------
{
    //par colonne !
    uint8 a0, b0, c0;
    uint8 a1, b1, c1;
    uint8 a2, b2, c2;
    uint8 ra, rb, rc;
    uint8 y;

    //deroulage de boucle externe
    int r = (i1-i0+1)%2;
    for(int i=i0; i<i1-r; i+=2){
        line_min3_ui8matrix_elu2_red(X, i ,j0,j1,Y);
    }
    if(r==1){
        //par colonne !
        uint8 a0, b0, c0;
        uint8 a1, b1, c1;
        uint8 a2, b2, c2;
        uint8 ra, rb, rc;
        uint8 y;
        
        int i = i1;

        for (int j=j0; j<=j1 ; j++) {
            //reduction colonne a //reduction colonne b     //reduction colonne c
            a0 = load2(X,i-1,j-1);       b0 = load2(X,i-1, j );       c0 = load2(X,i-1,j+1);
            a1 = load2(X, i ,j-1);       b1 = load2(X, i , j );       c1 = load2(X, i ,j+1);
            a2 = load2(X,i+1,j-1);       b2 = load2(X,i+1, j );       c2 = load2(X,i+1,j+1);
            min3( a0, a1, a2, ra);       min3( b0, b1, b2, rb);       min3( c0, c1, c2, rc);
            min3( ra, rb, rc, y);
            //operateur ligne
            store2(y,i,j,Y) ;
        }
    }    
}
// --------------------------------------------------------------------------------------
void min3_ui8matrix_elu2_red_factor(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// --------------------------------------------------------------------------------------
{
    //deroulage de boucle externe
    int r = (i1-i0+1)%2;
    for(int i=i0; i<i1-r; i+=2){
        line_min3_ui8matrix_elu2_red_factor(X, i ,j0,j1,Y);
    }
    if(r==1){
        int i = i1;

        //par colonne !
        uint8 a0, b0, c0;
        uint8 a1, b1, c1;
        uint8 a2, b2, c2;
        uint8 ra1, rb1, rc1;
        uint8 ra2, rb2, rc2;
        uint8 y;
        // pour i
        a0 = load2(X,i-1,j0-1);      b0 = load2(X,i-1, j0 );
        a1 = load2(X, i ,j0-1);      b1 = load2(X, i , j0 );
        a2 = load2(X,i+1,j0-1);      b2 = load2(X,i+1, j0 );
        min3( a0, a1, a2, ra1);      min3( b0, b1, b2, rb1);
        
        for (int j=j0; j<=j1 ; j++) {
            
            //reduction de la colonne c
            c0 = load2(X,i-1,j+1);
            c1 = load2(X, i ,j+1);
            c2 = load2(X,i+1,j+1);
            min3( c0, c1, c2, rc1);
            min3(ra1, rb1, rc1, y);
            //operateur ligne
            store2(y,i,j,Y) ;
            ra1 = rb1 ;
            rb1 = rc1 ;
        }
    }    
}
// ------------------------------------------------------------------------------------
void min3_ui8matrix_ilu3_elu2_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// ------------------------------------------------------------------------------------
{
}
// -------------------------------------------------------------------------------------------
void min3_ui8matrix_ilu3_elu2_red_factor(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// -------------------------------------------------------------------------------------------
{
    //deroulage de boucle externe
    int r = (i1-i0+1)%2;
    for(int i=i0; i<i1-r; i+=2){
        line_min3_ui8matrix_ilu3_elu2_red_factor(X, i  ,j0,j1,Y);
    }
    if(r==1){
        int i = i1;
        float a0, b0, c0, d0, e0;
        float a1, b1, c1, d1, e1;
        float a2, b2, c2, d2, e2;
        float ra1, rb1, rc1, rd1, re1;
        float y0, y1, y2;
        // pour i
        a0 = load2(X,i-1,j0-1);      b0 = load2(X,i-1, j0 );
        a1 = load2(X, i ,j0-1);      b1 = load2(X, i , j0 );
        a2 = load2(X,i+1,j0-1);      b2 = load2(X,i+1, j0 );
        min3( a0, a1, a2, ra1);      min3( b0, b1, b2, rb1);

        int r1 = (j1-j0+1)%3;

        for (int j=j0; j<j1-r1 ; j+=3) {
            //reduction de la colonne c
            c0 = load2(X,i-1,j+1);          d0 = load2(X,i-1,j+2);          e0 = load2(X,i-1,j+3);
            c1 = load2(X, i ,j+1);          d1 = load2(X, i ,j+2);          e1 = load2(X, i ,j+3);
            c2 = load2(X,i+1,j+1);          d2 = load2(X,i+1,j+2);          e2 = load2(X,i+1,j+3);
            min3( c0, c1, c2,rc1);          min3( d0, d1, d2,rd1);          min3( e0, e1, e2,re1);

            min3(ra1,rb1,rc1,y0);
            store2(y0,i,j,Y) ;

            min3(rb1,rc1,rd1,y1);
            store2(y1,i,j+1,Y) ;


            min3(rc1,rd1,re1,y2);
            store2(y2,i,j+2,Y) ;

            ra1 = rd1;
            rb1 = re1;
        }
        if(r1==1){
            c0 = load2(X,i-1,j1+1);
            c1 = load2(X, i ,j1+1);
            c2 = load2(X,i+1,j1+1);
            min3( c0, c1, c2, rc1);
            min3(ra1,rb1,rc1, y0 );
            //operateur ligne
            store2(y0,i,j1,Y) ;
        }else{
            if(r1==2){
                c0 = load2(X,i-1, j1 );          d0 = load2(X,i-1,j1+1);
                c1 = load2(X, i , j1 );          d1 = load2(X, i ,j1+1);
                c2 = load2(X,i+1, j1 );          d2 = load2(X,i+1,j1+1);
                min3( c0, c1, c2, rc1);          min3( d0, d1, d2, rd1);
                        
                min3(ra1,rb1,rc1, y0 );
                store2(y0,i,j1-1,Y) ;

                ra1 = rb1;
                rb1 = rc1;
                rc1 = rd1;

                min3(ra1,rb1,rc1,y1);
                store2(y1,i,j1,Y) ;
            }
        }
    }
}

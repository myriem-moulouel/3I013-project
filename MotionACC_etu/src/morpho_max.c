/* -------------------- */
/* --- morpho_max.c --- */
/* -------------------- */

/*
 * Copyright (c) 2004 - 2013, Lionel Lacassagne, All rights reserved
 * University of Paris Sud, Laboratoire de Recherche en Informatique
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#include "nrtype.h"
#include "nrdef.h"
#include "nrutil.h"
#include "macro.h"


// ------------------------------------------------------------------------
void line_max3_ui8matrix_basic(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ------------------------------------------------------------------------
{
	uint8 y;
	for(int j=j0;j<=j1;j++){
	    max9(load2(X,i-1,j-1), load2(X,i-1,j+0), load2(X,i-1,j+1),
		     load2(X,i+0,j-1), load2(X,i+0,j+0), load2(X,i+0,j+1),
		     load2(X,i+1,j-1), load2(X,i+1,j+0), load2(X,i+1,j+1), y);
		store2(Y,i,j,y);
	}
}
// ----------------------------------------------------------------------
void line_max3_ui8matrix_reg(uint8 **X, int i, int j0, int j1, uint8 **Y)
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
        max9(x11, x12, x13,
             x21, x22, x23,
             x31, x32, x33, y);
        store2(Y, i, j, y);
    }
}
// ----------------------------------------------------------------------
void line_max3_ui8matrix_rot(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ----------------------------------------------------------------------
{
    //avec scalarisation et rotation
    uint8 x11, x12, x13;
    uint8 x21, x22, x23;
    uint8 x31, x32, x33;
    uint8 y;

    x11 = load2(X,i-1,j0-1); x12 = load2(X,i-1,j0);
    x21 = load2(X, i ,j0-1); x22 = load2(X, i ,j0);
    x31 = load2(X,i+1,j0-1); x32 = load2(X,i+1,j0);

    for(int j=j0; j<=j1; j++){
        x13 = load2(X,i-1,j+1);
        x23 = load2(X, i ,j+1);
        x33 = load2(X,i+1,j+1);

        max9(x11, x12, x13,
             x21, x22, x23,
             x31, x32, x33, y);
        store2(Y,i,j,y);
        
        x11 = x12; x12 = x13; 
        x21 = x22; x22 = x23; 
        x31 = x32; x32 = x33; 
        
    }
}
// ----------------------------------------------------------------------
void line_max3_ui8matrix_red(uint8 **X, int i, int j0, int j1, uint8 **Y)
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
    max3(a0,a1,a2,ra);
    //reduction de la colonne b
    b0 = load2(X,i-1, j0 );
    b1 = load2(X, i , j0 );
    b2 = load2(X,i+1, j0 );
    max3(b0,b1,b2,rb);
    
    for (int j=j0; j<=j1 ; j++) {
        //reduction de la colonne c
        c0 = load2(X,i-1,j+1);
        c1 = load2(X, i ,j+1);
        c2 = load2(X,i+1,j+1);
        max3(c0,c1,c2,rc);
        max3(ra,rb,rc,y);
        //operateur ligne
        store2(Y,i,j,y) ;
        ra = rb ;
        rb = rc ;
    }
}
// -----------------------------------------------------------------------
void line_max3_ui8matrix_ilu3(uint8 **X, int i, int j0, int j1, uint8 **Y)
// -----------------------------------------------------------------------
{
    //avec scalarisation et rotation
    uint8 *x1 = load1(X,i-1); uint8 x11, x12, x13;
    uint8 *x2 = load1(X,  i); uint8 x21, x22, x23;
    uint8 *x3 = load1(X,i+1); uint8 x31, x32, x33;
    uint8 y;
    x11 = load1(x1,j0-1); x12 = load1(x1,j0);
    x21 = load1(x2,j0-1); x22 = load1(x2,j0);
    x31 = load1(x3,j0-1); x32 = load1(x3,j0);

    int r = (j1-j0+1)%3;
    for(int j=j0; j<j1-r; j+=3){
        x13 = load1(x1,j+1);
        x23 = load1(x2,j+1);
        x33 = load1(x3,j+1);
        max9(x11, x12, x13,
             x21, x22, x23,
             x31, x32, x33, y);
        store2(Y,i,j,y);

        x11 = x12; x12 = x13; x13 = load1(x1,j+2);
        x21 = x22; x22 = x23; x23 = load1(x2,j+2);
        x31 = x32; x32 = x33; x33 = load1(x3,j+2);
        max9(x11, x12, x13,
             x21, x22, x23,
             x31, x32, x33, y);
        store2(Y,i,j+1,y);

        x11 = x12; x12 = x13; x13 = load1(x1,j+3);
        x21 = x22; x22 = x23; x23 = load1(x2,j+3);
        x31 = x32; x32 = x33; x33 = load1(x3,j+3);
        max9(x11, x12, x13,
             x21, x22, x23,
             x31, x32, x33, y);
        store2(Y,i,j+2,y);

        x11 = x12; x12 = x13;
        x21 = x22; x22 = x23;
        x31 = x32; x32 = x33;
    }
    
    if(r==1){
        x13 = load1(x1,j1+1);
        x23 = load1(x2,j1+1);
        x33 = load1(x3,j1+1);
        max9(x11, x12, x13,
             x21, x22, x23,
             x31, x32, x33, y);
        store2(Y,i,j1,y);
    }else{
        if(r==2){
            x13 = load1(x1, j1 );
            x23 = load1(x2, j1 );
            x33 = load1(x3, j1 );
            max9(x11, x12, x13,
                 x21, x22, x23,
                 x31, x32, x33, y);
            store2(Y,i,j1-1,y);

            x11 = x12; x12 = x13; x13 = load1(x1,j1+1);
            x21 = x22; x22 = x23; x23 = load1(x2,j1+1);
            x31 = x32; x32 = x33; x33 = load1(x3,j1+1);
            max9(x11, x12, x13,
                 x21, x22, x23,
                 x31, x32, x33, y);
            store2(Y,i,j1,y);
        }
    }
}
// ---------------------------------------------------------------------------
void line_max3_ui8matrix_ilu3_red(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ---------------------------------------------------------------------------
{
    //par colonne !
    uint8 a0, b0, c0;
    uint8 a1, b1, c1;
    uint8 a2, b2, c2;
    uint8 ra, rb, rc;
    uint8 y;
    
    a0 = load2(X,i-1,j0-1);       b0 = load2(X,i-1, j0 );
    a1 = load2(X, i ,j0-1);       b1 = load2(X, i , j0 );
    a2 = load2(X,i+1,j0-1);       b2 = load2(X,i+1, j0 );
    max3( a0, a1, a2, ra );       max3( b0, b1, b2, rb );
    
    int r = (j1-j0+1)%3;

    for (int j=j0; j<j1-r ; j+=3) {

        //reduction colonne a //reduction colonne b     //reduction colonne c
        c0 = load2(X,i-1,j+1);
        c1 = load2(X, i ,j+1);
        c2 = load2(X,i+1,j+1);
        max3( c0, c1, c2, rc);

        max3( ra, rb, rc, y);
        //operateur ligne
        store2(Y,i,j,y) ;
        
        //reduction colonne a //reduction colonne b     //reduction colonne c
        ra = rb;                rb = rc;
        c0 = load2(X,i-1,j+2);
        c1 = load2(X, i ,j+2);
        c2 = load2(X,i+1,j+2);
        max3( c0, c1, c2, rc);

        max3( ra, rb, rc, y);
        //operateur ligne
        store2(Y,i,j+1,y) ;
        
        //reduction colonne a //reduction colonne b     //reduction colonne c
        ra = rb;                rb = rc;
        c0 = load2(X,i-1,j+3);
        c1 = load2(X, i ,j+3);
        c2 = load2(X,i+1,j+3);
        max3( c0, c1, c2, rc);

        max3( ra, rb, rc, y);
        //operateur ligne
        store2(Y,i,j+2,y) ;

        ra = rb;                rb = rc;
    }
    if(r==1){
        //reduction colonne a //reduction colonne b     //reduction colonne c
        c0 = load2(X,i-1,j1+1);
        c1 = load2(X, i ,j1+1);
        c2 = load2(X,i+1,j1+1);
        max3( c0, c1, c2, rc);

        max3( ra, rb, rc, y);
        //operateur ligne
        store2(Y,i,j1,y) ;
    }else{
        if(r==2){
            //reduction colonne a //reduction colonne b     //reduction colonne c
            c0 = load2(X,i-1, j1 );
            c1 = load2(X, i , j1 );
            c2 = load2(X,i+1, j1 );
            max3( c0, c1, c2, rc);

            max3( ra, rb, rc, y);
            //operateur ligne
            store2(Y,i,j1-1,y) ;

            //reduction colonne a //reduction colonne b     //reduction colonne c
            ra = rb;                rb = rc;
            c0 = load2(X,i-1,j1+1);
            c1 = load2(X, i ,j1+1);
            c2 = load2(X,i+1,j1+1);
            max3( c0, c1, c2, rc);

            max3( ra, rb, rc, y);
            //operateur ligne
            store2(Y,i,j1,y) ;
        }
    }
}
// ---------------------------------------------------------------------------
void line_max3_ui8matrix_elu2_red(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ---------------------------------------------------------------------------
{
    //version lu3_red avec rotation de registres
    //par colonne !
    uint8 a0, b0, c0;
    uint8 a1, b1, c1;
    uint8 a2, b2, c2;
    uint8 a3, b3, c3;
    uint8 ra1, rb1, rc1;
    uint8 ra2, rb2, rc2;
    uint8 y;
    // pour i
    a0 = load2(X,i-1,j0-1);      b0 = load2(X,i-1, j0 );
    a1 = load2(X, i ,j0-1);      b1 = load2(X, i , j0 );
    a2 = load2(X,i+1,j0-1);      b2 = load2(X,i+1, j0 );
    max3( a0, a1, a2, ra1);      max3( b0, b1, b2, rb1);
    // pour i+1
    a3 = load2(X,i+2,j0-1);      b3 = load2(X,i+2, j0 );
    max3( a1, a2, a3, ra2);      max3( b1, b2, b3, rb2);

    for (int j=j0; j<=j1 ; j++) {
        
        //reduction de la colonne c
        c0 = load2(X,i-1,j+1);
        c1 = load2(X, i ,j+1);
        c2 = load2(X,i+1,j+1);
        max3( c0, c1, c2, rc1);
        max3(ra1, rb1, rc1, y);
        //operateur ligne
        store2(Y,i,j,y) ;
        ra1 = rb1 ;
        rb1 = rc1 ;
        
        //reduction de la colonne c
        c3 = load2(X,i+2,j+1);
        max3( c1, c2, c3, rc2);
        max3(ra2, rb2, rc2, y);
        //operateur ligne
        store2(Y,i+1,j,y) ;
        ra2 = rb2 ;
        rb2 = rc2 ;
    }
}
// ----------------------------------------------------------------------------------
void line_max3_ui8matrix_elu2_red_factor(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ----------------------------------------------------------------------------------
{
    //par colonne !
    uint8 a0, b0, c0;
    uint8 a1, b1, c1;
    uint8 a2, b2, c2;
    uint8 a3, b3, c3;
    uint8 ra0,rb0,rc0;
    uint8 ra1,rb1,rc1;
    uint8 factor_a,factor_b,factor_c;
    uint8 y0, y1;
    //2 * max3 -> 4 comparison max //avec red factor 3*max ---> 3 comparaison// on gagne 1 comparaison par colonne
    a0       = load2(X,i-1,j0-1);          b0 = load2(X,i-1,j0);
    a1       = load2(X,i+0,j0-1);          b1 = load2(X,i+0,j0);
    a2       = load2(X,i+1,j0-1);          b2 = load2(X,i+1,j0);
    a3       = load2(X,i+2,j0-1);          b3 = load2(X,i+2,j0);
    //-----------------------------------------------------------
    factor_a = max2(a1, a2);              factor_b = max2(b1,b2);
    //----------------------------------------------------------
    ra0      = max2(a0, factor_a);        rb0 = max2(b0, factor_b);
    ra1      = max2(a3, factor_a);        rb1 = max2(b3, factor_b);
    for (int j=j0; j<=j1 ; j++) {
        //reduction colonne a //reduction colonne b     //reduction colonne c
        c0  = load2(X,i-1,j+1);
        c1  = load2(X,i+0,j+1);
        c2  = load2(X,i+1,j+1);
        c3  = load2(X,i+2,j+1);
        //---------------------------
        factor_c = max2(c1,c2);
        //---------------------------
        rc0      = max2(c0, factor_c);
        rc1      = max2(c3, factor_c);
        //--------------------------
        max3( ra0, rb0, rc0, y0);
        max3( ra1, rb1, rc1, y1);
        //-------------------------
        store2(Y, i ,j,y0);
        store2(Y,i+1,j,y1);

        ra0 = rb0;                  rb0 = rc0;
        ra1 = rb1;                  rb1 = rc1;
    }
}
// --------------------------------------------------------------------------------
void line_max3_ui8matrix_ilu3_elu2_red(uint8 **X, int i, int j0, int j1, uint8 **Y)
// --------------------------------------------------------------------------------
{
    //version lu3_red avec rotation de registres
    //par colonne !
    uint8 a0, b0, c0, d0, e0;
    uint8 a1, b1, c1, d1, e1;
    uint8 a2, b2, c2, d2, e2;
    uint8 a3, b3, c3, d3, e3;
    uint8 ra1, rb1, rc1, rd1, re1;
    uint8 ra2, rb2, rc2, rd2, re2;
    uint8 y0, y1, y2;

    // pour i
    a0 = load2(X,i-1,j0-1);      b0 = load2(X,i-1, j0 );
    a1 = load2(X, i ,j0-1);      b1 = load2(X, i , j0 );
    a2 = load2(X,i+1,j0-1);      b2 = load2(X,i+1, j0 );
    max3( a0, a1, a2, ra1);      max3( b0, b1, b2, rb1);
    // pour i+1
    a3 = load2(X,i+2,j0-1);      b3 = load2(X,i+2, j0 );
    max3( a1, a2, a3, ra2);      max3( b1, b2, b3, rb2);

    int r = (j1-j0+1)%3;

    for (int j=j0; j<j1-r; j+=3) {
        
        //reduction de la colonne c
        c0 = load2(X,i-1,j+1);      d0 = load2(X,i-1,j+2);      e0 = load2(X,i-1,j+3);
        c1 = load2(X, i ,j+1);      d1 = load2(X, i ,j+2);      e1 = load2(X, i ,j+3);
        c2 = load2(X,i+1,j+1);      d2 = load2(X,i+1,j+2);      e2 = load2(X,i+1,j+3);
        max3(c0, c1, c2, rc1);      max3(d0, d1, d2, rd1);      max3(e0, e1, e2, re1);
        max3(ra1,rb1,rc1, y0);      max3(rb1,rc1,rd1, y1);      max3(rc1,rd1,re1, y2);
        //operateur ligne
        store2(Y,i,j,y0) ;          store2(Y,i,j+1,y1) ;           store2(Y,i,j+2,y2) ;
        ra1 = rd1 ;                 rb1 = re1 ;
        
        //reduction de la colonne c
        c3 = load2(X,i+2,j+1);      d3 = load2(X,i+2,j+2);      e3 = load2(X,i+2,j+3);
        max3(c1, c2, c3, rc2);      max3(d1, d2, d3, rd2);      max3(e1, e2, e3, re2);
        max3(ra2,rb2,rc2, y0);      max3(rb2,rc2,rd2, y1);      max3(rc2,rd2,re2, y2);
        //operateur ligne
        store2(Y,i+1,j,y0) ;        store2(Y,i+1,j+1,y1) ;      store2(Y,i+1,j+2,y2) ;
        ra2 = rd2 ;                 rb2 = re2 ;
        
    }
    if(r==1){
        //reduction de la colonne c
        c0 = load2(X,i-1,j1+1);
        c1 = load2(X, i ,j1+1);
        c2 = load2(X,i+1,j1+1);
        max3(c0, c1, c2, rc1);
        max3(ra1,rb1,rc1, y0);
        //operateur ligne
        store2(Y,i,j1,y0) ;
        
        //reduction de la colonne c
        c3 = load2(X,i+2,j1+1);
        max3( c1, c2, c3, rc2);
        max3( ra2,rb2,rc2, y0);
        //operateur ligne
        store2(Y,i+1,j1,y0) ;
    }else{
        if(r==2){
            //reduction de la colonne c
            c0 = load2(X,i-1, j1 );     d0 = load2(X,i-1,j1+1);
            c1 = load2(X, i , j1 );     d1 = load2(X, i ,j1+1);
            c2 = load2(X,i+1, j1 );     d2 = load2(X,i+1,j1+1);
            max3( c0, c1, c2, rc1);     max3( d0, d1, d2, rd1);
            max3( ra1,rb1,rc1, y0);     max3( rb1,rc1,rd1, y1);
            //operateur ligne
            store2(Y,i,j1-1,y0) ;         store2(Y,i, j1 ,y1) ;
            
            //reduction de la colonne c
            c3 = load2(X,i+2, j1 );     d3 = load2(X,i+2,j1+1);
            max3( c1, c2, c3, rc2);     max3( d1, d2, d3, rd2);
            max3( ra2,rb2,rc2, y0);     max3( rb2,rc2,rd2, y1);
            //operateur ligne
            store2(Y,i+1,j1-1,y0) ;       store2(Y,i+1, j1 ,y1) ;
        }
    }
}
// ---------------------------------------------------------------------------------------
void line_max3_ui8matrix_ilu3_elu2_red_factor(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ---------------------------------------------------------------------------------------
{
	//par colonne !
    uint8 a0, b0, c0, d0, e0;
    uint8 a1, b1, c1, d1, e1;
    uint8 a2, b2, c2, d2, e2;
    uint8 a3, b3, c3, d3, e3;
    uint8 ra1,rb1,rc1,rd1,re1;
    uint8 ra2,rb2,rc2,rd2,re2;
    uint8 factor_a,factor_b,factor_c,factor_d,factor_e;
    uint8 y0, y1, y2, y3, y4, y5;

    // pour i
    a0 = load2(X,i-1,j0-1);      b0 = load2(X,i-1, j0 );
    a1 = load2(X, i ,j0-1);      b1 = load2(X, i , j0 );
    a2 = load2(X,i+1,j0-1);      b2 = load2(X,i+1, j0 );
    a3 = load2(X,i+2,j0-1);      b3 = load2(X,i+2, j0 );
    //-----------------------------------------------
    factor_a = max2(a1, a2);        factor_b = max2(b1,b2);
    //-----------------------------------------------
    ra1 = max2(a0, factor_a);       rb1 = max2(b0, factor_b);
    ra2 = max2(a3, factor_a);       rb2 = max2(b3, factor_b);
    //-----------------------------------------------
    int r = (j1-j0+1)%3;

    for (int j=j0; j<j1-r ; j+=3) {
        //Load
        c0 = load2(X,i-1,j+1);          d0 = load2(X,i-1,j+2);      e0 = load2(X,i-1,j+3);
        c1 = load2(X, i ,j+1);          d1 = load2(X, i ,j+2);      e1 = load2(X, i ,j+3);
        c2 = load2(X,i+1,j+1);          d2 = load2(X,i+1,j+2);      e2 = load2(X,i+1,j+3);
        c3 = load2(X,i+2,j+1);          d3 = load2(X,i+2,j+2);      e3 = load2(X,i+2,j+3);
        //---------------------------------------------------------------------------
        factor_c = max2(c1, c2);        factor_d = max2(d1,d2);     factor_e = max2(e1,e2);
        //---------------------------------------------------------------------------
        //calc
        rc1 = max2(c0, factor_c);       rd1 = max2(d0, factor_d);   re1 = max2(e0, factor_e);
        rc2 = max2(c3, factor_c);       rd2 = max2(d3, factor_d);   re2 = max2(e3, factor_e);
        //-----------------------------------------------------------------------------------
        max3(ra1, rb1, rc1, y0);        max3(rb1,rc1,rd1, y2);      max3(rc1,rd1,re1, y4);
        max3(ra2,rb2,rc2, y1);          max3(rb2,rc2,rd2, y3);      max3(rc2,rd2,re2, y5);
        //-----------------------------------------------------------------------------------
        //Store
        store2(Y,i,j,y0) ;              store2(Y,i,j+1,y2) ;        store2(Y,i,j+2,y4) ;
        store2(Y,i+1,j,y1) ;            store2(Y,i+1,j+1,y3) ;      store2(Y,i+1,j+2,y5) ;
        //----------------------------------
        //RR
        ra1 = rd1 ;                 rb1 = re1 ;
        ra2 = rd2 ;                 rb2 = re2 ;
    }
    if(r==1){
        c0 = load2(X,i-1,j1+1);
        c1 = load2(X, i ,j1+1);
        c2 = load2(X,i+1,j1+1);
        c3 = load2(X,i+2,j1+1);
        //--------------------
        factor_c = max2(c1,c2);
        //--------------------
        //calc
        rc1 = max2(c0, factor_c);
        rc2 = max2(c3, factor_c);
        //---------------------
        max3(ra1,rb1,rc1, y0);
        max3(ra2,rb2,rc2, y1);
        //---------------------
        //Store
        store2(Y,i,j1,y0) ;
        store2(Y,i+1,j1,y1) ;
    }else{
        if(r==2){
            //Load
            c0 = load2(X,i-1, j1 );      d0 = load2(X,i-1,j1+1);
            c1 = load2(X, i , j1 );      d1 = load2(X, i ,j1+1);
            c2 = load2(X,i+1, j1 );      d2 = load2(X,i+1,j1+1);
            c3 = load2(X,i+2, j1 );      d3 = load2(X,i+2,j1+1);
            //-----------------------------------------------
            //Facteur reduction
            factor_c = max2(c1,c2);      factor_d = max2(d1,d2);
            //-----------------------------------------------
            //Calc
            rc1 = max2(c0, factor_c);      rd1 = max2(d0, factor_d);
            rc2 = max2(c3, factor_c);      rd2 = max2(d3, factor_d);
            //-------------------------------------------------
            max3( ra1,rb1,rc1, y0);        max3( rb1,rc1,rd1, y2);
            max3( ra2,rb2,rc2, y1);        max3( rb2,rc2,rd2, y3);
            //------------------------------------------------
            store2(Y,i,j1-1,y0) ;          store2(Y,i, j1 ,y2);
            store2(Y,i+1,j1-1,y1) ;        store2(Y,i+1, j1 ,y3);
        }
    }
}
// ----------------------------------------------------------------------------
void max3_ui8matrix_basic(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// ----------------------------------------------------------------------------
{
    for(int i=i0; i<=i1; i++) {
        line_max3_ui8matrix_basic(X, i, j0, j1, Y);
    }
}
// --------------------------------------------------------------------------
void max3_ui8matrix_reg(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// --------------------------------------------------------------------------
{
    for(int i=i0; i<=i1; i++) {
        line_max3_ui8matrix_reg(X, i, j0, j1, Y);
    }
}
// --------------------------------------------------------------------------
void max3_ui8matrix_rot(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// --------------------------------------------------------------------------
{
    for(int i=i0; i<=i1; i++) {
        line_max3_ui8matrix_rot(X, i, j0, j1, Y);
    }
}
// --------------------------------------------------------------------------
void max3_ui8matrix_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// --------------------------------------------------------------------------
{
    for(int i=i0; i<=i1; i++) {
        line_max3_ui8matrix_red(X, i, j0, j1, Y);
    }
}
// ---------------------------------------------------------------------------
void max3_ui8matrix_ilu3(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// ---------------------------------------------------------------------------
{
    for(int i=i0; i<=i1; i++) {
        line_max3_ui8matrix_ilu3(X, i, j0, j1, Y);
    }
}
// -------------------------------------------------------------------------------
void max3_ui8matrix_ilu3_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// -------------------------------------------------------------------------------
{
    for(int i=i0; i<=i1; i++) {
        line_max3_ui8matrix_ilu3_red(X, i, j0, j1, Y);
    }
}
// -------------------------------------------------------------------------------
void max3_ui8matrix_elu2_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// -------------------------------------------------------------------------------
{
    //deroulage de boucle externe
    int r = (i1-i0+1)%2;
    for(int i=i0; i<i1-r; i+=2){
        line_max3_ui8matrix_elu2_red(X, i ,j0,j1,Y);
    }
    if(r==1){
        line_max3_ui8matrix_red(X,i1,j0,j1,Y);
    }    
}
// --------------------------------------------------------------------------------------
void max3_ui8matrix_elu2_red_factor(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// --------------------------------------------------------------------------------------
{
    //deroulage de boucle externe
    int r = (i1-i0+1)%2;
    for(int i=i0; i<i1-r; i+=2){
        line_max3_ui8matrix_elu2_red_factor(X, i ,j0,j1,Y);
    }
    if(r==1){
        line_max3_ui8matrix_red(X,i1,j0,j1,Y);
    }
}
// ------------------------------------------------------------------------------------
void max3_ui8matrix_ilu3_elu2_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// ------------------------------------------------------------------------------------
{
    //deroulage de boucle externe
    int r = (i1-i0+1)%2;
    for(int i=i0; i<i1-r; i+=2){
        line_max3_ui8matrix_ilu3_elu2_red(X, i ,j0,j1,Y);
    }
    if(r==1){
        line_max3_ui8matrix_ilu3_red(X,i1,j0,j1,Y);
    }
}
// -------------------------------------------------------------------------------------------
void max3_ui8matrix_ilu3_elu2_red_factor(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// -------------------------------------------------------------------------------------------
{
    //deroulage de boucle externe
    int r = (i1-i0+1)%2;
    for(int i=i0; i<i1-r; i+=2){
        line_max3_ui8matrix_ilu3_elu2_red_factor(X, i ,j0,j1,Y);
    }
    if(r==1){
        line_max3_ui8matrix_ilu3_red(X,i1,j0,j1,Y);
    }
}

/* ---------------- */
/* --- nrutil.h --- */
/* ---------------- */

#ifndef __MACRO_H__
#define __MACRO_H__

//#pragma message("  include  nrutil.h")

#ifdef __cplusplus
#pragma message ("C++")
extern "C" {
#endif

#define load1(X,i) X[i]
#define store1(X,i,valeur) X[i]=valeur

#define load2(X,i,j) X[i][j]

#define store2(Y,i,j,valeur) do{ Y[i][j]=valeur; }while(0)

//==========================================================

#define min2(a,b) (a<b)?a:b

#define min3(a,b,c,min) do{ \
    if(a<b){                \
        if(c<a){            \
            min=c;          \
        }else{              \
            min=a;          \
        }                   \
    }else{                  \
        if(c<b){            \
            min=c;          \
        }else{              \
            min=b;          \
        }                   \
    }                       \
}while(0)                   \

#define min9(a,b,c,d,e,f,g,h,i,min) do{             \
    uint8 m1, m2, m3;                               \
    min3(a,b,c,m1); min3(d,e,f,m2); min3(g,h,i,m3); \
    min3(m1,m2,m3,min);                             \
}while(0)                                           \

//==========================================================

#define max2(a,b) (a>b)?a:b

#define max3(a,b,c,max) do{ \
    if(a>b){                \
        if(c>a){            \
            max=c;          \
        }else{              \
            max=a;          \
        }                   \
    }else{                  \
        if(c>b){            \
            max=c;          \
        }else{              \
            max=b;          \
        }                   \
    }                       \
}while(0)                   \

#define max9(a,b,c,d,e,f,g,h,i,max) do{               \
    uint8 m1, m2, m3;                               \
    max3(a,b,c,m1); max3(d,e,f,m2); max3(g,h,i,m3); \
    max3(m1,m2,m3,max);                               \
}while(0)                                           \

#ifdef __cplusplus
}
#endif

#endif // __MACRO_H__

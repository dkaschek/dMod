#include <R.h>
 #include <math.h>
 void P_Boehm_JProteomeRes2014_m6b2lkv5 ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = 207.6*p[0] ;
y[1+i**l] = p[1] ;
y[2+i**l] = 207.6-207.6*p[0] ;
y[3+i**l] = p[2] ;
y[4+i**l] = p[3] ;
y[5+i**l] = p[4] ;
y[6+i**l] = p[5] ;
y[7+i**l] = p[6] ;
y[8+i**l] = pow(10.0,(p[7])) ;
y[9+i**l] = pow(10.0,(p[8])) ;
y[10+i**l] = pow(10.0,(p[9])) ;
y[11+i**l] = pow(10.0,(p[10])) ;
y[12+i**l] = pow(10.0,(p[11])) ;
y[13+i**l] = pow(10.0,(p[12])) ;
y[14+i**l] = p[13] ;
y[15+i**l] = p[14] ;
y[16+i**l] = p[15] ;
y[17+i**l] = pow(10.0,(p[16])) ;
y[18+i**l] = pow(10.0,(p[17])) ;
y[19+i**l] = pow(10.0,(p[18])) ; 
}
}
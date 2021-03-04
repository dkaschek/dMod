#include <R.h>
 #include <math.h>
 void P_Boehm_JProteomeRes2014_deriv_ajnki4d6 ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = 207.6 ;
y[2+i**l] = -207.6 ;
y[21+i**l] = 1.0 ;
y[43+i**l] = 1.0 ;
y[64+i**l] = 1.0 ;
y[85+i**l] = 1.0 ;
y[106+i**l] = 1.0 ;
y[127+i**l] = 1.0 ;
y[148+i**l] = pow(10.0,(p[7]))*log(10.0) ;
y[169+i**l] = pow(10.0,(p[8]))*log(10.0) ;
y[190+i**l] = pow(10.0,(p[9]))*log(10.0) ;
y[211+i**l] = pow(10.0,(p[10]))*log(10.0) ;
y[232+i**l] = pow(10.0,(p[11]))*log(10.0) ;
y[253+i**l] = pow(10.0,(p[12]))*log(10.0) ;
y[274+i**l] = 1.0 ;
y[295+i**l] = 1.0 ;
y[316+i**l] = 1.0 ;
y[337+i**l] = pow(10.0,(p[16]))*log(10.0) ;
y[358+i**l] = pow(10.0,(p[17]))*log(10.0) ;
y[379+i**l] = pow(10.0,(p[18]))*log(10.0) ; 
}
}
#include <R.h>
 #include <math.h>
 void P_Fujita_SciSignal2010_pwu0l07j ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = p[0] ;
y[1+i**l] = p[1] ;
y[2+i**l] = p[2] ;
y[3+i**l] = p[3] ;
y[4+i**l] = p[4] ;
y[5+i**l] = p[5] ;
y[6+i**l] = p[6] ;
y[7+i**l] = p[7] ;
y[8+i**l] = p[8] ;
y[9+i**l] = p[9] ;
y[10+i**l] = p[10] ;
y[11+i**l] = p[11] ;
y[12+i**l] = pow(10.0,(p[12])) ;
y[13+i**l] = pow(10.0,(p[13])) ;
y[14+i**l] = pow(10.0,(p[14])) ;
y[15+i**l] = pow(10.0,(p[15])) ;
y[16+i**l] = pow(10.0,(p[16])) ;
y[17+i**l] = pow(10.0,(p[17])) ;
y[18+i**l] = pow(10.0,(p[18])) ;
y[19+i**l] = pow(10.0,(p[19])) ;
y[20+i**l] = pow(10.0,(p[20])) ;
y[21+i**l] = pow(10.0,(p[21])) ;
y[22+i**l] = pow(10.0,(p[22])) ;
y[23+i**l] = pow(10.0,(p[23])) ;
y[24+i**l] = pow(10.0,(p[24])) ;
y[25+i**l] = p[25] ;
y[26+i**l] = p[26] ;
y[27+i**l] = p[27] ;
y[28+i**l] = p[28] ; 
}
}
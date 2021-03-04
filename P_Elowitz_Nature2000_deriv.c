#include <R.h>
 #include <math.h>
 void P_Elowitz_Nature2000_deriv_w0iwowb9 ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = pow(10.0,(p[0]))*log(10.0) ;
y[23+i**l] = pow(10.0,(p[1]))*log(10.0) ;
y[46+i**l] = pow(10.0,(p[2]))*log(10.0) ;
y[69+i**l] = pow(10.0,(p[3]))*log(10.0) ;
y[92+i**l] = pow(10.0,(p[4]))*log(10.0) ;
y[115+i**l] = pow(10.0,(p[5]))*log(10.0) ;
y[138+i**l] = pow(10.0,(p[6]))*log(10.0) ;
y[161+i**l] = pow(10.0,(p[7]))*log(10.0) ;
y[184+i**l] = pow(10.0,(p[8]))*log(10.0) ;
y[207+i**l] = pow(10.0,(p[9]))*log(10.0) ;
y[230+i**l] = pow(10.0,(p[10]))*log(10.0) ;
y[253+i**l] = pow(10.0,(p[11]))*log(10.0) ;
y[276+i**l] = pow(10.0,(p[12]))*log(10.0) ;
y[299+i**l] = pow(10.0,(p[13]))*log(10.0) ;
y[322+i**l] = pow(10.0,(p[14]))*log(10.0) ;
y[345+i**l] = pow(10.0,(p[15]))*log(10.0) ;
y[368+i**l] = pow(10.0,(p[16]))*log(10.0) ;
y[391+i**l] = pow(10.0,(p[17]))*log(10.0) ;
y[414+i**l] = 1.0 ;
y[437+i**l] = pow(10.0,(p[19]))*log(10.0) ;
y[460+i**l] = pow(10.0,(p[20]))*log(10.0) ;
y[483+i**l] = pow(10.0,(p[21]))*log(10.0) ; 
}
}
/** Code auto-generated by cOde 1.1.0 **/
#include <R.h> 
 #include <math.h> 

static double parms[13];
static double forc[0];
static double cons[0];
static double range[2];

#define nGridpoints 2 
#define nSplines 0 
#define precision 1e-05 

#define kb1 parms[0] 
 #define kc1 parms[1] 
 #define kc2 parms[2] 
 #define kb2 parms[3] 
 #define kc4 parms[4] 
 #define k5 parms[5] 
 #define y0_0 parms[6] 
 #define y1_0 parms[7] 
 #define y2_0 parms[8] 
 #define y3_0 parms[9] 
 #define y4_0 parms[10] 
 #define y5_0 parms[11] 
 #define y6_0 parms[12] 
#define tmin range[0]
#define tmax range[1]


void ccd4_initmod(void (* odeparms)(int *, double *)) {
	 int N=13;
	 odeparms(&N, parms);
}

void ccd4_initforc(void (* odeforcs)(int *, double *)) {
	 int N=0;
	 odeforcs(&N, forc);
}

/** Derivatives (ODE system) **/
void ccd4_derivs (int *n, double *t, double *y, double *ydot, double *RPAR, int *IPAR) {

	 double time = *t;

	 ydot[0] = -1.0*(kb1*y[0]);
 	 ydot[1] = -1.0*(kc1*y[1])-1.0*(kc2*y[1]);
 	 ydot[2] = 1.0*(kb1*y[0])-1.0*(kb2*y[2])+1.0*(kc1*y[1]);
 	 ydot[3] = 1.0*(kb1*y[0])+1.0*(kb2*y[2])+1.0*(kc2*y[1]);
 	 ydot[4] = 1.0*(kc2*y[1])-1.0*(kc4*y[4])+1.0*(k5*y[6]);
 	 ydot[5] = 1.0*(kc1*y[1])+1.0*(kc4*y[4])+1.0*(k5*y[6]);
 	 ydot[6] = -1.0*(k5*y[6]);

}

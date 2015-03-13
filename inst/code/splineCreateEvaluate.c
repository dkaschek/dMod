static NUBspline_1d_d* dataSpline[nSplines];


/** spline creation method **/
void createSplines(double *points, double *values)
{
	int i;
	int j;
	
	BCtype_d xBC = {NATURAL, NATURAL , 0.,0.};
	
	double myp[nSplines][nGridpoints];
	double myv[nSplines][nGridpoints];
	
	for(i = 0; i < nSplines; i++) {
		
		for(j=0; j < nGridpoints; j++) {
		
			myp[i][j] = points[i*nGridpoints + j];
			myv[i][j] = values[i*nGridpoints + j];
			
		}

		NUgrid* mygrid = create_general_grid(myp[i], nGridpoints);
		dataSpline[i] = create_NUBspline_1d_d(mygrid, xBC, myv[i]);
		
	}
  
  range[0] = points[0];
  range[1] = points[nGridpoints-1];

}





/** sample function for spline evaluation **/
void evaluateSplines(double* time, double* valueOut, double* gradOut) {
 	
	int i;
	double t = *time;
  double dt = 0.;
  if(t >= tmax) {
    dt = t-tmax+precision;
    t = tmax - precision;
  }
  if(t <= tmin) {
    dt = t-tmin-precision;
    t = tmin+precision; 
  }
  
	
	for(i=0; i<nSplines; i++) {
		
		double value = 0.;
		double grad = 0.;
	
		eval_NUBspline_1d_d_vg(dataSpline[i], t, &value, &grad);

		gradOut[i] = grad;
    valueOut[i] = value + dt*grad;
		
	}

}



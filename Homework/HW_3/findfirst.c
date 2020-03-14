/* C mex file to find index of first or last non-zero in an array */
/* As requested by John D'Errico, */
/* Program written by Richard C.A. Hindmarsh, 29/10/2002 */

/*
%%%%% © Natural Environment Research Council, (2002)
%%%%% Written by R.C.A. Hindmarsh.
*/

#include <math.h>
#include "mex.h"
#include <string.h>

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     { 
    double *a, *b, *c;
    char *buf;
    int m, n; 
    int i, iopt;
    long long ishibbs;
    
    /* Check for proper number of arguments */
    
    if (nrhs > 2) { 
   	mexErrMsgTxt("Maximum of two input arguments allowed"); 
    } else if (!nrhs ) { 
   	mexErrMsgTxt("Minimum of one input arguments required"); 
    } else if (nlhs > 2) {
	   mexErrMsgTxt("Maximum of two output arguments allowed"); 
    } 
    
    /* Find number of elements*/ 
    
    mxGetM(prhs[0]);

    m = mxGetNumberOfElements(prhs[0]);

    
    /* Assign pointers to the various parameters */ 
    a      = mxGetPr(prhs[0]); 

    /* Create a matrix for the return argument */ 
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL); 

    /* Assign pointer for the output parameter */ 
    b = mxGetPr(plhs[0]);
    b[0] = 0;
    if (nlhs == 2) {
       plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL); 
       c = mxGetPr(plhs[1]);
       c[0] = 0;
    }

    iopt = 0; /* Default to forward search */

    if (nrhs == 2) {
       if (mxGetClassID(prhs[1]) != mxCHAR_CLASS ){
          mexErrMsgTxt("Second input argument must be character string");
       }
       n = mxGetNumberOfElements(prhs[1]);
       buf=mxCalloc(n, sizeof(char));
       mxGetString(prhs[1], buf, n+1);
       if (!strcmp(buf,"first")){
          iopt = 0;
       }
       else if (!strcmp(buf,"last")){
          iopt = 1;
       }
       else if (!strcmp(buf,"why")){
          mexPrintf("Because John D'Errico told me to.");
          iopt = 2;
       }
       else{
          mexErrMsgTxt("Invalid value for second input argument");
       }
    }

   /* Start of search */
    if (!iopt) { /* Forward search */
       for (i = 0; i < m;i++){
          if (a[i] != 0){
             b[0] = a[i];
             if (nlhs == 2) {
                c[0] = i+1; 
             }
             break;
          }
       }
    }
    else if (iopt == 1) { /* Backward search */
       for (i = m-1; i >= 0;i--){
          if (a[i] != 0){
             b[0] = a[i];
             if (nlhs == 2) {
                c[0] = i+1; 
             }
             break;
          }
       }
    }    

    if (!b[0]) {  /* Return empty matrices if input matrix all zeros*/
       plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL); 
       if (nlhs == 2) {
          plhs[1] = mxCreateDoubleMatrix(0, 0, mxREAL); 
       }
}
}

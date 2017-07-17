#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  mwSize i, j, nobs, bins, nfac, *dim;
  double *eik, *mask, *outdata;

  eik = mxGetPr(prhs[0]);
  mask = mxGetPr(prhs[1]);
  bins = *mxGetPr(prhs[2]);
  dim = mxGetDimensions(prhs[0]);
  nobs = dim[0];
  nfac = dim[1];
  /* plhs[0] is first output */
  plhs[0] = mxCreateDoubleMatrix(bins, nfac, mxREAL);
  outdata = mxGetPr(plhs[0]);
  
  for(i=0;i<nobs;i++)
  {
      int w = (int)(mask[i]);
	  for(j=0;j<nfac;j++)
      {
		outdata[w + j*bins-1] += eik[i+j*nobs];
      }
  }
}


/*fftri2gridfast.c Interpolates from 2d triangular mesh to 2d rectangular grid
 *
 * Author: Chloros2 <chloros2@gmx.de>
 * Created: 2018-05-13
 *
 *   [varargout] = fftri2gridfast (tridata, X, Y) interpolates data given on
 *   a triangular mesh to a rectangular grid defined by X and Y. The columns
 *   tridata(:,1) and tridata(:,2) must contain the triangular mesh node
 *   coordinates. The following columns must contain the scalar values at
 *   the node points that need to be interpolated.
 *   The return value is the interpolation at the grid points X, Y. Returns
 *   NaN's if an interpolation point is outside the triangle mesh.
 *   Runtime can be considered approx x33 faster than fftri2gri.m.
 *
 *   Octave users compile the mex file with the command:
 *
 *       mkoctfile --mex -Wall fftri2gridfast.c
 *
 *
 *   Hint: We evaluate the PDE solution only on the grid vertices although the
 *   underlying FE space may have a higher order (P2 element, etc.).
 *   Therefore, there is a small loss of accuracy except P1 elements are used.
 *
 * Copyright (C) 2018 Chloros2 <chloros2@gmx.de>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see
 * <https://www.gnu.org/licenses/>.
 *
 */

#include "mex.h"

static inline
double max3(double i, double j, double k){
  return i > j? (i > k? i: k): (j > k? j: k);
}

static inline
double min3(double i, double j, double k){
  return i < j? (i < k? i: k): (j < k? j: k);
}

void
fftri2gridfast(double *M, double *X, double *Y, double **out,
               int nOuts, mwSize nNodes, mwSize nX,  mwSize nY){

  mwSize nTri=nNodes/3;
  double *invA0=(double *)mxMalloc(nTri*sizeof(double));
  double init=mxGetNaN( );
  //mexPrintf("nodes: %i\ntriangles: %i\n",nNodes,nTri);

  double *x, *y;
  /*Access first and second column as a vector */
  x=&M[0];
  y=&M[nNodes];

  /*Calculates the areas of all triangles and stores their reciprocal value */
  mwSize j=0;
  for (mwSize i=0; i<nTri; i++) {
     invA0[i]=1.0/((y[j+1]-y[j+2])*(x[j]-x[j+2])+(x[j+2]-x[j+1])*(y[j]-y[j+2]));
     j=j+3;
  }

  /*For all grid points (nX,nY) of the grid */
  for (mwSize mx=0; mx<nX; mx++){
     for (mwSize my=0; my<nY; my++){
        mwSize ofs = (mx*nY)+my;
        /*Variable number of output matrices */
        for (mwSize ncols=0; ncols<nOuts; ncols++){
          *(out[ncols]+ofs) = init;
        }
        mwSize i=0, j=0;
        /*Performs a quick search through all TRIs with a pre-criteria
          to improve speed */
        bool doSearchTri=true;
        while (doSearchTri && (i<nNodes)){
           /*If the point (X,Y) is outside a square defined by the TRI
             vertices, the point can not be within the TRI */
           bool presel=((X[mx]<=max3(x[i],x[i+1],x[i+2])) &&
                        (X[mx]>=min3(x[i],x[i+1],x[i+2])) &&
                        (Y[my]<=max3(y[i],y[i+1],y[i+2])) &&
                        (Y[my]>=min3(y[i],y[i+1],y[i+2])));
           /*Potential candiate - calculate Barycentric Coordinates */
           if (presel) {
              /* Sub-triangle areas */
              double Aa=((y[i+1]-y[i+2])*(X[mx]-x[i+2])+
                         (x[i+2]-x[i+1])*(Y[my]-y[i+2]))*invA0[j];
              double Ab=((y[i+2]-y[i])*(X[mx]-x[i+2])+
                         (x[i]-x[i+2])*(Y[my]-y[i+2]))*invA0[j];
              double Ac=1.0-Aa-Ab;
              /*If point is inside the triangle */
              if ((Aa>=0) && (Ab>=0) && (Ac>=0)){
                 /*Interpolates and stores to the output matrices */
                 for (mwSize ncols=0; ncols<nOuts; ncols++){
                    mwSize colofs=(2+ncols)*nNodes;
                    *(out[ncols]+ofs)=Aa*M[colofs+i]+
                                      Ab*M[colofs+i+1]+
                                      Ac*M[colofs+i+2];
                 }
                 doSearchTri=false;
              }
           }
           i=i+3;
           j=j+1;
        }
     }
  }
  mxFree(invA0);
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]){

  double *inMatrix0, *inMatrix1, *inMatrix2;
  double *outMatrix[100];
  mwSize mrows0, mrows1, mrows2;
  mwSize ncols0, ncols1, ncols2;

  //mexPrintf("nrhs: %i\nnlhs: %i\n",nrhs,nlhs);

  if (nrhs!=3) {
     mexErrMsgTxt("3 input arguments required");
  }
  if ((nlhs<1) || (nlhs>99)){
     mexErrMsgTxt("Between 1 and up to 99 output arguments");
  }

  inMatrix0 = mxGetPr(prhs[0]); //trigrid with cols (x_tri,y_tri,c1,c2,c3 ...)
  ncols0 = (mwSize)mxGetN(prhs[0]);
  mrows0 = (mwSize)mxGetM(prhs[0]);
  inMatrix1 = mxGetPr(prhs[1]); //rectgrid X
  ncols1 = (mwSize)mxGetN(prhs[1]);
  mrows1 = (mwSize)mxGetM(prhs[1]);
  inMatrix2 = mxGetPr(prhs[2]); //rectgrid Y
  ncols2 = (mwSize)mxGetN(prhs[2]);
  mrows2 = (mwSize)mxGetM(prhs[2]);

  //mexPrintf("ncols0: %i\nncols1: %i\nncols2: %i\n",ncols0,ncols1,ncols2);
  //mexPrintf("mrows0: %i\nmrows1: %i\nmrows2: %i\n",mrows0,mrows1,mrows2);

  mwSize nNodes, nX, nY, tmp;
  if (mrows1>ncols1){
    nX=mrows1;
    tmp=ncols1;
  }
  else{
    nX=ncols1;
    tmp=mrows1;
  }
  if (tmp!=1){
    mexErrMsgTxt("Input 2: must be a vector");
  }
  if (mrows2>ncols2){ 
    nY=mrows2;
    tmp=ncols2;
  }
  else{
    nY=ncols2;
    tmp=mrows2;
  }
  if (tmp!=1){
    mexErrMsgTxt("Input 3: must be a vector");
  }

  nNodes=mrows0;
  if (nNodes%3){
     mexErrMsgTxt("Input 1: rows must be a multiple of 3");
  }
  if (ncols0<3){
     mexErrMsgTxt("Input 1: ncols must be >= 3");
  }
  if (ncols0-2<nlhs){
     mexErrMsgTxt("Input 1: less columns for scalars than output vars");
  }

  //mexPrintf("nNodes: %i\nnX: %i\nnY: %i\nnOuts: %i\n",nNodes, nX, nY, nlhs);

  /*Create output matrices */
  for (mwSize i=0; i<nlhs; i++){
    plhs[i]=mxCreateDoubleMatrix(nY, nX, mxREAL);
    outMatrix[i]=mxGetPr(plhs[i]);
  }

  fftri2gridfast(inMatrix0, inMatrix1, inMatrix2, outMatrix, nlhs, nNodes, nX, nY);
}
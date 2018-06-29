/*plottri2grid.c Interpolates from 2D triangular mesh to 2D rectangular grid
 *
 * Author: Chloros2 <chloros2@gmx.de>
 * Created: 2018-05-13
 *
 *   [u,[v]] = plottri2grid (X, Y, tx, ty, tu[, tv]) interpolates the data [tu,tv]
 *   which is given on a triangular mesh defined by tx, ty onto a rectangular grid defined
 *   by X and Y. tx, ty, tu, tv must have a size of 3xnTriangle. tv, v is optional.
 *   The return value [u,v] is the interpolation at the grid points
 *   X, Y. Returns NaN's if an interpolation point is outside the triangle mesh.
 *   Runtime can be considered approx x33 faster than the Matlab version.
 *
 *   Octave users compile the mex file with the command:
 *
 *       mkoctfile --mex -Wall plottri2grid.c
 *
 *   Windows users compile the mex file with the command:
 *
 *       mex plottri2grid.cpp -largeArrayDims
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

#define min1(a,b) ( ((a)<(b)) ? (a):(b) )
#define max1(a,b) ( ((a)>(b)) ? (a):(b) )
#define min3(a,b,c) ( ((a)<(b)) ? (((a)<(c)) ? (a):(c)) : (((b)<(c)) ? (b):(c)) )
#define max3(a,b,c) ( ((a)>(b)) ? (((a)>(c)) ? (a):(c)) : (((b)>(c)) ? (b):(c)) )

void
plottri2grid(double *X, double *Y, double *x, double *y, double *u, double *v,
             double **out,int nOuts, mwSize nTri, mwSize nX,  mwSize nY){

  double *invA0=(double *)mxMalloc(nTri*sizeof(double));
  double init=mxGetNaN( );

  /*Areas of all triangles */
  mwSize j=0;
  mwSize nNodes=3*nTri;
  for (mwSize i=0; i<nTri; i++) {
     invA0[i]=1.0/((y[j+1]-y[j+2])*(x[j]-x[j+2])+(x[j+2]-x[j+1])*(y[j]-y[j+2]));
     j=j+3;
  }

  /*For all grid points of the grid */
  for (mwSize mx=0; mx<nX; mx++){
     for (mwSize my=0; my<nY; my++){
        mwSize ofs = (mx*nY)+my;
        for (mwSize ncols=0; ncols<nOuts; ncols++){
          *(out[ncols]+ofs) = init;
        }
        mwSize i=0, j=0;
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
                 /*Interpolates */
                 *(out[0]+ofs) = Aa*u[i]+Ab*u[i+1]+Ac*u[i+2];
                 if (nOuts==2){
                    *(out[1]+ofs) = Aa*v[i]+Ab*v[i+1]+Ac*v[i+2];
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

  double *inMatrix[6], *outMatrix[2];
  mwSize mrows[6], ncols[6];

  switch(nrhs) {
    //[u,v] = tri2gridfst (X, Y, tx, ty, tu, tv)
    case 6:
      if (nlhs!=2){
         mexErrMsgTxt("2 output arguments required");
      }
    break;
    //c = tri2gridfst (X, Y, tx, ty, tu) 
    case 5:
      if (nlhs!=1){
         mexErrMsgTxt("1 output arguments required");
      }
      inMatrix[5] = NULL;
    break;
    default: 
      mexErrMsgTxt("wrong number of input arguments");
    break;
  }

  for (int i=0; i<nrhs; i++){
      inMatrix[i] = mxGetPr(prhs[i]);
      ncols[i] = (mwSize)mxGetN(prhs[i]);
      mrows[i] = (mwSize)mxGetM(prhs[i]);
      //mexPrintf("ncols%i:%i mrows%i:%i\n",i,ncols[i],i,mrows[i]);
  }

  if (!((min1(mrows[0],ncols[0])==1) && (min1(mrows[1],ncols[1])==1))){
     mexErrMsgTxt("Input 1, 2: must be a vector");  
  }
  mwSize nX=max1(mrows[0],ncols[0]);
  mwSize nY=max1(mrows[1],ncols[1]);

  if (!((ncols[2]==ncols[3]) && (ncols[2]==ncols[4]))){
      mexErrMsgTxt("Arguments 3,4,5,6 must have same number of columns");
  }
  if (!((mrows[2]==3) && (mrows[3]==3) && (mrows[4]==3))){
      mexErrMsgTxt("Arguments 3,4,5,6 must have 3 rows");
  }

  if (nrhs==6){
     if (!((mrows[5]==3) && (ncols[5]==ncols[2]))){
        mexErrMsgTxt("6th argument incompatible to previous ones");
     }
  }

  mwSize nTri=ncols[2];

  for (mwSize i=0; i<nlhs; i++){
    plhs[i]=mxCreateDoubleMatrix(nY, nX, mxREAL);
    outMatrix[i]=mxGetPr(plhs[i]);
  }

  plottri2grid (inMatrix[0], inMatrix[1], inMatrix[2], inMatrix[3], inMatrix[4],
                inMatrix[5], outMatrix, nlhs, nTri, nX, nY);
}

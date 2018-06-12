/*fftet2gridfast.c Interpolates from 3D tetrahedral mesh to a rectangular grid
 *
 * Author: Chloros2 <chloros2@gmx.de>
 * Created: 2018-05-13
 *
 *   [varargout] = fftet2gridfast(tdata,X,Y,Z) interpolates data given on
 *   a tetrahedral mesh to a rectangular grid defined by X, Y and Z. The first
 *   three columns in tdata must contain the tetrahedral mesh node
 *   coordinates. The following columns must contain the scalar values at
 *   the four tetrahedra node points that need to be interpolated.
 *   The return value is the interpolation at the grid points defined by X, Y
 *   and Z. Returns NaN's if an interpolation point is outside the tetrahedral mesh.
 *   Runtime can be considered approx x440 faster than fftet2grid.m.
 *
 *   Octave users compile the mex file with the command:
 *
 *       mkoctfile --mex -Wall fftet2gridfast.c
 *
 *   Windows users compile the mex file with the command:
 *
 *       mex fftet2gridfast.cpp -largeArrayDims
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


#include <math.h>
#include "mex.h"

#define min1(a,b) ( ((a)<(b)) ? (a):(b) )
#define max1(a,b) ( ((a)>(b)) ? (a):(b) )

static inline
double max4(double i, double j, double k, double l){
  return max1(max1(i,j), max1(k,l));
}

static inline
double min4(double i, double j, double k, double l){
  return min1(min1(i,j), min1(k,l));
}

typedef struct{
  double i,j,k;
}Vector;

static inline
double dotProduct(Vector a, Vector b){
  return a.i*b.i+a.j*b.j+a.k*b.k;
}

static inline
Vector crossProduct(Vector a, Vector b){
  Vector c={a.j*b.k-a.k*b.j, a.k*b.i-a.i*b.k, a.i*b.j-a.j*b.i};
  return c;
}

static inline
Vector minus(Vector a, Vector b){
  Vector c={a.i-b.i, a.j-b.j, a.k-b.k};
  return c;
}

void fftet2gridfast(double *T, double *X, double *Y, double *Z, double **out,
                    int nOuts, mwSize N, mwSize M, mwSize nNodes){

  mwSize nTet=nNodes/3;
  double *invV0=(double *)mxMalloc(nTet*sizeof(double));
  double init=mxGetNaN( );

  double *x, *y, *z, *col;
  /*Access the 4 columns as a vector */
  x=&T[0];
  y=&T[nNodes];
  z=&T[2*nNodes];
  col=&T[3*nNodes];

  /*Volumes of all tetrahedrons */
  mwSize j=0;
  for (mwSize i=0; i<nTet; i++) {
     Vector ap={x[j], y[j], z[j]}, bp={x[j+1], y[j+1], z[j+1]},
            cp={x[j+2], y[j+2], z[j+2]}, dp={x[j+3], y[j+3], z[j+3]};
     invV0[i]=fabs(1.0/(dotProduct(crossProduct(minus(cp,ap),minus(bp,ap)),minus(dp,ap))));
     j=j+4;
  }

  /*For all grid points of the grid */
  for (mwSize mx=0; mx<N; mx++){
     for (mwSize my=0; my<M; my++){
        mwSize ofs=(mx*M)+my;
        /*Variable number of output matrices */
        for (mwSize ncols=0; ncols<nOuts; ncols++){
           *(out[ncols]+ofs) = init;
        }
        mwSize i=0, j=0;
        bool doSearchTet=true;
        while (doSearchTet && (i<nNodes)){
           mwSize idx=my+M*mx;
           /*If the grid point is outside a cube defined by the TET
             vertices, the point can not be within the TET */
           bool presel=((X[idx]<=max4(x[i],x[i+1],x[i+2],x[i+3])) &&
                        (X[idx]>=min4(x[i],x[i+1],x[i+2],x[i+3])) &&
                        (Y[idx]<=max4(y[i],y[i+1],y[i+2],y[i+3])) &&
                        (Y[idx]>=min4(y[i],y[i+1],y[i+2],y[i+3])) &&
                        (Z[idx]<=max4(z[i],z[i+1],z[i+2],z[i+3])) &&
                        (Z[idx]>=min4(z[i],z[i+1],z[i+2],z[i+3])) );
           /*Potential candiate - calculate Barycentric Coordinates */
           if (presel) {
              Vector ap={x[i], y[i], z[i]}, bp={x[i+1], y[i+1], z[i+1]},
                     cp={x[i+2], y[i+2], z[i+2]}, dp={x[i+3], y[i+3], z[i+3]},
                     xp={X[idx],Y[idx],Z[idx]};
              /* Sub-tet volumes */
              double Va=dotProduct(crossProduct(minus(dp,bp), minus(cp,bp)),
                                                minus(xp,bp));
              double Vb=dotProduct(crossProduct(minus(cp,ap), minus(dp,ap)),
                                                minus(xp,ap));
              double Vc=dotProduct(crossProduct(minus(dp,ap), minus(bp,ap)),
                                                minus(xp,ap));
              double Vd=dotProduct(crossProduct(minus(bp,ap), minus(cp,ap)),
                                                minus(xp,ap));
              /*If point is inside the tetrahedron */
              if ((Va>=0) && (Vb>=0) && (Vc>=0) && (Vd>=0)){
                 /*Interpolates and stores to the output matrices */
                 for (mwSize ncols=0; ncols<nOuts; ncols++){
                    mwSize colofs=ncols*nNodes;
                    *(out[ncols]+ofs)=invV0[j]*(Va*col[colofs+i]+Vb*col[colofs+i+1]+
                                                Vc*col[colofs+i+2]+Vd*col[colofs+i+3]);
                 }
                 doSearchTet=false;
              }
           }
           i=i+4;
           j=j+1;
        }
     }
  }
  mxFree(invV0);
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]){

  double *inMatrix0, *inMatrix1, *inMatrix2, *inMatrix3;
  double *outMatrix[100];
  mwSize ncols0, ncols1, ncols2, ncols3;
  mwSize mrows0, mrows1, mrows2, mrows3;

  /*mexPrintf("nrhs: %i\nnlhs: %i\n",nrhs,nlhs);*/
  if (nrhs!=4) {
     mexErrMsgTxt("4 input arguments required");
  }
  if ((nlhs<1) || (nlhs>99)){
     mexErrMsgTxt("Between 1 and up to 99 output arguments");
  }

  inMatrix0 = mxGetPr(prhs[0]); //mesh data first 3 columns + scalar 4th column
  ncols0 = (mwSize)mxGetN(prhs[0]);
  mrows0 = (mwSize)mxGetM(prhs[0]);
  inMatrix1 = mxGetPr(prhs[1]); //X
  ncols1 = (mwSize)mxGetN(prhs[1]);
  mrows1 = (mwSize)mxGetM(prhs[1]);
  inMatrix2 = mxGetPr(prhs[2]); //Y
  ncols2 = (mwSize)mxGetN(prhs[2]);
  mrows2 = (mwSize)mxGetM(prhs[2]);
  inMatrix3 = mxGetPr(prhs[3]); //Z
  ncols3 = (mwSize)mxGetN(prhs[3]);
  mrows3 = (mwSize)mxGetM(prhs[3]);

  if (ncols0<4){
     mexErrMsgTxt("Input 1: ncols must be > 3");
  }
  if (mrows0%4){
     mexErrMsgTxt("Input 1: rows must be a multiple of 4");
  }
  if (!((ncols1==ncols2) && (ncols1==ncols3) &&
        (mrows1==mrows2) && (mrows1==mrows3))){
     mexErrMsgTxt("Grid arguments must have same length and width");
  }

  for (mwSize i=0; i<nlhs; i++){
    plhs[i]=mxCreateDoubleMatrix(mrows1, ncols1, mxREAL);
    outMatrix[i]=mxGetPr(plhs[i]);
  }

 fftet2gridfast(inMatrix0, inMatrix1, inMatrix2, inMatrix3,
                outMatrix, nlhs, ncols1, mrows1, mrows0);
}

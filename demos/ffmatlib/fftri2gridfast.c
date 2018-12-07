/*fftri2gridfast.c Interpolates from 2D triangular mesh to a 2D cartesian or curved grid
 *
 * Author: Chloros2 <chloros2@gmx.de>
 * Created: 2018-11-13
 *
 *   This file is part of the ffmatlib which is hosted at
 *   https://github.com/samplemaker/freefem_matlab_octave_plot
 *
 *   [w1, [w2]] = fftri2gridfast (x, y, tx, ty, tu1, [tu2])
 *
 *   interpolates the real or complex data tu1, tu2 given on a triangular mesh
 *   defined by the two arguments tx, ty onto a meshgrid defined by the
 *   variables x, y. The mesh can be cartesian or curved. The argument
 *   tu1, tu2 must have a size of nTriangle columns and 3 rows. The return
 *   value w1, w2 is the interpolation of tu1, tu2 at the grid points defined
 *   by x, y. The result w1, w2 is real if tu1, tu2 is real or complex if
 *   tu1, tu2 is complex. fftri2gridfast.c uses barycentric coordinates to
 *   interpolate. Returns NaN's if an interpolation point is outside the
 *   triangle mesh. The argument tu2 is optional.
 *   fftri2gridfast.c is the mex implementation of the function
 *   fftri2grid.c.
 *
 *   Octave users on Linux with gcc compile the mex file with the command:
 *
 *       mkoctfile --mex -Wall fftri2gridfast.c
 *
 *   Matlab users on Windows with microsoft visual studio compile the mex
 *   file with the command:
 *
 *       mex fftri2gridfast.cpp -largeArrayDims
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
fftri2meshgrid(double *x, double *y, double *tx, double *ty,
               double *tu1Re, double *tu1Im, double *w1Re, double *w1Im,
               double *tu2Re, double *tu2Im, double *w2Re, double *w2Im,
               int nOuts, mwSize nTri, mwSize nx,  mwSize ny){

  double *invA0=(double *)mxMalloc(nTri*sizeof(double));
  double init=mxGetNaN( );

  /*Areas of all triangles */
  mwSize j=0;
  mwSize nNodes=3*nTri;
  for (mwSize i=0; i<nTri; i++) {
     invA0[i]=1.0/((ty[j+1]-ty[j+2])*(tx[j]-tx[j+2])+(tx[j+2]-tx[j+1])*(ty[j]-ty[j+2]));
     j=j+3;
  }
  /*For all grid points of the meshgrid */
  for (mwSize mx=0; mx<nx; mx++){
     for (mwSize my=0; my<ny; my++){
        mwSize ofs = (mx*ny)+my;
        *(w1Re+ofs) = init;
        if (tu1Im != NULL){
           *(w1Im+ofs) = init;
        }
        if (nOuts == 2){
           *(w2Re+ofs) = init;
           if (tu2Im != NULL){
              *(w2Im+ofs) = init;
           }
        }
        mwSize i=0, j=0;
        bool doSearchTri=true;
        while (doSearchTri && (i<nNodes)){
           /*If the point (X,Y) is outside a square defined by the TRI
             vertices, the point can not be within the TRI */
           bool presel=((x[ofs]<=max3(tx[i],tx[i+1],tx[i+2])) &&
                        (x[ofs]>=min3(tx[i],tx[i+1],tx[i+2])) &&
                        (y[ofs]<=max3(ty[i],ty[i+1],ty[i+2])) &&
                        (y[ofs]>=min3(ty[i],ty[i+1],ty[i+2])));
           /*Potential candiate - calculate Barycentric Coordinates */
           if (presel) {
              /* Sub-triangle areas */
              double Aa=((ty[i+1]-ty[i+2])*(x[ofs]-tx[i+2])+
                         (tx[i+2]-tx[i+1])*(y[ofs]-ty[i+2]))*invA0[j];
              double Ab=((ty[i+2]-ty[i])*(x[ofs]-tx[i+2])+
                         (tx[i]-tx[i+2])*(y[ofs]-ty[i+2]))*invA0[j];
              double Ac=1.0-Aa-Ab;
              /*If point is inside the triangle */
              /*Set a negative threshold due to numerical error in Aa, Ab ...
                if the interpolation point is on the triangle edge */
              if ((Aa>=-1e-13) && (Ab>=-1e-13) && (Ac>=-1e-13)){
                 /*Interpolates real and imaginary part */
                 *(w1Re+ofs) = Aa*tu1Re[i]+Ab*tu1Re[i+1]+Ac*tu1Re[i+2];
                 /*If input data tu is real no imaginary part of w is written */
                 if (tu1Im != NULL){
                    *(w1Im+ofs) = Aa*tu1Im[i]+Ab*tu1Im[i+1]+Ac*tu1Im[i+2];
                 }
                 if (nOuts==2){
                    *(w2Re+ofs) = Aa*tu2Re[i]+Ab*tu2Re[i+1]+Ac*tu2Re[i+2];
                    /*If input data tu is real no imaginary part of w is written */
                    if (tu2Im != NULL){
                       *(w2Im+ofs) = Aa*tu2Im[i]+Ab*tu2Im[i+1]+Ac*tu2Im[i+2];
                    }
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

  double *x, *y, *tx, *ty;
  double *tu1Re, *tu1Im, *w1Re, *w1Im;
  double *tu2Re, *tu2Im, *w2Re, *w2Im;
  mwSize nColX, mRowX, nColY, mRowY;
  mwSize nColTX, mRowTX, nColTY, mRowTY;
  mwSize nColTU1, mRowTU1, nColTU2, mRowTU2;

  tu1Re = NULL; tu1Im = NULL;
  w1Re = NULL; w1Im = NULL;
  tu2Re = NULL; tu2Im = NULL;
  w2Re = NULL; w2Im = NULL;

  switch(nrhs) {
    //[w1,w2] = tri2gridfst (x, y, tx, ty, tu1, tu2)
    case 6:
      if (nlhs!=2){
         mexErrMsgTxt("2 output arguments required");
      }
    break;
    //w = fftri2meshgrid (x, y, tx, ty, tu1)
    case 5:
      if (nlhs != 1){
         mexErrMsgTxt("1 output arguments required");
      }
    break;
    default:
      mexErrMsgTxt("wrong number of arguments");
    break;
  }

  /*x and y meshgrid */
  x = mxGetPr(prhs[0]);
  nColX = (mwSize)mxGetN(prhs[0]);
  mRowX = (mwSize)mxGetM(prhs[0]);

  y = mxGetPr(prhs[1]);
  nColY = (mwSize)mxGetN(prhs[1]);
  mRowY = (mwSize)mxGetM(prhs[1]);

  /*Triangle Mesh: tx */
  tx = mxGetPr(prhs[2]);
  nColTX = (mwSize)mxGetN(prhs[2]);
  mRowTX = (mwSize)mxGetM(prhs[2]);

  /*Triangle Mesh: ty */
  ty = mxGetPr(prhs[3]);
  nColTY = (mwSize)mxGetN(prhs[3]);
  mRowTY = (mwSize)mxGetM(prhs[3]);

  /*Get pointers to real and imaginary parts of tu1 */
  tu1Re = mxGetPr(prhs[4]);
  nColTU1 = (mwSize)mxGetN(prhs[4]);
  mRowTU1 = (mwSize)mxGetM(prhs[4]);
  /*Check if tu1 is complex. If no then output w1 must be real as well */
  bool tu1IsComplex = false;
  if (mxIsComplex(prhs[4])){
     tu1Im = mxGetPi(prhs[4]);
     tu1IsComplex = true;
  }

  if (!((nColX==nColY) && (mRowX==mRowY))){
      mexErrMsgTxt("Arguments 1,2 must have same dimensions");
  }
  if (!((nColTX==nColTY) && (nColTX==nColTU1))){
      mexErrMsgTxt("Arguments 3,4,5 must have same number of columns");
  }
  if (!((mRowTX==3) && (mRowTY==3) && (mRowTU1==3))){
      mexErrMsgTxt("Arguments 3,4,5 must have 3 rows");
  }

  bool tu2IsComplex = false;
  if (nrhs==6){
      /*Get pointers to real and imaginary parts of tu2 */
      tu2Re = mxGetPr(prhs[5]);
      nColTU2 = (mwSize)mxGetN(prhs[5]);
      mRowTU2 = (mwSize)mxGetM(prhs[5]);
     /*Check if tu2 is complex. If no then output w2 must be real as well */
     if (mxIsComplex(prhs[5])){
        tu2Im = mxGetPi(prhs[5]);
        tu2IsComplex = true;
     }
     if (!((mRowTU2==3) && (nColTU1==nColTU2))){
        mexErrMsgTxt("6th argument incompatible to previous ones");
     }
  }

  mwSize nX = nColX;
  mwSize nY = mRowX;
  mwSize nTri = nColTX;
  //mexPrintf("nXcols:%i nYrows:%i nTri:%i\n",nX,nY,nTri);

  /*If the input is complex the output must be complex as well */
  if (tu1IsComplex){
     /*Create two real matrices and set the output pointer to it */
     plhs[0] = mxCreateDoubleMatrix(nY, nX, mxCOMPLEX);
     w1Re = mxGetPr(plhs[0]);
     w1Im = mxGetPi(plhs[0]);
     //mexPrintf("data1 is complex\n");
  }else{
     /*Create one real matrix and set the output pointer to it */
     plhs[0] = mxCreateDoubleMatrix(nY, nX, mxREAL);
     w1Re = mxGetPr(plhs[0]);
     //mexPrintf("data1 is real\n");
  }

  if (nrhs==6){
     /*If the input is complex the output must be complex as well */
     if (tu2IsComplex){
        /*Create two real matrices and set the output pointer to it */
        plhs[1] = mxCreateDoubleMatrix(nY, nX, mxCOMPLEX);
        w2Re = mxGetPr(plhs[1]);
        w2Im = mxGetPi(plhs[1]);
        //mexPrintf("data2 is complex\n");
     }else{
        /*Create one real matrix and set the output pointer to it */
        plhs[1] = mxCreateDoubleMatrix(nY, nX, mxREAL);
        w2Re = mxGetPr(plhs[1]);
        //mexPrintf("data2 is real\n");
     }
  }

  fftri2meshgrid (x, y, tx, ty, tu1Re, tu1Im, w1Re, w1Im, tu2Re, tu2Im,
                  w2Re, w2Im, nlhs, nTri, nX, nY);
}

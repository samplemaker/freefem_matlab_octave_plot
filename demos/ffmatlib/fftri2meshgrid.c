/*fftri2meshgrid.c Interpolates from 2D triangular mesh to a 2D curved grid
 *
 * Author: Chloros2 <chloros2@gmx.de>
 * Created: 2018-11-13
 *
 *   [w] = fftri2meshgrid (x, y, tx, ty, tu) interpolates the real or complex
 *   data tu given on a triangle mesh defined by the two arguments tx, ty onto
 *   a meshgrid defined by the variables x, y. The mesh can be cartesian or
 *   curved. The argument tu must have a size of nTriangle columns and 3 rows.
 *   The return value w is the interpolation of tu at the grid points defined
 *   by x, y and is real if tu is real or complex if tu is complex. fftri2mesgrid
 *   uses barycentric coordinates to interpolate.
 *   Returns NaN's if an interpolation point is outside the triangle mesh.
 *
 *   Octave users compile the mex file with the command:
 *
 *       mkoctfile --mex -Wall fftri2meshgrid.c
 *
 *   Matlab users on Windows compile the mex file with the command:
 *
 *       mex fftri2meshgrid.cpp -largeArrayDims
 *
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
fftri2meshgrid(double *x, double *y, double *tx, double *ty, double *tu, double *tv,
               double *wu, double *wv, mwSize nTri, mwSize nx,  mwSize ny){

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
        *(wu+ofs) = init;
        if (tv != NULL){
           *(wv+ofs) = init;
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
                 *(wu+ofs) = Aa*tu[i]+Ab*tu[i+1]+Ac*tu[i+2];
                 /*If input data tu is real no imaginary part of w is written */
                 if (tv != NULL){
                    *(wv+ofs) = Aa*tv[i]+Ab*tv[i+1]+Ac*tv[i+2];
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

  double *x, *y, *tx, *ty, *tuRe, *tuIm, *wRe, *wIm;
  mwSize nColX, mRowX, nColY, mRowY;
  mwSize nColTX, mRowTX, nColTY, mRowTY, nColTU, mRowTU;

  switch(nrhs) {
    //w = fftri2meshgrid (x, y, tx, ty, tu)
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

  /*Get pointers to real and imaginary parts of tu */
  tuRe = mxGetPr(prhs[4]);
  nColTU = (mwSize)mxGetN(prhs[4]);
  mRowTU = (mwSize)mxGetM(prhs[4]);

  /*Check if tu is complex. If no output w shall be real as well */
  bool tuIsComplex = false;
  if (mxIsComplex(prhs[4])){
     tuIm = mxGetPi(prhs[4]);
     tuIsComplex = true;
  }else{
     tuIm = NULL;
  }

  if (!((nColX==nColY) && (mRowX==mRowY))){
      mexErrMsgTxt("Arguments 1,2 must have same dimensions");
  }
  if (!((nColTX==nColTY) && (nColTX==nColTU))){
      mexErrMsgTxt("Arguments 3,4,5 must have same number of columns");
  }
  if (!((mRowTX==3) && (mRowTY==3) && (mRowTU==3))){
      mexErrMsgTxt("Arguments 3,4,5 must have 3 rows");
  }

  mwSize nX = nColX;
  mwSize nY = mRowX;
  mwSize nTri = nColTX;
  mexPrintf("nXcols:%i nYrows:%i nTri:%i\n",nX,nY,nTri);

  /*If the input is complex the output must be complex as well */
  if (tuIsComplex){
     /*Create two real matrices and set the output pointer to it */
     plhs[0] = mxCreateDoubleMatrix(nY, nX, mxCOMPLEX);
     wRe = mxGetPr(plhs[0]);
     wIm = mxGetPi(plhs[0]);
     mexPrintf("data is complex\n");
  }else{
     /*Create one real matrix and set the output pointer to it */
     plhs[0] = mxCreateDoubleMatrix(nY, nX, mxREAL);
     wRe = mxGetPr(plhs[0]);
     wIm = NULL;
     mexPrintf("data is real\n");
  }

  fftri2meshgrid (x, y, tx, ty, tuRe, tuIm, wRe, wIm, nTri, nX, nY);
}

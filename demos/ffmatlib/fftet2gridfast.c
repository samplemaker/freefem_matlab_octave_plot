/*fftet2gridfast.c Interpolates from 3d tetrahedral mesh to a rectangular grid
 *
 * Author: Chloros2 <chloros2@gmx.de>
 * Created: 2018-05-13
 *
 *   [varargout] = fftet2gridfast(a,b,c,d,u,X,Y,Z) interpolates data given on
 *   a tetrahedral mesh to a rectangular grid defined by X, Y and Z. The column
 *   matrices a,b,c and d must contain the tetrahedral mesh node
 *   coordinates. The matrix u must contain the scalar values at
 *   the four tetrahedra node points that need to be interpolated.
 *   The return value is the interpolation at the grid points defined by X, Y
 *   and Z. Returns NaN's if an interpolation point is outside the tetrahedral mesh.
 *   Runtime can be considered approx x480 faster than fftet2grid.m.
 *
 *   Octave users compile the mex file with the command:
 *
 *       mkoctfile --mex -Wall fftet2gridfast.c
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

void fftet2gridfast(double *a, double *b, double *c, double *d, double *col,
                    double *out, double *X, double *Y, double *Z,
                    mwSize N, mwSize M, mwSize nTet){

  double *invV0=mxMalloc(nTet*sizeof(double));
  double init=mxGetNaN( );

  /*Calculates the volumes of all tetrahedrons and stores their reciprocal value */
  for (mwSize i=0; i<nTet; i++) {
     Vector ap={a[i], a[i+nTet], a[i+2*nTet]}, bp={b[i],b[i+nTet],b[i+2*nTet]},
            cp={c[i], c[i+nTet], c[i+2*nTet]}, dp={d[i],d[i+nTet],d[i+2*nTet]};

     invV0[i]=fabs(1.0/(dotProduct(crossProduct(minus(cp,ap),minus(bp,ap)),minus(dp,ap))));
     //mexPrintf("%f %f %f %f\n",a[i],a[i+nTet],a[i+2*nTet],invV0[i]);
  }

  /*For all grid points (NxM) of the grid */
  for (mwSize mx=0; mx<N; mx++){
     for (mwSize my=0; my<M; my++){
        mwSize ofs=(mx*M)+my;
        *(out+ofs)=init;
        mwSize i=0;
        /*Performs a quick search through all TETs with a pre-criteria
          to improve speed */
        bool doSearchTet=true;
        while (doSearchTet && (i<nTet)){
           mwSize idx=my+M*mx;
           /*If the point (x, y, z) is outside a cube defined by the TET
             vertices, the point can not be within the TET */
           bool presel=((X[idx]<=max4(a[i],b[i],c[i],d[i])) &&
                        (X[idx]>=min4(a[i],b[i],c[i],d[i])) &&
                        (Y[idx]<=max4(a[i+nTet],b[i+nTet],c[i+nTet],d[i+nTet])) &&
                        (Y[idx]>=min4(a[i+nTet],b[i+nTet],c[i+nTet],d[i+nTet])) &&
                        (Z[idx]<=max4(a[i+2*nTet],b[i+2*nTet],c[i+2*nTet],d[i+2*nTet])) &&
                        (Z[idx]>=min4(a[i+2*nTet],b[i+2*nTet],c[i+2*nTet],d[i+2*nTet])) );
           /*Potential candiate - calculate Barycentric Coordinates */
           if (presel) {
              Vector ap={a[i],a[i+nTet],a[i+2*nTet]}, bp={b[i],b[i+nTet],b[i+2*nTet]},
                     cp={c[i],c[i+nTet],c[i+2*nTet]}, dp={d[i],d[i+nTet],d[i+2*nTet]},
                     xp={X[idx],Y[idx],Z[idx]};
              /* Partial volumes */
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
                 /*Interpolate */
                 *(out+ofs)=invV0[i]*(Va*col[i]+Vb*col[i+nTet]+
                                      Vc*col[i+2*nTet]+Vd*col[i+3*nTet]);
                 doSearchTet=false;
              }
           }
           i=i+1;
        }
     }
  }
  mxFree(invV0);
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]){

  double *inMatrix0, *inMatrix1, *inMatrix2, *inMatrix3,
         *inMatrix4, *inMatrix5, *inMatrix6, *inMatrix7;
  double *outMatrix;
  mwSize ncols0, ncols1, ncols2, ncols3, ncols4, ncols5, ncols6, ncols7;
  mwSize mrows0, mrows1, mrows2, mrows3, mrows4, mrows5, mrows6, mrows7;

  //mexPrintf("nrhs: %i\nnlhs: %i\n",nrhs,nlhs);

  inMatrix0 = mxGetPr(prhs[0]); //tet point a(x,y,z)
  ncols0 = mxGetN(prhs[0]);
  mrows0 = mxGetM(prhs[0]);
  inMatrix1 = mxGetPr(prhs[1]); //tet point b(x,y,z)
  ncols1 = mxGetN(prhs[1]);
  mrows1 = mxGetM(prhs[1]);
  inMatrix2 = mxGetPr(prhs[2]); //tet point c(x,y,z)
  ncols2 = mxGetN(prhs[2]);
  mrows2 = mxGetM(prhs[2]);
  inMatrix3 = mxGetPr(prhs[3]); //tet point d(x,y,z)
  ncols3 = mxGetN(prhs[3]);
  mrows3 = mxGetM(prhs[3]);
  inMatrix4 = mxGetPr(prhs[4]); //tet color for (a,b,c,d)
  ncols4 = mxGetN(prhs[4]);
  mrows4 = mxGetM(prhs[4]);

  inMatrix5 = mxGetPr(prhs[5]); //X
  ncols5 = mxGetN(prhs[5]);
  mrows5 = mxGetM(prhs[5]);
  inMatrix6 = mxGetPr(prhs[6]); //Y
  ncols6 = mxGetN(prhs[6]);
  mrows6 = mxGetM(prhs[6]);
  inMatrix7 = mxGetPr(prhs[7]); //Z
  ncols7 = mxGetN(prhs[7]);
  mrows7 = mxGetM(prhs[7]);

/*  mexPrintf("ncols0: %i\nncols1: %i\nncols2: %i\nncols3: %i\nncols4 %i\nncols5: %i\nncols6: %i\nncols7: %i\n",
              ncols0,ncols1,ncols2,ncols3,ncols4,ncols5,ncols6,ncols7);
    mexPrintf("mrows0: %i\nmrows1: %i\nmrows2: %i\nmrows3: %i\nmrows4: %i\nnmrows5: %i\nmrows6: %i\nmrows7: %i\n",
              mrows0,mrows1,mrows2,mrows3,mrows4,mrows5,mrows6,mrows7); */

  if (!( (ncols0==ncols1) && (ncols0==ncols2) && (ncols0==ncols3) &&
         (mrows0==mrows1) && (mrows0==mrows2) && (mrows0==mrows3) && (mrows0==mrows4))){
     mexErrMsgTxt("TET-coordinates / colors have wrong length and width");
  }

  if (!((ncols0==3) && (ncols4==4))){
     mexErrMsgTxt("Coordinates must have 3 and colors 4 columns");
  }

  if (!((ncols5==ncols6) && (ncols5==ncols7) &&
        (mrows5==mrows6) && (mrows5==mrows7))){
     mexErrMsgTxt("Grid must have same length and width");
  }

  plhs[0]=mxCreateDoubleMatrix(mrows5, ncols5, mxREAL);
  outMatrix=mxGetPr(plhs[0]);

  fftet2gridfast(inMatrix0, inMatrix1, inMatrix2, inMatrix3, inMatrix4,
                 outMatrix, inMatrix5, inMatrix6, inMatrix7,
                 ncols5, mrows5, mrows0);
}

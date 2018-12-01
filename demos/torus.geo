//torus.geo Creates the mesh of a torus embedded into a cuboid (GMSH)
//
// Author: Chloros2 <chloros2@gmx.de>
// Created: 2018-11-30
//
// Copyright (C) 2018 Chloros2 <chloros2@gmx.de>
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see
// <https://www.gnu.org/licenses/>.
//

/// DEFINITIONS ///

//Volume of the TORUS
PHYSTORUSVOLUME = 101;
//Volume of CUBOID minus TORUS
PHYSCUBOIDTORUSVOLUME = 102;
//Surface of the TORUS
PHYSTORUSSUFACE = 201;
//Surface of CUBOID minus TORUS
PHYSCUBOIDTORUSSURFACE = 102;
//Outer surface of the CUBOID (e.g. for boundary conditions)
PHYSCUBOIDSURFACE = 103;

//characteristic mesh length
clen = 0.3;
tlen = 0.15;

TSURFACE[] = {};
TVOLUME[] = {};
CSURFACE[] = {};

/// CUBOID ///

//Cuboid dimensions
CXY = 1.3;
P1 = newp; Point(P1) = {-CXY, CXY, 0, clen};
P2 = newp; Point(P2) = {CXY, CXY, 0, clen};
P3 = newp; Point(P3) = {CXY, -CXY, 0, clen};
P4 = newp; Point(P4) = {-CXY, -CXY, 0, clen};
L1 = newl; Line(L1) = {P1, P2};
L2 = newl; Line(L2) = {P2, P3};
L3 = newl; Line(L3) = {P3, P4};
L4 = newl; Line(L4) = {P4, P1};
LL1 = newll; Line Loop(LL1) = {L4, L1, L2, L3};
S1 = news; Plane Surface(S1) = {LL1};
num[] = Extrude {0, 0, 2.0} {
 Surface{S1};
 //Layers{5};
 //Recombine;
};

//Collect necessary outer surfaces
CSURFACE[0]=num[2];
CSURFACE[1]=num[3];
CSURFACE[2]=num[4];
CSURFACE[3]=num[5];
CSURFACE[4]=num[0];
CSURFACE[5]=S1;

/// TORUS ///

TR = 0.2; //torus small radius
DY = 0.7; //torus big radius
DZ = 1.0; //shifts
P11 = newp; Point(P11) = {0,DY,DZ,clen};
P12 = newp; Point(P12) = {0,DY,DZ+TR,tlen};
P13 = newp; Point(P13) = {0,DY+TR,DZ,tlen};
P14 = newp; Point(P14) = {0,DY-TR,DZ,tlen};
P15 = newp; Point(P15) = {0,DY,DZ-TR,tlen};
C11=newc; Circle(C11) = {P12,P11,P13};
C12=newc; Circle(C12) = {P13,P11,P15};
C13=newc; Circle(C13) = {P15,P11,P14};
C14=newc; Circle(C14) = {P14,P11,P12};
LL11 = newll; Line Loop(LL11) = {C11,C12,C13,C14};
S11 = news; Plane Surface(S11) = {LL11};
numT1[] = Extrude {{0, 0, 1}, {0, 0, 0}, 2*Pi/3} {
 Surface{S11};
};
numT2[] = Extrude {{0, 0, 1}, {0, 0, 0}, 2*Pi/3} {
 Surface{numT1[0]};
};
numT3[] = Extrude {{0, 0, 1}, {0, 0, 0}, 2*Pi/3} {
 Surface{numT2[0]};
};

TVOLUME[0]=numT1[1];
TVOLUME[1]=numT2[1];
TVOLUME[2]=numT3[1];

//Collect necessary torus surface in order to be joined by loops
TSURFACE[0]=numT1[2];
TSURFACE[1]=numT1[3];
TSURFACE[2]=numT1[4];
TSURFACE[3]=numT1[5];
TSURFACE[4]=numT2[2];
TSURFACE[5]=numT2[3];
TSURFACE[6]=numT2[4];
TSURFACE[7]=numT2[5];
TSURFACE[8]=numT3[2];
TSURFACE[9]=numT3[3];
TSURFACE[10]=numT3[4];
TSURFACE[11]=numT3[5];

//start and end surface from torus crossection
//TSURFACE[4]=numT1[0];
//TSURFACE[5]=S11;

/// VOLUME AND SURFACE ///

SL11 = newsl; Surface Loop(SL11) = {TSURFACE[]};
V11 = newv; Volume(V11) = {SL11};

//exclude the torus from the cuboid --> the cuboid surfaces are the
//set of the rectangle + torus surfaces (all boundaries!)
SL1 = newsl; Surface Loop(SL1) = {TSURFACE[],CSURFACE[]};
V1 = newv; Volume(V1) = {SL1};

//the cuboid rectangle outer surface to apply boundary conditions
SL2 = newsl; Surface Loop(SL2) = {CSURFACE[]};

Physical Volume(PHYSTORUSVOLUME) = V11;
//torus volume already excluded by the surface loop!
Physical Volume(PHYSCUBOIDTORUSVOLUME) = V1;

Physical Surface(PHYSTORUSSUFACE) = {TSURFACE[]};
//exclude the torus from the cuboid --> the cuboid surfaces are the
//set of the rectangle + torus surfaces (all boundaries!)
Physical Surface(PHYSCUBOIDTORUSSURFACE) = {CSURFACE[],TSURFACE[]};
//the cuboid rectangle outer surface to apply boundary conditions
Physical Surface(PHYSCUBOIDSURFACE) = {CSURFACE[]};

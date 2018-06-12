%demo8_slice3d_2dgrid_vectors.m Slicinging a 3D mesh and plotting the
%                               vector field at the cross section
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-05-19
%
% Copyright (C) 2018 Chloros2 <chloros2@gmx.de>
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see
% <https://www.gnu.org/licenses/>.
%

clear all;

addpath('ffmatlib');

S1=[-0 -0 0.2];
S2=[1 0.1 0];
S3=[0 1.6 0.5]/2.4;

[M] = ffreadfile('File1','temp_demo8_tetdata3d_box.txt', ...
                 'Delimiter',';','Format','%f %f %f %f %f %f %f');
tFluxData=M(:,1:6);
tTData=[M(:,1:3), M(:,7)];
[sliceTData] = slicetet2data(tTData,S1,S2,S3);

%Plots the vector field
N=20;
M=20;
[X,Y,Z] = gridplane3d(S1,S2,S3,N,M);
%[qx,qy,qz] = fftet2gridfast(tFluxData,X,Y,Z);
[qx,qy,qz] = fftet2grid(tFluxData,X,Y,Z);
figure;
quiver3(X,Y,Z,qx,qy,qz,2.0);
hold on;

%Plots the temperature
N=100;
M=100;
[X,Y,Z] = gridplane3d(S1,S2,S3,N,M);
%[C] = fftet2gridfast(sliceTData,X,Y,Z);
[C] = fftet2grid(sliceTData,X,Y,Z);
surf(X,Y,Z,C,'EdgeColor','none');
colormap(jet(192));
caxis([min(min(C)) max(max(C))]);
hcb=colorbar;
title(hcb,'dT');
hold on;
zlabel('z');
ylabel('y');
xlabel('x');
view(3);
title('Heat Flux and Temperature (Quiver3)');
daspect([1 1 1]);

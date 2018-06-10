%demo8_slice3d_2dgrid_vectors.m Slice a 3d mesh, interpolate cross section to
%                               rectangular grid and show heat flux vector field
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

%%%%%% Interpolation of a cross section on a rectangular grid

%Slicing plane definition - three points defining a plane
S1=[-0 -0 0.2];
S2=[1 0.1 0];
S3=[0 1.6 0.5]/2.4;

%Creates the rectangular grid
N=50;
M=50;
[X,Y,Z] = gridplane3d(S1,S2,S3,N,M);

[tdata] = ffreadfile('File1','temp_demo8_tetdata3d_box.txt', ...
                     'Delimiter',';','Format','%f %f %f %f %f %f %f');

[qx qy qz C] = fftet2gridfast(tdata,X,Y,Z);

figure;
%Shows heat flux vector field at cross section
quiver3(X,Y,Z,qx,qy,qz,3.0);
hold on;
%Shows temperature at cross section
surf(X,Y,Z,C);
colormap(jet(192));
caxis([min(min(C)) max(max(C))]);
hcb=colorbar;
title(hcb,'Temp');
hold on;
zlabel('z');
ylabel('y');
xlabel('x');
view(3);
title('Heat Flux Vector Field at Cross Section');
daspect([1 1 1]);

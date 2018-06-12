%demo7_slice3d_2dgrid.m Slicing a 3D mesh and interpolating the cross section
%                       on a rectangular grid
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

S1=[-0.01 -0.01 -0.01];
S2=[1.01 0.3 -0.01];
S3=[-0.01 0.82 0.42];

[bdata,tdata] = ffreadfile('File1','temp_demo7_bddata3d_box.txt', ...
                           'File2','temp_demo7_tetdata3d_box.txt', ...
                           'Delimiter',';','Format','%f %f %f %f');
[sliceTData] = slicetet2data(tdata,S1,S2,S3);

%Creates a rectangular grid
N=100;
M=100;
[X,Y,Z] = gridplane3d(S1,S2,S3,N,M);
%Do the interpolation
tic;
[C] = fftet2grid(sliceTData,X,Y,Z);
toc;
%Consider changing fftet2grid() into fftet2gridfast() which is approx.
%x50 (Matlab) and x450 (Octave) faster (MEX implementation)
%[C] = fftet2gridfast(sliceTData,X,Y,Z);

%%%%%% Shows the interpolation

figure();
surf(X,Y,Z,C,'EdgeColor','none');
colormap(jet(192));
caxis([min(min(C)) max(max(C))]);
hcb=colorbar;
title(hcb,'dT[K]');
hold on;
plot3([S1(1) S2(1)],[S1(2) S2(2)],[S1(3) S2(3)], ...
      '-m','LineWidth',2);
plot3([S1(1) S3(1)],[S1(2) S3(2)],[S1(3) S3(3)], ...
      '-m','LineWidth',2);
plot3([S2(1) (S2(1)+(S3(1)-S1(1)))],[S2(2) (S2(2)+(S3(2)-S1(2)))], ...
      [S2(3) (S2(3)+(S3(3)-S1(3)))],'-m','LineWidth',2);
plot3([S3(1) (S2(1)+(S3(1)-S1(1)))],[S3(2) (S2(2)+(S3(2)-S1(2)))], ...
      [S3(3) (S2(3)+(S3(3)-S1(3)))],'-m','LineWidth',2);
zlabel('z');
ylabel('y');
xlabel('x');
view(3);
title('Grid Interpolation');
daspect([1 1 1]);

%%%%%% Shows the boundary data and the slicing plane

figure();
ax=axes();
[sz1,sz2]=size(Z);
ZZ=0.5*ones(sz1,sz2);
surf(X,Y,Z,ZZ,'EdgeColor','none');
%alpha(0.7);
hold on;
[BX,BY,BZ,BC] = slicebd2patch(bdata,S1,S2,S3);
[~,sz2]=size(BC);
patch(BX,BY,BZ,repmat([0;1;0],1,sz2),'EdgeColor','none');
text(S1(1),S1(2),S1(3),'S1','HorizontalAlignment','center', ...
     'FontSize',15,'FontWeight','bold','Color','m');
text(S2(1),S2(2),S2(3),'S2','HorizontalAlignment','center', ...
     'FontSize',15,'FontWeight','bold','Color','m');
text(S3(1),S3(2),S3(3),'S3','HorizontalAlignment','center', ...
     'FontSize',15,'FontWeight','bold','Color','m');
zlabel('z');
ylabel('y');
xlabel('x');
title('Plane Definition');
view(3);
daspect([1 1 1]);
%viridis(3)
map=[0.2670040   0.0048743   0.3294152
     0.1281485   0.5651070   0.5508924
     0.9932479   0.9061566   0.1439362];
colormap(map);
props = {'CameraViewAngle','DataAspectRatio','PlotBoxAspectRatio'};
set(ax,props,get(ax,props));
set(gcf,'color',[0.9 0.9 0.9]);

%%%%%% Plots a 2d projection

[UN,UM] = gridplane2d(S1,S2,S3,N,M);
figure();
surf(UN,UM,C,'EdgeColor','none');
colormap(jet(192));
caxis([min(min(C)) max(max(C))]);
title('Plane Projection');
ylabel('M');
xlabel('N');
view(2);
axis tight equal;
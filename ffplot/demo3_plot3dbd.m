%demo3_plot3dbd.m Plot FreeFem++ 3d simulation results.
%                 Plot boundary surface.
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

[X,Y,Z,C] = ffread2patch('temp_demo3_bddata3d_box.txt','Delimiter',';','Format','auto');

%%%%%% 3D surf plot. Surface colored by value.

figure;
patch(X,Y,Z,C,'EdgeColor',[0 0 0],'LineWidth',1);
colormap(jet(250));
caxis([min(min(C)) max(max(C))]);
hcb=colorbar;
title(hcb,'dT[K]');
zlabel('z');
ylabel('y');
xlabel('x');
title('A 3d plot of a surface boundary');
view(3);
daspect([1 1 1]);
%Alpha transparency can be added if this uncommented:
%alpha(0.7);

%%%%%% Again 3D surf plot this time but with fixed size during rotation

figure;
ax=axes();
patch(X,Y,Z,C,'EdgeColor',[0 0 0],'LineWidth',1);
colormap(jet(250));
caxis([min(min(C)) max(max(C))]);
hcb=colorbar;
title(hcb,'dT[K]');
zlabel('z');
ylabel('y');
xlabel('x');
view(3);
daspect([1 1 1]);
props = {'CameraViewAngle','DataAspectRatio','PlotBoxAspectRatio'};
set(ax,props,get(ax,props));
set(gcf,'color',[0.9 0.9 0.9]);
set(ax, 'Visible','off');


[X,Y,Z] = ffread2patch('temp_demo3_bddata3d_cyl.txt','Delimiter',';','Format','auto');

%%%%%% 3D surface plot of a mesh

figure;
patch(X,Y,Z,[1 1 1],'EdgeColor',[0 0 1],'LineWidth',1);
zlabel('z');
ylabel('y');
xlabel('x');
title('A 3d plot of a surface boundary');
view(3);
daspect([1 1 1]);
%alpha(0.7);

%demo1_getstarted2.m Plot FreeFem++ 2d simulation results
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-05-13
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

%Set path where to find ffread2patch.m
addpath('ffmatlib');

%Read FreeFem++ file and convert data in order to plot the facets
[X,Y,C]=ffread2patch('temp_demo1_getstarted.txt');

%Draw triangles. X, Y contain vertice coordinates and C the color information.
%Set a black edge color in order to show the mesh.
patch(X,Y,C,C,'EdgeColor',[0 0 0],'LineWidth',1);
%Create a colormap with jet array
colormap(jet(192));
caxis([min(min(C)) max(max(C))]);
%Create a colorbar
colorbar;
%Viewpoint specification: 3D plot.
view(3);
%Set aspect ratio for x,y,z
daspect([1 1 1*(max(max(C))-min(min(C)))]);

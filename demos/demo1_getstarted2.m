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

%Specifies the path where the ffmatlib can be found
addpath('ffmatlib');

%Reads the file content and converts to patch() - plot data
[X,Y,C]=ffread2patch('temp_demo1_getstarted.txt');

%%%%%%% 3D surf plot

%Plots the facets and display the mesh
patch(X,Y,C,C,'EdgeColor',[0 0 0],'LineWidth',1);
%Creates a colormap with jet array
colormap(jet(192));
caxis([min(min(C)) max(max(C))]);
%Creates a colorbar
colorbar;
%Sets the view point specification to 3d
view(3);
%Sets 1:1:1 aspect ratio
daspect([1 1 1*(max(max(C))-min(min(C)))]);
grid;

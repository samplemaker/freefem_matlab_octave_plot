%demo2_plot2d.m Plot FreeFem++ 2d simulation results
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

clear all;
addpath('ffmatlib');

[X,Y,C] = ffread2patch('temp_demo2_tridata2d.txt','Delimiter',';','Format','auto');

%%%%%% 2D density plot

figure;
%No Edge Color (hide the mesh)
patch(X,Y,C,'EdgeColor','none');
colormap(jet(250));
caxis([min(min(C)) max(max(C))]);
hcb=colorbar;
title(hcb,'dT[K]');
ylabel('y');
xlabel('x');
title('A 2d density plot');
view(2);
axis tight equal;

%%%%%% 2D mesh plot

figure;
%Plots white facets with blue edge color in order to show the mesh
patch(X,Y,[1 1 1],'EdgeColor',[0 0 1],'LineWidth',1);
ylabel('y');
xlabel('x');
title('A 2d mesh plot');
view(2);
axis tight equal;

%%%%%%% 3D surf plot

figure;
%Same works for 3D plots
patch(X,Y,C,C,'EdgeColor','none', ...
      'FaceLighting','gouraud','AmbientStrength', 0.3);
colormap(jet(250));
caxis([min(min(C)) max(max(C))]);
hcb=colorbar;
title(hcb,'dT[K]');
camlight('left');
zlabel('u');
ylabel('y');
xlabel('x');
title('A 3d surf plot');
view(3);
daspect([1 1 6*(max(max(C))-min(min(C)))]);

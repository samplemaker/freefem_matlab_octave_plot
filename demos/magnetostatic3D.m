%magnetostatic3D.m Magnetostatic problem of a toroidal current
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-11-30
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

%see torus.geo
PHYSTORUSSUFACE = 201;
PHYSCUBOIDSURFACE = 103;

[p,b,t,nv,nbe,nt,labels]=ffreadmesh('torus.mesh');
[Ex,Ey,Ez]=ffreaddata('torus.txt');

figure;
ffpdeplot3D(p,b,t,'XYZStyle','monochrome','BDLabels',[PHYSTORUSSUFACE]);
hold on;
ffpdeplot3D(p,b,t,'XYZStyle','noface','BDLabels',[PHYSCUBOIDSURFACE]);
axis tight;
str={sprintf('nVertex: %i',nv);
     sprintf('nTets:   %i',nt);
     sprintf('nTris:   %i',nbe)};
annotation('textbox',[0.05 0.05 0.2 0.15],'String',str, ...
           'FitBoxToText','on','FontName','Courier','EdgeColor',[1 1 1]);
ylabel('y');
xlabel('x');
zlabel('z');

figure;
ffpdeplot3D(p,b,t,'FlowData',[Ex,Ey,Ez],'FGridParam3D',[9,9,9],'Boundary','on','BDLabels',PHYSTORUSSUFACE,'XYZStyle','monochrome');
ylabel('y');
xlabel('x');
zlabel('z');

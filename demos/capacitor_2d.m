%capacitor_2d.m The parallel plate capacitor problem
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

%Reads a FreeFem++ mesh created with the savemesh(Th,"mesh.msh"); command
[p,b,t,nv,nbe,nt,labels]=ffreadmesh('capacitorp1.msh');
%Reads the PDE data
[u]=ffreaddata('capacitor_potential_p1only.txt');
[Ex,Ey]=ffreaddata('capacitor_field_p1only.txt');

%%%%%% 2D Patch (density map) Plot

handles=ffpdeplot(p,b,t, ...
                  'XYData',u, ...
                  'Mesh','on', ...
                  'Boundary','on', ...
                  'XLim',[-2 2],'YLim',[-2 2], ...
                  'CBTitle','U[V]', ...
                  'Title','2D Patch Plot (Electrostatic Potential)');

ylabel('y');
xlabel('x');

%%%%%% Mesh Plot

figure;

handles=ffpdeplot(p,b,t, ...
                  'Mesh','on', ...
                  'Boundary','on', ...
                  'Title','Boundary/Edge (Capacitor Electrodes) and 2D Mesh');

ylabel('y');
xlabel('x');

%%%%%% 3D Surf Plot

figure;

handles=ffpdeplot(p,b,t, ...
                  'XYData',u, ...
                  'ZStyle','continuous', ...
                  'Mesh','on', ...
                  'CBTitle','U[V]', ...
                  'Title','3D Patch Plot (Electrostatic Potential)');
ylabel('y');
xlabel('x');
zlabel('u');
grid;

%%%%%% Combine Quiver and Contour

figure;

handles=ffpdeplot(p,b,t, ...
                  'XYData',u, ...
                  'Mesh','off', ...
                  'Boundary','on', ...
                  'XLim',[-2 2],'YLim',[-2 2], ...
                  'Contour','on', ...
                  'CStyle','monochrome', ...
                  'CColor','b', ...
                  'XYStyle','off', ...
                  'CGridParam',[150, 150], ...
                  'ColorBar','off', ...
                  'FlowData',[Ex,Ey], ...
                  'FGridParam',[25, 25], ...
                  'Title','Quiver+Contour Interpolation Plot');

ylabel('y');
xlabel('x');

%%%%%% Show Labels

figure;

handles=ffpdeplot(p,b,t, ...
                  'Mesh','off', ...
                  'Boundary','on', ...
                  'BDLabels',labels, ...
                  'Title','Boundary Labels');

%%%%%% 3D Surf Plot Gouraud lighting

figure;

handles=ffpdeplot(p,b,t, ...
                  'XYData',u, ...
                  'ZStyle','continuous', ...
                  'Mesh','off', ...
                  'CBTitle','U[V]', ...
                  'Title','Gouraud');
ylabel('y');
xlabel('x');
zlabel('u');
grid;

lighting gouraud;
view([-47,24]);
camlight('headlight');

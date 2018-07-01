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
[nv,nbe,nt,points,boundary,triangles]=ffreadmesh('capacitorp1.msh');

fid=fopen('capacitor_potential_p1only.txt','r');
data=textscan(fid,'%f','Delimiter','\n');
fclose(fid);
u=cell2mat(data);
[sz1,sz2]=size(u);
fprintf('Size of data (nDof): %i %i\n', sz1,sz2);

fid=fopen('capacitor_field_p1only.txt','r');
data=textscan(fid,'%f %f','Delimiter','\n');
fclose(fid);
v=cell2mat(data)';


%%%%%% 2D Patch (density map) Plot

handles=pdeplot2dff(points,boundary,triangles, ...
                    'XYData',u, ...
                    'Mesh','on', ...
                    'Edge','on', ...
                    'XLim',[-2 2],'YLim',[-2 2], ...
                    'Title','2D Patch Plot (Electrostatic Potential)');

title(handles(2),'U[V]');
ylabel('y');
xlabel('x');

%%%%%% Mesh Plot

figure;

handles=pdeplot2dff(points,boundary,triangles, ...
                    'Mesh','on', ...
                    'Edge','on', ...
                    'Title','Boundary/Edge (Capacitor Electrodes) and 2D Mesh');

ylabel('y');
xlabel('x');

%%%%%% 3D Surf Plot

figure;

handles=pdeplot2dff(points,boundary,triangles, ...
                    'XYData',u, ...
                    'ZStyle','continuous', ...
                    'Mesh','on', ...
                    'Title','3D Patch Plot (Electrostatic Potential)');
ylabel('y');
xlabel('x');
zlabel('u');
title(handles(2),'U[V]');
grid;

%%%%%% Combine Quiver and Contour

figure;

handles=pdeplot2dff(points,boundary,triangles, ...
                    'XYData',u, ...
                    'Mesh','off', ...
                    'Edge','on', ...
                    'XLim',[-2 2],'YLim',[-2 2], ...
                    'Contour','on', ...
                    'CColor',[0 0 1], ...
                    'XYStyle','off', ...
                    'CGridParam',[150, 150], ...
                    'ColorBar','off', ...
                    'FlowData',v, ...
                    'FGridParam',[60, 60], ...
                    'Title','Quiver+Contour Interpolation Plot');

ylabel('y');
xlabel('x');

%%%%%% 3D Surf Plot Gouraud lighting

figure;

handles=pdeplot2dff(points,boundary,triangles, ...
                    'XYData',u, ...
                    'ZStyle','continuous', ...
                    'Mesh','off', ...
                    'Title','');
ylabel('y');
xlabel('x');
zlabel('u');
title(handles(2),'U[V]');
grid;

lighting gouraud;
view([-47,24]);
camlight('headlight');
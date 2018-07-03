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
[points,boundary,triangles,nv,nbe,nt,labels]=ffreadmesh('capacitorp1.msh');
%Reads the PDE data
[u]=ffreaddata('capacitor_potential_p1only.txt');
[Ex,Ey]=ffreaddata('capacitor_field_p1only.txt');

%%%%%% 2D Patch (density map) Plot

handles=ffpdeplot(points,boundary,triangles, ...
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

handles=ffpdeplot(points,boundary,triangles, ...
                  'Mesh','on', ...
                  'Edge','on', ...
                  'Title','Boundary/Edge (Capacitor Electrodes) and 2D Mesh');

ylabel('y');
xlabel('x');

%%%%%% 3D Surf Plot

figure;

handles=ffpdeplot(points,boundary,triangles, ...
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

handles=ffpdeplot(points,boundary,triangles, ...
                  'XYData',u, ...
                  'Mesh','off', ...
                  'Edge','on', ...
                  'XLim',[-2 2],'YLim',[-2 2], ...
                  'Contour','on', ...
                  'CColor',[0 0 1], ...
                  'XYStyle','off', ...
                  'CGridParam',[150, 150], ...
                  'ColorBar','off', ...
                  'FlowData',[Ex,Ey], ...
                  'FGridParam',[60, 60], ...
                  'Title','Quiver+Contour Interpolation Plot');

ylabel('y');
xlabel('x');

%%%%%% Show Labels

figure;

handles=ffpdeplot(points,boundary,triangles, ...
                  'Mesh','off', ...
                  'Edge','on', ...
                  'Title','Boundary Labels');

for i=1:numel(labels)
    xpts=points(1,:);
    ypts=points(2,:);
    labpts=points(3,:);
    pos=find((labpts==labels(i)),1,'first');
    text(xpts(pos),ypts(pos),[num2str(labels(i))]);
end

%%%%%% 3D Surf Plot Gouraud lighting

figure;

handles=ffpdeplot(points,boundary,triangles, ...
                  'XYData',u, ...
                  'ZStyle','continuous', ...
                  'Mesh','off', ...
                  'Title','Gouraud');
ylabel('y');
xlabel('x');
zlabel('u');
title(handles(2),'U[V]');
grid;

lighting gouraud;
view([-47,24]);
camlight('headlight');

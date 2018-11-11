%periodic_bc.m Plot a periodic boundary condition problem
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-10-31
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
[p,b,t,nv,nbe,nt,labels]=ffreadmesh('periodic.msh');
%Reads the PDE data
[u]=ffreaddata('periodic.txt');

%%%%%% Mesh Plot

figure;

handles=ffpdeplot(p,b,t, ...
                  'Mesh','on', ...
                  'Boundary','on', ...
                  'Title','Mesh only');

ylabel('y');
xlabel('x');

%%%%%% Show Labels

figure;

handles=ffpdeplot(p,b,t, ...
                  'Mesh','off', ...
                  'Boundary','on', ...
                  'BDLabels',labels, ...
                  'BDShowText','on', ...
                  'Title','Boundary Labels');

%%%%%% 2D Patch (density map) Plot + Contour

figure;

handles=ffpdeplot(p,b,t, ...
                  'XYData',u, ...
                  'Mesh','off', ...
                  'Boundary','on', ...
                  'CBTitle','u', ...
                  'Contour','on', ...
                  'Title','Patch Plot + Contour');

ylabel('y');
xlabel('x');
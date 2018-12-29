%testsubplot.m Check if the subplot command is working
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


figure('position', [0, 0, 800, 300])
subplot(1,2,1);
handles=ffpdeplot(p,b,t, ...
                  'XYData',u, ...
                  'Mesh','off', ...
                  'XLim',[-2 2],'YLim',[-2 2], ...
                  'CBTitle','U[V]', ...
                  'Title','2D Patch Plot (Electrostatic Potential)');
hold on;
ffpdeplot(p,b,t, ...
                  'Mesh','off', ...
                  'Boundary','on', ...
                  'BDLabels',[4,3], ...
                  'BDShowtext', 'on', ...
                  'Title','Boundary Labels');

subplot(1,2,2);
handles=ffpdeplot(p,b,t, ...
                  'Mesh','on', ...
                  'Boundary','on', ...
                  'Title','Boundary/Edge (Capacitor Electrodes) and 2D Mesh');

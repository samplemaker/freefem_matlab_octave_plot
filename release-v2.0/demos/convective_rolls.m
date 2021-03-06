%vortex_rolls.m Free convection problem between to flat plates
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
[p,b,t]=ffreadmesh('convective_rolls.msh');
[vh]=ffreaddata('convective_rolls_vh.txt');
[qh]=ffreaddata('convective_rolls_qh.txt');
[Tfluid]=ffreaddata('convective_rolls_temperature.txt');
[psi]=ffreaddata('convective_rolls_psi.txt');
[vx,vy]=ffreaddata('convective_rolls_ux_uy.txt');

figure('Resize','off','ToolBar','none','MenuBar','none','Position',[50 50 1000 460]);
subplot(2,1,1);

%%%%%% Additive plotting: Combine Patch and Contour with different data
%Scalar data overlayed with extra contour data from a different origin
%like e.g. a weather forecast plot employing a temperature patch plot
%with overlayed isolevel pressure values

%Patch plot of a temperature
ffpdeplot(p,b,t, ...
          'vhseq',qh, ...
          'XYData',Tfluid, ...
          'ColorMap','jet', ...
          'CBTitle','T[degC]', ...
          'Title','Stream lines + Temperature');

hold on;

%Add Stream lines of the velocity field
ffpdeplot(p,b,t, ...
          'vhseq',qh, ...
          'XYData',psi, ...
          'XYStyle','off', ...
          'Contour','on', ...
          'CGridParam',[150, 150], ...
          'CStyle','dashedneg', ...
          'ColorMap','off', ...
          'ColorRange','off', ...
          'ColorBar','off');

ylabel('y');
xlabel('x');
axis tight;

%%%%%% Combine Quiver and Patch (Combine P1 and P2 by incrementel plotting)

subplot(2,1,2);

ffpdeplot(p,b,t, ...
          'vhseq',qh, ...
          'XYData',Tfluid, ...
          'ColorMap','jet', ...
          'CBTitle','T[degC]', ...
          'Title','Velocity + Temperature');

hold on;

ffpdeplot(p,b,t, ...
          'vhseq',vh, ...
          'FlowData',[vx, vy], ...
          'FGridParam',[60, 15]);

ylabel('y');
xlabel('x');
axis tight;


%%%%%% Default Contour plot

figure('Resize','off','ToolBar','none','MenuBar','none','Position',[50 50 1000 460]);
subplot(2,1,1);

ffpdeplot(p,b,t, ...
          'vhseq',qh, ...
          'XYData',psi, ...
          'Contour','on', ...
          'XYStyle','off', ...
          'CGridParam',[150, 150], ...
          'ColorMap','jet', ...
          'CColor','flat', ...
          'CBTitle','psi [(m^3/s)/m]', ...
          'Title','Streamlines');

ylabel('y');
xlabel('x');
axis tight;

subplot(2,1,2);

% Isotherms

ffpdeplot(p,b,t, ...
          'vhseq',qh, ...
          'XYData',Tfluid, ...
          'Contour','on', ...
          'ColorMap','jet', ...
          'CStyle','solid', ...
          'CColor','k', ...
          'CGridParam',[150, 150], ...
          'CBTitle','dT [K]', ...
          'Title','Isotherms');

ylabel('y');
xlabel('x');
axis tight;


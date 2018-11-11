%complex_pde_v4.m Conformal plot of a complex PDE problem
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

[p,b,t,nv,nbe,nt,labels]=ffreadmesh('capacitorp1.msh');
[uC]=ffreaddata('capacitor_potential_p1only.txt');
xpts=p(1,:);
ypts=p(2,:);
xdata=[xpts(t(1,:)); xpts(t(2,:)); xpts(t(3,:))];
ydata=[ypts(t(1,:)); ypts(t(2,:)); ypts(t(3,:))];
udata=[uC(t(1,:)), uC(t(2,:)), uC(t(3,:))].';

N = 100;
s = linspace(0,2*pi(),N);
Z = 3.5*(cos(s)+1i*sin(s)).*sin(0.5*s);
W = fftri2gridcplx(Z,xdata,ydata,udata);

figure('position', [0, 0, 800, 300])

subplot(1,2,1);
hold on;
ffpdeplot(p,b,t, ...
          'XYData',uC, ...
          'ZStyle','continuous', ...
          'Mesh','off', ...
          'ColorBar','off', ...
          'Title','Single curve Interpolation');

plot3(real(Z),imag(Z),real(W),'g');
ylabel('y');
xlabel('x');
zlabel('u');
grid;
lighting gouraud;
view([-47,24]);
camlight('headlight');

subplot(1,2,2);
plot(s,real(W),'b');
grid;
xlabel('s');
ylabel('u');
title('Interpolation Values');
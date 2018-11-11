%complex_pde_v2.m Conformal plot of a complex PDE problem
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

[g_phi, g_r] = ffcplxmesh([1.5,0], [-1,2*pi()], [10,10], [10,10]);
%Map to polar coordinates
p_phi = exp(g_phi);
p_r = exp(g_r);
w_phi = fftri2gridcplx(p_phi,xdata,ydata,udata);
w_r = fftri2gridcplx(p_r,xdata,ydata,udata);

figure;
hold on;
handles=ffpdeplot(p,b,t, ...
                  'XYData',uC, ...
                  'ZStyle','continuous', ...
                  'Mesh','off', ...
                  'CBTitle','U[V]', ...
                  'Title','Curved Interpolation');

plot3(real(p_phi),imag(p_phi),real(w_phi),'g');
plot3(real(p_r),imag(p_r),real(w_r),'g');
% plot3(real(p_phi),imag(p_phi),real(w_phi),'*r');
% plot3(real(p_r),imag(p_r),real(w_r),'*r');

ylabel('y');
xlabel('x');
zlabel('u');
grid;

lighting gouraud;
view([-47,24]);
camlight('headlight');


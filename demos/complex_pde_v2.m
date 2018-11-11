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

[p,b,t,nv,nbe,nt,labels]=ffreadmesh('complex_test.msh');
[uC]=ffreaddata('complex_test.txt');
xpts=p(1,:);
ypts=p(2,:);
xdata=[xpts(t(1,:)); xpts(t(2,:)); xpts(t(3,:))];
ydata=[ypts(t(1,:)); ypts(t(2,:)); ypts(t(3,:))];
udata=[uC(t(1,:)), uC(t(2,:)), uC(t(3,:))].';

figure;
hold on;
[Z1, Z2] = ffcplxmesh([0,0], [2*pi(),2*pi()], [3,3], [5,5]);
if ~isempty(Z1) %numx = 0
    W1 = fftri2gridcplx(Z1,xdata,ydata,udata);
    plot(real(W1),imag(W1),'b','LineWidth',2);
end
if ~isempty(Z2) %numy = 0
    W2 = fftri2gridcplx(Z2,xdata,ydata,udata);
    plot(real(W2),imag(W2),'r','LineWidth',2);
end
grid;
daspect([1,1,1]);

figure('position', [0, 0, 800, 300]);
subplot(1,2,1);
ffpdeplot(p,b,t,'XYData',real(uC),'Boundary','on','CBTitle','Im(u)','Mesh','on');
subplot(1,2,2);
ffpdeplot(p,b,t,'XYData',imag(uC),'Boundary','on','CBTitle','Re(u)','Mesh','on');


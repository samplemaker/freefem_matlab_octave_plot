%runplotdemo3d.m Plot some 3d FreeFem++ simulation results
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-05-19
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

%rearrange all columns in order to plot the facets using the patch() command
[XX,YY,ZZ,CC] = ff2patch('tridata3dbox.txt','Delimiter',';','Format','auto');

%3D surf plot with scalar value
figure;
%patch(XX,YY,ZZ,CC,'FaceColor','interp');
patch(XX,YY,ZZ,CC,'EdgeColor',[0 0 0],'LineWidth',1);
colormap(jet(250));
caxis([min(min(CC)) max(max(CC))]);
colorbar;
zlabel('z');
ylabel('y');
xlabel('x');
view(3);
daspect([1 1 1]);
%alpha(0.7);

%3D surf plot of a mesh only
[XX,YY,ZZ] = ff2patch('tridata3dcyl.txt','Delimiter',';','Format','auto');
figure;
patch(XX,YY,ZZ,[1 1 1],'EdgeColor',[0 0 1],'LineWidth',1);
zlabel('z');
ylabel('y');
xlabel('x');
view(3);
daspect([1 1 1]);
%alpha(0.7);

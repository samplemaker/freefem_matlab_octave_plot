%runplotdemo2d.m Plot some FreeFem++ simulation results
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

%rearrange all columns in order to plot the facets using the patch() command
[XX,YY,CC] = ff2patch('tridata2ddisc.txt','Delimiter',';','Format','auto');

%3D surf plot
figure;
patch(XX,YY,CC,CC, ...
      'EdgeColor','none', ...
      'FaceLighting','gouraud', ...
      'AmbientStrength', 0.3);
colormap(jet(250));
caxis([min(min(CC)) max(max(CC))]);
colorbar;
camlight('left');
zlabel('u');
ylabel('y');
xlabel('x');
view(3);
daspect([1 1 6*(max(max(CC))-min(min(CC)))]);

%2D density plot
figure;
patch(XX,YY,CC,'EdgeColor','none');
colormap(jet(250));
caxis([min(min(CC)) max(max(CC))]);
colorbar;
ylabel('y');
xlabel('x');
view(2);
axis equal;

%2D mesh plot
figure;
patch(XX,YY,[1 1 1],'EdgeColor',[0 0 1],'LineWidth',1);
ylabel('y');
xlabel('x');
view(2);
axis equal;

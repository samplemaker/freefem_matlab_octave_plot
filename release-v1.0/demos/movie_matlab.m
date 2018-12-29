%movie_matlab.m Creates a movie
%
% Note: This Code is not compatible with Octave!
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
%read simulation data
addpath('ffmatlib');
[p,b,t,nv,nbe,nt,labels] = ffreadmesh('movie.msh');
n=250;
for j=1:n
    name = sprintf('movie_temp_%i.txt', j+10000);
    u(:,j) = ffreaddata(name);
    name = sprintf('movie_psi_%i.txt', j+10000);
    psi(:,j) = ffreaddata(name);
end
%capture frames from the simulation data
figure('Resize','off','ToolBar','none','MenuBar','none','Position',[50 50 1000 230]);
set(gcf,'Renderer','OpenGL');
opengl('software'); %a getframe issue on older Matlab versions
ffpdeplot(p,b,t,'XYData',u(:,1), ...
          'Contour','on', ...
          'CXYData',psi(:,1), ...
          'CStyle','patchdashedneg', ...
          'CGridParam',[250, 100], ...
          'Boundary','off', ...
          'Title','Frame 0','CBTitle','dT[K]', ...
          'ColorMap',cool(150),'ColorRange', [0 50]);
xlabel('x');
ylabel('y');
axis tight manual;
hax = gca;
set(hax, 'NextPlot', 'replaceChildren');
m = moviein(n);
for j = 1:n
    titlestr = sprintf('Free Convection Cavity [%0.1fsec]', (j-1)*0.04);
    ffpdeplot(p,b,t,'XYData',u(:,j), ...
              'Contour','on', ...
              'CXYData',psi(:,j), ...
              'CStyle','patchdashedneg', ...
              'CGridParam',[250, 100], ...
              'Boundary','off', ...
              'Title',titlestr,'CBTitle','dT[K]', ...
              'ColorMap',cool(150),'ColorRange', [0 50]);
    drawnow;
    m(:,j) = getframe(gcf);
end
%create AVI file from animation
v = VideoWriter('movie.avi');
set(v, 'FrameRate', 15);
set(v, 'Quality', 95);
open(v);
for j = 1:n
    writeVideo(v, m(:,j));
end
close(v);
close all;

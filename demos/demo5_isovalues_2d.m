%demo5_isovalues_2d.m Plot isolines of FreeFem ++ 2d simulation results
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

[tridata] = ffreadfile('File1','temp_demo5_isovalues.txt', ...
                       'Delimiter',';','Format','%f %f %f');

%%%%%% Interpolation on a rectangular grid

%Choose a grid resolution
N=200;
x=linspace(-1,1,N);
y=linspace(-1,1,N);
tic;
C=fftri2grid(tridata,x,y);
toc;
[X,Y] = meshgrid(x,y);

%%%%%% Isovalue plot

figure();
[c,h]=contour(X,Y,1000*C,8);
texth=clabel(c,h,'fontsize', 8);
for i=1:size(texth)
    textstr=get(texth(i),'String');
    textnum=str2double(textstr);
    textstrnew=sprintf('%0.1f', textnum);
    set(texth(i),'String',textstrnew);
end
title('isovalues x1000');
axis tight equal;
%grid;

%%%%%% Surf plot

figure();
surf(X,Y,C,'EdgeColor','none');
colormap(jet(192));
caxis([min(min(C)) max(max(C))]);
colorbar;
ylabel('y');
xlabel('x');
view(3);
daspect([1 1 1*(max(max(C))-min(min(C)))]);

%%%%%% Patch plot

figure;
[PX,PY,PC]=fftri2patch(tridata);
patch(PX,PY,PC,PC);
colormap(jet(250));
caxis([min(min(C)) max(max(C))]);
hcb=colorbar;
ylabel('y');
xlabel('x');
view(3);
daspect([1 1 1*(max(max(C))-min(min(C)))]);
grid;

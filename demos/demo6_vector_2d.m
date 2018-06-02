%demo2_plot2d_isovalues.m Plot a 2d vector field from FreeFem ++ 2d simulation
%                         results
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

[fdata] = ffreadfile('File1','temp_demo6_vector.txt', ...
                     'Delimiter',';','Format','%f %f %f %f');

%%%%%% Interpolation on a rectangular grid

%Choose grid resolution
N=10;
M=14;
x=linspace(0,1,N);
y=linspace(0,1,M);
tic;
[UH,VH]=fftri2grid(fdata,x,y);
toc;
[X,Y] = meshgrid(x,y);

%%%%%% 2d vector field plot

figure();
quiver(X,Y,UH,VH);
ylabel('y');
xlabel('x');
axis tight equal;
grid;


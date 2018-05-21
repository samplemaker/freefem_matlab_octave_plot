%rundemoslice3d.m Slice a 3d plot
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

%define the slicing plane; three points S1 .. S3 containing the point coordinates
%x,y,z respectively
S1=[0 0 0]';
S2=[1 0.5 0]';
S3=[1 1.3 1]';

%S1=[0 0.75 0]';
%S2=[1 0.75 0]';
%S3=[1 0.75 1]';

tic;
[XX YY ZZ CC]=ffslice3d('bdtridata.txt', ...
                        'tetrahedrondata.txt', ...
                        S1,S2,S3, ...
                        'Delimiter',';','Format','%f %f %f %f');
toc;

%plot the sliced object
figure;
patch(XX,YY,ZZ,CC);
colormap(jet(250));
caxis([min(min(CC)) max(max(CC))]);
colorbar;
zlabel('z');
ylabel('y');
xlabel('x');
view(3);
daspect([1 1 1]);

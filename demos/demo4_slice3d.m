%demo4_slice3d.m Slices a 3D plot
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

addpath('ffmatlib');

[bdata,tdata] = ffreadfile('File1','temp_demo4_bddata3d_box.txt', ...
                           'File2','temp_demo4_tetdata3d_box.txt', ...
                           'Delimiter',';','Format','%f %f %f %f');

%Slicing plane definition - 3 points (x,y,z)
S1=[0 0 0];
S2=[1 0.5 0];
S3=[1 1.3 1];

[BX,BY,BZ,BC] = slicebd2patch(bdata,S1,S2,S3);
[SX,SY,SZ,SC] = slicetet2patch(tdata,S1,S2,S3);

%%%%%% Plots cross section

figure;
patch(SX,SY,SZ,SC);
colormap(jet(192));
caxis([min(min(SC)) max(max(SC))]);
hcb=colorbar;
title(hcb,'dT[K]');
zlabel('z');
ylabel('y');
xlabel('x');
view(3);
title('Cross Section');
daspect([1 1 1]);

%%%%%% Plots boundary

figure;
patch(BX,BY,BZ,BC);
colormap(jet(192));
caxis([min(min(BC)) max(max(BC))]);
hcb=colorbar;
title(hcb,'dT[K]');
zlabel('z');
ylabel('y');
xlabel('x');
view(3);
title('Boundary');
daspect([1 1 1]);

%%%%%% Combine and plot cross section + boundary

figure;
patch([SX BX],[SY BY],[SZ BZ],[SC BC]);
colormap(jet(192));
caxis([min(min([SC BC])) max(max([SC BC]))]);
hcb=colorbar;
title(hcb,'dT[K]');
zlabel('z');
ylabel('y');
xlabel('x');
view(3);
title('Cross Section + Boundary');
daspect([1 1 1]);

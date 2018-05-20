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
fileID = fopen('tetrahedrondata.txt','r');
fdata = textscan(fileID,'%f %f %f %f %d','Delimiter',';');
[x,y,z,u,no] = fdata{:};

tx = vec2mat(x,4);
ty = vec2mat(y,4);
tz = vec2mat(z,4);
tc = vec2mat(u,4);

[sz1 ntetrahedron] = size(tx);
%split each tetrahedron into 4 triangles containing 3 nodes (vertices) each
XX=zeros(3,4*ntetrahedron);
YY=zeros(3,4*ntetrahedron);
ZZ=zeros(3,4*ntetrahedron);
CC=zeros(3,4*ntetrahedron);
DD=zeros(3,4*ntetrahedron);

i=1;
for j=1:ntetrahedron
  [XX(:,i) XX(:,i+1) XX(:,i+2) XX(:,i+3) ...
   YY(:,i) YY(:,i+1) YY(:,i+2) YY(:,i+3) ...
   ZZ(:,i) ZZ(:,i+1) ZZ(:,i+2) ZZ(:,i+3) ...
   CC(:,i) CC(:,i+1) CC(:,i+2) CC(:,i+3)]=tet2tri(tx(:,j),ty(:,j),tz(:,j),tc(:,j));
  %colorize each triangle in another color
  DD(:,i:i+3)=i;
  i=i+4;
end

figure;
%patch(XX,YY,ZZ,CC,'FaceColor','interp');
patch(XX,YY,ZZ,CC);
colormap(jet(250));
caxis([min(min(CC)) max(max(CC))]);
colorbar;
zlabel('z');
ylabel('y');
xlabel('x');
view(3);
daspect([1 1 1]);
%alpha(0.3);

figure;
patch(XX,YY,ZZ,DD);
colormap(jet(250));
caxis([min(min(DD)) max(max(DD))]);
colorbar;
zlabel('z');
ylabel('y');
xlabel('x');
view(3);
daspect([1 1 1]);
%alpha(0.3);

%slicing test
%identify all tetrahedrons affected by the slice
%maybe sort out all tetrahedrons in the slice plane and create new meshplot???
%erase marker
ID=false(1,4*ntetrahedron);
x0=0.0; y0=1.0;z0=0.0;
i=1;
for j=1:ntetrahedron
  if ~all((XX(:,i:i+3) >= x0) & ...
          (YY(:,i:i+3) >= y0) & ...
          (ZZ(:,i:i+3) >= z0) )
     ID(i:i+3)=[true true true true];
  end
  i=i+4;
end
%remove
XX(:,ID)=[];
YY(:,ID)=[];
ZZ(:,ID)=[];
DD(:,ID)=[];

figure;
patch(XX,YY,ZZ,DD);
colormap(jet(250));
caxis([min(min(DD)) max(max(DD))]);
colorbar;
zlabel('z');
ylabel('y');
xlabel('x');
view(3);
daspect([1 1 1]);

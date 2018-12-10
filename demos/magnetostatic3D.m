%magnetostatic3D.m Magnetostatic problem of a toroidal current
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-11-30
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

%see torus.geo
PHYSTORUSSUFACE = 201;
PHYSCUBOIDSURFACE = 103;

[p,b,t,nv,nbe,nt,labels]=ffreadmesh('torus.mesh');
[Bx,By,Bz]=ffreaddata('torus_induction.txt');
[Jx,Jy,Jz]=ffreaddata('torus_current.txt');
[u]=ffreaddata('torus_p.txt');

figure;
ffpdeplot3D(p,b,t,'XYZStyle','monochrome','BDLabels',[PHYSTORUSSUFACE]);
hold on;
ffpdeplot3D(p,b,t,'XYZStyle','noface','BDLabels',[PHYSCUBOIDSURFACE]);
axis tight;
str={sprintf('nVertex: %i',nv);
     sprintf('nTets:   %i',nt);
     sprintf('nTris:   %i',nbe)};
annotation('textbox',[0.05 0.05 0.2 0.15],'String',str, ...
           'FitBoxToText','on','FontName','Courier','EdgeColor',[1 1 1]);
ylabel('y');
xlabel('x');
zlabel('z');
lighting gouraud;
view([-47,24]);
camlight('headlight');

figure;
ffpdeplot3D(p,b,t,'FlowData',[Bx,By,Bz],'FGridParam3D',[15,15,15],'Boundary','on', ...
            'BDLabels',PHYSTORUSSUFACE,'XYZStyle','monochrome','FMode3D','random');
ylabel('y');
xlabel('x');
zlabel('z');
grid;
axis tight;
lighting gouraud;
view([-47,24]);
camlight('headlight');

S1=[0.0 -1.3 0.0];
S2=[0.0 -1.3 2.0];
S3=[0.0 1.3 0.0];
figure;
ffpdeplot3D(p,b,t,'FlowData',[Bx,By,Bz],'Slice',S1,S2,S3,'BoundingBox','on', ...
            'BDLabels',PHYSTORUSSUFACE,'XYZStyle','monochrome','FGridParam',[25,25]);
ylabel('y');
xlabel('x');
zlabel('z');
axis tight off;
lighting gouraud;
view([-47,24]);
camlight('headlight');

S1=[0.0 -1.3 0.0];
S2=[0.0 -1.3 2.0];
S3=[0.0 1.3 0.0];
figure;
ffpdeplot3D(p,b,t,'FlowData',[Jx,Jy,Jz],'Slice',S1,S2,S3,'BoundingBox','on', ...
            'Boundary','off','FGridParam',[50,50]);

S1=[0.0 -1.3 0.0; -1.3 1.3 1];
S2=[0.0 -1.3 2.0; -1.3 -1.3 1];
S3=[0.0 1.3 0.0; 1.3 1.3 1];
figure;
ffpdeplot3D(p,b,t,'XYZData',u,'Slice',S1,S2,S3,'SGridParam',[50,50], ...
            'Boundary','off','ColorMap',jet(150),'ColorRange','centered', ...
            'Colorbar','on','BoundingBox','on');
ylabel('y');
xlabel('x');
zlabel('z');
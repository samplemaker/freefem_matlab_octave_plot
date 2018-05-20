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

lcol=[153 101 21] ./ 255;

figure;
fileID = fopen('tetrahedrondata.txt','r');
fdata = textscan(fileID,'%f %f %f %f %d','Delimiter',';');
[x,y,z,c,no] = fdata{:};
%scatter3(x,y,z,10,'filled');
plot3(x,y,z,'.','MarkerSize',25);

hold on;
i=1;
for j=1:(numel(no)/4)
  plot3([x(i) x(i+1)],[y(i) y(i+1)],[z(i) z(i+1)],'Color',lcol);
  plot3([x(i+1) x(i+2)],[y(i+1) y(i+2)],[z(i+1) z(i+2)],'Color',lcol);
  plot3([x(i+2) x(i+3)],[y(i+2) y(i+3)],[z(i+2) z(i+3)],'Color',lcol);
  plot3([x(i+3) x(i)],[y(i+3) y(i)],[z(i+3) z(i)],'Color',lcol);
  i=i+4;
end

for i=1:numel(no)
  text(x(i)+0.1,y(i),z(i), ...
       [num2str(no(i))], ...
       'HorizontalAlignment','left','FontSize',11);
end
title('all nodes + numbering');
daspect([1 1 1]);

figure;
fileID = fopen('boundarynodes.txt','r');
fdata = textscan(fileID,'%f %f %f %f %d','Delimiter',';');
[x,y,z,c,no] = fdata{:};
%scatter3(x,y,z,10,'filled');
plot3(x,y,z,'.','MarkerSize',25);

hold on;
i=1;
for j=1:(numel(no)/3)
  plot3([x(i) x(i+1)],[y(i) y(i+1)],[z(i) z(i+1)],'Color',lcol);
  plot3([x(i+1) x(i+2)],[y(i+1) y(i+2)],[z(i+1) z(i+2)],'Color',lcol);
  plot3([x(i+2) x(i)],[y(i+2) y(i)],[z(i+2) z(i)],'Color',lcol);
  i=i+3;
end

for i=1:numel(no)
  text(x(i)+0.1,y(i),z(i), ...
       [num2str(no(i))], ...
       'HorizontalAlignment','left','FontSize',11);
end
title('boundary nodes only');
daspect([1 1 1]);

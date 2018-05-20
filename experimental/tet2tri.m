%convert a tetrahedron into four triangles
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

function [X1 X2 X3 X4 Y1 Y2 Y3 Y4 Z1 Z2 Z3 Z4 C1 C2 C3 C4] = tet2tri (tx,ty,tz,tc)
  %first triangle
  X1=[tx(1);tx(2);tx(3)];
  Y1=[ty(1);ty(2);ty(3)];
  Z1=[tz(1);tz(2);tz(3)];
  C1=[tc(1);tc(2);tc(3)];

  %2nd triangle
  X2=[tx(1);tx(2);tx(4)];
  Y2=[ty(1);ty(2);ty(4)];
  Z2=[tz(1);tz(2);tz(4)];
  C2=[tc(1);tc(2);tc(4)];

  %3nd triangle
  X3=[tx(1);tx(3);tx(4)];
  Y3=[ty(1);ty(3);ty(4)];
  Z3=[tz(1);tz(3);tz(4)];
  C3=[tc(1);tc(3);tc(4)];

  %4th triangle
  X4=[tx(2);tx(3);tx(4)];
  Y4=[ty(2);ty(3);ty(4)];
  Z4=[tz(2);tz(3);tz(4)];
  C4=[tc(2);tc(3);tc(4)];
end

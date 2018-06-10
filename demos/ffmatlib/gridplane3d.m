%gridplane3d.c Returns a rectangular 2D grid plane in 3D space defined by three points.
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-05-13
%
%   [X,Y,Z] = gridplane3d (SO,SN,SM,N,M) creates a rectangular grided
%   parallelogram which is defined by the corners SO,SN,SM. The grid size
%   is N,M. The distance between SN and SO is divided into N grids.
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

function [X,Y,Z] = gridplane3d (SO,SN,SM,N,M)
    SO=colvec(SO);
    SN=colvec(SN);
    SM=colvec(SM);
    X=NaN(M,N);
    Y=NaN(M,N);
    Z=NaN(M,N);
    for i=1:N
          v1=(i-1)/(N-1);
          %[x,y,z]=(SO+vN*(SN-SO)+vM*(SM-SO));
          X(:,i)=linspace(SO(1)+v1*(SN(1)-SO(1)),SO(1)+v1*(SN(1)-SO(1))+(SM(1)-SO(1)),M);
          Y(:,i)=linspace(SO(2)+v1*(SN(2)-SO(2)),SO(2)+v1*(SN(2)-SO(2))+(SM(2)-SO(2)),M);
          Z(:,i)=linspace(SO(3)+v1*(SN(3)-SO(3)),SO(3)+v1*(SN(3)-SO(3))+(SM(3)-SO(3)),M);
    end
end

function [S] = colvec(S)
    if size(S,2)>1
        S=S(:);
    end
end
%fftri2gridcplx.m Interpolates from 2D triangular mesh to 2D complex grid
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-11-11
%
%   [u] = ffplottri2grid (z, tx, ty, tu) interpolates the complex data [tu]
%   which is given on a triangular mesh defined by tx, ty onto a complex
%   curved meshgrid defined by z. tx, ty, tu must have a size of 3xnTriangle. 
%   The return value [u] is the interpolation at the grid points
%   z. Returns NaN's if an interpolation point is outside the triangle mesh.
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
function [u] = fftri2gridcplx(z, tx, ty, tu)
[ny,nx]=size(z);
x=real(z);
y=imag(z);
ax=tx(1,:);
ay=ty(1,:);
bx=tx(2,:);
by=ty(2,:);
cx=tx(3,:);
cy=ty(3,:);
invA0=(1.0)./((by-cy).*(ax-cx)+(cx-bx).*(ay-cy));
u=NaN(ny,nx);
for mx=1:nx
    for my=1:ny
        px=x(my,mx);
        py=y(my,mx);
        Aa=((by-cy).*(px-cx)+(cx-bx).*(py-cy)).*invA0;
        Ab=((cy-ay).*(px-cx)+(ax-cx).*(py-cy)).*invA0;
        Ac=1.0-Aa-Ab;
        pos=find(((Aa>=-1e-13) & (Ab>=-1e-13) & (Ac>=-1e-13)),1,'first');
        if ~isempty(pos)
            u(my,mx)=Aa(pos).*tu(1,pos)+ ...
                     Ab(pos).*tu(2,pos)+ ...
                     Ac(pos).*tu(3,pos);
        end
    end
end
end

%ffinterpolate.m Interpolates PDE simulation data on a meshgrid
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-11-11
%
%   Wrapper function
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
function [w] = ffinterpolate(p, b, t, x, y, u)
    xpts=p(1,:);
    ypts=p(2,:);
    xdata=[xpts(t(1,:)); xpts(t(2,:)); xpts(t(3,:))];
    ydata=[ypts(t(1,:)); ypts(t(2,:)); ypts(t(3,:))];
    udata=[u(t(1,:)), u(t(2,:)), u(t(3,:))].';
    if exist('fftri2meshgrid','file')
        w = fftri2meshgrid(x,y,xdata,ydata,udata);
    else
        fprintf('Note: To improve runtime build MEX function fftri2meshgrid() from fftri2meshgrid.c\n');
        w = fftri2meshgridint(x,y,xdata,ydata,udata);
    end
 end

%Interpolates from 2D triangular mesh to 2D mesh grid
%[u] = fftri2meshgridint (x, y, tx, ty, tu) interpolates the complex data [tu]
%which is given on a triangular mesh defined by tx, ty onto a
%meshgrid defined by x, y. tx, ty, tu must have a size of 3xnTriangle. 
%The return value [u] is the interpolation at the grid points
%x,y. Returns NaN's if an interpolation point is outside the triangle mesh.
function [u] = fftri2meshgridint(x, y, tx, ty, tu)
    if ~isequal(size(x), size(y))
        error('meshgrid sizes must be equal');
    end
    [ny,nx]=size(x);
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

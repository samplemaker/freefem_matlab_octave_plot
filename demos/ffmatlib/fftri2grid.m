%fftri2grid.m Interpolates from 2D triangular mesh to a 2D cartesian or curved grid
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-11-11
%
%   This file is part of the ffmatlib which is hosted at
%   https://github.com/samplemaker/freefem_matlab_octave_plot
%
%   [w1, [w2]] = fftri2grid (x, y, tx, ty, tu1, [tu2])
%
%   interpolates the real or complex data tu1, tu2 given on a triangular mesh
%   defined by the two arguments tx, ty onto a meshgrid defined by the
%   variables x, y. The mesh can be cartesian or curved. The argument
%   tu1, tu2 must have a size of nTriangle columns and 3 rows. The return
%   value w1, w2 is the interpolation of tu1, tu2 at the grid points defined
%   by x, y. The result w1, w2 is real if tu1, tu2 is real or complex if
%   tu1, tu2 is complex. fftri2gridfast.c uses barycentric coordinates to
%   interpolate. Returns NaN's if an interpolation point is outside the
%   triangle mesh. The argument tu2 is optional.
%   fftri2gridfast.c is the mex implementation of the function
%   fftri2grid.c.
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

function [u1, u2] = fftri2grid(x, y, tx, ty, tu1, tu2)
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
    if (nargin == 5)
       u1=NaN(ny,nx);
       for mx=1:nx
           for my=1:ny
               px=x(my,mx);
               py=y(my,mx);
               Aa=((by-cy).*(px-cx)+(cx-bx).*(py-cy)).*invA0;
               Ab=((cy-ay).*(px-cx)+(ax-cx).*(py-cy)).*invA0;
               Ac=1.0-Aa-Ab;
               pos=find(((Aa>=-1e-13) & (Ab>=-1e-13) & (Ac>=-1e-13)),1,'first');
               if ~isempty(pos)
                   u1(my,mx)=Aa(pos).*tu1(1,pos)+ ...
                             Ab(pos).*tu1(2,pos)+ ...
                             Ac(pos).*tu1(3,pos);
               end
           end
       end
     else
       u1=NaN(ny,nx);
       u2=NaN(ny,nx);
       for mx=1:nx
           for my=1:ny
               px=x(my,mx);
               py=y(my,mx);
               Aa=((by-cy).*(px-cx)+(cx-bx).*(py-cy)).*invA0;
               Ab=((cy-ay).*(px-cx)+(ax-cx).*(py-cy)).*invA0;
               Ac=1.0-Aa-Ab;
               pos=find(((Aa>=-1e-13) & (Ab>=-1e-13) & (Ac>=-1e-13)),1,'first');
               if ~isempty(pos)
                   u1(my,mx)=Aa(pos).*tu1(1,pos)+ ...
                             Ab(pos).*tu1(2,pos)+ ...
                             Ac(pos).*tu1(3,pos);
                   u2(my,mx)=Aa(pos).*tu2(1,pos)+ ...
                             Ab(pos).*tu2(2,pos)+ ...
                             Ac(pos).*tu2(3,pos);
               end
           end
       end
     end
end
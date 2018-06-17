%fftri2grid.m Interpolates from 2D triangular mesh to 2D rectangular grid
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-05-13
%
%   [varargout] = fftri2grid (tridata, X, Y) interpolates data given on
%   a triangular mesh to a rectangular grid defined by X and Y. The columns
%   tridata(:,1) and tridata(:,2) must contain the triangular mesh node
%   coordinates. The following columns must contain the scalar values at
%   the node points that need to be interpolated.
%   The return value is the interpolation at the grid points X, Y. Returns
%   NaN's if an interpolation point is outside the triangle mesh.
%
%   Hint: We evaluate the PDE solution only on the grid vertices although the
%   underlying FE space may have a higher order (P2 element, etc.).
%   Therefore, there is a small loss of accuracy except P1 elements are used.
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
function [varargout] = fftri2grid(tridata, X, Y)
    switch nargin
        case {3}
        otherwise
            printhelp();
            error('wrong number arguments');
    end
    [npts,nvars]=size(tridata);
    if mod(npts,3) ~= 0
        printhelp();
        error('number of triangles not a multiple of 3');
    end
    if nvars < 3
        printhelp();
        error('wrong number of columns - must be >= 3');
    end
    %Splitting into triangles
    tx=arrangecols(tridata(:,1),3);
    ty=arrangecols(tridata(:,2),3);
    uvarin=cell(1,nvars-2);
    varargout=cell(1,nvars-2);
    for i=1:nvars-2
        uvarin{i}=arrangecols(tridata(:,i+2),3);
        varargout{i}=NaN(numel(Y),numel(X));
    end
    %Making copies saves 50% of running time instead of using
    %tx(1,:) directly
    ax=tx(1,:); %x values of the first triangle point for all triangles
    ay=ty(1,:);
    bx=tx(2,:); %x values of the second triangle point for all triangles
    by=ty(2,:);
    cx=tx(3,:);
    cy=ty(3,:);
    invA0=(1.0)./((by-cy).*(ax-cx)+(cx-bx).*(ay-cy));
    for mx=1:numel(X)
        for my=1:numel(Y)
            px=X(mx);
            py=Y(my);
            Aa=((by-cy).*(px-cx)+(cx-bx).*(py-cy)).*invA0;
            Ab=((cy-ay).*(px-cx)+(ax-cx).*(py-cy)).*invA0;
            Ac=1.0-Aa-Ab;
            %Is inside which triangle?
            %Can possibly sit on a border (multiple output)
            pos=find(((Aa>=0) & (Ab>=0) & (Ac>=0)),1,'first');
            %Out of triangle: No else because varargout contains already NaN's
            if ~isempty(pos)
                for i=1:nvars-2
                    varargout{i}(my,mx)=Aa(pos).*uvarin{i}(1,pos)+ ...
                                        Ab(pos).*uvarin{i}(2,pos)+ ...
                                        Ac(pos).*uvarin{i}(3,pos);
                end
            end
        end
    end
end

function [M] = arrangecols(V,c)
    r = length(V)/c;
    M = reshape(V,c,r);
end

function printhelp()
    fprintf('%s\n\n','Invalid call to fftri2grid.  Correct usage is:');
    fprintf('%s\n',' -- [A,B, ...] = fftri2grid (tridata, X, Y)');
end

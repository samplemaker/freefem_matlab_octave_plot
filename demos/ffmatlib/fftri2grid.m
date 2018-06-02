%fftri2grid.m Interpolate from 2d triangular mesh to 2d rectangular grid
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-05-13
%
%   [varargout] = fftri2grid (tridata, X, Y) interpolates data given on
%   a triangular mesh to a rectangular mesh defined by X and Y. The columns
%   tridata(:,1) and tridata(:,2) must contain the triangular mesh node
%   coordinates. The following columns must contain the scalar values ​​at
%   the node points that need to be interpolated.
%   The return value is the interpolation at the grid points X, Y. Returns
%   NaN's if an interpolation point is outside the triangle mesh.
%
%   Hint: We use the PDE solution only on the grid vertices, although the
%   underlying FE space may have a higher order (P2 element, etc.).
%   Therefore, there is a small loss of accuracy, except P1 elements are used
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
        error('wrong number of columns - must be > 3');
    end
    %Splitting into triangles
    tx=arrangecols(tridata(:,1),3);
    ty=arrangecols(tridata(:,2),3);
    ntriangles=npts/3;
    tu=zeros(3*nvars,ntriangles);
    j=0;
    for i=3:nvars
        %tu=arrangecols(tridata(:,3),3);
        %ua=tu(1,:);
        %ub=tu(2,:);
        %uc=tu(3,:);
        tu(1+j:3+j,:)=arrangecols(tridata(:,i),3);
        j=j+3;
        %C=NaN(numel(Y),numel(X));        
        strc=sprintf('C%1.0f',i);
        vars.(strc)=NaN(numel(Y),numel(X));
    end
    %Making copies saves 50% of running time instead of using
    %tx(1,:) directly
    ax=tx(1,:); %x values ​​of the first triangle point for all triangles
    ay=ty(1,:);
    bx=tx(2,:); %x values ​​of the second triangle point for all triangles
    by=ty(2,:);
    cx=tx(3,:);
    cy=ty(3,:);
    for mx=1:numel(X)
        for my=1:numel(Y)
            px=X(mx);
            py=Y(my);
            %Calculates barycentric coordinates
            fac=(1.0)./((by-cy).*(ax-cx)+(cx-bx).*(ay-cy));
            wa=((by-cy).*(px-cx)+(cx-bx).*(py-cy)).*fac;
            wb=((cy-ay).*(px-cx)+(ax-cx).*(py-cy)).*fac;
            wc=1.0-wa-wb;
            %Is inside which triangle?
            %Can possibly sit on a border (multiple output)
            pos=find(((wa>=0) & (wb>=0) & (wc>=0)),1,'first');
            if ~isempty(pos)
                j=0;
                for i=3:nvars
                    %C(my,mx)=wa(pos).*ua(pos)+wb(pos).*ub(pos)+wc(pos).*uc(pos);
                    strc=sprintf('C%1.0f',i);
                    vars.(strc)(my,mx)=wa(pos).*tu(1+j,pos)+ ...
                                       wb(pos).*tu(2+j,pos)+ ...
                                       wc(pos).*tu(3+j,pos);
                    j=j+3;
                end
            end
        end
    end
    varargout=cell(1,nvars-2);
    for i = 1:(nvars-2)
        strc=sprintf('C%1.0f',i+2);
        varargout{i} = vars.(strc);
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

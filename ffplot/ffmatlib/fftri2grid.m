%fftri2grid.m Interpolate from 2d triangular mesh to 2d rectangular grid
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-05-13
%
%   [varargout] = fftri2grid (tridata, X, Y) interpolates on a
%   rectangular grid defined by X and Y. tridata(:,1) and tridata(:,2)
%   specify the triangle coordinates and following columns scalar value
%   on the vertices respectively.
%   The return value is the interpolation at X,Y.
%   Returns NaN's if a interpolation point lies outside
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
    tx=arrangecols(tridata(:,1),3); %want to calculate triangle and point wise
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
    %saves 50% runtime to create copies instead of using indexing tx(1,:)
    %again and again
    ax=tx(1,:); %x werte des ersten punktes für alle dreiecke
    ay=ty(1,:);
    bx=tx(2,:); %x werte des zweiten punktes für alle dreiecke
    by=ty(2,:);
    cx=tx(3,:);
    cy=ty(3,:);
    for mx=1:numel(X)
        for my=1:numel(Y)
            px=X(mx);
            py=Y(my);
            %calculate barycentric coordinates
            fac=(1.0)./((by-cy).*(ax-cx)+(cx-bx).*(ay-cy));
            wa=((by-cy).*(px-cx)+(cx-bx).*(py-cy)).*fac;
            wb=((cy-ay).*(px-cx)+(ax-cx).*(py-cy)).*fac;
            wc=1.0-wa-wb;
            %is inside which triangle?
            pos=find(((wa>=0) & (wb>=0) & (wc>=0)),1,'first');%can possibly sit on a border
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


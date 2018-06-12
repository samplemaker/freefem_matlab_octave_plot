%fftet2grid.c Interpolates from 3d tetrahedral mesh to a rectangular grid
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-05-13
%
%   [varargout] = fftet2grid(tetdata,X,Y,Z) interpolates data given on
%   a tetrahedral mesh to a rectangular grid defined by X, Y and Z. The first
%   three columns in tdata must contain the tetrahedral mesh node
%   coordinates. The following columns must contain the scalar values at
%   the four tetrahedra node points that need to be interpolated.
%   The return value is the interpolation at the grid points defined by X, Y
%   and Z. Returns NaN's if an interpolation point is outside the tetrahedral mesh.
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

function [varargout] = fftet2grid(tetdata,X,Y,Z)
    switch nargin
        case {4}
        otherwise
            printhelp();
            error('4 input arguments required');
    end
    [npts,nvars]=size(tetdata);
    if mod(npts,4) ~= 0
        printhelp();
        error('number of tet points not a multiple of 4');
    end
    if nvars < 4
        printhelp();
        error('wrong number of columns - must be > 3');
    end

     mesh=tetdata(:,1:3);
     a = mesh(1:4:end,:);
     b = mesh(2:4:end,:);
     c = mesh(3:4:end,:);
     d = mesh(4:4:end,:);

    [M,N]=size(X);
    [ntet,~]=size(a);

    varargout=cell(1,nvars-3);
    u=cell(1,nvars-3);
    for i=1:nvars-3
        tmp=tetdata(:,i+3);
        u{i}=reshape([tmp(1:4:end,:);tmp(2:4:end,:);tmp(3:4:end,:);tmp(4:4:end,:)],ntet,4);
        varargout{i}=NaN(M,N);
    end

    invVTET=abs(1.0./(dot(cross(c-a,b-a),d-a,2)));

    for i=1:N
        for j=1:M
            xp=[X(j,i), Y(j,i), Z(j,i)];
            Xp=repmat(xp,ntet,1);
            Va=dot(cross(d-b,c-b),Xp-b,2);
            Vb=dot(cross(c-a,d-a),Xp-a,2);
            Vc=dot(cross(d-a,b-a),Xp-a,2);
            Vd=dot(cross(b-a,c-a),Xp-a,2);
            pos=find(((Va>=0) & (Vb>=0) & (Vc>=0) & (Vd>=0)),1,'first');
            if ~isempty(pos)
                for k=1:nvars-3
                    varargout{k}(j,i)=((u{k}(pos,1).*Va(pos)+ ...
                                        u{k}(pos,2).*Vb(pos)+ ...
                                        u{k}(pos,3).*Vc(pos)+ ...
                                        u{k}(pos,4).*Vd(pos))).*invVTET(pos);
                end
            end
        end
    end
end

function printhelp()
    fprintf('%s\n\n','Invalid call to fftet2grid. Correct usage is:');
    fprintf('%s\n',' -- [X,Y,Z,C, ...] = fftet2grid(tetdata,X,Y,Z)');
end

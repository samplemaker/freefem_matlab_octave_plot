%slicetet2tet.m Cuts 3D mesh elements (tetrahedrons) and returns
%               all affected tetrahedrons
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-05-13
%
%   [varargout] = slicetet2tet (tdata,S1,S2,S3) The mesh data is cut
%   by a cutting plane defined by the three points S1, S2, S3. tdata must
%   contain the coordinates in the first three columns and in the following
%   colums scalar values.
%
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
function [varargout] = slicetet2tet(fdata,S1,S2,S3)
    switch nargin
        case {4}
        otherwise
            printhelp();
            error('4 input arguments required');
    end
    [npts,nvars]=size(fdata);
    if mod(npts,4) ~= 0
        printhelp();
        error('number of tet points not a multiple of 4');
    end
    if nvars < 4
        printhelp();
        error('wrong number of columns - must be > 3');
    end

    S1=colvec(S1);
    S2=colvec(S2);
    S3=colvec(S3);

% Theory
%
% let Xn:=cross((S2-S1),(S3-S1)) be perpendicular to the slicing plane
% and Xp a point in the plane. it turns out that for any point X,
% i.)   in the plane          --> dot(N,(X-Xp)) == 0
% ii.)  in front of the plane --> dot(N,(X-Xp)) > 0
% iii.) behind of the plane   --> dot(N,(X-Xp)) < 0
%

    %Nodal coordinates
    M=fdata(:,1:3);
    [npts,~]=size(M);
    Xn=cross((S2-S1),(S3-S1));
    Xn0=repmat(Xn',npts,1);
    S10=repmat(S1',npts,1);
    pos=arrangecols(dot(Xn0,(M-S10),2),4);
    %Tetrahedrons that have one or more points in front of and one or more
    %points behind the cutting plane must be identified. Tetrahedra with one
    %or more points exactly in the plane do count as well
    front=any((pos>0));
    behind=any((pos<0));
    isin=any((pos==0));
    %Boolarray indicating sliced or touched tetrahedra
    keep=(((behind==1) & (front==1)) | (isin==1));
    %Extracts all affected triangles
    varargout=cell(1,nvars);
    for i=1:nvars
      tmp=arrangecols(fdata(:,i),4);
      varargout{i}=tmp(:,keep);
    end
end

function [M] = arrangecols(V,c)
    r = length(V)/c;
    M = reshape(V,c,r);
end

function [S] = colvec(S)
    if size(S,2)>1
        S=S(:);
    end
end

function printhelp()
    fprintf('%s\n\n','Invalid call to slicetet2tet. Correct usage is:');
    fprintf('%s\n',' -- [X, Y, Z, C] = slicetet2tet (tdata, S1, S2, S3)');
end

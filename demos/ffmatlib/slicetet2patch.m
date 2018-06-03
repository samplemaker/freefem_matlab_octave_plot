%slicetet2patch.m Cuts 3D mesh elements (tetrahedrons) and converts the
%                 cross section into patch plot data.
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-05-13
%
%   [X,Y,Z,C] = slicetet2patch (tdata,S1,S2,S3) The mesh data is cut
%   by a cutting plane defined by the three points S1, S2, S3. tdata must
%   contain the coordinates in the first three columns and in the forth
%   column a scalar value. 
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
function [SX,SY,SZ,SC] = slicetet2patch(fdata,S1,S2,S3)
    switch nargin
        case {4}
        otherwise
            printhelp();
            error('wrong number arguments');
    end

    S1=colvec(S1);
    S2=colvec(S2);
    S3=colvec(S3);

%%%%%%% Theory
%
% let Xn:=cross((S2-S1),(S3-S1)) be perpendicular to the slicing plane
% and Xp a point in the plane. it turns out that for any point X,
% i.)   in the plane          --> dot(N,(X-Xp)) == 0
% ii.)  in front of the plane --> dot(N,(X-Xp)) > 0
% iii.) behind of the plane   --> dot(N,(X-Xp)) < 0
%

    %Finds tetrahedra cut or touched by the cutting plane
    x=fdata(:,1);y=fdata(:,2);z=fdata(:,3);c=fdata(:,4);
    M=[x y z];
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
    %Extracts all affected triangles and removes the rest
    tmpX=arrangecols(x,4);
    TX=tmpX(:,keep);
    tmpY=arrangecols(y,4);
    TY=tmpY(:,keep);
    tmpZ=arrangecols(z,4);
    TZ=tmpZ(:,keep);
    tmpC=arrangecols(c,4);
    TC=tmpC(:,keep);
    %Extracts from each tetrahedron 4 triangles, don't bother about duplicates
    SX=tet2tri(TX);
    SY=tet2tri(TY);
    SZ=tet2tri(TZ);
    SC=tet2tri(TC);
    %eliminate duplicates
    %T=[SX;SY;SZ;SC]'; %Combine in order to find identical triangles
    %[C,ia,ic]=unique(T(:,1:9),'rows'); %Find the unique rows of A based on the data in the first 9 cols
    %TU=T(ia,:); %Use ia to index into T and retrieve the rows
    %SX=TU(:,1:3)';
    %SY=TU(:,4:6)';
    %SZ=TU(:,7:9)';
    %SC=TU(:,10:12)';
end

%Converts tetrahedrons into vertex/triangles
function [Y] = tet2tri(X)
    [~,sz2] = size(X);
    M=[[X(1,:); X(2,:); X(3,:)]; ...
      [X(1,:); X(2,:); X(4,:)]; ...
      [X(1,:); X(3,:); X(4,:)]; ...
      [X(2,:); X(3,:); X(4,:)]];
    Y=reshape(M,3,4*sz2);
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
    fprintf('%s\n\n','Invalid call to slicetet2patch. Correct usage is:');
    fprintf('%s\n',' -- [X, Y, Z, C] = slicetet2patch (tdata, S1, S2, S3)');
end

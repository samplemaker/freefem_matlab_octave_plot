%slicetet2patch.m Slice 3d mesh elements (tetrahedra) and convert to
%                 patch plot data. The mesh is sliced by a slicing
%                 plane defined by the three points S1,S2,S3.
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-05-13
%
%   [X,Y,Z,C] = slicetet2patch (tdata,S1,S2,S3)
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%   Theory   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% let Xn:=cross((S2-S1),(S3-S1)) be perpendicular to the slicing plane
% and Xp a point in the plane. it turns out that for any point X,
% i.)   in the plane          --> dot(N,(X-Xp)) == 0
% ii.)  in front of the plane --> dot(N,(X-Xp)) > 0
% iii.) behind of the plane   --> dot(N,(X-Xp)) < 0
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %find out tetrahedra which are cut or touched by the slicing plane
    %[x,y,z,c]=fdata{:};
    x=fdata(:,1);y=fdata(:,2);z=fdata(:,3);c=fdata(:,4);
    M=[x y z];
    [npts,~]=size(M);
    Xn=cross((S2-S1),(S3-S1));
    Xn0=repmat(Xn',npts,1);
    S10=repmat(S1',npts,1);
    pos=arrangecols(dot(Xn0,(M-S10),2),4);
    %tetrahedra which do have one or more points in front of, as well
    %as one or more points behind the slicing plane have to be identified.
    %tetrahedra with a point in the plane do count as well
    front=any((pos>0));
    behind=any((pos<0));
    isin=any((pos==0));
    %boolarray identifing sliced or touched tetrahedra
    keep=(((behind==1) & (front==1)) | (isin==1));
    %extract all affected tetrahedra and remove the chunck
    tmpX=arrangecols(x,4);
    TX=tmpX(:,keep);
    tmpY=arrangecols(y,4);
    TY=tmpY(:,keep);
    tmpZ=arrangecols(z,4);
    TZ=tmpZ(:,keep);
    tmpC=arrangecols(c,4);
    TC=tmpC(:,keep);
    %extract from each tetrahedron 4 triangles, don't bother about duplicates
    SX=tet2tri(TX);
    SY=tet2tri(TY);
    SZ=tet2tri(TZ);
    SC=tet2tri(TC);
end

%convert tetrahedrons into vertex/triangles
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

function printhelp()
    fprintf('%s\n\n','Invalid call to slicetet2patch. Correct usage is:');
    fprintf('%s\n',' -- [X, Y, Z, C] = slicetet2patch (tdata, S1, S2, S3)');
end

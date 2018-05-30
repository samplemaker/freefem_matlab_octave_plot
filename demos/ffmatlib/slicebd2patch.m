%slicebd2patch.m Slice 3d boundary (triangle) data and convert to
%                patch plot data. The data is sliced by a slicing
%                plane defined by the three points S1,S2,S3.
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-05-13
%
%   [X,Y,Z,C] = slicebd2patch (data, S1, S2, S3)
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
function [BX,BY,BZ,BC] = slicebd2patch(fdata,S1,S2,S3)
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

    %proceed with the boundary triangles
    %[x,y,z,c]=fdata{:};
    x=fdata(:,1);y=fdata(:,2);z=fdata(:,3);c=fdata(:,4);
    M=[x y z];
    [npts,~]=size(M);
    Xn=cross((S2-S1),(S3-S1));
    Xn0=repmat(Xn',npts,1);
    S10=repmat(S1',npts,1);
    %used to check which points are in front of the plane
    pos=dot(Xn0,(M-S10),2);
    %bool array indicating affected points
    keep=(any(arrangecols((pos<=0),3))==1);
    %extract all affected triangle and remove the chunck
    tmpX=arrangecols(x,3);
    BX=tmpX(:,keep);
    tmpY=arrangecols(y,3);
    BY=tmpY(:,keep);
    tmpZ=arrangecols(z,3);
    BZ=tmpZ(:,keep);
    tmpC=arrangecols(c,3);
    BC=tmpC(:,keep);
end

function [M] = arrangecols(V,c)
    r = length(V)/c;
    M = reshape(V,c,r);
end

function printhelp()
    fprintf('%s\n\n','Invalid call to slicebd2patch. Correct usage is:');
    fprintf('%s\n',' -- [X, Y, Z, C] = slicebd2patch (bdata, S1, S2, S3)');
end

%ffslicebd3d.m Returns all boundary triangles affected if a 3d mesh is sliced with
%            a slicing plane defined by the three points S1,S2,S3
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-05-13
%
%   [XX YY ZZ CC] = ffslicebd3d (bdfilename,S1,S2,S3,varargin)
%
%   [XX YY ZZ CC] = ffslicebd3d (...,'PARAM1',val1,'PARAM2',val2,...) specifies
%   parameter name/value pairs to control the input file format
%
%       Parameter       Value
%      'Delimiter'     Specifies the delimiter of the input file
%                      (default='\t')
%      'Format'        Specifies the column format of the input file
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
function [BXX BYY BZZ BCC] = ffslicebd3d(bdfile,S1,S2,S3,varargin)
  switch nargin
    case {4,6,8}
    otherwise
      printhelp();
      error('wrong number arguments');
  end
  numvarargs = length(varargin);
  optsnames = {'Delimiter','Format'};
  varargout = {'\t'};
  switch numvarargs
    case 0
    case {2,4}
      optargs(1:numvarargs) = varargin;
      for j = 1:length(optsnames)
        for i = 1:length(varargin)
          if strcmp(optsnames(j),optargs(i))
            varargout(j) = optargs(i+1);
          end
        end
      end
    otherwise
      printhelp()
      error('wrong number arguments');
  end
  [delimiter, format] = varargout{:};

  bdfid = fopen(bdfile,'r');
  if bdfid < 0
    error('cannot open boundary file');
  end

      %%%%%%%%%%%%%%%%%%%%%%%%%%   Theory   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% let Xn:=cross((S2-S1),(S3-S1)) be perpendicular to the slicing plane
% and Xp a point in the plane. it turns out that for any point X,
% i.)   in the plane          --> dot(N,(X-Xp)) == 0
% ii.)  in front of the plane --> dot(N,(X-Xp)) > 0
% iii.) behind of the plane   --> dot(N,(X-Xp)) < 0
  
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %proceed with the boundary triangles
  fdata=textscan(bdfid,format,'Delimiter',delimiter);
  fclose(bdfid);
  [x,y,z,c]=fdata{:};
  M=[x y z];
  [npts sz2]=size(M);
  Xn=cross((S2-S1),(S3-S1));
  Xn0=repmat(Xn',npts,1);
  S10=repmat(S1',npts,1);
  %used to check which points are in front of the plane
  pos=dot(Xn0,(M-S10),2);
  %bool array indicating affected points
  keep=(any(arrangecols((pos<=0),3))==1);
  nkeep=sum(keep==1);
  %extract all affected triangle and remove the chunck
  tmpX=arrangecols(x,3);
  BXX=tmpX(:,keep);
  tmpY=arrangecols(y,3);
  BYY=tmpY(:,keep);
  tmpZ=arrangecols(z,3);
  BZZ=tmpZ(:,keep);
  tmpC=arrangecols(c,3);
  BCC=tmpC(:,keep);
end

%convert tetrahedrons into vertex/triangles
function [Y] = tet2tri(X)
  [sz1 sz2] = size(X);
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
  fprintf('%s\n\n','Invalid call to ffslicebd3d. Correct usage is:');
  fprintf('%s\n',' -- [X, Y, Z, C] = ffslicebd3d (bdfilename, S1, S2, S3, ''Delimiter'','';'',''Format'',''%f %f %f %f'')');
end
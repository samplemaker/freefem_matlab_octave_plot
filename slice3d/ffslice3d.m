%ffslice3d.m Returns all triangles affected if a 3d mesh is sliced with
%            a slicing plane defined by the three points S1,S2,S3
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-05-13
%
%   [XX YY ZZ CC] = ffslice3d (bdfilename,tetfilename,S1,S2,S3,varargin)
%
%   [XX YY ZZ CC] = ffslice3d (...,'PARAM1',val1,'PARAM2',val2,...) specifies
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
function [XX YY ZZ CC] = ffslice3d(bdfile,tetfile,S1,S2,S3,varargin)
  switch nargin
    case {5,7,9}
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

  tetfid = fopen(tetfile,'r');
  if tetfid < 0
    error('cannot open tetfile');
  end

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

  %step 1.)
  %find out tetrahedra which are cut or touched by the slicing plane
  fdata=textscan(tetfid,format,'Delimiter',delimiter);
  fclose(tetfid);
  [x,y,z,c]=fdata{:};
  M=[x y z];
  [npts sz2]=size(M);
  Xn=cross((S2-S1),(S3-S1));
  Xn0=repmat(Xn',npts,1);
  S10=repmat(S1',npts,1);
  pos=dot(Xn0,(M-S10),2);
  %all tetrahedra which do have one or more points in front of, as well
  %as one more points behind the slicing plane have to be identified.
  %tetrahedra with point in the point do count as well
  front=any(arrangecols((pos>0),4));
  behind=any(arrangecols((pos<0),4));
  ison=any(arrangecols((pos==0),4));
  %boolarray identifing sliced or touched tetrahedra
  keep=(((behind==1) & (front==1)) | (ison==1));
  nkeep=sum(keep==1);
  %extract all affected tetrahedra and remove the chunck
  tmpX=arrangecols(x,4);
  TXX=tmpX(:,keep);
  tmpY=arrangecols(y,4);
  TYY=tmpY(:,keep);
  tmpZ=arrangecols(z,4);
  TZZ=tmpZ(:,keep);
  tmpC=arrangecols(c,4);
  TCC=tmpC(:,keep);
  %we will store the triangles here
  SXX=zeros(3,4*nkeep);
  SYY=zeros(3,4*nkeep);
  SZZ=zeros(3,4*nkeep);
  SCC=zeros(3,4*nkeep);
  %extract from each tetrahedron 4 triangles, don't bother about duplicates
  k=1;
  for j=1:nkeep
      [SXX(:,k:k+3) SYY(:,k:k+3) ...
       SZZ(:,k:k+3) SCC(:,k:k+3)]= tet2tri(TXX(:,j),TYY(:,j), ...
                                           TZZ(:,j),TCC(:,j));  
      k=k+4;
  end
  
  %step 2.)
  %proceed with the boundary triangles
  fdata=textscan(bdfid,format,'Delimiter',delimiter);
  fclose(bdfid);
  [x,y,z,c]=fdata{:};
  M=[x y z];
  [npts sz2]=size(M);
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
  
  %append to the other stuff and return
  XX=[SXX BXX];
  YY=[SYY BYY];
  ZZ=[SZZ BZZ];
  CC=[SCC BCC];
end

function [M] = arrangecols(V,c)
  r = length(V)/c;
  M = reshape(V,c,r);
end

%by brut force. any suggestions??!?!
%return the vertex/triangle points from a tetrahedron
function [X Y Z C] = tet2tri (x,y,z,c)
  X=[x(1) x(1) x(1) x(2);
     x(2) x(2) x(3) x(3);
     x(3) x(4) x(4) x(4)];

  Y=[y(1) y(1) y(1) y(2);
     y(2) y(2) y(3) y(3);
     y(3) y(4) y(4) y(4)];
 
  Z=[z(1) z(1) z(1) z(2);
     z(2) z(2) z(3) z(3);
     z(3) z(4) z(4) z(4)];
     
  C=[c(1) c(1) c(1) c(2);
     c(2) c(2) c(3) c(3);
     c(3) c(4) c(4) c(4)];
end

function printhelp()
  fprintf('%s\n\n','Invalid call to ffslice3d. Correct usage is:');
  fprintf('%s\n',' -- [X, Y, Z, C] = ffslice3d (bdfilename, tetfilename, S1, S2, S3, ''Delimiter'','';'',''Format'',''%f %f %f %f'')');
end

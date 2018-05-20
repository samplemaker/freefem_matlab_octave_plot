%ffslice3d.m Returns all triangles affected if a 3d mesh is sliced with
%            a slicing plane defined by the three points Sx,Sy,Sz
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-05-13
%
%   [XX YY ZZ CC] = ffslice3d (bdfilename,tetfilename,Sx,Sy,Sz,varargin)
%
%   [XX YY ZZ CC] = ffslice3d (...,'PARAM1',val1,'PARAM2',val2,...) specifies parameter
%   name/value pairs to control the input file format
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
function [XX YY ZZ CC] = ffslice3d(bdfile,tetfile,Sx,Sy,Sz,varargin)
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
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %we starting to determine affected tetrahedra which are cut by the slicing plane
  
  fid = fopen(tetfile,'r');
  if fid < 0
    error('cannot open file');
  end
  fdata = textscan(fid,format,'Delimiter',delimiter);
  fclose(fid);
  [x,y,z,c] = fdata{:};
  M=[x y z];
  [sz1 sz2]=size(M);
  
  Ns=cross((Sy-Sx),(Sz-Sx));
  Nd=repmat(Ns',sz1,1);
  Sxd=repmat(Sx',sz1,1);
  %used to check which points are in front of the plane
  location=dot(Nd,(M-Sxd),2);
  %all tetrahedra which do have one more point in front of, as well as behind
  %the slicing plane have to be sorted out
  ptfront=(location>0);
  tfront=arrangecols(ptfront,4);
  tfront=any(tfront);
  %bool array indicating points
  pbehind=(location<0);
  tbehind=arrangecols(pbehind,4);
  %boolarray indication for tetrahedra
  tbehind=any(tbehind);
  %point on the plane
  pison=(location==0);
  tison=arrangecols(pison,4);
  tison=any(tison);
  %boolarray identifing sliced tetrahedra
  keep=(((tbehind==1) & (tfront==1)) | (tison==1));
  ntetrahedron=sum(keep==1);
  tx=arrangecols(x,4);
  ty=arrangecols(y,4);
  tz=arrangecols(z,4);
  tc=arrangecols(c,4);
  %split each tetrahedron into 4 triangles containing
  %3 nodes (vertices) eachBXX(:,k:k+3)=x(j:j+2)
  XX=zeros(3,4*ntetrahedron);
  YY=zeros(3,4*ntetrahedron);
  ZZ=zeros(3,4*ntetrahedron);
  CC=zeros(3,4*ntetrahedron);
  %remove the chunck
  [sz1 sz2] = size(tx);
  k=1;
  i=1;
  for j=1:sz2
    if keep(j)
      [XX(:,k:k+3) YY(:,k:k+3) ZZ(:,k:k+3) CC(:,k:k+3)]=tet2tri(tx(:,j),ty(:,j),tz(:,j),tc(:,j));  
      k=k+4;
    end
    i=i+4;
  end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %now proceed further with the boundary triangles
  
  fid = fopen(bdfile,'r');
  if fid < 0
    error('cannot open file');
  end
  fdata = textscan(fid,format,'Delimiter',delimiter);
  fclose(fid);
  [x,y,z,c] = fdata{:};
  M=[x y z];
  [sz1 sz2]=size(M);
  Nd=repmat(Ns',sz1,1);
  Sxd=repmat(Sx',sz1,1);
  %used to check which points are in front of the plane
  location=dot(Nd,(M-Sxd),2);
  %bool array indicating affected points
  pa=(location<=0);
  ta=arrangecols(pa,3);
  %boolarray indicating affected triangles
  ta=any(ta);
  keep=(ta==1);
  ntriangle=sum(keep==1);
  %remove the chunck
  BXX=zeros(3,ntriangle);
  BYY=zeros(3,ntriangle);
  BZZ=zeros(3,ntriangle);
  BCC=zeros(3,ntriangle);
  k=1;
  j=1;
  for i=1:(sz1/3)
    if keep(i)
      BXX(:,k)=x(j:j+2);
      BYY(:,k)=y(j:j+2);
      BZZ(:,k)=z(j:j+2);
      BCC(:,k)=c(j:j+2);
      k=k+3;
    end
    j=j+3;
  end
  
  %append to the other stuff
  XX=[XX BXX];
  YY=[YY BYY];
  ZZ=[ZZ BZZ];
  CC=[CC BCC];
end

function [M] = arrangecols(V,c)
  r = length(V)/c;
  M = reshape(V,c,r);
end

%return the vertex/triangle points from a tetrahedron
function [X Y Z C] = tet2tri (tx,ty,tz,tc)
  %first triangle
  X1=[tx(1);tx(2);tx(3)];
  Y1=[ty(1);ty(2);ty(3)];
  Z1=[tz(1);tz(2);tz(3)];
  C1=[tc(1);tc(2);tc(3)];
  %2nd triangle
  X2=[tx(1);tx(2);tx(4)];
  Y2=[ty(1);ty(2);ty(4)];
  Z2=[tz(1);tz(2);tz(4)];
  C2=[tc(1);tc(2);tc(4)];
  %3nd triangle
  X3=[tx(1);tx(3);tx(4)];
  Y3=[ty(1);ty(3);ty(4)];
  Z3=[tz(1);tz(3);tz(4)];
  C3=[tc(1);tc(3);tc(4)];
  %4th triangle
  X4=[tx(2);tx(3);tx(4)];
  Y4=[ty(2);ty(3);ty(4)];
  Z4=[tz(2);tz(3);tz(4)];
  C4=[tc(2);tc(3);tc(4)];
  
  X=[X1 X2 X3 X4];
  Y=[Y1 Y2 Y3 Y4];
  Z=[Z1 Z2 Z3 Z4];
  C=[C1 C2 C3 C4];
end

function printhelp()
  fprintf('%s\n\n','Invalid call to ffslice3d. Correct usage is:');
  fprintf('%s\n',' -- [X, Y, Z, C] = ffslice3d (bdfilename, tetfilename, ''Delimiter'','';'',''Format'',''%f %f %f %f'')');
end

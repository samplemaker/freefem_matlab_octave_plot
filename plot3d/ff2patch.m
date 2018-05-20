%ff2patch.m Convert FreeFem++ vertex/triangle data to patch plot data
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-05-13
%
%   [varargout] = ff2patch (filename,varargin) rearranges vertex coordinates and 
%   color data in such an order that a patch command can be invoked. Each column of
%   the input file is processed and returned as separate variable.
%
%   [varargout] = ff2patch (...,'PARAM1',val1,'PARAM2',val2,...) specifies parameter
%   name/value pairs to control the input file format
%
%       Parameter       Value
%      'Delimiter'     Specifies the delimiter of the input file
%                      (default='\t')
%      'Format'        Specifies the column format of the input file
%                      (default='auto' which means 'autodetect')
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
function [varargout] = ff2patch(filename,varargin)
  switch nargin
    case {1,3,5}
    otherwise
      printhelp();
      error('wrong number arguments');
  end
  numvarargs = length(varargin);
  optsnames = {'Delimiter','Format'};
  varargout = {'\t','auto'};
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
  fid = fopen(filename,'r');
  if fid < 0
    error('cannot open file');
  end
  if strcmp(format,'auto')
    firstLine = fgetl(fid);
    frewind(fid);
    numCols = numel(strsplit(firstLine,delimiter));
    fdata = textscan(fid,repmat('%f ',[1, numCols]),'Delimiter',delimiter);
  else
    fdata = textscan(fid,format,'Delimiter',delimiter);
  end
  fclose(fid);
  M = cell2mat(fdata);
  [sz1 sz2] = size(M);
  if mod(sz1,3) ~= 0
    error('number of triangles not a multiple of 3');
  end
  for i = 1:sz2
    varargout{i} = arrangecols(M(:,i),3);
  end
end

function [M] = arrangecols(V,c)
  r = length(V)/c;
  M = reshape(V,c,r);
end

function printhelp()
  fprintf('%s\n\n','Invalid call to ff2patch.  Correct usage is:');
  fprintf('%s\n',' -- [varargout] = ff2patch (filename, varargin)');
  fprintf('%s\n',' -- [X, Y, ...] = ff2patch (filename, ...)');
  fprintf('%s\n',' -- [X, Y, C] = ff2patch (filename, ''Delimiter'','';'',''Format'',''%f %f %f'')');
  fprintf('%s\n',' -- [X, Y, C] = ff2patch (filename)');
end

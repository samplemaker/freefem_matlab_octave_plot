%ffreadfile.m Reads one or two FreeFem++ simulation result files.
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-05-13
%
%   [varargout] = ffreadfile (varargin)
%
%   [varargout] = ffreadfile (...,'PARAM1',val1,'PARAM2',val2,...) specifies
%   parameter name/value pairs to control the input file format
%
%       Parameter       Value
%      'File1'         The file name of the first file
%      'File2'         The file name of the second file
%      'Delimiter'     Specifies the delimiter of the input files
%                      (default=';')
%      'Format'        Specifies the column format of the input files
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
function [varargout] = ffreadfile(varargin)
    switch nargin
        case {2,4,6,8}
        otherwise
            printhelp();
            error('wrong number arguments');
    end
    numvarargs = length(varargin);
    optsnames = {'File1','File2','Delimiter','Format'};
    vararginval = {[],[],';','%f %f %f %f'};
    switch numvarargs
        case {2,4,6,8}
            optargs(1:numvarargs) = varargin;
            for j = 1:length(optsnames)
                for i = 1:length(varargin)
                    if strcmp(optsnames(j),optargs(i))
                        vararginval(j) = optargs(i+1);
                    end
                end
            end
        otherwise
        printhelp();
        error('wrong number arguments');
    end
    [file1, file2, delimiter, format] = vararginval{:};

    idx=1;
    if ~isempty(file1)
        file1fid = fopen(file1,'r');
        if file1fid < 0
            error('cannot open file1');
        end
        M=textscan(file1fid,format,'Delimiter',delimiter);
        file1data = cell2mat(M);
        varargout{idx} = file1data;
        idx = idx + 1;
        fclose(file1fid);
    end
    if ~isempty(file2)
        tetfid = fopen(file2,'r');
        if tetfid < 0
            error('cannot open file2');
        end
        M=textscan(tetfid,format,'Delimiter',delimiter);
        file2data = cell2mat(M);
        varargout{idx} = file2data;
        idx = idx + 1;
        fclose(tetfid);
    end
    if (idx == 1)
        printhelp();
        error('must at least specify one file');
    end
end

function printhelp()
    fprintf('%s\n\n','Invalid call to ffreadfile.  Correct usage is:');
    fprintf('%s\n',' -- [file1data, file2data] = ffreadfile (''File1'',name1,''File2'',name2,''Delimiter'',...,''Format'',...)');
    fprintf('%s\n',' -- [file1data] = ffreadfile (''File1'',name1)');
end

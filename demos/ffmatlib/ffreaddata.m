%ffreaddata.m Reads FreeFem++ Data Files into Matlab / Octave
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-05-13
%
%   This file is a part of the ffmatlib which is hosted at
%   https://github.com/samplemaker/freefem_matlab_octave_plot
%
%   [varargout] = ffreadmesh(filename) reads a FreeFem++
%   data file created by the FreeFem++
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
%

function [varargout] = ffreaddata(filename)

    verbose=1;

    if (nargin ~= 1)
        printhelp();
        error('wrong number arguments');
    end
    fid = fopen(filename,'r');
    if fid < 0
        error('cannot open file %s', filename);
    end
    firstLine = fgetl(fid);
    frewind(fid);
    numCols = numel(strsplit(strtrim(firstLine),' '));
    fdata = textscan(fid,repmat('%f ',[1, numCols]),'Delimiter','\n');
    fclose(fid);
    M = cell2mat(fdata);
    [sz1,sz2]=size(M);
    if verbose
        fprintf('Size of data: (nDof) %i %i\n', sz1,sz2);
    end
    varargout=cell(1,sz2);
    for i = 1:sz2
        varargout{i} = M(:,i);
    end
end

function printhelp()
    fprintf('%s\n\n','Invalid call to ffreaddata. Correct usage is:');
    fprintf('%s\n',' -- [varargout] = ffreaddata (filename)');
end

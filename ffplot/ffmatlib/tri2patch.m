%tri2patch.m Convert FreeFem++ vertex/triangle data to patch plot data
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-05-13
%
%   [varargout] = tri2patch (X) rearranges triangle/vertex coordinates and 
%   color data in such an order that a patch command can be invoked.
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
function [varargout] = tri2patch(X)
    switch nargin
        case {1}
        otherwise
            printhelp();
            error('wrong number arguments');
    end
    [sz1,sz2] = size(X);
    if mod(sz1,3) ~= 0
        error('number of triangles not a multiple of 3');
    end
    varargout=cell(1,sz2);
    for i = 1:sz2
        varargout{i} = arrangecols(X(:,i),3);
    end
end

function [M] = arrangecols(V,c)
    r = length(V)/c;
    M = reshape(V,c,r);
end

function printhelp()
    fprintf('%s\n\n','Invalid call to tri2patch.  Correct usage is:');
    fprintf('%s\n',' -- [varargout] = tri2patch (X)');
end

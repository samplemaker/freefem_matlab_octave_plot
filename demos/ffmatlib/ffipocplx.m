%ffipocplx.m Invokes fftri2gridclpx
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-11-11
%
%   Wrapper function
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
function [W1,W2] = ffipocplx(p,b,t,uC,Z1,Z2)
    xpts=p(1,:);
    ypts=p(2,:);
    xdata=[xpts(t(1,:)); xpts(t(2,:)); xpts(t(3,:))];
    ydata=[ypts(t(1,:)); ypts(t(2,:)); ypts(t(3,:))];
    udata=[uC(t(1,:)), uC(t(2,:)), uC(t(3,:))].';
    W1 = fftri2gridcplx(Z1,xdata,ydata,udata);
    if (nargin == 6)
        W2 = fftri2gridcplx(Z2,xdata,ydata,udata);
    end
end


%ffinterpolate.m Interpolates PDE simulation data on a cartesian or curved meshgrid
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-11-11
%
%   This file is part of the ffmatlib which is hosted at
%   https://github.com/samplemaker/freefem_matlab_octave_plot
%
%   [w1, w2] = ffinterpolate(p, b, t, x, y, u1, u2);
%
%   interpolates the real valued or complex data u1, u2 given on a
%   triangular mesh defined by the points p, triangle t and boundary b
%   arguments onto a cartesian or curved meshgrid defined by the arguments
%   x, y. The return values are real if u1, u2 is real or complex if u1, u2
%   is complex. This is a wrapper function invoking fftri2grid() or if
%   built fftri2gridfast.c
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

function [w1, w2] = ffinterpolate(p, ~, t, x, y, u1, u2)
    xpts=p(1,:);
    ypts=p(2,:);
    xdata=[xpts(t(1,:)); xpts(t(2,:)); xpts(t(3,:))];
    ydata=[ypts(t(1,:)); ypts(t(2,:)); ypts(t(3,:))];
    u1data=[u1(t(1,:)), u1(t(2,:)), u1(t(3,:))].';
    if (nargin == 6)
       if exist('fftri2gridfast','file')
           w1 = fftri2gridfast(x,y,xdata,ydata,u1data);
       else
           fprintf('Note: To improve runtime build MEX function fftri2gridfast() from fftri2gridfast.c\n');
           w1 = fftri2grid(x,y,xdata,ydata,u1data);
       end
    else
       u2data=[u2(t(1,:)), u2(t(2,:)), u2(t(3,:))].';
       if exist('fftri2gridfast','file')
           [w1,w2] = fftri2gridfast(x,y,xdata,ydata,u1data,u2data);
       else
           fprintf('Note: To improve runtime build MEX function fftri2gridfast() from fftri2gridfast.c\n');
           [w1,w2] = fftri2grid(x,y,xdata,ydata,u1data,u2data);
       end
    end
 end


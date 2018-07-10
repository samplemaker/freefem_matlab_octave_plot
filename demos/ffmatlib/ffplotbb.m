%ffplotbb.m Plot bounding box of slicing planes
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-06-15
%
%   This file is a part of the ffmatlib which is hosted at
%   https://github.com/samplemaker/freefem_matlab_octave_plot
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
function [] = ffplotbb(slice1,slice2,slice3)
    [nslices, sz2]=size(slice1);
    if (~isequal(size(slice1), size(slice2), size(slice3)) | (sz2~=3))
        error('dimension check failed for slicing data');
    end
    for i=1:nslices
        S1=slice1(i,:);
        S2=slice2(i,:);
        S3=slice3(i,:);
        plot3([S1(1) S2(1)],[S1(2) S2(2)],[S1(3) S2(3)], ...
              '-m','LineWidth',2);
        plot3([S1(1) S3(1)],[S1(2) S3(2)],[S1(3) S3(3)], ...
              '-m','LineWidth',2);
        plot3([S2(1) (S2(1)+(S3(1)-S1(1)))],[S2(2) (S2(2)+(S3(2)-S1(2)))], ...
              [S2(3) (S2(3)+(S3(3)-S1(3)))],'-m','LineWidth',2);
        plot3([S3(1) (S2(1)+(S3(1)-S1(1)))],[S3(2) (S2(2)+(S3(2)-S1(2)))], ...
              [S3(3) (S2(3)+(S3(3)-S1(3)))],'-m','LineWidth',2);
        if nslices>1
            str1=sprintf('S1,%i',i);
            str2=sprintf('S2,%i',i);
            str3=sprintf('S3,%i',i);
        else
            str1='S1';
            str2='S2';
            str3='S3';
        end
        text(S1(1),S1(2),S1(3),str1,'HorizontalAlignment','center', ...
             'FontSize',15,'FontWeight','bold','Color','m');
        text(S2(1),S2(2),S2(3),str2,'HorizontalAlignment','center', ...
             'FontSize',15,'FontWeight','bold','Color','m');
        text(S3(1),S3(2),S3(3),str3,'HorizontalAlignment','center', ...
             'FontSize',15,'FontWeight','bold','Color','m');
    end
end

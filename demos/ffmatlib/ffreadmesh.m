%ffreadmesh.m Reads FreeFem++ Mesh File
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-05-13
%
%   [nv,nt,ns,points,tri,bd]=readffmesh(filename) reads a FreeFem++
%   mesh file created with the FreeFem++ savemesh(Th,"mesh.msh") command.
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
function [nv,nt,ns,points,tri,bd]=ffreadmesh(filename)
    fid=fopen(filename,'r');
    if fid < 0
      error('cannot open file');
    end
    headerline=textscan(fid,'%f %f %f',1,'Delimiter','\n');
    %n vertex
    nv=headerline{1};
    %n triangle
    nt=headerline{2};
    %n edges
    ns=headerline{3};
    tmp=textscan(fid,'%f %f %f',nv,'Delimiter','\n');
    %vertex coordinates [x,y] and boundary label
    points=cell2mat(tmp)';
    %triangle definition - vertex numbers (counter clock wise) and region label
    tmp=textscan(fid,'%f %f %f %f',nt,'Delimiter','\n');
    tri=cell2mat(tmp)';
    %boundary definition
    tmp=textscan(fid,'%f %f %f',ns,'Delimiter','\n');
    bd=cell2mat(tmp)';
    fclose(fid);
end

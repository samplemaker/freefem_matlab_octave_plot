%ffreadmesh.m Reads FreeFem++ Mesh File
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-05-13
%
%   [nv,nt,ns,points,tri,bd]=ffreadmesh(filename) reads a FreeFem++
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

%savemesh(Th,"2d.msh"); FreeFem++ format
%savemesh(Th,"2d.mesh"); writes two files
%savemesh(Th3d,"3d.mesh"); GMSH format

function [nv,nt,ns,p,t,b]=ffreadmesh(filename)

    mesh_format_FF=1;
    mesh_format_GMSH=2;
    mesh_format_NONE=-1;
    meshformat=mesh_format_NONE;

    fid=fopen(filename,'r');
    if fid < 0
      error('cannot open file');
    end

    fline=fgetl(fid);
    tests=regexpi(fline,'^(MeshVersionFormatted)\s*\d?\s*$','tokens');
    if ~isempty(tests)
        meshformat=mesh_format_GMSH;
    else
        tests=regexpi(fline,'^(\d+)\s+(\d+)\s+(\d+)\s*$','tokens');
        if ~isempty(tests)
            testf=str2double(tests{:});
            if ~any(isnan(testf))
                if numel(testf==3)
                    meshformat=mesh_format_FF;
                end
            end
        end
    end

    switch meshformat
        case mesh_format_FF
          
              fline = fgetl(fid);
              dimension=numel(strsplit(strtrim(fline),' '))-1;
              fprintf('mesh_format_FF; dimension=%i\n',dimension);
              frewind(fid);
              if ~(dimension==2)
                   error('only supported dimension is 2');
              end
              %start over
              headerline=textscan(fid,'%f %f %f',1,'Delimiter','\n');
              %n vertex
              nv=headerline{1};
              %n triangle
              nt=headerline{2};
              %n edges
              ns=headerline{3};
              tmp=textscan(fid,repmat('%f ',[1, 3]),nv,'Delimiter','\n');
              %vertex coordinates [x,y] and boundary label
              p=cell2mat(tmp)';
              %triangle definition - vertex numbers (counter clock wise) and region label
              tmp=textscan(fid,repmat('%f ',[1, 4]),nt,'Delimiter','\n');
              t=cell2mat(tmp)';
              %boundary definition
              tmp=textscan(fid,repmat('%f ',[1, 3]),ns,'Delimiter','\n');
              b=cell2mat(tmp)';
              fclose(fid);
              fprintf('nvertex:%i ntriangle:%i nboundary:%i\n',nv,nt,ns);
              fprintf('NaNs %i %i %i\n',any(any(isnan(p))),any(any(isnan(t))),any(any(isnan(b))));
              fprintf('sizes %ix%i %ix%i %ix%i\n',size(p),size(t),size(b));

        case mesh_format_GMSH

              i=1;
              while (isempty(regexpi(fline,'^(Vertices)\s*$','tokens')) && ~feof(fid))
                 fline = fgetl(fid);
                 i=i+1;
              end
              nv=str2double(fgetl(fid));
              fline=fgetl(fid);
              dimension=numel(strsplit(strtrim(fline),' '))-1;
              fprintf('mesh_format_GMSH; dimension %i\n',dimension);
              if ~(dimension==3)
                  error('only supported dimension is 3');
              end
              %(q1_x, q1_y, q1_z, Blabel1)
              tmp=textscan(fid,repmat('%f ',[1, 4]),nv-1,'Delimiter','\n');
              p=[str2double(strsplit(strtrim(fline),' ')); cell2mat(tmp)]';

              i=1;
              while (isempty(regexpi(fline,'^(Tetrahedra)\s*$','tokens')) && ~feof(fid))
                 fline = fgetl(fid);
                 i=i+1;
              end
              nt=str2double(fgetl(fid));
              %(1_1, 1_2, 1_3, 1_4, Rlabel1), (2_1, 2_2, 2_3, 2_4, Rlabel2) ...
              tmp=textscan(fid,repmat('%f ',[1, 5]),nt,'Delimiter','\n');
              t=cell2mat(tmp)'; 

              i=1;
              while (isempty(regexpi(fline,'^(Triangles)\s*$','tokens')) && ~feof(fid))
                 fline = fgetl(fid);
                 i=i+1;
              end
              ns=str2double(fgetl(fid));
              %(1_1, 1_2, 1_3, Blabel1), (2_1, 2_2, 2_3, Blabel2) ...
              tmp=textscan(fid,repmat('%f ',[1, 4]),ns,'Delimiter','\n');
              b=cell2mat(tmp)';
              fprintf('vertex nbvx:%i tetrahedra nbtet:%i ntriangle boundary nbtri:%i\n',nv,nt,ns);
              fprintf('NaNs %i %i %i\n',any(any(isnan(p))),any(any(isnan(t))),any(any(isnan(b))));
              fprintf('sizes %ix%i %ix%i %ix%i\n',size(p),size(t),size(b));
              fclose(fid);

        otherwise
            error('unknown mesh file format');
    end

end

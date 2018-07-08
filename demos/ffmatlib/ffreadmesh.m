%ffreadmesh.m Reads FreeFem++ Mesh Files into Matlab / Octave
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-05-13
%
%   This file is a part of the ffmatlib which is hosted at
%   https://github.com/samplemaker/freefem_matlab_octave_plot
%
%   [p,b,t,nv,nbe,nt,labels] = ffreadmesh(filename) reads a FreeFem++
%   mesh file created by the FreeFem++ savemesh(Th,"2dmesh.msh") or
%   savemesh(Th3d,"3dmesh.mesh") command.
%   2D FreeFem++ Format:
%      p: Matrix containing the nodal points
%      b: Matrix containing the boundary edges
%      t: Matrix containing the triangles
%      nv:  Number of points/vertices (Th.nv) in the Mesh
%      nt:  Number of triangles (Th.nt) in the Mesh
%      nbe: Number of (boundary) edges (Th.nbe)
%      labels: Labels found in the mesh file
%   3D GMSH Format:
%      p: Matrix containing the nodal points
%      b: Matrix containing the boundary triangles
%      t: Matrix containing the tetrahedra
%      nv:  Number of points/vertices (nbvx, Th.nv) in the Mesh
%      nt:  Number of tetrahedra (nbtet, Th.nt) in the Mesh
%      nbe: Number of (boundary) triangles (nbtri, Th.nbe)
%      labels: Labels found in the mesh file
%
%   savemesh(Th,"2d.msh"); FreeFem++ format
%   savemesh(Th,"2d.mesh"); writes two files ??
%   savemesh(Th3d,"3d.mesh"); GMSH format
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

function [p,b,t,nv,nbe,nt,labels]=ffreadmesh(filename)

    verbose=1;

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
              frewind(fid);
              if ~(dimension==2)
                   error('only supported dimension is 2');
              end
              %start over
              headerline=textscan(fid,'%f %f %f',1,'Delimiter','\n');
              %points/vertices
              nv=headerline{1};
              %triangles
              nt=headerline{2};
              %boundary/edges
              nbe=headerline{3};
              tmp=textscan(fid,repmat('%f ',[1, 3]),nv,'Delimiter','\n');
              %vertex coordinates [x,y] and boundary label
              p=cell2mat(tmp)';
              %triangle definition - vertex numbers (counter clock wise) and region label
              tmp=textscan(fid,repmat('%f ',[1, 4]),nt,'Delimiter','\n');
              t=cell2mat(tmp)';
              %boundary definition (a set of edges)
              tmp=textscan(fid,repmat('%f ',[1, 3]),nbe,'Delimiter','\n');
              b=cell2mat(tmp)';
              fclose(fid);
              labels=unique(b(3,b(3,:)~=0));
              nlabels=numel(labels);
              if verbose
                  fprintf('FreeFem++ .msh; dimension=%i\n',dimension);
                  fprintf('[Vertices nv:%i; Triangles nt:%i; Edge (Boundary) nbe:%i]\n',nv,nt,nbe);
                  fprintf('NaNs: %i %i %i\n',any(any(isnan(p))),any(any(isnan(t))),any(any(isnan(b))));
                  fprintf('Sizes: %ix%i %ix%i %ix%i\n',size(p),size(t),size(b));
                  fprintf('Labels found: %i\n' ,nlabels);
                  if nlabels<10
                      fprintf(['They are: ' repmat('%i ',1,size(labels,2)) '\n'],labels);
                  end
              end

        case mesh_format_GMSH

              i=1;
              while (isempty(regexpi(fline,'^(Vertices)\s*$','tokens')) && ~feof(fid))
                 fline = fgetl(fid);
                 i=i+1;
              end
              %points/vertices (nbvx, Th.nv)
              nv=str2double(fgetl(fid));
              fline=fgetl(fid);
              dimension=numel(strsplit(strtrim(fline),' '))-1;
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
              %ntetrahedra (nbtet, Th.nt)
              nt=str2double(fgetl(fid));
              %(1_1, 1_2, 1_3, 1_4, Rlabel1), (2_1, 2_2, 2_3, 2_4, Rlabel2) ...
              tmp=textscan(fid,repmat('%f ',[1, 5]),nt,'Delimiter','\n');
              t=cell2mat(tmp)';

              i=1;
              while (isempty(regexpi(fline,'^(Triangles)\s*$','tokens')) && ~feof(fid))
                 fline = fgetl(fid);
                 i=i+1;
              end
              %nboundary/ntriangle (nbtri, Th.nbe)
              nbe=str2double(fgetl(fid));
              %(1_1, 1_2, 1_3, Blabel1), (2_1, 2_2, 2_3, Blabel2) ...
              tmp=textscan(fid,repmat('%f ',[1, 4]),nbe,'Delimiter','\n');
              b=cell2mat(tmp)';
              fclose(fid);
              labels=unique(b(4,b(4,:)~=0));
              nlabels=numel(labels);
              if verbose
                  fprintf('GMSH .mesh; dimension=%i\n',dimension);
                  fprintf('[Vertices nv:%i; Tetrahedras nt:%i; Triangles (Boundary) nbe:%i]\n',nv,nt,nbe);
                  fprintf('NaNs: %i %i %i\n',any(any(isnan(p))),any(any(isnan(t))),any(any(isnan(b))));
                  fprintf('Sizes: %ix%i %ix%i %ix%i\n',size(p),size(t),size(b));
                  fprintf('Labels found: %i\n' ,nlabels);
                  if nlabels<10
                      fprintf(['They are: ' repmat('%i ',1,size(labels,2)) '\n'],labels);
                  end
              end

        otherwise
            error('unknown mesh file format');
    end

end

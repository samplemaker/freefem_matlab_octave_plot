%ffpdeplot3D.m Plot FreeFem ++ 3D mesh- and 3d simulation data
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-06-15
%
%   This file is a part of the ffmatlib which is hosted at
%   https://github.com/samplemaker/freefem_matlab_octave_plot
%
%   [handles,varargout] = ffpdeplot3D (points, triangles, tetrahedra, varargin)
%
%   is a function specially tailored to FreeFem++ simulation data. The FEM-Mesh
%   is entered by vertex coordinates, the boundary values, and the tetrahedron
%   definition as provided by the savemesh(Th3D, "mesh_file.mesh") command.
%
%   [handles,varargout] = ffpdeplot3D (...,'PARAM1',val1,'PARAM2',val2,...)
%   specifies parameter name/value pairs to control the input file format
%
%       Parameter       Value
%      'XYData'      Data in order to colorize the plot
%                       FreeFem++ data
%      'XYZStyle'    Plot style for boundary
%                       'interp' (default) | 'noface' | 'monochrome' 
%      'Boundary'    Shows the domain boundary / edges
%                       'on' (default) | 'off'
%      'BoundingBox' Shows the bounding box of a slice
%                       'on' | 'off' (default)
%      'BDLabels'    Draws boundary / edges with a specific label
%                       [] (default) | [label1,label2,...]
%      'Slice'       3 point slicing plane definition
%                       [] | [S1'; S2'; S3']
%      'SGridParam'  Number of grid points used for the slice
%                       'auto' (default) | [N,M]
%      'ColorMap'    ColorMap value or matrix of such values
%                       'cool' (default) | colormap name | three-column matrix of RGB triplets
%      'ColorBar'    Indicator in order to include a colorbar
%                       'on' (default) | 'off'
%      'ColorRange'  Range of values to adjust the color thresholds
%                       'minmax' (default) | [min,max]
%      'FlowData'    Data for quiver3 plot
%                       FreeFem++ point data
%      'FGridParam'  Number of grid points used for quiver3 plot
%                       'auto' (default) | [N,M]
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
function [hh,varargout] = ffpdeplot3D(points,triangles,tetrahedra,varargin)

    if (~mod(nargin,2) || (nargin<3))
        printhelp();
        error('wrong number arguments');
    end

    numvarargs = length(varargin);

    %ATTENTION: No parameters after 'Slice' otherwise pos gets screwed
    optsnames = {'XYZData', 'XYZStyle', 'Boundary', ...
                 'ColorMap', 'ColorBar', 'ColorRange', ...
                 'BDLabels', 'SGridParam', ...
                 'FlowData', 'FGridParam', 'BoundingBox', 'Slice'};

    vararginval = {[], 'interp', 'on', ...
                   'cool', 'off', 'minmax', ...
                   [], [75,75], ...
                   [], [20,20], 'off', [], [], []};

    if (numvarargs>0)
        if (~mod(numvarargs,2))
            i=1;
            while (i<numvarargs)
                pos=find(strcmpi(optsnames,varargin(i)));
                if ~isempty(pos)
                    if strcmpi('Slice',varargin(i))
                        vararginval(pos)=varargin(i+1);
                        vararginval(pos+1)=varargin(i+2);
                        vararginval(pos+2)=varargin(i+3);
                        i=i+4;
                    else
                        vararginval(pos)=varargin(i+1);
                        i=i+2;
                    end
                else
                    printhelp();
                    fprintf('%s\n',char(varargin(i)));
                    error('unknown input parameter');
                end
            end
        else
            printhelp();
            error('wrong number arguments');
        end
    end

    [xyzrawdata, xyzstyle, showboundary, ...
     setcolormap, showcolbar, colorrange, ...
     bdlabels, sgridparam, ...
     flowdata, fgridparam, boundingbox, slice1, slice2, slice3] = vararginval{:};

    hax=newplot;
    fig=get(hax,'Parent');
    oldnextplotval{1}=get(hax,'nextplot');
    set(hax,'nextplot','add');
    oldnextplotval{2}=get(fig,'nextplot');
    set(fig,'nextplot','add');

    useslicecolormap=false;
    usesliceboundingbox=false;
    points=rowvec(points);
    triangles=rowvec(triangles);
    tetrahedra=rowvec(tetrahedra);
    xpts=points(1,:);
    ypts=points(2,:);
    zpts=points(3,:);
    %Prepare spatial boundary coordinates data if necessary
    %Transform boundary to real coordinates
    if strcmpi(showboundary,'on')
        if ~isempty(bdlabels)
            keep=(triangles(4,:)==bdlabels(1));
            for i=2:numel(bdlabels)
                keep=(keep | (triangles(4,:)==bdlabels(i)));
            end
            xbddata=[xpts(triangles(1,keep)); xpts(triangles(2,keep)); xpts(triangles(3,keep))];
            ybddata=[ypts(triangles(1,keep)); ypts(triangles(2,keep)); ypts(triangles(3,keep))];
            zbddata=[zpts(triangles(1,keep)); zpts(triangles(2,keep)); zpts(triangles(3,keep))];
        else
            xbddata=[xpts(triangles(1,:)); xpts(triangles(2,:)); xpts(triangles(3,:))];
            ybddata=[ypts(triangles(1,:)); ypts(triangles(2,:)); ypts(triangles(3,:))];
            zbddata=[zpts(triangles(1,:)); zpts(triangles(2,:)); zpts(triangles(3,:))];
        end
    end
    if ~isempty(xyzrawdata)
        if (~isempty(slice1) && ~isempty(slice2) && ~isempty(slice3))
            tdata=preparetetdata(points,tetrahedra,xyzrawdata);
            if (~strcmpi(sgridparam,'auto') && isnumeric(sgridparam))
                N=sgridparam(1);
                M=sgridparam(2);
            else
                [~,nt]=size(tetrahedra);
                N=(nt)^0.33;
                M=N;
            end
            [nslices, sz2]=size(slice1);
            if (~isequal(size(slice1), size(slice2), size(slice3)) || (sz2~=3))
                error('dimension check failed for slicing data');
            end
            for i=1:nslices
                S1=slice1(i,:);
                S2=slice2(i,:);
                S3=slice3(i,:);
                [sliceTData] = slicetet2data(tdata,S1,S2,S3);
                [X,Y,Z] = gridplane3d(S1,S2,S3,N,M);
                if exist('fftet2gridfast','file')
                    [C] = fftet2gridfast(sliceTData,X,Y,Z);
%                    tx=sliceTData(:,1);
%                    ty=sliceTData(:,2);
%                    tz=sliceTData(:,3);
%                    V=sliceTData(:,4);
%                    [C] = fftet2gridfast(X,Y,Z,tx,ty,tz,V);
                else
                    fprintf('Note: To improve runtime compile MEX function fftet2gridfast()\n');
                    [C] = fftet2gridint(sliceTData,X,Y,Z);
                end
                surf(X,Y,Z,C,'EdgeColor','none');
                %Debug - show cross section
                %[SX,SY,SZ] = slicetet2patch(sliceTData(:,1:3));
                %patch(SX,SY,SZ,[0 1 1],'EdgeColor',[0 0 1],'LineWidth',1,'FaceColor','none');
                if ~strcmpi(boundingbox,'off')
                    plotboundingbox(slice1,slice2,slice3);
                    usesliceboundingbox=true;
                end
            end
            useslicecolormap=true;
            colormap(setcolormap);
            if (isnumeric(colorrange))
                caxis(colorrange);
            else
                caxis([min(min(C)) max(max(C))]);
            end
            if strcmpi(showcolbar,'on')
                hcb=colorbar;
            end
        end
        if strcmpi(showboundary,'on')
            if strcmpi(xyzstyle,'interp')
                [~,nv]=size(points);
                [cdata]=preparebddata(nv,triangles,xyzrawdata,bdlabels);
                patch(xbddata,ybddata,zbddata,cdata,'EdgeColor',[0 0 0],'LineWidth',1);
                %slicing colormap rules boundary colormaps
                if ~useslicecolormap
                    colormap(setcolormap);
                    if (isnumeric(colorrange))
                        caxis(colorrange);
                    else
                        caxis([min(min(cdata)) max(max(cdata))]);
                    end
                    if strcmpi(showcolbar,'on')
                        hcb=colorbar;
                    end
                end
            end
        end
    end
    if strcmpi(showboundary,'on')
        switch xyzstyle
            case('noface')
                patch(xbddata,ybddata,zbddata,[0 1 1],'EdgeColor',[0 0 0],'LineWidth',1,'FaceColor','none');
            case('monochrome')
                patch(xbddata,ybddata,zbddata,[0 1 1],'EdgeColor',[0 0 0],'LineWidth',1);
            otherwise
                %assume interp which is drawn elsewhere
        end
    end

    if ~isempty(flowdata)
        if ~((~isempty(slice1) && ~isempty(slice2) && ~isempty(slice3)))
            error('a quiver3 plot needs slicing data');
        end
        fdata=preparetetdata(points,tetrahedra,flowdata);
        %TODO: Only one slicing plane allowed currently
        %Need to gather all data and then call quiver3 on times
        %otherwise the arrow length is wrong
        S1=slice1(1,:);
        S2=slice2(1,:);
        S3=slice3(1,:);
        [sliceTData] = slicetet2data(fdata,S1,S2,S3);
        if (~strcmpi(fgridparam,'auto') && isnumeric(fgridparam))
            N=fgridparam(1);
            M=fgridparam(2);
        else
            N=15;
            M=15;
        end
        [X,Y,Z] = gridplane3d(S1,S2,S3,N,M);
        if exist('fftet2gridfast','file')
            [Ex,Ey,Ez] = fftet2gridfast(sliceTData,X,Y,Z);
        else
            fprintf('Note: To improve runtime build MEX function fftet2gridfast()\n');
            [Ex,Ey,Ez] = fftet2gridint(sliceTData,X,Y,Z);
        end
        quiver3(X,Y,Z,Ex,Ey,Ez,1.0);
        if (~strcmpi(boundingbox,'off') && ~usesliceboundingbox)
            plotboundingbox(slice1,slice2,slice3);
        end
    end
    daspect([1 1 1]);
    view(3);
    props = {'CameraViewAngle','DataAspectRatio','PlotBoxAspectRatio'};
    set(hax,props,get(hax,props));
    set(fig,'color',[1 1 1]);
end

function [udata] = preparetetdata(points,tetrahedra,data)
    xpts=points(1,:);
    ypts=points(2,:);
    zpts=points(3,:);
    M=rowvec(data);
    [ndim,ndof]=size(M);
    [~,nv]=size(points);
    if (ndof==nv)
        [~,nt]=size(tetrahedra);
        %Transform all tetrahedra vertices to spatial coordinates
        xtetdata=[xpts(tetrahedra(1,:)); xpts(tetrahedra(2,:)); xpts(tetrahedra(3,:)); xpts(tetrahedra(4,:))];
        ytetdata=[ypts(tetrahedra(1,:)); ypts(tetrahedra(2,:)); ypts(tetrahedra(3,:)); ypts(tetrahedra(4,:))];
        ztetdata=[zpts(tetrahedra(1,:)); zpts(tetrahedra(2,:)); zpts(tetrahedra(3,:)); zpts(tetrahedra(4,:))];
        %Continue with the FE space function (could be a vector field)
        udata=zeros(4*nt,ndim+3);
        udata(:,1:3)=[reshape(xtetdata,4*nt,1), reshape(ytetdata,4*nt,1),reshape(ztetdata,4*nt,1)];
        for i=1:ndim
            cols=M(i,:);
            u=[cols(tetrahedra(1,:)); cols(tetrahedra(2,:)); cols(tetrahedra(3,:)); cols(tetrahedra(4,:))];
            udata(:,i+3)=reshape(u,4*nt,1);
        end
    else
        error('domain: unable to recognize input data format');
    end
end

function [varargout] = preparebddata(nv,triangles,data,bdlabels)
    M=rowvec(data);
    [ndim,ndof]=size(M);
    varargout=cell(1,ndim);
    if (ndof==nv)
        %Data in points/vertex format (P1-Simulation): Works for P1 FE-space only
        if ~isempty(bdlabels)
            keep=(triangles(4,:)==bdlabels(1));
            for i=2:numel(bdlabels)
                keep=(keep | (triangles(4,:)==bdlabels(i)));
            end
            for i=1:ndim
                cols=M(i,:);
                varargout{i}=[cols(triangles(1,keep)); cols(triangles(2,keep)); cols(triangles(3,keep))];
            end
        else
            for i=1:ndim
                cols=M(i,:);
                varargout{i}=[cols(triangles(1,:)); cols(triangles(2,:)); cols(triangles(3,:))];
            end
        end
    else
        error('boundary: unable to recognize input data format');
    end
end

%Creates the spatial grid plane
%[x,y,z]=(SO+u*(SN-SO)+v*(SM-SO))
%which has NxM grid points
function [X,Y,Z] = gridplane3d(SO,SN,SM,N,M)

    n=(0:1/(N-1):1);
    m=(0:1/(M-1):1);
    [u,v] = meshgrid(n,m);
    X=(SO(1)+u*(SN(1)-SO(1))+v*(SM(1)-SO(1)));
    Y=(SO(2)+u*(SN(2)-SO(2))+v*(SM(2)-SO(2)));
    Z=(SO(3)+u*(SN(3)-SO(3))+v*(SM(3)-SO(3)));

    %\TODO: Create a local [R,S] gridplane for a view(2) plot
end

function [fdataout] = slicetet2data(fdata,S1,S2,S3)

    [npts,nvars]=size(fdata);

    S1=colvec(S1);
    S2=colvec(S2);
    S3=colvec(S3);

    % let Xn:=cross((S2-S1),(S3-S1)) be perpendicular to the slicing plane
    % and Xp a point in the plane. it turns out that for any point X,
    % i.)   in the plane          --> dot(N,(X-Xp)) == 0
    % ii.)  in front of the plane --> dot(N,(X-Xp)) > 0
    % iii.) behind of the plane   --> dot(N,(X-Xp)) < 0

    %Nodal coordinates
    M=fdata(:,1:3);
    [npts,~]=size(M);
    Xn=cross((S2-S1),(S3-S1));
    Xn0=repmat(Xn',npts,1);
    S10=repmat(S1',npts,1);
    pos=arrangecols(dot(Xn0,(M-S10),2),4);
    %Tetrahedrons that have one or more points in front of and one or more
    %points behind the cutting plane must be identified. Tetrahedra with one
    %or more points exactly in the plane do count as well
    front=any((pos>0));
    behind=any((pos<0));
    isin=any((pos==0));
    %Boolarray indicating sliced or touched tetrahedra
    keep=(((behind==1) & (front==1)) | (isin==1));
    %Extracts all affected triangles
    nTotTets=sum(keep);
    fdataout=zeros(4*nTotTets,nvars);
    for i=1:nvars
      tmp=arrangecols(fdata(:,i),4);
      fdataout(:,i)=reshape(tmp(:,keep),4*nTotTets,1);
    end
end

function [varargout] = fftet2gridint(tetdata,X,Y,Z)

    [npts,nvars]=size(tetdata);
    mesh=tetdata(:,1:3);
    a = mesh(1:4:end,:);
    b = mesh(2:4:end,:);
    c = mesh(3:4:end,:);
    d = mesh(4:4:end,:);

    [M,N]=size(X);
    [ntet,~]=size(a);

    varargout=cell(1,nvars-3);
    u=cell(1,nvars-3);
    for i=1:nvars-3
        tmp=tetdata(:,i+3);
        u{i}=reshape([tmp(1:4:end,:);tmp(2:4:end,:);tmp(3:4:end,:);tmp(4:4:end,:)],ntet,4);
        varargout{i}=NaN(M,N);
    end

    invVTET=abs(1.0./(dot(cross(c-a,b-a),d-a,2)));

    for i=1:N
        for j=1:M
            xp=[X(j,i), Y(j,i), Z(j,i)];
            Xp=repmat(xp,ntet,1);
            Va=dot(cross(d-b,c-b),Xp-b,2);
            Vb=dot(cross(c-a,d-a),Xp-a,2);
            Vc=dot(cross(d-a,b-a),Xp-a,2);
            Vd=dot(cross(b-a,c-a),Xp-a,2);
            pos=find(((Va>=0) & (Vb>=0) & (Vc>=0) & (Vd>=0)),1,'first');
            if ~isempty(pos)
                for k=1:nvars-3
                    varargout{k}(j,i)=((u{k}(pos,1).*Va(pos)+ ...
                                        u{k}(pos,2).*Vb(pos)+ ...
                                        u{k}(pos,3).*Vc(pos)+ ...
                                        u{k}(pos,4).*Vd(pos))).*invVTET(pos);
                end
            end
        end
    end
end

function [M] = arrangecols(V,c)
    r = length(V)/c;
    M = reshape(V,c,r);
end

function [] = plotboundingbox(slice1,slice2,slice3)
    [nslices, sz2]=size(slice1);
    if (~isequal(size(slice1), size(slice2), size(slice3)) || (sz2~=3))
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

%Converts tetrahedrons into vertex/triangles
function [varargout] = slicetet2patch(slicetetdata)
      [sz1,sz2] = size(slicetetdata);
      for i=1:sz2
          X=reshape(slicetetdata(:,i),4,sz1/4);
          M=[[X(1,:); X(2,:); X(3,:)]; ...
             [X(1,:); X(2,:); X(4,:)]; ...
             [X(1,:); X(3,:); X(4,:)]; ...
             [X(2,:); X(3,:); X(4,:)]];
          varargout{i}=reshape(M,3,sz1);
       end
end

function [S] = rowvec(S)
    [sz1,sz2]=size(S);
    if sz1>sz2
        S=S';
    end
end

function [S] = colvec(S)
    [sz1,sz2]=size(S);
    if sz1<sz2
        S=S';
    end
end


function printhelp()
    fprintf('%s\n\n','Invalid call to ffpdeplot3D. Correct usage is:');
    fprintf('%s\n',' -- [handles,varargout] = ffpdeplot3D (points, triangles, tetrahedra, varargin)');
    fprintf('%s\n',' -- ffpdeplot3D(p,b,t,''XYZData'',u)');
    fprintf('%s\n',' -- ffpdeplot3D(p,b,t,''BDLabels'',[30,31],''XYZStyle'',''monochrome'')');
    fprintf('%s\n',' -- ffpdeplot3D(p,b,t,''XYZData'',u,''Slice'',S1,S2,S3,''Boundary'',''off'',''ColorMap'',''jet'')');
    fprintf('%s\n',' -- ffpdeplot3D(p,b,t,''XYZData'',u,''Slice'',S1,S2,S3,''XYZStyle'',''noface'',''ColorMap'',''jet'')');
    fprintf('%s\n',' -- ffpdeplot3D(p,b,t,''FlowData'',[Ex,Ey,Ez],''Slice'',S1,S2,S3,''Boundary'',''off'',''BoundingBox'',''on'')');
    fprintf('\n');
    fprintf('''XYData''       Data in order to colorize the plot\n');
    fprintf('''XYZStyle''     Plot style for boundary (default=''interp'')\n');
    fprintf('''Boundary''     Shows the domain boundary / edges (default=''on'')\n');
    fprintf('''BoundingBox''  Shows the bounding box of a slice (default=''off'')\n');
    fprintf('''BDLabels''     Draws boundary / edges with a specific label\n');
    fprintf('''Slice''        3 point slicing plane definition\n');
    fprintf('''SGridParam''   Number of grid points used for the slice (default=''auto'')\n');
    fprintf('''ColorMap''     ColorMap value or matrix of such values (default=''cool'')\n');
    fprintf('''ColorBar''     Indicator in order to include a colorbar (default=''on'')\n');
    fprintf('''ColorRange''   Range of values to adjust the color thresholds (default=''minmax'')\n');
    fprintf('''FlowData''     Data for quiver3 plot \n');
    fprintf('''FGridParam''   Number of grid points used for quiver3 plot (default=''auto'')\n');
    fprintf('\n');
end

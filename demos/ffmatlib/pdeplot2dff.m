%pdeplot2dff.m Plot FreeFem ++ 2D mesh- and 2d simulation data
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-06-15
%
%   This file is a part of the ffmatlib which is hosted at
%   https://github.com/samplemaker/freefem_matlab_octave_plot
%
%   [varargout] = pdeplot2dff (points, boundary, triangles, varargin)
%
%   is a function specially tailored to FreeFem++ simulation data that
%   offers most of the features of the classic pdeplot() command. The FEM-Mesh
%   is entered by vertex coordinates, the boundary values, and the triangle
%   definition as provided by the savemesh(Th, "mesh_file.msh") command.
%   The simulation data can be entered either as point data (native support
%   for P1 simulation data) or as interpolation at the nodes (workaround for
%   P1, P2 and other FEM simulation results).
%
%   [varargout] = pdeplot2dff (...,'PARAM1',val1,'PARAM2',val2,...)
%   specifies parameter name/value pairs to control the input file format
%
%       Parameter       Value
%      'XYData'      Data in order to colorize the plot
%                       FreeFem++ point data | FreeFem++ triangle data
%      'XYStyle'     Coloring choice
%                       'interp' (default) | 'off'
%      'ZStyle'      Draws 3D surface plot instead of flat 2D Map plot
%                       'continuous' | 'off' (default)
%      'ColorMap'    ColorMap value or matrix of such values
%                       'cool' (default) | colormap name | three-column matrix of RGB triplets
%      'ColorBar'    Indicator in order to include a colorbar
%                       'on' (default) | 'off'
%      'ColorRange'  Range of values to adjust the color thresholds
%                       'minmax' (default) | [min,max]
%      'Mesh'        Switches the mesh off / on
%                       'on' | 'off' (default)
%      'Edge'        Shows the domain boundary / edges
%                       'on' | 'off' (default)
%      'ELabs'       Draws boundary / edges with a specific label
%                       [] (default) | [label1,label2,...]
%      'Contour'     Isovalue plot
%                       'off' (default) | 'on'
%      'CColor'      Isovalue color
%                       [0,0,0] (default) | 'auto' | RGB triplet | 'r' | 'g' | 'b' |
%      'CXYData'     Use extra (overlay) data to draw the contour plot
%                       FreeFem++ points | FreeFem++ triangle data
%      'CStyle'      Contour line style
%                       'plain' (default) | 'dashed'
%      'CLevels'     Number of isovalues used in the contour plot
%                       (default=10)
%      'CGridParam'  Number of grid points used for the contour plot
%                       'auto' (default) | [N,M]
%      'Title'       Title
%                       (default=[])
%      'XLim'        Range for the x-axis
%                       'minmax' (default) | [min,max]
%      'YLim'        Range for the y-axis
%                       'minmax' (default) | [min,max]
%      'ZLim'        Range for the z-axis
%                       'minmax' (default) | [min,max]
%      'DAspect'     Data unit length of the xy- and z-axes
%                       'off' | 'xyequal' (default) | [ux,uy,uz]
%      'FlowData'    Data for quiver plot
%                       FreeFem++ point data | FreeFem++ triangle data
%      'FGridParam'  Number of grid points used for quiver plot
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
function [varargout] = pdeplot2dff(points,boundary,triangles,varargin)

    if (~mod(nargin,2) || (nargin<3))
        printhelp();
        error('wrong number arguments');
    end

    numvarargs = length(varargin);
    optsnames = {'XYData', 'XYStyle', 'Mesh', 'Edge', 'ELabs', 'ZStyle', ...
                 'Contour', 'CGridParam', 'CColor', 'CXYData', 'CLevels', 'CStyle', ...
                 'ColorMap', 'ColorBar', 'ColorRange', ...
                 'Title', 'XLim', 'YLim', 'ZLim', 'DAspect', ...
                 'FlowData', 'FGridParam'};

    vararginval = {[], 'interp', 'off', 'off', [], 'off', ...
                   'off', 'auto', [0,0,0], [], 10, 'plain', ...
                   'cool', 'on', 'minmax', ...
                   [], 'minmax', 'minmax', 'minmax', 'xyequal', ...
                   [], 'auto'};

    if (numvarargs>0)
        if (~mod(numvarargs,2))
            for i=1:2:(numvarargs-1)
                pos=find(strcmpi(optsnames,varargin(i)));
                if ~isempty(pos)
                    vararginval(pos)=varargin(i+1);
                else
                    printhelp();
                    error('''%s'' unknown input parameter',char(varargin(i)));
                end
            end
        else
            printhelp();
            error('wrong number arguments');
        end
    end

    [xyrawdata, xystyle, showmesh, showedge, edgelabel, zstyle, ...
     contourplt, cgridparam, ccolor, contourrawdata, isolevels, contourstyle, ...
     setcolormap, showcolbar, colorrange, ...
     plottitle, plotxlim, plotylim, plotzlim, axisaspect, ...
     flowdata, fgridparam] = vararginval{:};

    hh=[];
    clab=[];
    is2dmode=true;
    points=rowvec(points);
    triangles=rowvec(triangles);
    xpts=points(1,:);
    ypts=points(2,:);
    xdata=[xpts(triangles(1,:)); xpts(triangles(2,:)); xpts(triangles(3,:))];
    ydata=[ypts(triangles(1,:)); ypts(triangles(2,:)); ypts(triangles(3,:))];
    %All types of 'XYData'-plots (all but quiver)
    if ~isempty(xyrawdata)
        [cdata]=preparedata(points,triangles,xyrawdata);
        %3D and 2D Plots invoking patch() command
        if strcmpi(contourplt,'off')
            %2D Map/Density Plot
            if strcmpi(zstyle,'off')
                %A colored Plot
                if (strcmpi(xystyle,'interp'))
                    %With Mesh
                    if ~strcmpi(showmesh,'off')
                        hh=patch(xdata,ydata,cdata,'EdgeColor',[0 0 0]);
                    else
                        %Without Mesh
                        hh=patch(xdata,ydata,cdata,'EdgeColor','none');
                    end
                %Uncolored Plot (Mesh only). Case is dublicate.
                %But it is explicitely recommended not to color the plot (XYStyle)
                else
                    hh=patch(xdata,ydata,[1 1 1],'FaceColor','none', ...
                             'EdgeColor',[0 0 1]);
                end
                %Once xyrawdata is given there is no way to create an empty plot
                view(2);
            %3D Surface Plots
            else
                is2dmode=false;
                %A colored Plot
                if (strcmpi(xystyle,'interp'))
                    %With Mesh
                    if ~strcmpi(showmesh,'off')
                        hh=patch(xdata,ydata,cdata,cdata,'EdgeColor',[0 0 0]);
                    else
                        %Without Mesh
                        hh=patch(xdata,ydata,cdata,cdata,'EdgeColor','none');
                    end
                %Uncolored Plot (Mesh only)
                else
                    hh=patch(xdata,ydata,cdata,[1 1 1],'FaceColor','none', ...
                             'EdgeColor',[0 0 1]);
                end
                %Once xyrawdata is given there is no way to create an empty plot
                if (~strcmpi(plotzlim,'minmax'))
                     zlim(plotzlim);
                     zf=1.3*(max(plotzlim)-min(plotzlim));
                else
                     zf=1.3*(max(max(cdata))-min(min(cdata)));
                end
                if strcmpi(axisaspect,'xyequal')
                    yf=(max(max(ydata))-min(min(ydata)));
                    xf=(max(max(xdata))-min(min(xdata)));
                    daspect([max(xf,yf) max(xf,yf) zf]);
                else
                    if isnumeric(axisaspect)
                        daspect(axisaspect);
                    end
                end
                view(3);
            end
        %Contour Plot
        else
            if (~strcmpi(cgridparam,'auto') && isnumeric(cgridparam))
                N=cgridparam(1);
                M=cgridparam(2);
            else
                [~,nt]=size(triangles);
                N=sqrt(nt);
                M=N;
            end
            ymin=min(min(ydata));
            ymax=max(max(ydata));
            xmin=min(min(xdata));
            xmax=max(max(xdata));
            x=linspace(xmin,xmax,N);
            y=linspace(ymin,ymax,M);
            [X,Y]=meshgrid(x,y);
            %e.g. weather forecast (temperature patch and isolevel pressure)
            if ~isempty(contourrawdata)
                ccdata=preparedata(points,triangles,contourrawdata);
            else
                ccdata=cdata;
            end
            if exist('plottri2grid','file')
                C=plottri2grid(x,y,xdata,ydata,ccdata);
            else
                if (N>100)
                    fprintf('Note: To improve runtime compile MEX function plottri2grid()\n');
                end
                C=plottri2grid1int(x,y,xdata,ydata,ccdata);
            end
            %Contour + Patch
            if strcmpi(xystyle,'interp')
                hh=patch(xdata,ydata,cdata,'EdgeColor','none');
                hold on;
                %[clab,hhc]=plotcontour(X,Y,C,isolevels,'LineColor',[0 0 0]);
                [clab,hhc]=plotcontour(X,Y,C,isolevels,ccolor,contourstyle);
                hh=[hh; hhc];
                hold off;
            %Isovalues without Patch
            else
                [clab,hhc]=plotcontour(X,Y,C,isolevels,ccolor,contourstyle);
                hh=hhc;
            end
            view(2);
        end
        colormap(setcolormap);
        if (isnumeric(colorrange))
            caxis(colorrange);
        else
            caxis([min(min(cdata)) max(max(cdata))]);
        end
        if strcmpi(showcolbar,'on')
            hcb=colorbar;
            hh=[hh; hcb];
        end
    %Uncolored, flat 2D mesh plot (no xyrawdata)
    else
        if ~strcmpi(showmesh,'off')
            hh=patch(xdata,ydata,[1 1 1],'FaceColor','none','EdgeColor',[0 0 1]);
            view(2);
        end
    end

    %Finally, the Quiver- and Border Plot is executed. These can potentially
    %be additive, therefore ...
    hold on;
    if ~isempty(flowdata)
        [udata,vdata]=preparedata(points,triangles,flowdata);
        if (~strcmpi(fgridparam,'auto') && isnumeric(fgridparam))
            N=fgridparam(1);
            M=fgridparam(2);
        else
            N=15;
            M=15;
        end
        ymin=min(min(ydata));
        ymax=max(max(ydata));
        xmin=min(min(xdata));
        xmax=max(max(xdata));
        x=linspace(xmin,xmax,N);
        y=linspace(ymin,ymax,M);
        [X,Y]=meshgrid(x,y);
        if exist('plottri2grid','file')
            [U,V]=plottri2grid(x,y,xdata,ydata,udata,vdata);
        else
            [U,V]=plottri2grid2int(x,y,xdata,ydata,udata,vdata);
        end
        idx=(~isnan(U)) & (~isnan(V));
        hq=quiver(X(idx),Y(idx),U(idx),V(idx));
        if ~isempty(hh)
            hh=[hh; hq];
        else
            hh=hq;
        end
    end

    %Adds the domain boundary (border) to all plot types
    if ~strcmpi(showedge,'off')
        boundary=rowvec(boundary);
        if ~isempty(edgelabel)
            for i=1:numel(edgelabel)
                keep=(boundary(3,:)==edgelabel(i));
                line([xpts(boundary(1,keep));xpts(boundary(2,keep))], ...
                     [ypts(boundary(1,keep));ypts(boundary(2,keep))],'Color','red','LineWidth',2);
            end
        else
            line([xpts(boundary(1,:));xpts(boundary(2,:))], ...
                 [ypts(boundary(1,:));ypts(boundary(2,:))],'Color','red','LineWidth',2);
        end
    end

    if ~isempty(plottitle)
        title(plottitle);
    end
    if ~strcmpi(plotxlim,'minmax')
        xlim(plotxlim);
    end
    if ~strcmpi(plotylim,'minmax')
        ylim(plotylim);
    end
    if (is2dmode) && (~strcmpi(axisaspect,'off'))
        if strcmpi(axisaspect,'xyequal')
            daspect([1 1 1]);
        else
            if isnumeric(axisaspect)
                daspect(axisaspect);
            end
        end
    end
    varargout{1}=hh;
    varargout{2}=clab;
end

function [S] = rowvec(S)
    [sz1,sz2]=size(S);
    if sz1>sz2
        S=S';
    end
end

%To improve the runtime, the interpolation functions for contour() and
%quiver() are doubled, rather than merging everything into one function
function [u,v] = plottri2grid2int(x, y, tx, ty, tu, tv)
    u=NaN(numel(y),numel(x));
    v=NaN(numel(y),numel(x));
    ax=tx(1,:);
    ay=ty(1,:);
    bx=tx(2,:);
    by=ty(2,:);
    cx=tx(3,:);
    cy=ty(3,:);
    invA0=(1.0)./((by-cy).*(ax-cx)+(cx-bx).*(ay-cy));
    for mx=1:numel(x)
        for my=1:numel(y)
            px=x(mx);
            py=y(my);
            Aa=((by-cy).*(px-cx)+(cx-bx).*(py-cy)).*invA0;
            Ab=((cy-ay).*(px-cx)+(ax-cx).*(py-cy)).*invA0;
            Ac=1.0-Aa-Ab;
            pos=find(((Aa>=0) & (Ab>=0) & (Ac>=0)),1,'first');
            if ~isempty(pos)
                u(my,mx)=Aa(pos).*tu(1,pos)+ ...
                         Ab(pos).*tu(2,pos)+ ...
                         Ac(pos).*tu(3,pos);
                v(my,mx)=Aa(pos).*tv(1,pos)+ ...
                         Ab(pos).*tv(2,pos)+ ...
                         Ac(pos).*tv(3,pos);
            end
        end
    end
end

function [u] = plottri2grid1int(x, y, tx, ty, tc)
    u=NaN(numel(y),numel(x));
    ax=tx(1,:);
    ay=ty(1,:);
    bx=tx(2,:);
    by=ty(2,:);
    cx=tx(3,:);
    cy=ty(3,:);
    invA0=(1.0)./((by-cy).*(ax-cx)+(cx-bx).*(ay-cy));
    for mx=1:numel(x)
        for my=1:numel(y)
            px=x(mx);
            py=y(my);
            Aa=((by-cy).*(px-cx)+(cx-bx).*(py-cy)).*invA0;
            Ab=((cy-ay).*(px-cx)+(ax-cx).*(py-cy)).*invA0;
            Ac=1.0-Aa-Ab;
            pos=find(((Aa>=0) & (Ab>=0) & (Ac>=0)),1,'first');
            if ~isempty(pos)
                u(my,mx)=Aa(pos).*tc(1,pos)+ ...
                         Ab(pos).*tc(2,pos)+ ...
                         Ac(pos).*tc(3,pos);
            end
        end
    end
end

function [varargout] = preparedata(points,triangles,data)
    M=rowvec(data);
    [ndim,ndof]=size(M);
    [~,nv]=size(points);
    varargout=cell(1,ndim);
    if (ndof==nv)
        %Data in points/vertex format (P1-Simulation): Works for P1 FE-space only
        for i=1:ndim
            cols=M(i,:);
            varargout{i}=[cols(triangles(1,:)); cols(triangles(2,:)); cols(triangles(3,:))];
        end
    else
        %Data in triangle format (work around): Works for P1 as well as for other Elements
        [~,nt]=size(triangles);
        if ~((nt*3)==ndof)
            error('unable to recognize input data format');
        end
        for i=1:ndim
            varargout{i}=reshape(M(i,:),3,nt);
        end
    end
end

function [varargout] = plotcontour(X,Y,C,isolevels,ccolor,contourstyle)
    switch (contourstyle)
        %case ('dashednegative')
        %case ('plain')
        case ('dashed')
            if strcmpi(ccolor,'auto')
                [clab,hhc]=contour(X,Y,C,isolevels,'--');
            else
                [clab,hhc]=contour(X,Y,C,isolevels,'--','LineColor',ccolor);
            end
        otherwise
            if strcmpi(ccolor,'auto')
                [clab,hhc]=contour(X,Y,C,isolevels);
            else
                [clab,hhc]=contour(X,Y,C,isolevels,'LineColor',ccolor);
            end
    end
    varargout{2}=hhc;
    varargout{1}=clab;
end

function printhelp()
    fprintf('%s\n\n','Invalid call to pdeplot2dff. Correct usage is:');
    fprintf('%s\n',' -- [varargout] = pdeplot2dff (points,boundary,triangles,varargin)');
    fprintf('%s\n',' -- [varargout] = pdeplot2dff (points,boundary,triangles,''Edge'',''on'')');
    fprintf('%s\n',' -- [varargout] = pdeplot2dff (points,boundary,triangles,''Edge'',''on'',''Mesh'',''on'')');
    fprintf('%s\n',' -- [varargout] = pdeplot2dff (points,boundary,triangles,''XYData'',u)');
    fprintf('%s\n',' -- [varargout] = pdeplot2dff (points,boundary,triangles,''XYData'',u,''ZStyle'',''continuous'')');
    fprintf('%s\n',' -- [varargout] = pdeplot2dff (points,boundary,triangles,''XYData'',u,''Contour'',''on'')');
    fprintf('%s\n',' -- [varargout] = pdeplot2dff (points,boundary,triangles,''FlowData'',v,''Edge'',''on'')');
    fprintf('\n');
    fprintf('''XYData''      Data in order to colorize the plot\n');
    fprintf('''XYStyle''     Coloring choice (default=''interp'')\n');
    fprintf('''ZStyle''      Draws 3D surface plot instead of flat 2D Map plot (default=''off'')\n');
    fprintf('''ColorMap''    ColorMap value or matrix of such values (default=''on'')\n');
    fprintf('''ColorBar''    Indicator in order to include a colorbar\n');
    fprintf('''ColorRange''  Range of values to adjust the color thresholds (default=''minmax'')\n');
    fprintf('''Mesh''        Switches the mesh off / on (default=''off'')\n');
    fprintf('''Edge''        Shows the domain boundary / edges (default=''off'')\n');
    fprintf('''ELabs''       Draws boundary / edges with a specific label (default=[])\n');
    fprintf('''Contour''     Isovalue plot (default=''off'')\n');
    fprintf('''CColor''      Isovalue color (default=[0,0,0])\n');
    fprintf('''CXYData''     Use extra (overlay) data to draw the contour plot\n');
    fprintf('''CStyle''      Contour line style (default=''plain'')\n');
    fprintf('''CLevels''     Number of isovalues used in the contour plot (default=10)\n');
    fprintf('''CGridParam''  Number of grid points used for the contour plot (default=''off'')\n');
    fprintf('''Title''       Title (default=[])\n');
    fprintf('''XLim''        Range for the x-axis (default=''minmax'')\n');
    fprintf('''YLim''        Range for the y-axis (default=''minmax'')\n');
    fprintf('''ZLim''        Range for the z-axis (default=''minmax'')\n');
    fprintf('''DAspect''     Data unit length of the xy- and z-axes (default=''xyequal'')\n');
    fprintf('''FlowData''    Data for quiver plot\n');
    fprintf('''FGridParam''  Number of grid points used for quiver plot (default=''off'')\n');
    fprintf('\n');
end

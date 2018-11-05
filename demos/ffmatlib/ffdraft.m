%ffpdeplot.m Plot FreeFem++ 2D mesh and space functions
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-06-15
%
%   This file is a part of the ffmatlib which is hosted at
%   https://github.com/samplemaker/freefem_matlab_octave_plot
%
%   [handles,varargout] = ffpdeplot (points, boundary, triangles, varargin)
%
%   is a function specially tailored to FreeFem++ simulation data that
%   offers most of the features of the classic pdeplot() command. The FEM-Mesh
%   is entered by vertex coordinates, the boundary values, and the triangle
%   definition as provided by the savemesh(Th, "mesh_file.msh") command.
%   The simulation data can be entered either as point data (native support
%   for P1 simulation data) or as interpolation at the nodes (workaround for
%   P1, P2 and other FEM simulation results).
%
%   [handles,varargout] = ffpdeplot (...,'PARAM1',val1,'PARAM2',val2,...)
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
%                       'on' (default) | 'off' | 'northoutside' ...
%      'CBTitle'     Colorbar Title
%                       (default=[])
%      'ColorRange'  Range of values to adjust the color thresholds
%                       'minmax' (default) | 'centered' | 'cropminmax' | 'cropcentered' | [min,max]
%      'Mesh'        Switches the mesh off / on
%                       'on' | 'off' (default)
%      'Boundary'    Shows the domain boundary / edges
%                       'on' | 'off' (default)
%      'BDLabels'    Draws boundary / edges with a specific label
%                       [] (default) | [label1,label2,...]
%      'BDColors'    Colorize boundary / edges with color (linked to 'BDLabels')
%                       'r' (default) | ['r','g', ... ]
%      'BDShowText'  Shows the labelnumber on the boundary / edges
%                       'on' | 'off' (default)
%      'Contour'     Isovalue plot
%                       'off' (default) | 'on'
%      'CXYData'     Use extra (overlay) data to draw the contour plot
%                       FreeFem++ points | FreeFem++ triangle data
%      'CStyle'      Contour plot style
%                       'patch' (default) | 'patchdashed' | 'patchdashedneg' | 'monochrome' | 'colormap'
%      'CColor'      Isovalue color
%                       [0,0,0] (default) | RGB triplet | 'r' | 'g' | 'b' |
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
function [hh,varargout] = ffpdeplot(points,boundary,triangles,varargin)

    if (~mod(nargin,2) || (nargin<3))
        printhelp();
        error('wrong number arguments');
    end

    numvarargs = length(varargin);
    optsnames = {'XYData', 'XYStyle', 'Mesh', 'Boundary', 'BDLabels', ...
                 'BDColors', 'BDShowText', 'ZStyle', 'Contour', ...
                 'CGridParam', 'CColor', 'CXYData', 'CLevels', 'CStyle', ...
                 'ColorMap', 'ColorBar', 'CBTitle', 'ColorRange', ...
                 'Title', 'XLim', 'YLim', 'ZLim', 'DAspect', ...
                 'FlowData', 'FGridParam'};

    vararginval = {[], 'interp', 'off', 'off', [], ...
                   'r', 'off', 'off', 'off', ...
                   'auto', [0,0,0], [], 10, 'patch', ...
                   'cool', 'on', [], 'minmax', ...
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
                    fprintf('%s\n',char(varargin(i)));
                    error('unknown input parameter');
                end
            end
        else
            printhelp();
            error('wrong number arguments');
        end
    end

    [xyrawdata, xystyle, showmesh, showboundary, bdlabels, ...
     bdcolors, bdshowtext, zstyle, contourplt, ...
     cgridparam, ccolor, contourrawdata, isolevels, contourstyle, ...
     setcolormap, showcolbar, colorbartitle, colorrange, ...
     plottitle, plotxlim, plotylim, plotzlim, axisaspect, ...
     flowdata, fgridparam] = vararginval{:};

    %"newplot" checks the values of the "NextPlot" properties and takes the
    %appropriate action based on these values.
    %Note: If there is no current figure, newplot creates one
    hax=newplot;
    fig=get(hax,'Parent');
    %Note: "hold on" sets the figure and axes "NextPlot-properties" to "add"
    oldnextplotval{1}=get(hax,'nextplot');
    set(hax,'nextplot','add');
    oldnextplotval{2}=get(fig,'nextplot');
    set(fig,'nextplot','add');
    %fprintf('nextplot hax: %s, fig: %s\n',oldnextplotval{1},oldnextplotval{2});
    %set(hax,'FontName','Times');
    %set(fig,'MenuBar','none');

    %Used for various output handels: ColorBar, Patch, Contour, Quiver
    hh=[];
    %Used for Contour labels, if there are any
    varargout{1}=[];
    varargout{2}=[];
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
                %Uncolored Plot (Mesh only).
                %Note: Case is duplicate. In this case coloring is explicetly
                %disabled by the XYStyle parameter.
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

            %Set the grid resolution depending on the cropped area, not the entire mesh
            %ymin=min(min(ydata));
            %ymax=max(max(ydata));
            %xmin=min(min(xdata));
            %xmax=max(max(xdata));
            if strcmpi(plotylim,'minmax')
                ymin = min(min(ydata));
                ymax = max(max(ydata));
            else
                ymin = min(plotylim);
                ymax = max(plotylim);
            end
            if strcmpi(plotxlim,'minmax')
                xmin = min(min(xdata));
                xmax = max(max(xdata));
            else
                xmin = min(plotxlim);
                xmax = max(plotxlim);
            end

            x=linspace(xmin,xmax,N);
            y=linspace(ymin,ymax,M);
            [X,Y]=meshgrid(x,y);
            %Scalar data overlayed with some extra contour data
            %like e.g. a weather forecast plot employing a temperature patch plot
            %with overlayed isolevel pressure values
            if ~isempty(contourrawdata)
                if strcmpi(contourstyle,'colormap') || strcmpi(contourstyle,'monochrome')
                    error('CStyle is ''monochrome'' or ''colormap'' but CXYData is given - use XYData instead!');
                end
                ccdata=preparedata(points,triangles,contourrawdata);
            else
                ccdata=cdata;
            end
            if exist('ffplottri2grid','file')
                C=ffplottri2grid(x,y,xdata,ydata,ccdata);
            else
                if (N>100)
                    fprintf('Note: To improve runtime build MEX function ffplottri2grid() from ffplottri2grid.c\n');
                end
                C=ffplottri2gridint(x,y,xdata,ydata,ccdata);
            end
            switch (contourstyle)
                %Filled iso value plot
                case ('colormap')
                    [clab,hhc]=contour(X,Y,C,isolevels);
                    hh=hhc;
                    varargout{1}=clab;
                %Single color iso value
                case ('monochrome')
                    [clab,hhc]=contour(X,Y,C,isolevels,'LineColor',ccolor);
                    hh=hhc;
                    varargout{1}=clab;
                %Patch with overlayed dashed style iso value
                case ('patchdashed')
                    hh=patch(xdata,ydata,cdata,'EdgeColor','none');
                    [clab,hhc]=contour(X,Y,C,isolevels,'--','LineColor',ccolor);
                    hh=[hh; hhc];
                    varargout{1}=clab;
                %Same as before but positive values with solid lines and negative
                %with dashed lines
                case ('patchdashedneg')
                    hh=patch(xdata,ydata,cdata,'EdgeColor','none');
                    if (length(isolevels)==1)
                        step = (max(max(C))-min(min(C)))/(isolevels+1);
                        isolevels = min(min(C))+step:step:max(max(C))-step;
                    end
                    isopos = isolevels(isolevels>=0);
                    if (length(isopos)==1)
                        isopos = [isopos isopos];
                    end;
                    isoneg = isolevels(isolevels<0);
                    if(length(isoneg)==1)
                        isoneg = [isoneg isoneg];
                    end;
                    [clabneg,hhcneg]=contour(X,Y,C,isoneg,'--','LineColor',ccolor);
                    [clabpos,hhcpos]=contour(X,Y,C,isopos,'LineColor',ccolor);
                    hh=[hh;hhcneg;hhcpos];
                    varargout{1}=clabneg;
                    varargout{2}=clabpos;
                %Map all others to the 'patch' case
                otherwise
                    hh=patch(xdata,ydata,cdata,'EdgeColor','none');
                    [clab,hhc]=contour(X,Y,C,isolevels,'LineColor',ccolor);
                    hh=[hh; hhc];
                    varargout{1}=clab;
            end
            view(2);
        end
        colormap(setcolormap);
        if (isnumeric(colorrange))
            if (min(colorrange) == max(colorrange))
               error('''ColorRange'': Must be a numeric 2-element vector where LIM1 < LIM2');
            end
            caxis(colorrange);
        else
            if strcmpi(colorrange,'minmax') || strcmpi(colorrange,'centered')
               %Set to [min max] of the whole mesh
               if strcmpi(colorrange,'minmax')
                   caxis([min(min(cdata)) max(max(cdata))]);
               else
                   %Colormap symmetric around zero ('centered')
                   caxis([-max(max(abs(cdata))) max(max(abs(cdata)))]);
               end
            else
               %Set to [min max] of cropped area (auto ranging)
               [~,sz2]=size(cdata);
               keep = true(1,sz2);
               if ~strcmpi(plotxlim,'minmax')
                  keepx = (xdata > min(plotxlim)) & (xdata < max(plotxlim));
                  keep = keep & (keepx(1,:) & keepx(2,:) & keepx(3,:));
               end
               if ~strcmpi(plotylim,'minmax')
                  keepy = (ydata > min(plotylim)) & (ydata < max(plotylim));
                  keep = keep & (keepy(1,:) & keepy(2,:) & keepy(3,:));
               end
               crng = [min(min(cdata(:,keep))) max(max(cdata(:,keep)))];
               if (isempty(crng) || (min(crng) == max(crng)))
                  %E.g. too much magnification
                  error('''ColorRange'': No color spread in the cropped section / try ''minmax''');
               end
               if strcmpi(colorrange,'cropminmax')
                   caxis(crng);
               else
                   %Colormap symmetric around zero ('cropcentered')
                   caxis([-max(abs(crng)) max(abs(crng))]);
               end
            end
        end
        if ~(strcmpi(showcolbar,'off'))
            if strcmpi(showcolbar,'on')
                hcb=colorbar;
            else
                hcb=colorbar(showcolbar);
            end
            hh=[hcb; hh];
            if ~isempty(colorbartitle)
                title(hcb,colorbartitle);
            end
        end
    %Uncolored, flat 2D mesh plot (no xyrawdata is given)
    else
        if ~strcmpi(showmesh,'off')
            hh=patch(xdata,ydata,[1 1 1],'FaceColor','none','EdgeColor',[0 0 1]);
            view(2);
        end
        %If "mesh" is switched of we can exit at this point without any plot
    end

    %Quiver
    if ~isempty(flowdata)
        [udata,vdata]=preparedata(points,triangles,flowdata);
        if (~strcmpi(fgridparam,'auto') && isnumeric(fgridparam))
            N=fgridparam(1);
            M=fgridparam(2);
        else
            N=20;
            M=20;
        end

        %Set the grid resolution depending on the cropped area, not the entire mesh
        %ymin=min(min(ydata));
        %ymax=max(max(ydata));
        %xmin=min(min(xdata));
        %xmax=max(max(xdata));
        if strcmpi(plotylim,'minmax')
            ymin = min(min(ydata));
            ymax = max(max(ydata));
        else
            ymin = min(plotylim);
            ymax = max(plotylim);
        end
        if strcmpi(plotxlim,'minmax')
            xmin = min(min(xdata));
            xmax = max(max(xdata));
        else
            xmin = min(plotxlim);
            xmax = max(plotxlim);
        end

        x=linspace(xmin,xmax,N);
        y=linspace(ymin,ymax,M);
        [X,Y]=meshgrid(x,y);
        if exist('ffplottri2grid','file')
            [U,V]=ffplottri2grid(x,y,xdata,ydata,udata,vdata);
        else
            [U,V]=ffplottri2gridint(x,y,xdata,ydata,udata,vdata);
        end
        idx=(~isnan(U)) & (~isnan(V));
        hq=quiver(X(idx),Y(idx),U(idx),V(idx));
        if ~isempty(hh)
            hh=[hh; hq];
        else
            hh=hq;
        end
    end

    %Adds the domain boundary (border) to all plot types, or creates a border plot
    %if nothing has been drawn yet
    if ~strcmpi(showboundary,'off')
        boundary=rowvec(boundary);
        if ~isempty(bdlabels)
            if (numel(bdlabels) ~= numel(bdcolors)) && (numel(bdcolors) > 1)
                error('''BDColors'': Numel of BDColors not equal to numel of BDLabels');
            end
            for i=1:numel(bdlabels)
                keep=(boundary(3,:)==bdlabels(i));
                if any(keep)
                    if (numel(bdcolors) > 1)
                        line([xpts(boundary(1,keep));xpts(boundary(2,keep))], ...
                             [ypts(boundary(1,keep));ypts(boundary(2,keep))],'Color',char(bdcolors(i)),'LineWidth',2);
                    else
                        line([xpts(boundary(1,keep));xpts(boundary(2,keep))], ...
                             [ypts(boundary(1,keep));ypts(boundary(2,keep))],'Color',char(bdcolors),'LineWidth',2);
                    end
                    if strcmpi(bdshowtext,'on')
                        textpos=find(keep,1,'first');
                        text(xpts(boundary(1,textpos)),ypts(boundary(1,textpos)),num2str(bdlabels(i)));
                    end
                end
            end
        else
            %No BDLabels specified: All labels are displayed without text
            if (numel(bdcolors) == 1)
                line([xpts(boundary(1,:));xpts(boundary(2,:))], ...
                     [ypts(boundary(1,:));ypts(boundary(2,:))],'Color',char(bdcolors),'LineWidth',2);
            else
                error('''BDColors'': Multple BDColors only with BDLabel specifier possible');
            end
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
    %Restore axes and figure handles
    set(hax,'nextplot',oldnextplotval{1});
    set(fig,'nextplot',oldnextplotval{2});
end

function [S] = rowvec(S)
    [sz1,sz2]=size(S);
    if sz1>sz2
        S=S';
    end
end

%Triangle to rectangular grid interpolation. Employs barycentric interpolation.
%Note: In order to improve runtime keep code for scalar and vector field
%problems separated. Cases can be distinguished by the number of calling arguments
function [u,v] = ffplottri2gridint(x, y, tx, ty, tu, tv)
    ax=tx(1,:);
    ay=ty(1,:);
    bx=tx(2,:);
    by=ty(2,:);
    cx=tx(3,:);
    cy=ty(3,:);
    invA0=(1.0)./((by-cy).*(ax-cx)+(cx-bx).*(ay-cy));
    if (nargin == 6)
        u=NaN(numel(y),numel(x));
        v=NaN(numel(y),numel(x));
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
    else
        u=NaN(numel(y),numel(x));
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
                end
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

function printhelp()
    fprintf('%s\n\n','Invalid call to ffpdeplot. Correct usage is:');
    fprintf('%s\n',' -- [handles,varargout] = ffpdeplot (points,boundary,triangles,varargin)');
    fprintf('%s\n',' -- [handles,varargout] = ffpdeplot (points,boundary,triangles,''Boundary'',''on'')');
    fprintf('%s\n',' -- [handles,varargout] = ffpdeplot (points,boundary,triangles,''Boundary'',''on'',''Mesh'',''on'')');
    fprintf('%s\n',' -- [handles,varargout] = ffpdeplot (points,boundary,triangles,''XYData'',u)');
    fprintf('%s\n',' -- [handles,varargout] = ffpdeplot (points,boundary,triangles,''XYData'',u,''ZStyle'',''continuous'')');
    fprintf('%s\n',' -- [handles,varargout] = ffpdeplot (points,boundary,triangles,''XYData'',u,''Contour'',''on'')');
    fprintf('%s\n',' -- [handles,varargout] = ffpdeplot (points,boundary,triangles,''FlowData'',v,''Boundary'',''on'')');
    fprintf('\n');
    fprintf('''XYData''      Data in order to colorize the plot\n');
    fprintf('''XYStyle''     Coloring choice (default=''interp'')\n');
    fprintf('''ZStyle''      Draws 3D surface plot instead of flat 2D Map plot (default=''off'')\n');
    fprintf('''ColorMap''    ColorMap value or matrix of such values (default=''on'')\n');
    fprintf('''ColorBar''    Indicator in order to include a colorbar\n');
    fprintf('''CBTitle''     Colorbar Title (default=[])\n');
    fprintf('''ColorRange''  Range of values to adjust the color thresholds (default=''minmax'')\n');
    fprintf('''Mesh''        Switches the mesh off / on (default=''off'')\n');
    fprintf('''Boundary''    Shows the domain boundary / edges (default=''off'')\n');
    fprintf('''BDLabels''    Draws boundary / edges with a specific label (default=[])\n');
    fprintf('''BDColors''    Colorize boundary / edges with color (default=''r'')\n');
    fprintf('''BDShowText''  Shows the labelnumber on the boundary / edges (default=''off'')\n');
    fprintf('''Contour''     Isovalue plot (default=''off'')\n');
    fprintf('''CColor''      Isovalue color (default=[0,0,0])\n');
    fprintf('''CXYData''     Use extra (overlay) data to draw the contour plot\n');
    fprintf('''CStyle''      Contour line style (default=''patch'')\n');
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

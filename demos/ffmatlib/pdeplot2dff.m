%pdeplot2dff.m Wrapper function reading and plotting FreeFem++
%              2D mesh and FreeFem++ data
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-05-13
%
%   [handles] = pdeplot2dff (points,boundary,triangles,varargin)
%   is a customized FreeFem++ wrapper function to implement some
%   of the classic pdeplot() features. The input arguments
%   include the vertex coordinates, the triangles and the boundary
%   definition as provided by the FreeFem++ savemesh(Th,"mesh.msh")
%   command.
%
%   [handles] = pdeplot2dff (...,'PARAM1',val1,'PARAM2',val2,...)
%   specifies parameter name/value pairs to control the input file format
%
%       Parameter       Value
%      'XYData'     Scalar value in order to colorize the plot
%      'ZStyle'     3D surface instead of 2D Map plot
%                      (default='off')
%      'ColorMap'   Specifies the colormap
%                      (default='jet')
%      'ColorBar'   Indicator in order to include a colorbar
%                      (default='on')
%      'Mesh'       Switches the mesh off/on
%                      (default='off')
%      'Edge'       Shows the PDE boundary
%                      (default='off')
%      'Contour'    Isovalue plot
%                      (default='off' ; accepted values : 'on', 'off' or 'value')
%      'Levels'     Number of isovalues for contour plot
%                      (default=10)
%      'ColorRange' Range of values to ajust the colormap
%                      (default='minmax', or specify [min,max])
%      'Title'      Title
%                      (default=[])
%      'XLim'       Range for the x-axis
%                      (default='minmax')
%      'YLim'       Range for the y-axis
%                      (default='minmax')
%      'ZLim'       Range for the z-axis
%                      (default='minmax')
%      'ZAspect'    3D Plot aspect ratio
%                      (default='on' ; accepted values : 'on' or [a1,a2,a3]))
%      'FlowData'   Data for quiver plot
%      'GridParam'  Number of grid points used for contour() and quiver() plots 
%                      (default='off' ; accepted values : 'off' or [N,M])

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
function [hh] = pdeplot2dff(points,boundary,triangles,varargin)

    if (~mod(nargin,2) || (nargin<3))
        printhelp();
        error('wrong number arguments');
    end

    numvarargs = length(varargin);
    optsnames = {'XYData','ColorMap','Mesh','Edge','ZStyle','Contour','Levels', ...
                 'ColorBar','ColorRange','Title','XLim','YLim','ZLim','ZAspect','FlowData','GridParam'};
    vararginval = {[],'jet','off','off','off','off',10,'on','minmax',[], ...
                   'minmax','minmax','minmax','on',[],'off'};

    if (numvarargs>0)
        if (~mod(numvarargs,2))
            optargs(1:numvarargs) = varargin;
            i=1;
            while (i<numvarargs)
                found=false;
                for j=1:length(optsnames)
                    if strcmpi(optargs(i),optsnames(j))
                        vararginval(j)=optargs(i+1);
                        found=true;
                    end
                end
                if (~found)
                    printhelp();
                    error('''%s'' unknown parameter',char(optargs(i)));
                end
                i=i+2;
            end
        else
            printhelp();
            error('wrong number arguments');
        end
    end

    [xydata,setcolormap,showmesh,showedge,zstyle, contourplt,isolevels, ...
     showcolbar,colorrange,plottitle,plotxlim,plotylim,plotzlim,zaspect, ...
     flowdata,gridparam] = vararginval{:};

    points=rowvec(points);
    triangles=rowvec(triangles);
    xpts=points(1,:);
    ypts=points(2,:);
    dataX=[xpts(triangles(1,:)); xpts(triangles(2,:)); xpts(triangles(3,:))];
    dataY=[ypts(triangles(1,:)); ypts(triangles(2,:)); ypts(triangles(3,:))];
    %All types of xydata-plots
    if (~isempty(xydata))
        cpts=rowvec(xydata);
        dataC=[cpts(triangles(1,:)); cpts(triangles(2,:)); cpts(triangles(3,:))];
        if strcmpi(contourplt,'off')
            %2D Map/Density Plots with Color
            if strcmpi(zstyle,'off')
                if strcmpi(showmesh,'off')
                    hh=patch(dataX,dataY,dataC,'EdgeColor','none');
                else
                    hh=patch(dataX,dataY,dataC,'EdgeColor',[0 0 0],'LineWidth',1);
                end
                view(2);
            %3D Surface Plots with Color
            else
                if strcmpi(showmesh,'off')
                    hh=patch(dataX,dataY,dataC,dataC,'EdgeColor','none');
                else
                    hh=patch(dataX,dataY,dataC,dataC,'EdgeColor',[0 0 0],'LineWidth',1);
                end
                zfactor=2.0;
                if (~strcmpi(plotzlim,'minmax'))
                     zlim(plotzlim);
                     zfactor=zfactor*(max(plotzlim)-min(plotzlim));
                else
                     zfactor=zfactor*(max(max(dataC))-min(min(dataC)));
                end
                if (~strcmpi(zaspect,'on'))
                     daspect(zaspect);
                else
                     daspect([1 1 zfactor]);
                end
                view(3);
            end
        %Contour Plot with Color
        else
            hh=patch(dataX,dataY,dataC,'EdgeColor','none');
            view(2);
            hold on;
            [~,sz2]=size(dataX);
            if (strcmpi(gridparam,'off'))
                N=sqrt(sz2);
                M=N;
                if (N > 150)
                   N=150;
                   M=150;
                end
            else
                N=gridparam(1);
                M=gridparam(2);
            end
            ymin=min(min(dataY));
            ymax=max(max(dataY));
            xmin=min(min(dataX));
            xmax=max(max(dataX));
            x=linspace(xmin,xmax,N);
            y=linspace(ymin,ymax,M);
            [X,Y]=meshgrid(x,y);
            C=scalartri2grid(dataX,dataY,dataC,x,y);
            [lab,hhc]=contour(X,Y,C,isolevels,'LineColor',[0 0 0]);
            if strcmpi(contourplt,'value')
                 clabel(lab,hhc);
            end
            hh=[hh; hhc];
            hold off;
        end
        colormap(setcolormap);
        if(isnumeric(colorrange))
            caxis(colorrange);
        else
            caxis([min(min(dataC)) max(max(dataC))]);
        end
        if strcmpi(showcolbar,'on')
            hcb=colorbar;
            hh=[hh; hcb];
        end
    %Uncolored, flat 2D mesh plot (no xydata)
    else
        if ~strcmpi(showmesh,'off')
            hh=patch(dataX,dataY,[1 1 1],'EdgeColor',[0 0 1],'LineWidth',1);
            view(2);
        end
    end

    %Finally, the Quiver- and Border Plot is executed. These can potentially
    %be additive, therefore ...
    hold on;
    if (~isempty(flowdata))
        cpts=rowvec(flowdata);
        upts=cpts(1,:); vpts=cpts(2,:);
        dataU=[upts(triangles(1,:)); upts(triangles(2,:)); upts(triangles(3,:))];
        dataV=[vpts(triangles(1,:)); vpts(triangles(2,:)); vpts(triangles(3,:))];
        if (strcmpi(gridparam,'off'))
            N=15;
            M=15;
        else
            N=gridparam(1);
            M=gridparam(2);
        end
        ymin=min(min(dataY));
        ymax=max(max(dataY));
        xmin=min(min(dataX));
        xmax=max(max(dataX));
        x=linspace(xmin,xmax,N);
        y=linspace(ymin,ymax,M);
        [X,Y]=meshgrid(x,y);
        [U,V]=flowtri2grid(dataX,dataY,dataU,dataV,x,y);
        idx = ~isnan(U) & ~isnan(V);
        quiver(X(idx),Y(idx),U(idx),V(idx));
    end

    %Adds the domain boundary (border) to all plot types
    if ~strcmpi(showedge,'off')
        boundary=rowvec(boundary);
        line([xpts(boundary(1,:));xpts(boundary(2,:))], ...
             [ypts(boundary(1,:));ypts(boundary(2,:))],'Color','red','LineWidth',2);
    end

    if (~isempty(plottitle))
        title(plottitle);
    end
    if (~strcmpi(plotxlim,'minmax'))
        xlim(plotxlim);
    end
    if (~strcmpi(plotylim,'minmax'))
        ylim(plotylim);
    end

end

function [S] = rowvec(S)
    [sz1,sz2]=size(S);
    if sz1>sz2
        S=S';
    end
end

%To improve the runtime, the interpolation functions for contour() and
%quiver() are doubled, rather than merging everything into one function
function [u,v] = flowtri2grid(tx, ty, tu, tv, X, Y)
    u=NaN(numel(Y),numel(X));
    v=NaN(numel(Y),numel(X));
    ax=tx(1,:);
    ay=ty(1,:);
    bx=tx(2,:);
    by=ty(2,:);
    cx=tx(3,:);
    cy=ty(3,:);
    invA0=(1.0)./((by-cy).*(ax-cx)+(cx-bx).*(ay-cy));
    for mx=1:numel(X)
        for my=1:numel(Y)
            px=X(mx);
            py=Y(my);
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

function [u] = scalartri2grid(tx, ty, tc, X, Y)
    u=NaN(numel(Y),numel(X));
    ax=tx(1,:);
    ay=ty(1,:);
    bx=tx(2,:);
    by=ty(2,:);
    cx=tx(3,:);
    cy=ty(3,:);
    invA0=(1.0)./((by-cy).*(ax-cx)+(cx-bx).*(ay-cy));
    for mx=1:numel(X)
        for my=1:numel(Y)
            px=X(mx);
            py=Y(my);
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

function printhelp()
    fprintf('%s\n\n','Invalid call to pdeplot2dff.  Correct usage is:');
    fprintf('%s\n',' -- [handles] = pdeplot2dff (points,boundary,triangles,varargin)');
    fprintf('%s\n',' -- [handles] = pdeplot2dff (points,boundary,triangles,''XYData'',u,''ColorMap'',''jet'',''Mesh'',''on'')');
    fprintf('%s\n',' -- [handles] = pdeplot2dff (points,boundary,triangles,''XYData'',u,''ZStyle'',''on'')');
    fprintf('%s\n',' -- [handles] = pdeplot2dff (points,boundary,triangles,''XYData'',u,''Edge'',''on'',''Contour'',''on'',''Levels'',10)');
    fprintf('\n');
    fprintf('''XYData''     Scalar value in order to colorize the plot\n');
    fprintf('''ZStyle''     3D surface instead of 2D Map plot (default=''off'')\n');
    fprintf('''ColorMap''   Specifies the colormap (default=''jet'')\n');
    fprintf('''ColorBar''   Indicator in order to include a colorbar (default=''on'')\n');
    fprintf('''Mesh''       Switches the mesh off/on (default=''off'')\n');
    fprintf('''Edge''       Shows the PDE boundary (default=''off'')\n');
    fprintf('''Contour''    Isovalue plot (default=''off'' ; accepted values : ''on'', ''off'' or ''value'')\n');
    fprintf('''Levels''     Number of isovalues for contour plot (default=10)\n');
    fprintf('''ColorRange'' Range of values to ajust the colormap (default=''minmax'', or specify [min,max])\n');
    fprintf('''Title''      Title (default=[])\n');
    fprintf('''XLim''       Range for the x-axis (default=''minmax'')\n');
    fprintf('''YLim''       Range for the y-axis (default=''minmax'')\n');
    fprintf('''ZLim''       Range for the z-axis (default=''minmax'')\n');
    fprintf('''ZAspect''    3D Plot aspect ratio (default=''on'' ; accepted values : ''on'' or [a1,a2,a3])\n');
    fprintf('''FlowData''   Data for quiver plot\n');
    fprintf('''GridParam''  Number of grid points used for contour() and quiver() plots (''off'' or [N,M])\n');
    fprintf('\n');
end

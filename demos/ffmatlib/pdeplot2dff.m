%pdeplot2dff.m FreeFem++ - Wrapper Function for PDEplot()
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-05-13
%
%   [handles] = pdeplot2dff (points,triangles,boundary,varargin)
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
%                      (default='off')
%      'ZStyle'     Appear as 3D surface instead of 2D Map plot
%                      (default='off')
%      'ColorMap'   Specifies the colormap
%                      (default='jet')
%      'ColorBar'   Indicator to include a colorbar
%                      (default='on')
%      'Mesh'       Switch mesh off/on 
%                      (default='off')
%      'Edge'       Show PDE boundary
%                      (default='off')
%      'Contour'    Isovalue plot
%                      (default='off')
%      'Levels'     Number of isovalues for contour plot
%                      (default=10)
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
function [hh] = pdeplot2dff(points,triangles,boundary,varargin)
    switch nargin
        case {3,5,7,9,11,13}
        otherwise
            printhelp();
            error('wrong number arguments');
    end
    numvarargs = length(varargin);
    optsnames = {'XYData','ColorMap','Mesh','Edge','ZStyle','Contour','Levels','ColorBar'};
    vararginval = {[],'jet','off','off','off','off',10,'on'};
    switch numvarargs
        case 0
        case {2,4,6,8,10}
            optargs(1:numvarargs) = varargin;
            for j = 1:length(optsnames)
                for i = 1:length(varargin)
                    if strcmp(optsnames(j),optargs(i))
                        vararginval(j) = optargs(i+1);
                    end
                end
            end
        otherwise
            printhelp();
            error('wrong number arguments');
    end
    [xydata,setcolormap,showmesh,showedge,zstyle,contourplt,isolevels,showcolbar] = vararginval{:};
    xpts=points(1,:);
    ypts=points(2,:);
    triDataX=[xpts(triangles(1,:)); xpts(triangles(2,:)); xpts(triangles(3,:))];
    triDataY=[ypts(triangles(1,:)); ypts(triangles(2,:)); ypts(triangles(3,:))];
    if (~isempty(xydata))
        cpts=rowvec(xydata);
        colData=[cpts(triangles(1,:)); cpts(triangles(2,:)); cpts(triangles(3,:))];
        if strcmp(contourplt,'off')
            if strcmp(showmesh,'off')
                if strcmp(zstyle,'off')
                    hh=patch(triDataX,triDataY,colData,'EdgeColor','none');
                    view(2);
                else
                    hh=patch(triDataX,triDataY,colData,colData,'EdgeColor','none');
                    daspect([1 1 3.0*(max(max(colData))-min(min(colData)))]);
                    view(3);
                end
            else
                if strcmp(zstyle,'off')
                    hh=patch(triDataX,triDataY,colData,'EdgeColor',[0 0 0],'LineWidth',1);
                    view(2);
                else
                    hh=patch(triDataX,triDataY,colData,colData,'EdgeColor',[0 0 0],'LineWidth',1);
                    daspect([1 1 3.0*(max(max(colData))-min(min(colData)))]);
                    view(3);
                end
            end
        else
            hh=patch(triDataX,triDataY,colData,'EdgeColor','none');
            view(2);
            hold on;
            [~,sz2]=size(triDataX);
            N=sqrt(sz2);
            if (N > 150)
                N=150;
            end
            ymin=min(min(triDataY));
            ymax=max(max(triDataY));
            xmin=min(min(triDataX));
            xmax=max(max(triDataX));
            x=linspace(xmin,xmax,N);
            y=linspace(ymin,ymax,N);
            [X,Y]=meshgrid(x,y);
            C=tri2grid(triDataX,triDataY,colData,x,y);
            [~,hhc]=contour(X,Y,C,isolevels,'LineColor',[0 0 0]);
            hh=[hh; hhc];
            hold off;
        end
        colormap(setcolormap);
        caxis([min(min(colData)) max(max(colData))]);
        if strcmp(showcolbar,'on')
            hcb=colorbar;
            hh=[hh; hcb];
        end
    else
        if ~strcmp(showmesh,'off')
            hh=patch(triDataX,triDataY,[1 1 1],'EdgeColor',[0 0 1],'LineWidth',1);
            view(2);
        end
    end
    if ~strcmp(showedge,'off')
        line([xpts(boundary(1,:));xpts(boundary(2,:))], ...
             [ypts(boundary(1,:)); ypts(boundary(2,:))],'Color','red','LineWidth',2);
    end
end

function [S] = rowvec(S)
    if size(S,1)>1
        S=S';
    end
end

function [colout] = tri2grid(tx, ty, colin, X, Y)
    colout=NaN(numel(Y),numel(X));
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
               colout(my,mx)=Aa(pos).*colin(1,pos)+ ...
                             Ab(pos).*colin(2,pos)+ ...
                             Ac(pos).*colin(3,pos);
            end
        end
    end
end

function printhelp()
    fprintf('%s\n\n','Invalid call to pdeplot2dff.  Correct usage is:');
    fprintf('%s\n',' -- [handles] = pdeplot2dff (points,triangles,boundary,varargin)');
    fprintf('%s\n',' -- [handles] = pdeplot2dff (points,triangles,boundary,''XYData'',u,''ColorMap'',''jet'',''Mesh'',''on'')');
    fprintf('%s\n',' -- [handles] = pdeplot2dff (points,triangles,boundary,''XYData'',u,''ZStyle'',''on'')');
    fprintf('%s\n',' -- [handles] = pdeplot2dff (points,triangles,boundary,''XYData'',u,''Edge'',''on'',''Contour'',''on'',''Levels'',10)');
    fprintf('\n');
end

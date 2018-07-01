clear all;

%Specifies the path where the ffmatlib can be found
addpath('ffmatlib');

%Reads a FreeFem++ mesh created with the savemesh(Th,"mesh.msh"); command
[nv,nbe,nt,points,boundary,triangles]=ffreadmesh('demo_meshp1.msh');

%Reads the PDE data
fid=fopen('demo_data_points_p1.txt','r');
data=textscan(fid,'%f','Delimiter','\n');
fclose(fid);
u=cell2mat(data);
[sz1,sz2]=size(u);
fprintf('size of data (nDof): %i %i\n', sz1,sz2);

fid=fopen('demo_flowdata_points_p1.txt','r');
data=textscan(fid,'%f %f','Delimiter','\n');
fclose(fid);
v=cell2mat(data)';

%%%%%%% 2D Mesh Plot
figure();
pdeplot2dff(points,boundary,triangles, ...
            'Mesh','on', ...
            'Title','Mesh without boundary');

%%%%%%% Including Boundary
figure();
pdeplot2dff(points,boundary,triangles, ...
            'Mesh','on', ...
            'Edge','on', ...
            'Title','Mesh with boundary');

%%%%%%% Boundary only
figure();
pdeplot2dff(points,boundary,triangles, ...
            'Edge','on', ...
            'Title','Boundary only');

%%%%%%% 2D PDE Map Plot
figure();
pdeplot2dff(points,boundary,triangles, ...
            'XYData',u, ...
            'Title','2D Density Plot');

%%%%%%% 2D PDE ColorMap Plot
figure();
mymap='jet';
pdeplot2dff(points,boundary,triangles, ...
            'XYData',u,'ColorMap',mymap, ...
            'Title',strcat('2D Density Plot - ColMap: ',mymap));

%%%%%%% 2D PDE Map Plot - customized HSV Map
figure();
hsv=[4./6.,1,0.5; 4./6.,1,1; 5./6.,1,1; 1,1,1; 1,0.5,1];
[sz1,~]=size(hsv);
cm_data=interp1(linspace(0,1,sz1),hsv,linspace(0,1,64));
usedmap=hsv2rgb(cm_data);
pdeplot2dff(points,boundary,triangles, ...
            'XYData',u, ...
            'ColorMap',usedmap, ...
            'Title','2D Density Plot HSV custom MAP');

%%%%%%% 2D PDE Map Plot with Mesh
figure();
pdeplot2dff(points,boundary,triangles, ...
            'XYData',u,'Mesh','on', ...
            'Title','2D Density Plot with Mesh');

%%%%%%% 2D PDE Map Plot with Mesh
figure();
pdeplot2dff(points,boundary,triangles, ...
            'XYData',u,'Mesh','on', 'XYStyle', 'off', ...
            'Title','2D Density Plot with Mesh without color (XYStyle)');

%%%%%%% 2D PDE Map Plot with Mesh without Colorbar
figure();
pdeplot2dff(points,boundary,triangles, ...
            'XYData',u,'Mesh','on','ColorBar','off', ...
            'Title','2D Density Plot with Mesh without Colbar');

%%%%%%% 3D Surface
figure();
pdeplot2dff(points,boundary,triangles, ...
            'XYData',u,'ZStyle','continuous', ...
            'Title','Classic 3D Surf Plot');

%%%%%%% 3D Surface with Mesh
figure();
pdeplot2dff(points,boundary,triangles, ...
            'XYData',u,'ZStyle','continuous','Mesh','on', ...
            'Title','3D Plot with Mesh');

%%%%%%% 3D Surface with Mesh without color
figure();
pdeplot2dff(points,boundary,triangles, ...
            'XYData',u,'ZStyle','continuous','Mesh','on', 'XYStyle','off', ...
            'Title','3D Plot with Mesh withou color (XYStyle)');

%%%%%%% Contour Plot
figure();
pdeplot2dff(points,boundary,triangles, ...
            'XYData',u,'Edge','on','Contour','on', 'XYStyle','off', ...
            'ColorBar','off', ...
            'Title','Contour BW');

%%%%%%% Contour Plot with colors
figure();
pdeplot2dff(points,boundary,triangles, ...
            'XYData',u,'Edge','on','Contour','on','CColor','auto', 'XYStyle','off', ...
            'Title','Contour with colors');

%%%%%%% Contour Plot Mixed with 2D Patch
figure();
pdeplot2dff(points,boundary,triangles, ...
            'XYData',u,'Edge','on','Contour','on', ...
            'Title','Contour mixed with Patch Plot');

%%%%%%% Contour Plot Mixed with 2D Patch and dashed lines
figure();
pdeplot2dff(points,boundary,triangles, ...
            'XYData',u,'Edge','on','Contour','on', 'CStyle','dashed', ...
            'Title','Contour mixed with Patch Plot');

%%%%%%% Contour Plot
figure();
pdeplot2dff(points,boundary,triangles, ...
            'XYData',u,'Edge','on', ...
            'Contour','on','CGridParam',[10,10],'CLevels',15, ...
            'Title','Contour Plot with ugly set GridParam');

%%%%%%% Contour Plot with Labels
figure();
[handles,clab]=pdeplot2dff(points,boundary,triangles, ...
               'XYData',u,'Edge','on','Contour','on','CLevels',5, ...
               'Title','Contour mixed with Patch Plot');


texth=clabel(clab,handles(2),'fontsize', 8);
for i=1:size(texth)
    textstr=get(texth(i),'String');
    textnum=str2double(textstr);
    textstrnew=sprintf('%0.3f', textnum);
    set(texth(i),'String',textstrnew);
end

%%%%%%% 2D PDE Map Plot axis limits
figure();
pdeplot2dff(points,boundary,triangles, ...
            'XYData',u,'XLim',[-0.25 0.75],'YLim',[0.25 1.25], ...
            'Title','2D Map Patch Plot with Limits');

%%%%%%% 3D Surface axis limits
figure();
pdeplot2dff(points,boundary,triangles, ...
            'XYData',u,'ZStyle','continuous', ...
            'XLim',[-0.25 0.75],'YLim',[0.25 1.25],'ZLim',[-0.1 0.1], ...
            'Title','3D Plot with xy and z-Limits');

%%%%%%% 3D Surface aspect ratio
figure();
pdeplot2dff(points,boundary,triangles, ...
            'XYData',u,'ZStyle','continuous','DAspect',[1 1 0.02], ...
            'Title','3D Plot Aspect ratio');

%%%%%%% 2D PDE Map Plot ColorRange
figure();
pdeplot2dff(points,boundary,triangles, ...
            'XYData',u,'ColorRange',[0 0.1], ...
            'Title','3D Plot Color Range set Test');

%%%%%%% Quiver
figure();
pdeplot2dff(points,boundary,triangles, ...
            'FlowData',v,'Edge','on', ...
            'Title','Quiver default');

axis tight;

%%%%%%% Quiver
figure();
pdeplot2dff(points,boundary,triangles, ...
            'FlowData',v,'FGridParam',[26,30],'Edge','on', ...
            'Title','Quiver Plot with GridParam');

axis tight;

%%%%%%% Superposition Quiver + 2D
figure();
pdeplot2dff(points,boundary,triangles, ...
            'XYData',u,'FlowData',v,'FGridParam',[26,30],'Edge','on', ...
            'Title','Quiver Plot combined with 2D Map Patch Plot');

axis tight;

%%%%%%% 3D Surface Lighting - possibly available on Matlab only
figure();
pdeplot2dff(points,boundary,triangles, ...
            'XYData',u,'ZStyle','continuous', ...
            'Title','Classic 3D Surf Plot Gouraud');

lighting gouraud;
view([-166,40]);
camlight('left');
grid;

pause(10);
close all;

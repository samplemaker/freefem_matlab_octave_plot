clear all;

%Specifies the path where the ffmatlib can be found
addpath('ffmatlib');

%Reads a FreeFem++ mesh created with the savemesh(Th,"mesh.msh"); command
[nv,nt,ns,points,triangles,boundary]=ffreadmesh('demo_mesh.msh');

%Reads the PDE data
fid=fopen('demo_data.txt','r');
data=textscan(fid,'%f','Delimiter','\n');
fclose(fid);
u=cell2mat(data);

%%%%%%% 2D Mesh Plot
figure();
pdeplot2dff(points,boundary,triangles,'Mesh','on');
ylabel('y');
xlabel('x');
title('Mesh without boundary');
axis tight equal;

%%%%%%% Including Boundary
figure();
pdeplot2dff(points,boundary,triangles,'Edge','on','Mesh','on');
ylabel('y');
xlabel('x');
title('Mesh with boundary');
axis tight equal;

%%%%%%% Including Boundary
figure();
pdeplot2dff(points,boundary,triangles,'Edge','on');
ylabel('y');
xlabel('x');
title('Boundary');
axis equal;

%%%%%%% 2D PDE Map Plot
figure();
pdeplot2dff(points,boundary,triangles,'XYData',u);
ylabel('y');
xlabel('x');
title('2D Density Plot');
axis tight equal;

%%%%%%% 2D PDE Map Plot
figure();
pdeplot2dff(points,boundary,triangles,'XYData',u,'ColorMap','hot');
ylabel('y');
xlabel('x');
title('2D Density Plot ColMap Hot');
axis tight equal;

%%%%%%% 2D PDE Map Plot with Mesh
figure();
pdeplot2dff(points,boundary,triangles,'XYData',u,'Mesh','on');
ylabel('y');
xlabel('x');
title('2D Density Plot with Mesh');
axis tight equal;

%%%%%%% 2D PDE Map Plot with Mesh without Colorbar
figure();
pdeplot2dff(points,boundary,triangles,'XYData',u,'Mesh','on','ColorBar','off');
ylabel('y');
xlabel('x');
title('2D Density Plot with Mesh without Colbar');
axis tight equal;

%%%%%%% 3D Surface
figure();
pdeplot2dff(points,boundary,triangles,'XYData',u,'ZStyle','on');
ylabel('y');
xlabel('x');
zlabel('u');
title('3D Plot Default');

%%%%%%% 3D Surface with Mesh
figure();
pdeplot2dff(points,boundary,triangles,'XYData',u,'ZStyle','on','Mesh','on');
ylabel('y');
xlabel('x');
zlabel('u');
title('3D Plot with Mesh');

%%%%%%% Contour Plot
figure();
pdeplot2dff(points,boundary,triangles,'XYData',u,'Edge','on','Contour','on');
ylabel('y');
xlabel('x');
title('Default Contour Plot');
axis tight equal;

%%%%%%% Contour Plot
figure();
pdeplot2dff(points,boundary,triangles,'XYData',u,'Edge','on','Contour','on','GridParam',[10,10],'Levels',15);
ylabel('y');
xlabel('x');
title('Contour Plot with ugly GridParam');
axis tight equal;

%%%%%%% 2D PDE Map Plot axis limits
figure();
pdeplot2dff(points,boundary,triangles,'XYData',u,'Title','Map plot with Title and limits','XLim',[-0.25 0.75],'YLim',[0.25 1.25]);
ylabel('y');
xlabel('x');
axis tight equal;
title('2D Map Plot with Limits');

%%%%%%% 3D Surface axis limits
figure();
pdeplot2dff(points,boundary,triangles,'XYData',u,'ZStyle','on','XLim',[-0.25 0.75],'YLim',[0.25 1.25],'ZLim',[-0.1 0.1]);
ylabel('y');
xlabel('x');
zlabel('u');
title('3D Plot with Limits');

%%%%%%% 3D Surface aspect ratio
figure();
pdeplot2dff(points,boundary,triangles,'XYData',u,'ZStyle','on','ZAspect',[1 1 0.02]);
ylabel('y');
xlabel('x');
zlabel('u');
title('3D Plot Aspect ratio');

%%%%%%% 2D PDE Map Plot ColorRange
figure();
pdeplot2dff(points,boundary,triangles,'XYData',u,'ColorRange',[0 0.1]);
ylabel('y');
xlabel('x');
axis tight equal;
title('3D Plot Color Range Test');

%%%%%%% Contour Plot
figure();
pdeplot2dff(points,boundary,triangles,'XYData',u,'Edge','on','Contour','value','Levels',5);
ylabel('y');
xlabel('x');
title('Contour Plot including Values');
axis tight equal;


%%%%%%%%%%%%%%%%%%%%%%%%%%% for the overlaid plot variant and quiver
fid=fopen('demo_flowdata.txt','r');
data=textscan(fid,'%f %f','Delimiter','\n');
fclose(fid);
v=cell2mat(data)';

%%%%%%% Quiver
figure();
pdeplot2dff(points,boundary,triangles,'FlowData',v,'Edge','on');
ylabel('y');
xlabel('x');
title('Quiver Default');
axis equal;

%%%%%%% Quiver
figure();
pdeplot2dff(points,boundary,triangles,'FlowData',v,'GridParam',[26,30],'Edge','on');
ylabel('y');
xlabel('x');
title('Quiver Plot with GridParam');
axis equal;

%%%%%%% Superposition Quiver + 2D
figure();
pdeplot2dff(points,boundary,triangles,'XYData',u,'FlowData',v,'GridParam',[26,30],'Edge','on');
ylabel('y');
xlabel('x');
title('Quiver Plot and 2D Color');
axis tight equal;

%%%%%%% Modified Contour Plot
figure();
pdeplot2dff(points,boundary,triangles,'XYData',v(1,:),'Edge','on','Contour',v(2,:),'Levels',15);
ylabel('y');
xlabel('x');
title('Contour Plot');
axis tight equal;

pause(10);
close all;

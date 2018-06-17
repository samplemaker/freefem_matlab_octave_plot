clear all;

%Specifies the path where the ffmatlib can be found
addpath('ffmatlib');

%Read a FreeFem++ mesh created with the savemesh(Th,"mesh.msh"); command
[nv,nt,ns,points,triangles,boundary]=ffreadmesh('demo_mesh.msh');

%Read the PDE data
fid=fopen('demo_data.txt','r');
data=textscan(fid,'%f','Delimiter','\n');
fclose(fid);
u=cell2mat(data);

%2D Mesh Plot
figure();
pdeplot2dff(points,triangles,boundary,'Mesh','on');
ylabel('y');
xlabel('x');
title('Mesh without boundary');
axis tight equal;

%Including Boundary
figure();
pdeplot2dff(points,triangles,boundary,'Edge','on','Mesh','on');
ylabel('y');
xlabel('x');
title('Mesh with boundary');
axis tight equal;

%Including Boundary
figure();
pdeplot2dff(points,triangles,boundary,'Edge','on');
ylabel('y');
xlabel('x');
title('Boundary');
axis equal;

%2D PDE Map Plot
figure();
pdeplot2dff(points,triangles,boundary,'XYData',u);
ylabel('y');
xlabel('x');
title('2D Density Plot');
axis tight equal;

%2D PDE Map Plot
figure();
pdeplot2dff(points,triangles,boundary,'XYData',u,'ColorMap','hot');
ylabel('y');
xlabel('x');
title('2D Density Plot');
axis tight equal;

%2D PDE Map Plot with Mesh
figure();
pdeplot2dff(points,triangles,boundary,'XYData',u,'Mesh','on');
ylabel('y');
xlabel('x');
title('2D Density Plot with Mesh');
axis tight equal;

%2D PDE Map Plot with Mesh without Colorbar
figure();
pdeplot2dff(points,triangles,boundary,'XYData',u,'Mesh','on','ColorBar','off');
ylabel('y');
xlabel('x');
title('2D Density Plot with Mesh');
axis tight equal;

%3D Surface
figure();
pdeplot2dff(points,triangles,boundary,'XYData',u,'ZStyle','on');
ylabel('y');
xlabel('x');
zlabel('u');
title('3D Plot');

%3D Surface with Mesh
figure();
pdeplot2dff(points,triangles,boundary,'XYData',u,'ZStyle','on','Mesh','on');
ylabel('y');
xlabel('x');
zlabel('u');
title('3D Plot with Mesh');

%Contour Plot
figure();
pdeplot2dff(points,triangles,boundary,'XYData',u,'Edge','on','Contour','on','Levels',15);
ylabel('y');
xlabel('x');
title('Contour Plot');
axis tight equal;

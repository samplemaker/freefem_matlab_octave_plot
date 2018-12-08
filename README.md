# How to plot FreeFem++ Simulations in Matlab&copy; and Octave

Once you have successfully simulated a PDE problem using FreeFem++ you may want to have a look at the simulation results from within Matlab or Octave. `ffmatlib` provides useful commands in order to load FreeFem++ meshes and simulation data and to call the underlying Matlab/Octave plot routines like `contour()`, `quiver()` as well as `patch()`.

![](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/cap_3d_gouraud.png)

## Getting started

  * Click on the button `Clone or download` (above) and then on the button `Download ZIP`
  * Unzip and change to the directory `demos` and run all FreeFem++ *.epd scripts to create simulation data for plotting
  * Run the matlab `*.m` demo files with Matlab or Octave

Hint: The ffmatlib functions are stored in the folder `ffmatlib`. Use the `addpath(path to ffmatlib)` command if you are working in a different directory.

## Demos

<a name="2dcapacitorexample"></a>

### 2D-Parallel Plate Capacitor

[capacitor_2d.m](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/demos/capacitor_2d.m)  
[capacitor_2d_p1.edp](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/demos/capacitor_2d_p1.edp)  

[Screenshot: 3D Patch](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/cap_3d_gouraud.png)  
[Screenshot: Mesh](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/cap_2d_mesh.png)  
[Screenshot: Contour and Quiver](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/cap_2d_contour.png)  
[Screenshot: 2D Patch with Mesh](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/cap_2d_patch.png)  
[Screenshot: Boundary and Labels](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/cap_2d_labels.png)  

<a name="convectexample"></a>

### 2D-Horizontal Roll Vortices

[convective_rolls.m](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/demos/convective_rolls.m)  
[convective_rolls.edp](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/demos/convective_rolls.edp)  

[Screenshot: Patch + Contour](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/convect_temperature_streamlines.png)  
[Screenshot: Patch + Quiver](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/convect_temperature_velocity.png)  
[Screenshot: 2D Patch + Mesh](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/convect_mesh.png)  
[Screenshot: 2D Streamlines](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/convect_streamlines.png)  

<a name="pdeplotexample"></a>

### 2D-Various Examples

[demo_pdeplot.m](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/demos/demo_pdeplot.m)  
[demo_pdeplot_2d_p1.edp](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/demos/demo_pdeplot_2d_p1.edp)  

[Screenshot: 2D Patch](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/l_shape_patch_mesh.png)  
[Screenshot: Contour](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/l_shape_patch_contour.png)  
[Screenshot: Quiver](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/l_shape_patch_quiver.png)  

<a name="3dcapacitorexample"></a>

### 3D-Parallel Plate Capacitor

[capacitor_3d.m](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/demos/capacitor_3d.m)  
[capacitor_3d.edp](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/demos/capacitor_3d.edp)  

[Screenshot: 3D Slice](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/cap_3d_slices.png)  
[Screenshot: 3D Vector field](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/cap_3d_spatial_vectorfield.png)  

<a name="3dcoilexample"></a>

### 3D-Toroidal Current

[magnetostatic3D.m](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/demos/magnetostatic3D.m)  
[magnetostatic3D.edp](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/demos/magnetostatic3D.edp)  
[torus.geo](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/demos/torus.geo)  

[Screenshot: 3D Vector field](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/toroid_3d_spatial_vectorfield.png)  
[Screenshot: Vector field - Slice](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/toroid_3d_spatial_vectorfield2.png)  

## Function Reference

| Name | Description |
| --- | --- |
| [ffpdeplot()](#ffpdeplotfct) | Creates contour(), quiver() as well as patch() plots from FreeFem++ 2D simulation data |
| [ffinterpolate()](#ffinterpolatefct) | Interpolates from 2D triangular mesh to 2D cartesian or curved grid |
| [fftri2grid()](#fftri2gridfct) | Interpolates from 2D triangular mesh to 2D cartesian or curved grid |
| [ffpdeplot3D()](#ffpdeplot3Dfct) | Creates cross-sections, quiver3() as well as boundary plots from FreeFem++ 3D simulation data |
| [ffreadmesh()](#ffreadmeshfct) | Reads FreeFem++ Mesh Files into Matlab/Octave |
| [ffreaddata()](#ffreaddatafct) | Reads FreeFem++ Data Files into Matlab/Octave |

<a name="ffpdeplotfct"></a>

## ffpdeplot()

`ffpdeplot()` is a function specially tailored to FreeFem++ that offers most of the features of the classic Matlab `pdeplot()` command. `contour()` plots (= 2D iso value) and `patch()` plots (= 2D map data) can be created as well as combinations of both. In addition `quiver()` plots (= 2D vector field) can be created and domain border edges can be displayed. The display of the flow data as well as the border edges is additive and can be superimposed on the contour data as well as the patch data or their combinations.

The FEM mesh is entered by vertex coordinates, the boundary values, and the triangles in terms of connectivity as provided by the FreeFem++ `savemesh(Th, "mesh_file.msh")` command. The simulation data can be entered as values at the mesh nodes.


#### Synopsis

```Matlab
[handles,varargout] = ffpdeplot(p,b,t,varargin)
```

#### Description / Name-Value Pair Arguments

The contents of the points `p`, boundaries `b` and triangles `t` arguments are explained in the section [ffreadmesh()](#ffreadmeshfct). `ffpdeplot()` can be called with name-value pair arguments as per following table:

| Parameter | Value |
| --- | --- |
| 'XYData' |     Data in order to colorize the plot |
|           |       FreeFem++ point data \| FreeFem++ triangle data |
| 'XYStyle' |    Coloring choice |
|           |       'interp' (default) \| 'off' |
| 'ZStyle' |     Draws 3D surface plot instead of flat 2D Map plot |
|           |       'continuous' \| 'off' (default) |
| 'ColorMap' |   ColorMap value or matrix of such values |
|           |       'cool' (default) \| colormap name \| three-column matrix of RGB triplets |
| 'ColorBar' |   Indicator in order to include a colorbar |
|            |      'on' (default) \| 'off' |
| 'CBTitle'   |   Colorbar Title  |
|          |        (default=[]) |
| 'ColorRange' | Range of values to adjust the color thresholds |
|          |        'minmax' (default) \| [min,max] |
| 'Mesh' |       Switches the mesh off / on |
|         |         'on' \| 'off' (default) |
| 'Boundary' |   Shows the boundary / edges |
|          |        'on' \| 'off' (default) |
| 'BDLabels' |   Draws boundary / edges with a specific label |
|          |        [] (default) | [label1,label2,...] |
| 'BDColors' |   Colorize boundary / edges with color (linked to 'BDLabels') |
|         |         'r' (default) | three-column matrix of RGB triplets |
| 'BDShowText' |   Shows the labelnumber on the boundary / edges (linked to 'BDLabels') |
|         |         'on' | 'off' (default) |
| 'BDTextSize' |   Size of labelnumbers on the boundary / edges |
|         |         scalar value greater than zero |
| 'BDTextWeight' |   Character thickness of labelnumbers on the boundary / edges |
|         |         'normal' (default) | 'bold' |
| 'Contour' |    Isovalue plot |
|           |       'off' (default) \| 'on' |
| 'CColor' |     Isovalue color |
|           |       [0,0,0] (default) \| RGB triplet three-element row vector \| 'r' \| 'g' \| 'b' \| |
| 'CXYData' |    Use extra (overlay) data to draw the contour plot |
|           |       FreeFem++ points \| FreeFem++ triangle data |
| 'CStyle'  |    Contour line style |
|           |       'patch' (default) \| 'patchdashed' \| 'patchdashedneg' \| 'monochrome' \| 'colormap' |
| 'CLevels' |    Number of isovalues used in the contour plot |
|           |       (default=10) |
| 'CGridParam' | Number of grid points used for the contour plot |
|         |         'auto' (default) \| [N,M] |
| 'Title' |      Title |
|          |        (default=[]) |
| 'XLim' |       Range for the x-axis |
|        |          'minmax' (default) \| [min,max] |
| 'YLim' |       Range for the y-axis |
|        |         'minmax' (default) \| [min,max] |
| 'ZLim' |       Range for the z-axis |
|         |         'minmax' (default) \| [min,max] |
| 'DAspect' |    Data unit length of the xy- and z-axes |
|          |        'off' \| 'xyequal' (default) \| [ux,uy,uz] |
| 'FlowData' |   Data for quiver plot |
|            |      FreeFem++ point data \| FreeFem++ triangle data |
| 'FGridParam' | Number of grid points used for quiver plot |
|             |     'auto' (default) \| [N,M] |

The return value `handles` contains handles to the plot figures. The return value `varargout` contains references to the contour labels.

#### Examples

First load the mesh and the simulation data:

```Matlab
[p,b,t]=ffreadmesh('capacitorp1.msh');
[u]=ffreaddata('capacitor_potential_p1only.txt');
[Ex,Ey]=ffreaddata('capacitor_field_p1only.txt');
```

2D Patch Plot (2D map / density) without boundary:
```Matlab
ffpdeplot(p,[],t,'XYData',u);
```

Plot of the domain boundary:
```Matlab
ffpdeplot(p,b,t,'Boundary','on');
```

2D Patch (2D Map or Density) Plot with boundary:
```Matlab
ffpdeplot(p,b,t,'XYData',u,'Mesh','on','Boundary','on');
```

3D Surf Plot:
```Matlab
ffpdeplot(p,b,t,'XYData',u,'ZStyle','continuous');
```

Contour Plot (isovalues):
```Matlab
ffpdeplot(p,b,t,'XYData',u,'Contour','on','Boundary','on');
```

Quiver Plot (vector field):
```Matlab
ffpdeplot(p,b,t,'FlowData',[Ex, Ey],'Boundary','on');
```

<a name="fftri2gridfct"></a>

## fftri2grid() / fftri2gridfast()

Interpolates the real valued or complex data `tu1`, `tu2` given on a triangular mesh defined by `tx` and `ty` onto a cartesian- or curved meshgrid defined by `x` and `y`. The data `tu2` is optional and can be omitted.<br>

#### Synopsis

```Matlab
[w1,[w2]] = fftri2grid (x,y,tx,ty,tu1,[tu2])
[w1,[w2]] = fftri2gridfast (x,y,tx,ty,tu1,[tu2])
```

#### Description

`fftri2grid()` uses barycentric interpolation. `tx`, `ty` must contain the triangle vertice coordinates. The arguments `tx`, `ty`, `tu1` and `tu2` must have a size of nTriangle-columns and 3 rows. The returned data `w1`, `w2` is the interpolation of `tu1`, `tu2` at the grid points defined by `x`, `y` and is real if `tu1`, `tu2` is real or complex if `tu1`, `tu2` is complex. The function returns `NaN's` if an interpolation point is outside the triangle mesh. `ffri2gridfast.c` if the faster mex implementation of this function and need to be build before use. `ffri2grid.m` is the slower matlab code variant. For more information see also [Notes on MEX Compilation](#notesoncompilation).<br>
Note that this is a library function and should not be used directly. To interpolate data use the wrapper function `ffinterpolate.m` instead.

#### Examples

```Matlab
[p,b,t]=ffreadmesh('capacitorp1.msh');
u=ffreaddata('capacitor_potential_p1only.txt');
xpts=p(1,:);
ypts=p(2,:);
xdata=[xpts(t(1,:)); xpts(t(2,:)); xpts(t(3,:))];
ydata=[ypts(t(1,:)); ypts(t(2,:)); ypts(t(3,:))];
udata=[u(t(1,:)), u(t(2,:)), u(t(3,:))].';
x=linspace(-5,5,500);
y=linspace(-5,5,500);
[X,Y]=meshgrid(x,y);
U=fftri2grid(X,Y,xdata,ydata,udata);
surf(X,Y,U,'EdgeColor','none');
view(3);
```

<a name="ffinterpolatefct"></a>

## ffinterpolate()

Interpolates the real valued or complex data `u1`, `u2` given on a triangular mesh defined by `p`, `b` and `t` onto a cartesian- or curved meshgrid defined by `x` and `y`. The data `u2` is optional and can be omitted.<br>

#### Synopsis

```Matlab
[w1,[w2]] = ffinterpolate (p,b,t,x,y,u1,u2)
```

#### Description

`ffinterpolate()` uses barycentric interpolation. The function returns `NaN's` if an interpolation point is outside the triangle mesh. The contents of the `p`, `b` and `t` arguments are explained in the section [ffreadmesh()](#ffreadmeshfct). The content of `u` is described in the section [ffreaddata](#ffreaddatafct). To improve runtime there is a MEX implementation of the interpolation section. If Matlab/Octave finds an executable of `ffri2gridfast.c` within its search path the faster C-implementation is used instead of `ffri2grid.m`.

#### Examples

<a name="ffpdeplot3Dfct"></a>

## ffpdeplot3D()

The purpose of the library function `ffpdeplot3D()` is to create cross-sections, to selectively plot boundaries and to create quiver3() plots from 3D simulation data. See also the example [3D-Parallel Plate Capacitor](#3dcapacitorexample). Note: This function is still under construction.

#### Synopsis

```Matlab
[] = ffpdeplot3D(p,b,t,varargin)
```

#### Description / Name-Value Pair Arguments

The contents of the points `p`, boundaries `b` and triangles `t` arguments are explained in the section [ffreadmesh()](#ffreadmeshfct). Although the function can be run as pure Matlab code it is strongly recommended to build the MEX file of the underlying library function in `fftet2gridfast.c` to improve execution speed. See also chapter [Notes on MEX Compilation](#notesoncompilation).

`ffpdeplot3D()` can be called with name-value pair arguments as per following table:

| Parameter | Value |
| --- | --- |
| 'XYData' |     Data in order to colorize the plot |
|           |       FreeFem++ data |
| 'XYStyle' |    Plot style for boundary |
|           |       'interp' (default) \| 'noface' \| 'monochrome' |
| 'Boundary' |    Shows the domain boundary / edges |
|           |       'on' (default) \| 'off' |
| 'BoundingBox' |   Shows the bounding box of a slice |
|           |       'on' \| 'off' (default) |
| 'BDLabels' |   Draws boundary / edges with a specific label |
|            |      [] (default) \| [label1,label2,...] |
| 'Slice'   |   3 point slicing plane definition  |
|          |        [] \| three-column matrix of [x,y,z] triplets |
| 'SGridParam' | Number of grid points used for the slice |
|          |        'auto' (default) \| [N,M] |
| 'Project2D' | View cross section in 2D |
|           |       'on' \| 'off' (default) |
| 'ColorMap' |       ColorMap value or matrix of such values |
|         |         'cool' (default) \| colormap name \| three-column matrix of RGB triplets |
| 'ColorBar' |   Indicator in order to include a colorbar |
|          |        'on' (default) \| 'off' |
| 'CBTitle'   |   Colorbar Title  |
|          |        (default=[]) |
| 'ColorRange' |    Range of values to adjust the color thresholds |
|          |        'minmax' (default) \| [min,max] |
| 'FlowData' |    Data for quiver3 plot |
|           |       FreeFem++ point data |
| 'FGridParam' |    Number of grid points used for quiver3 plot at cross-section |
|           |       'auto' (default) \| [N,M] |
| 'FGridParam3D' |    Number of grid points used for a spatial quiver3 plot |
|           |       'auto' (default) \| [N,M,L] |
| 'FLim3D' |    Bounding box for a spatial quiver3 plot |
|           |       'auto' (default) \| [xmin,xmax;ymin,ymax;zmin,zmax] |
| 'FMode3D' |    Arrow distribution choice |
|           |       'cartesian' (default) | 'random' |


#### Examples

First load the mesh and the simulation data:

```Matlab
[p,b,t]=ffreadmesh('cap3d.mesh');
[u]=ffreaddata('cap3dpot.txt');
[Ex,Ey,Ez]=ffreaddata('cap3dvec.txt');
```

To plot the domain boundaries the entire mesh is drawn with a monochrome face color:
```Matlab
ffpdeplot3D(p,b,t,'XYZStyle','monochrome');
```

In three dimensions some boundaries can be buried inside the domain and therefore invisible. In a FreeFem++ mesh all boundaries are labeled. If a specific label is specified it is possible to make the corresponding domain boundary visible. The following example plots only the boundaries/borders labeled with the numbers 30 and 31:
```Matlab
ffpdeplot3D(p,b,t,'BDLabels',[30,31],'XYZStyle','monochrome');
```

To color the domain boundaries according the PDE solution `u` run:
```Matlab
ffpdeplot3D(p,b,t,'XYZData',u,'ColorMap','jet');
```

`ffpdeplot3D` with the argument `Slice` draws slices (cross-sections) for the volumetric data `u`. A slicing plane is defined by the parallelogram spanned by the three corner points S1, S2 and S3. Multiple slicing planes can be put together (series of cross-sections). As an example the following statement sequence creates and draws two orthogonal cross-sections of the volumetric data `u`:
```Matlab
S1=[-0 0.375 0.0; ...
    0.375 0 0.0];
S2=[0.0 0.375 0.5; ...
    0.375 0 0.5];
S3=[0.75 0.375 0.0; ...
    0.375 0.75 0.0];
ffpdeplot3D(p,b,t,'XYZData',u,'Slice',S1,S2,S3,'Boundary','off','ColorMap','jet(200)', ...
            'SGridParam',[300,300],'BoundingBox','on')
```

A cross-section can also be viewed in a 2D projection:
```Matlab
S1=[-0 0.375 0.0];
S2=[0.0 0.375 0.5];
S3=[0.75 0.375 0.0];
ffpdeplot3D(p,b,t,'XYZData',u,'Slice',S1,S2,S3,'SGridParam',[300,300], 'Project2D', 'on', ...
            'Boundary','off','ColorMap',jet(200),'ColorBar','on');
```

The following command plots a cross-section and additionally draws the complete mesh (the mesh facets are transparent):
```Matlab
ffpdeplot3D(p,b,t,'XYZData',u,'Slice',S1,S2,S3,'XYZStyle','noface','ColorMap','jet')
```

Three dimensional vector fields can be plotted over a cross-section. The following example creates a quiver3() plot for a cross-section defined by the three points S1, S2, and S3:
```Matlab
ffpdeplot3D(p,b,t,'FlowData',[Ex,Ey,Ez],'Slice',S1,S2,S3,'Boundary','on','BoundingBox','on', ...
            'BDLabels',[30,31],'XYZStyle','monochrome');
```

If the parameter specification for the slice is omitted the vector field is instead drawn on a rectangular grid defined by the FGridParam3D and FLim3D parameters:
```Matlab
ffpdeplot3D(p,b,t,'FlowData',[Ex,Ey,Ez],'FGridParam3D',[8,8,5],'FLim3D', ...
            [0.125,0.625;0.125,0.625;0.1,0.4],'BDLabels',[30,31],'XYZStyle','monochrome');
```

<a name="ffreadmeshfct"></a>

## ffreadmesh()

Reads a FreeFem++ mesh file created by the FreeFem++ `savemesh(Th,"2dmesh.msh")` or `savemesh(Th3d,"3dmesh.mesh")` command to the Matlab/Octave workspace.

#### Synopsis

```Matlab
[p,b,t,nv,nbe,nt,labels] = ffreadmesh(filename)
```

#### Description

A mesh file consists of three parts:  

1. a mesh point list containing the nodal coordinates  
2. a list of boundary elements including the boundary labels  
3. list of triangles or tetrahedra defining the mesh in terms of connectivity  

These three blocks are stored in the variables `p`, `b` and `t` respectively.

**2D FreeFem++(*.msh)**

| Parameter | Value |
| --- | --- |
| p | Matrix containing the nodal points |
| b | Matrix containing the boundary edges |
| t | Matrix containing the triangles |
| nv | Number of points/vertices in the Mesh (Th.nv) |
| nt | Number of triangles in the Mesh (Th.nt) |
| nbe | Number of (boundary) edges (Th.nbe) |
| labels | Labels found in the mesh file |

**3D INRIA Medit(*.mesh)**

| Parameter | Value |
| --- | --- |
| p | Matrix containing the nodal points |
| b | Matrix containing the boundary triangles |
| t | Matrix containing the tetrahedra |
| nv | Number of points/vertices in the Mesh (nbvx, Th.nv) |
| nt | Number of tetrahedra in the Mesh (nbtet, Th.nt) |
| nbe | Number of (boundary) triangles (nbtri, Th.nbe) |
| labels | Labels found in the mesh file |

#### Examples

Reads a mesh file into the Matlab/Octave workspace:
```Matlab
[p,b,t,nv,nbe,nt,labels]=ffreadmesh('capacitorp1.msh');
fprintf('[Vertices nv:%i; Triangles nt:%i; Boundary Edges nbe:%i]\n',nv,nt,nbe);
fprintf('NaNs: %i %i %i\n',any(any(isnan(p))),any(any(isnan(t))),any(any(isnan(b))));
fprintf('Sizes: %ix%i %ix%i %ix%i\n',size(p),size(t),size(b));
fprintf('Labels found: %i\n',nlabels);
fprintf(['They are: ' repmat('%i ',1,size(labels,2)) '\n'],labels);
```

<a name="ffreaddatafct"></a>

## ffreaddata()

Reads a FreeFem++ data file created by the FreeFem++ [scripts](#exportfromff) to the Matlab/Octave workspace.

```Matlab
[varargout] = ffreadmesh(filename)
```

Note: The data can be real or complex.

#### Examples

Reads scalar data and a two dimensional vector field to the Matlab/Octave workspace:

```Matlab
[u]=ffreaddata('capacitor_potential_p1only.txt');
[Ex,Ey]=ffreaddata('capacitor_field_p1only.txt');
```

<a name="exportfromff"></a>

## Exporting data from FreeFem++

In order to create a plot from a FreeFem++ PDE simulation with Matlab / Octave the mesh and the FE-Space function must be written into ASCII data files.  

A FreeFem++ mesh can easily be exported via the built-in `savemesh` command as per follows:

Saves a 2D Mesh:
```Matlab
savemesh(Th,"capacitorp1.msh");
```

Saves a 3D Mesh:
```Matlab
savemesh(Thn3d,"cap3d.mesh");
```

FE-Space functions must also be exported to ASCII files. However two different data formats can be read by the `ffmatlib`:  

1.) The FE-Space function is given on each mesh node (preferred method)  
2.) The FE-Space function is given as a triangle/vertice list  

In this sense P1-Element simulations can be written directly into an ASCII file because the degree of freedom (=Vh.ndof) of a P1-space function is equal to the number of nodes in the mesh. In order to export simulation data created from higher order FE-Elements the FE-Space data must first be converted to P1-Element data. For this purpose you can use the `=` operator in FreeFem++. I.e. you can simply write `Vh u=v` where `u` is from type P1 and `v` is from type P2.

Saves the P1-Element scalar function `u`:
```Matlab
ofstream file("capacitor_potential_p1only.txt"); 
for (int j=0; j<u[].n; j++)
   file << u[][j] << endl;
}
```

Saves the P1-Element 2D vector field `[Ex,Ey]`:
```Matlab
ofstream file("capacitor_field_p1only.txt");
for (int j=0; j<Ex[].n; j++)
   file << Ex[][j] << " " << Ey[][j] << endl;
}
```

Note: If periodic boundary conditions are set you have to cycle through all mesh vertices and write the interpolation of the FE-Space function for each vertex:

```Matlab
ofstream file("periodic.txt");
int nbvertices = Th.nv;
for (int i=0; i<nbvertices; i++){
   file << u(Th(i).x, Th(i).y) << "\n";
}
```

Finally, in order to import those files into Matlab/Octave see the sections [ffreadmesh](#ffreadmeshfct) and [ffreaddata](#ffreaddatafct).

<a name="notesoncompilation"></a>

## Notes on MEX Compilation

Octave/Linux:<br>
In Octave under a Linux system with gcc as compiler go into the folder `./ffmatlib/` and invoke the following two commands:

`mkoctfile --mex -Wall  fftri2gridfast.c`  
`mkoctfile --mex -Wall  fftet2gridfast.c`

Matlab/Windows:<br>
In Matlab under a Windows system with Microsoft Visual Studio as compiler invoke:

`mex  fftri2gridfast.c -v -largeArrayDims COMPFLAGS='$COMPFLAGS /Wall'`  
`mex  fftet2gridfast.c -v -largeArrayDims COMPFLAGS='$COMPFLAGS /Wall'`

Note: If your build fails with Microsoft Visual Studio 10 ensure to enable the C99-standard or you can try to change the file name into *.cpp, forcing MVSD to use a C++ compiler.<br>

Since release R2018 Matlab implements a new memory layout to handle complex numbers. Old MEX files are incompatible with the new Interleaved Complex API. New ffmatlib MEX-files are available but completely untested. If you want to be a pilot testing volunteer and if you want to segfault your computer feel free to send me an email.

## Notes on Hardware Acceleration

It should be emphasized that the responsiveness of the plots is highly dependent on the degree of freedom of the PDE problem and the capabilities of the graphics hardware. For larger problems (lots of thousand of vertices), a dedicated graphics card rather than on board graphics should be used. Make use of hardware acceleration extensively. Some notes on trouble shooting and tweaking:<br><br>
 If `get(gcf,'RendererMode')` is set to auto Matlab/Octave will decide on its own which renderer is the best for the current graphic task.

  * `get(figure_handle,'Renderer')` returns the current figure() renderer
  * `set(figure_handle,'Renderer','OpenGL')` forces a figure() to switch to OpenGL
  * `set(figure_handle,'Renderer','painters')` forces a figure() to switch to vector graphics

Generally OpenGL can be considered to be faster than painters. To get an OpenGL info you can type `opengl info` within Matlab. Ensure the line `Software` shows `false` otherwise OpenGL will run in software mode. If hardware-accelerated OpenGL is available on the system, you can manually change the modes using the `opengl software` and `opengl hardware` commands.

## Software

  * [FreeFem++][freefem]
  * [Octave][octave]
  * [Matlab][matlab]

[freefem]:    http://www.freefem.org/
             "FreeFem++ solver for partial differential equations"
[octave]:     https://www.gnu.org/software/octave/
             "GNU Octave scientific programming language"
[matlab]:     https://www.mathworks.com/
             "Matlab scientific programming language"

## Acknowledgments

Many thanks to David Fabre ([StabFEM](https://github.com/erbafdavid/StabFem)) for feature and implementation suggestions and code review.

## The License

GPLv3+

Have fun ...

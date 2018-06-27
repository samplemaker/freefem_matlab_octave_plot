# How to plot FreeFem++ Simulations in Matlab and Octave

Once you have successfully simulated a PDE problem using FreeFem++ you may want to have a look at your simulation results from within Matlab&copy; or Octave. In this repository you will find some code snippets showing how to make this wish come true.

## Basic Theory

The purpose of a FE mesh is to describe a spatial geometry - for example a CAD object. A FE mesh is built from mesh elements. Mesh elements have the shape of a triangle or rectangle, can vary in size and the nodes (vertices) are not necessarily bound to a rectangular grid. If we look at a meshed object the surface may be colored according to the solution of a PDE which - for the sake of simplicity - shall be given at the mesh nodes only. Finally parts of the FE mesh may be obscured depending on the point of view.<br>
At the other hand the `patch()` command which is built into Matlab&copy; and Octave renders a set of polygons (=facets, patches). The patch drawing primitives are defined by a color value and the spatial coordinates at the polygon vertices. We can associate these drawing primitives with FE mesh elements and looking towards a plot implementation all the prerequisites stated will be met.<br>
In the current implementation the mesh data (i.e. the PDE solution at the mesh nodes and the nodal coordinates) has to be written to a text file via the FreeFem++ script. In order to plot the problem with the `patch()` command this file must be parsed and processed by the `ffmatlib` library. Basically, the `ffmatlib` library splits and rearranges the continuous vertice data because `patch()` expects its input data to be clustered patch-wise. A detailed documentation of the `patch()` command can be found at [1](https://de.mathworks.com/help/matlab/ref/patch.html), [2](https://de.mathworks.com/help/matlab/visualize/introduction-to-patch-objects.html) and [3](https://de.mathworks.com/help/matlab/creating_plots/how-patch-data-relates-to-a-colormap.html).

![](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/3dsurf_10.png)

## Getting started

  * Click on the button "Clone or download" and then on the button `Download ZIP`
  * Unzip and change to the directory `demos` and run all FreeFem++ *.epd scripts to create simulation data for plotting
  * Run the matlab `*.m` demo files with Matlab / Octave

Hint: The ffmatlib functions are stored in the folder `ffmatlib`. Use the `addpath(path to ffmatlib)` command if you are working in a different directory.

## Overview

**Plot Examples**

| Plot type | Description |
| --- | --- |
| [pdeplot2dff()](#pdeplot2dffexample) | Wrapper function reading and plotting FreeFem++ 2D mesh and FreeFem++ data |
| [2D Map Plot](#2ddensity) | Creates a 2D color "density" plot of a function R<sup>2</sup> &rarr; R |
| [Surf Plot of 2D functions](#2dsurf) | Creates a 3D surface plot of a function R<sup>2</sup> &rarr; R |
| [2D Mesh Plot](#2dmesh) | Creates a 2D mesh plot |
| [2D Contour Plot](#2dcontour) | Interpolation to a rectangular grid and 2D isovalue plot |
| [2D Vector Plot](#2dvector) | Interpolation to a rectangular grid and 2D quiver plot of a 2D vector field |
| [3D Boundary Plot](#3dboundaryplot) | Creates a plot of a 3D domain boundary colored by a function R<sup>3</sup> &rarr; R |
| [3D Crop Boundary Plot](#3dslice) | Creates a plot of a 3D domain boundary cross section |
| [3D Boundary Mesh Plot](#3dboundaryplot) | Creates a mesh plot of a 3D domain boundary |
| [3D Slice Plot](#3dslice) | Creates the plot of a cross section of a 3D mesh colored by a function R<sup>3</sup> &rarr; R |
| [3D Slice (smoothed) Plot](#3dipocrosssection) | Interpolation of a 3D cross section to a 2D rectangular grid colored by a value |
| [3D Slice Vectorfield](#3dvectorcrosssection)  | Interpolation of a 3D cross section to a 2D rectangular grid including quiver3 vector plot |

**Functions**

| Name | Description |
| --- | --- |
| [pdeplot2dff()](#pdeplot2dfffct) | Wrapper function reading and plotting FreeFem++ 2D mesh and FreeFem++ data |
| [ffread2patch()](#general2d3dboundaryplots) | Reads FreeFem++ simulation results and converts the vertex data into patch plot data |
| [ffreadfile()](#readingffsimulationresults) | Reads one or two FreeFem++ simulation result files |
| [fftri2patch()](#general2d3dboundaryplots) | Converts vertex/triangle data into patch plot data |
| [slicetet2patch()](#cutting3dproblems) | Cuts 3D mesh elements (tetrahedrons) and converts the cross section into patch plot data |
| [slicetet2data()](#3dcrosssectioninterpolation) | Cuts 3D mesh elements (tetrahedrons) and returns all affected tetrahedrons |
| [fftet2grid()](#3dcrosssectioninterpolation) | Interpolates from 3d tetrahedral mesh to a rectangular grid |
| [fftet2gridfast()](#3dcrosssectioninterpolation) | MEX replacement for `fftet2grid.m` which is much faster |
| [slicebd2patch()](#cutting3dproblems) | Cuts the boundary data and converts remaining rest into patch plot data |
| [slicer_gui()](#slicergui) | Graphical user interface for cuttig FreeFem++ simulations |
| [fftri2grid()](#2dinterpolationcountourquiver) | Interpolates from 2D triangular mesh to 2D rectangular grid |
| [fftri2gridfast()](#2dinterpolationcountourquiver) | MEX replacement for `fftri2grid.m` which is much faster |

## Matlab / Octave Plot Examples

<a name="pdeplot2dffexample"></a>

### pdeplot2dff()

`pdeplot2dff()` is a customized FreeFem++ wrapper function to implement some of the classic `pdeplot()` features.

[demo_pdeplot.m](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/demos/demo_pdeplot.m)  

[Screenshot: PDEplot](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/pdeplot1.png)  

<a name="2ddensity"></a>

### 2D Map (Density) Plot of a 2D Function

The 2D plot examples focus on displaying functions of the type R<sup>2</sup> &rarr; R or 2D meshes.

[demo1_getstarted1.m](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/demos/demo1_getstarted1.m)  

[Screenshot: Density](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/2dgetstarted1.png)  

<a name="2dsurf"></a>

### Surf Plot (3D Plot of a 2D Function)

The 2D plot examples focus on displaying functions of the type R<sup>2</sup> &rarr; R or 2D meshes.

[demo1_getstarted2.m](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/demos/demo1_getstarted2.m)  

[Screenshot: Surf](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/2dgetstarted2.png)  

<a name="2dmesh"></a>

### Various 2D Plots and 2D Mesh Plots

The 2D plot examples focus on displaying functions of the type R<sup>2</sup> &rarr; R or 2D meshes.

[demo2_plot2d.m](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/demos/demo2_plot2d.m)  

[Screenshot: Density](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/2ddensity.png)  
[Screenshot: Surf](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/2dsurf.png)  
[Screenshot: Mesh](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/2dmesh.png)  

<a name="3dboundaryplot"></a>

### 3D Boundary Plots and 3D Mesh Plots

The 3D plot examples focus on displaying functions of the type R<sup>3</sup> &rarr; R (i.e. a 3D object boundary colored with a scalar value like a temperature) or 3D mesh surfaces.

[demo3_plot3dbd.m](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/demos/demo3_plot3dbd.m)  

[Screenshot: Surf Plot](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/3dsurf_2.png)  
[Screenshot: Mesh Surface](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/3dmesh.png)

<a name="3dslice"></a>

### 3D Cross Sections without Interpolation

To make the inside visible it is possible to cut a 3D FreeFem++ simulation along a slicing plane.

[demo4_slice3d.m](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/demos/demo4_slice3d.m)  
[demo4_start_slicer_gui.m](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/demos/demo4_start_slicer_gui.m)  

[Screenshot: Slice](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/3dsurf_slice7.png)  
[Screenshot: Boundary](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/3dsurf_slice8.png)  
[Screenshot: Cross Section](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/3dsurf_slice9.png)  
[Screenshot: GUI-Slicer](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/GUI_Slicer.png)  

<a name="2dcontour"></a>

### 2D Isovalue Plot (Contour)

The following example creates an Isolevel or Isovalue diagram. The actual plot is made with the Matlab / Octave command `contour()`.

[demo5_isovalues_2d.m](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/demos/demo5_isovalues_2d.m)  

[Screenshot: Isovalues](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/2dsurf_isovalues.png)  
[Screenshot: Surf](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/2dsurf_surf.png)  

Note: Consider using a Mex file for significant runtime improvement: `./ffmatlib/fftri2gridfast.c`

<a name="2dvector"></a>

### 2D Vector Fields (Quiver)

A 2D vector field plot shows arrows of a 2D vector field. Such a plot will be created in the following example. The actual plot is done with the Matlab / Octave command `quiver()`.

[demo6_vector_2d.m](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/demos/demo6_vector_2d.m)  

[Screenshot: Quiverplot](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/2dvectorfield.png)  

<a name="3dipocrosssection"></a>

### Interpolation of a 3D Mesh on a Cross Section

The cross section of a 3D mesh basically returns tetrahedra that are cut or touched by the slicing plane. However, the surface of such a cross section is rough. The following example interpolates the cross section on a rectangular gridded slicing plane and thus creates a smooth representation of the PDE solution on the slicing plane. The cross section can be plotted using the Matlab / Octave command `surf()`.

[demo7_slice3d_2dgrid.m](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/demos/demo7_slice3d_2dgrid.m)  

[Screenshot: Plane](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/3dplanedefinition.png)  
[Screenshot: Cross section](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/3dinterpolation.png)  
[Screenshot: Projection](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/3dcrossprojection.png)  

Note: Consider using a Mex file for significant runtime improvement: `./ffmatlib/fftet2gridfast.c`

<a name="3dvectorcrosssection"></a>

### 3D Vector Fields on Cross Sections (Quiver3)

A 3D vector field can be visualized with the Matlab / Octave command `quiver3()`. The following example creates a cross section of a 3D FreeFem++ simulation result. The vector field is displayed at the cross section with the `quiver3()` command.

[demo8_slice3d_2dgrid_vectors.m](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/demos/demo8_slice3d_2dgrid_vectors.m)  

[Screenshot: 3D Vector](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/3dvectorfield.png)  

Note: Consider using a Mex file for significant runtime improvement: `./ffmatlib/fftet2gridfast.c`

<a name="exportfromff"></a>

## Writing simulation results with FreeFem++ scripts

### 2D Problems

From within the FreeFem++ script write the mesh elements (triangles defined by its vertices) and the solution of the PDE at the nodes respectively:

```Matlab
ofstream tridata ("export_tri.txt");
for (int i=0; i<Th.nt; i++){
  for (int j=0; j<3; j++){
    tridata << Th[i][j].x << ";"
            << Th[i][j].y << ";"
            << u[][Vh(i,j)] << "\n";
  }
}
```

### 2D Vector Fields

Typically used in order to produce quiver plots:

```Matlab
ofstream file ("temp_demo6_vector.txt");
for (int i = 0; i < Th.nt; i++){
  for (int j = 0; j < 3; j++){
    file << Th[i][j].x << ";"
         << Th[i][j].y << ";"
         << uh[][Vh(i,j)] << ";"
         << vh[][Vh(i,j)] << "\n";
  }
}
```

### 3D Boundary Values

If the domain boundary (surface) is to be displayed write the boundary elements:

```Matlab
int idx;
int nbelement=Th3d.nbe;
ofstream bddata ("boundary_file.txt");
for (int k=0;k<nbelement;++k){
  for (int num=0;num<3;num++){
    idx=Th3d.be(k)[num];
    bddata << Th3d(idx).x << ";"
           << Th3d(idx).y << ";"
           << Th3d(idx).z << ";"
           <<  u(Th3d(idx).x,Th3d(idx).y,Th3d(idx).z) << "\n";
  }
}
```

### 3D Cross Sections (Slicing)

A cross section can be made from 3D boundary data as well as from 3D mesh data. Mesh elements (tetrahedra) can be written as:

```Matlab
ofstream tetdata ("tet_file.txt");
for (int i=0; i<Th3d.nt; i++){
  for (int j=0; j<4; j++){
    tetdata << Th3d[i][j].x << ";"
            << Th3d[i][j].y << ";"
            << Th3d[i][j].z << ";"
            << u[][Vh(i,j)] << "\n";

  }
}

```

### Higher order Vector Fields in 3D

If multiple components have to be exported they can be append at the end (here a temperature u and a heat flux vector field (qx,qy,qz) is exported):

```Matlab
ofstream tetdata ("tet_file.txt");
for (int i=0; i<Th3d.nt; i++){
  for (int j=0; j<4; j++){
    tetdata << Th3d[i][j].x << ";"
            << Th3d[i][j].y << ";"
            << Th3d[i][j].z << ";"
            << qx[][Vh(i,j)] << ";"
            << qy[][Vh(i,j)] << ";"
            << qz[][Vh(i,j)] << ";"
            << u[][Vh(i,j)] << "\n";
  }
}
```

## FFMATLIB Function Reference

<a name="pdeplot2dfffct"></a>

### pdeplot2dff()

`pdeplot2dff()` is a function specially tailored to FreeFem++ simulation data that offers most of the features of the classic `pdeplot()` command. The FEM-Mesh is entered by vertex coordinates, the boundary values, and the triangle definition as provided by the `savemesh(Th, "mesh_file.msh")` command. The simulation data can be entered either as point data (native support for P1 simulation data) or as interpolation at the nodes (workaround for
P1, P2 and other FEM simulation results).

```Matlab
[handles] = pdeplot2dff(points,boundary,triangles,varargin)
```

| Parameter | Value |
| --- | --- |
| 'XYData' |     Data in order to colorize the plot |
|           |       FreeFem++ point data \| FreeFem++ triangle data |
| 'XYStyle' |    Coloring choice |
|           |       'interp' \| (default) \| 'off' |
| 'ZStyle' |     Draws 3D surface plot instead of flat 2D Map plot |
|           |       'continuous' \| 'off' |(default) |
| 'ColorMap' |   ColorMap value or matrix of such values |
|           |       'cool' \| (default) \| colormap name \| three-column matrix of RGB triplets |
| 'ColorBar' |   Indicator in order to include a colorbar |
|            |      'on' \| (default) \| 'off' |
| 'ColorRange' | Range of values to ajust the color thresholds |
|          |        'minmax' \| (default) \| [min,max] |
| 'Mesh' |       Switches the mesh off / on |
|         |         'on' \| 'off' \| (default) |
| 'Edge' |       Shows the PDE boundary / edges |
|          |        'on' \| 'off' \| (default) |
| 'Contour' |    Isovalue plot |
|           |       'off' \| (default) \| 'on' |
| 'CColor' |     Color if level curves |
|           |       'off' \| (default) \| 'on' |
| 'CLevels' |    Number of isovalues used in the contour plot |
|           |       (default=10) |
| 'CGridParam' | Number of grid points used for the contour plot |
|         |         'auto' \| (default) \| [N,M] |
| 'Title' |      Title |
|          |        (default=[]) |
| 'XLim' |       Range for the x-axis |
|        |          'minmax' \| (default) \| [min,max] |
| 'YLim' |       Range for the y-axis |
|        |         'minmax' \| (default) \| [min,max] |
| 'ZLim' |       Range for the z-axis |
|         |         'minmax' \| (default) \| [min,max] |
| 'DAspect' |    Data unit length of the xy- and z-axes |
|          |        'off' \| 'xyequal' \| (default) \| [ux,uy,uz] |
| 'FlowData' |   Data for quiver plot |
|            |      FreeFem++ point data \| FreeFem++ triangle data |
| 'FGridParam' | Number of grid points used for quiver plot |
|             |     'auto' \| (default) \| [N,M] |


Examples:

In order to read a mesh file created by `savemesh(Th,"mesh.msh")` use the command:

```Matlab
[nv,nbe,nt,points,boundary,triangles] = ffreadmesh(filename)
```

A 2D Patch Plot:
```Matlab
[nv,nbe,nt,points,boundary,triangles]=ffreadmesh('demo_mesh.msh');
fid=fopen('demo_data.txt','r');
data=textscan(fid,'%f','Delimiter','\n');
fclose(fid);
u=cell2mat(data);
pdeplot2dff(points,boundary,triangles,'XYData',u,'Mesh','on','Edge','on');
```

A Plot of the Domain Boundary:
```Matlab
pdeplot2dff(points,boundary,triangles,'Edge','on');
```

A 3D Surf plot:
```Matlab
pdeplot2dff(points,boundary,triangles,'XYData',u,'ZStyle','continuous');
```

Contour Plot:
```Matlab
pdeplot2dff(points,boundary,triangles,'XYData',u,'Edge','on','Contour','on');
```

Quiver Plot:
```Matlab
pdeplot2dff(points,boundary,triangles,'FlowData',v,'Edge','on');
```

<a name="readingffsimulationresults"></a>

### Reading FreeFem++ Simulation Results

The `ffreadfile()` command can be used to read FreeFem++ ascii log files previously written in FreeFem++ simulation runs. It is possible to read one or two files. The `ffreadfile()` command does not interpret the input data, the data is merely passed through. Hence the data can be a 3D border, a 2D mesh or 3D mesh data. The output data is stored in one or two matrices, which can be processed with further ffmatlib commands. The number of columns to be read is determined by the format identifier. The following call reads a 3D mesh as well as 3D boundary data, both of the form (x,y,z,c):

```Matlab
[file1data,file2data] = ffreadfile ('File1','filename1.txt','File2','filename2.txt', ...
                                    'Delimiter',';','Format','%f %f %f %f')');
```

Read 2D mesh data of the form (x,y,c):

```Matlab
[filedata] = ffreadfile('File1','filename.txt','Delimiter',';','Format','%f %f %f');
```

<a name="general2d3dboundaryplots"></a>

### General 2D plotting and 3D Boundary Value Plots

The function `ffread2patch()` reads a FreeFem++ simulation log file and rearranges its content in such a way that diagrams can be created with the `patch()` command. It's arguments depend on the number of columns and on the separation character. `ffread2patch()` can process 2D mesh elements or 3D boundary data (triangles):

```Matlab
[X,Y,C, ...] = ffread2patch('filename.txt','Delimiter',';','Format','auto');
```
The number of columns can be set explicitly via the format identifier:

```Matlab
[X,Y,Z,C] = ffread2patch('filename.txt','Delimiter',';','Format','%f %f %f %f');
```
X, Y and C are matrices which can be supplied to `patch()`.

Instead of using the `ffread2patch()` command, the reading and conversion process can be divided into two different steps by the following commands:

```Matlab
[tridata] = ffreadfile('File1','filename.txt','Delimiter',';','Format','%f %f %f');
```

```Matlab
[X,Y,C] = fftri2patch(tridata);
```

<a name="cutting3dproblems"></a>

### Cutting 3D Problems without interpolation

It is possible to cut both 3D boundary and 3D mesh element data. When the boundary data of a problem is cropped, all triangles before the cutting plane are removed. Slicing of a 3D mesh data yields a set of tetrahedrons cut or touched by the slicing plane (cross section). The command `slicebd2patch()` cuts the boundary data and `slicetet2patch()` intersects the 3D mesh data. Both commands convert the output data into patch plot data.

```Matlab
[BX,BY,BZ,BC] = slicebd2patch (boundarydata,S1,S2,S3);
```

```Matlab
[SX,SY,SZ,SC] = slicetet2patch (meshdata,S1,S2,S3);
```

The input arguments `S1..S3` contain (x,y,z) coordinates of three points defining the slicing plane. The output data `(BX,BY,BZ,BC)` and `(SX,SY,SZ,SC)` or the superposition of both `([SX BX],[SY BY],[SZ BZ],[SC BC])` can be plot with the `patch()` command.<br>
The cross section, which consists of a set of tetrahedra, has a rough surface. Note that it is also possible to smooth the cross section with the `fftet2grid()` command.

<a name="slicergui"></a>

### Slicer_GUI

There's a graphical user interface that makes it easier to create cross sections of PDE problems. The function call is:

```Matlab
slicer_gui(bddatafile,tetdatafile);
```

<a name="2dinterpolationcountourquiver"></a>

### 2D Mesh Interpolation on a Rectangular Grid

In order to create 2D iso-level curve plots with the Matlab / Octave command `contour()` or vector field plots with the command `quiver()`, the simulation data must be interpolated on a rectangular grid. First load the data with the `ffreadfile()` command into a matrix - lets say `tridata`. Let X and Y define rectangular mesh grid. `tridata` can be interpolated on the grid with the command `fftri2grid()`.

```Matlab
C = fftri2grid(tridata,X,Y);
```

```Matlab
[U,V] = fftri2grid(tridata,X,Y);
```

If grid points are outside from any triangles the function will return NaN's. The method used is a Barycentric Coordinate interpolation.

Hint: To significantly speed up `fftri2grid()`, the interpolation routine is available as MEX implementation. The source code must be compiled before use. Take a look at the section [Notes on compilation](#notesoncompilation). The function is a 1:1 replacement and can be called via `fftri2gridfast(meshdata,X,Y,Z)`.

<a name="3dcrosssectioninterpolation"></a>

### 3D Mesh Interpolation on a Rectangular Grid

To create plots with a smooth cross section or 3D vector field plots using the `quiver3()` 3D command, the simulation data must be interpolated on a rectangular 3D grid in advance. Single value simulation results as well as higher order vector fields can be processed:

```Matlab
[C] = fftet2grid(meshdata,X,Y,Z);
```

If `meshdata` contains multiple data columns the function is to be called via:

```Matlab
[V1,V2, ...] = fftet2grid(meshdata,X,Y,Z);
```

The set (X, Y, Z) defines a rectangular 3D grid to which the PDE solution is to be interpolated. `meshdata` is a FreeFem++ simulation result that can be read with `ffreadfile()`. To improve runtime, it is better to call `slicetet2data()` before in order to reduce the mesh data. The method used is a Barycentric Coordinate interpolation. 

Hint: To significantly speed up `fftet2grid()`, the interpolation routine is available as MEX implementation. The source code must be compiled before use. Take a look at the section [Notes on compilation](#notesoncompilation). The function is a 1:1 replacement and can be called via `fftet2gridfast(meshdata,X,Y,Z)`.

For more information see section [Cutting 3D Problems](#cutting3dproblems)

<a name="notesoncompilation"></a>

## Notes on Compilation

Octave:<br>
In Octave seek to the folder `./ffmatlib/` and invoke the command 

`mkoctfile --mex -Wall fftet2gridfast.c`  
`mkoctfile --mex -Wall fftri2gridfast.c`

Windows:<br>
Under Windows with Microsoft Visual Studio invoke 

`mex fftet2gridfast.c -v -largeArrayDims COMPFLAGS='$COMPFLAGS /Wall'`  
`mex fftri2gridfast.c -v -largeArrayDims COMPFLAGS='$COMPFLAGS /Wall'`

It should be noted that the C99 standard must be used. If your build fails with Microsoft Visual Studio 10, you can try changing the file name to * .cpp, forcing MVSD to use a C ++ compiler.

## Notes on Hardware Acceleration

It should be emphasized that the responsiveness of the plots is highly dependent on the degree of freedom of the PDE problem and the capabilities of the graphics hardware. For larger problems (lots of thousand of vertices), a dedicated graphics card rather than on board graphics should be preferred. Make use of hardware acceleration extensively. Some notes on trouble shooting and tweaking:<br><br>
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

## The License

GPLv3+

Have fun ...

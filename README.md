# How to plot FreeFem++ Simulations in Matlab and Octave

Once you have successfully simulated a PDE problem using FreeFem++ you may want to have a look at your simulation results from within Matlab&copy; or Octave. `ffmatlib` provides some useful commands in order to load FreeFem++ meshes and simulation data and to call the underlying Matlab/Octave plot routines like `contour()`, `quiver()` as well as `patch()`.

![](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/cap_3d_gouraud.png)

## Getting started

  * Click on the button `Clone or download` (see above) and then on the button `Download ZIP`
  * Unzip and change to the directory `demos` and run all FreeFem++ *.epd scripts to create simulation data for plotting
  * Run the matlab `*.m` demo files with Matlab or Octave

Hint: The ffmatlib functions are stored in the folder `ffmatlib`. Use the `addpath(path to ffmatlib)` command if you are working in a different directory.

## Demos

<a name="capacitorexample"></a>

### 2D-Parallel Plate Capacitor

[capacitor_2d.m](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/demos/capacitor_2d.m)  
[capacitor_2d_p1.edp](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/demos/capacitor_2d_p1.edp)  

[Screenshot: 3D Patch](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/cap_3d_gouraud.png)  
[Screenshot: Mesh](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/cap_2d_mesh.png)  
[Screenshot: Contour and Quiver](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/cap_2d_contour.png)  
[Screenshot: 2D Patch with Mesh](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/cap_2d_patch.png)  
[Screenshot: Boundary and Labels](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/cap_2d_labels.png)  

<a name="pdeplotexample"></a>

### 2D-Various ffpdeplot() Examples

[demo_pdeplot.m](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/demos/demo_pdeplot.m)  
[demo_pdeplot_2d_p1.edp](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/demos/demo_pdeplot_2d_p1.edp)  

[Screenshot: 2D Patch](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/l_shape_patch_mesh.png)  
[Screenshot: Contour](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/l_shape_patch_contour.png)  
[Screenshot: Quiver](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/l_shape_patch_quiver.png)  

<a name="convectexample"></a>

### Horizontal Roll Vortices

[convective_rolls.m](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/demos/convective_rolls.m)  
[convective_rolls.edp](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/demos/convective_rolls.edp)  

[Screenshot: Patch + Contour](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/convect_temperature_streamlines.png)  
[Screenshot: Patch + Quiver](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/convect_temperature_velocity.png)  
[Screenshot: 2D Patch + Mesh](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/convect_mesh.png)  
[Screenshot: 2D Streamlines](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/convect_streamlines.png)  

## Function Reference

| Name | Description |
| --- | --- |
| [ffpdeplot()](#ffpdeplotfct) | Creates contour(), quiver() as well as patch() plots from FreeFem++ 2D simulation data |
| [ffreadmesh()](#ffreadmeshfct) | Reads FreeFem++ Mesh Files into Matlab/Octave |
| [ffreaddata()](#ffreaddatafct) | Reads FreeFem++ Data Files into Matlab/Octave |
| [ffplottri2grid()](#ffplottri2gridfct) | Interpolates from 2D triangular mesh to 2D rectangular grid |

<a name="ffpdeplotfct"></a>

## ffpdeplot()

`ffpdeplot()` is a function specially tailored to FreeFem++ that offers most of the features of the classic Matlab `pdeplot()` command. The FEM mesh is entered by vertex coordinates, the boundary values, and the triangles in terms of connectivity as provided by the FreeFem++ `savemesh(Th, "mesh_file.msh")` command. The simulation data can be entered either as point data (native support for P1 simulation data) or as interpolation at the triangle nodes.

#### Syntax

```Matlab
[varargout] = ffpdeplot(p,b,t,varargin)
```

#### Description / Name-Value Pair Arguments

The contents of the points `p`, boundary conditions `b` and triangles `t` arguments are explained in the section [ffreadmesh()](#ffreadmeshfct). `ffpdeplot()` can be called with name-value pair arguments as per following table:

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
| 'Contour' |    Isovalue plot |
|           |       'off' (default) \| 'on' |
| 'CColor' |     Isovalue color |
|           |       [0,0,0] (default) \| 'auto' \| RGB triplet \| 'r' \| 'g' \| 'b' \| |
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

To create a plot read the mesh and the simulation data:

```Matlab
[p,b,t]=ffreadmesh('capacitorp1.msh');
[u]=ffreaddata('capacitor_potential_p1only.txt');
[Ex,Ey]=ffreaddata('capacitor_field_p1only.txt');
```

2D Patch (2D Map or Density) Plot without Boundary:
```Matlab
ffpdeplot(p,[],t,'XYData',u);
```

Plot of the Domain Boundary:
```Matlab
ffpdeplot(p,b,t,'Boundary','on');
```

2D Patch (2D Map or Density) Plot with Boundary:
```Matlab
ffpdeplot(p,b,t,'XYData',u,'Mesh','on','Boundary','on');
```

3D Surf Plot:
```Matlab
ffpdeplot(p,b,t,'XYData',u,'ZStyle','continuous');
```

Contour Plot:
```Matlab
ffpdeplot(p,b,t,'XYData',u,'Contour','on','Boundary','on');
```

Quiver Plot:
```Matlab
ffpdeplot(p,b,t,'FlowData',[Ex, Ey],'Boundary','on');
```

<a name="ffreadmeshfct"></a>

## ffreadmesh()

Reads a FreeFem++ mesh file created by the FreeFem++ `savemesh(Th,"2dmesh.msh")` or `savemesh(Th3d,"3dmesh.mesh")` command to the Matlab/Octave workspace.

#### Syntax

```Matlab
[p,b,t,nv,nbe,nt,labels] = ffreadmesh(filename)
```

#### Description

A mesh file consists of three parts:  

1. the mesh points as nodal coordinates  
2. a list of boundary elements including the boundary labels  
3. list of triangles or tetrahedra defining the mesh in terms of connectivity  

These three blocks are stored in the variables `p`, `b` and `t`.

**2D FreeFem++ Format**

| Parameter | Value |
| --- | --- |
| p | Matrix containing the nodal points |
| b | Matrix containing the boundary edges |
| t | Matrix containing the triangles |
| nv | Number of points/vertices (Th.nv) in the Mesh |
| nt | Number of triangles (Th.nt) in the Mesh |
| nbe | Number of (boundary) edges (Th.nbe) |
| labels | Labels found in the mesh file |

**3D GMSH Format**

| Parameter | Value |
| --- | --- |
| p | Matrix containing the nodal points |
| b | Matrix containing the boundary triangles |
| t | Matrix containing the tetrahedra |
| nv | Number of points/vertices (nbvx, Th.nv) in the Mesh |
| nt | Number of tetrahedra (nbtet, Th.nt) in the Mesh |
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

#### Examples

Reads a mesh file into the Matlab/Octave workspace:

```Matlab
[u]=ffreaddata('capacitor_potential_p1only.txt');
[Ex,Ey]=ffreaddata('capacitor_field_p1only.txt');
```

<a name="ffplottri2gridfct"></a>

## ffplottri2grid()

interpolates the data `tu`, `tv` given on a triangular mesh defined by `tx` and `ty` on a rectangular grid defined by the two vectors `x` and `y`.<br>
To create contour or quiver plots `ffpdeplot()` has its own interpolation routine in the form of a vectorized Matlab/Octave code. However to improve runtime there is a MEX implementation of this code section. If Matlab/Octave finds an executable of `ffplottri2grid.c` within its search path the faster C-implementation is used instead of the internal interpolation routine.

#### Syntax

```Matlab
[u] = ffplottri2grid (x,y,tx,ty,tu)
```

```Matlab
[u,v] = ffplottri2grid (x,y,tx,ty,tu,tv)
```

#### Description

`ffplottri2grid()` uses barycentric interpolation. `tx`, `ty` are 3xnbTriangle matrices containing the triangle vertice coordinates. `tu`, `tv` must be the same size and contain the data at the triangle vertices. The return value `u`, `v` is the interpolation at the mesh grid `mesh(x,y)`. The function returns `NaN's` if an interpolation point is outside the triangle mesh. For more information see also [Notes on MEX Compilation](#notesoncompilation).

#### Examples

```Matlab
[p,b,t]=ffreadmesh('capacitorp1.msh');
u=ffreaddata('capacitor_potential_p1only.txt')';
x=linspace(-5,5,500);
y=linspace(-5,5,500);
xpts=p(1,:);
ypts=p(2,:);
xdata=[xpts(t(1,:)); xpts(t(2,:)); xpts(t(3,:))];
ydata=[ypts(t(1,:)); ypts(t(2,:)); ypts(t(3,:))];
udata=[u(t(1,:)); u(t(2,:)); u(t(3,:))];
U=ffplottri2grid(x,y,xdata,ydata,udata);
[X,Y]=meshgrid(x,y);
surf(X,Y,U,'EdgeColor','none');
view(3);
```

<a name="exportfromff"></a>

## Export FEspace Data from FreeFem++

There are two input data formats supported by `ffpdeplot`:  

1.) FEspace function is given on the mesh nodes  
2.) FEspace function is given at the triangle vertices   

The first method is the prefered one but works only for P1-Element simulations. In order to also be able to plot simulation data created with higher order FE-Elements the data must be converted into P1-Element data or interpolated on the triangle vertices (method 2).

From within the FreeFem++ script write the simulation data with the following statement sequence:

### 2D Problems (P1-Elements)

Saves the Mesh:

```Matlab
savemesh(Th,"capacitorp1.msh");
```
Saves the FEspace function:

```Matlab
ofstream file("capacitor_potential_p1only.txt"); 
for (int j=0; j<u[].n; j++)
   file << u[][j] << endl;
}
```

Saves a 2D vector field:

```Matlab
ofstream file("capacitor_field_p1only.txt");
for (int j=0; j<Ex[].n; j++)
   file << Ex[][j] << " " << Ey[][j] << endl;
}
```

<a name="notesoncompilation"></a>

## Notes on MEX Compilation

Octave:<br>
In Octave seek to the folder `./ffmatlib/` and invoke the command 

`mkoctfile --mex -Wall  ffplottri2grid.c`

Windows:<br>
Under Windows with Microsoft Visual Studio invoke 

`mex  ffplottri2grid.c -v -largeArrayDims COMPFLAGS='$COMPFLAGS /Wall'`

If your build fails with Microsoft Visual Studio 10 ensure to enable C99-standart or you can try to change the file name into *.cpp, forcing MVSD to use a C++ compiler.

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

## Acknowledgments

Special thanks to David Fabre ([StabFEM](https://github.com/erbafdavid/StabFem)) for suggestions and verification of the code

## The License

GPLv3+

Have fun ...

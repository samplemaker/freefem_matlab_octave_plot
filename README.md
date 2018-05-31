# How to plot FreeFem++ simulations in Matlab and Octave

Once you have successfully simulated a PDE problem using FreeFem++ you may want to have a look at your simulation results from within Matlab&copy; or Octave. In this repository you will find some code snippets showing how to make this wish come true.

## Basic theory

A FEM mesh is describing a CAD object or any other type of spatial geometry. It is built from mesh elements. Mesh elements have a shape like a triangle or rectangle, may vary in size and the nodes (vertices) are not necessarily bounded to a rectangular grid. As we look at such an object the surface may be colorized according to the solution of a PDE which is given at the nodal coordinates and parts of the object may be obscured depending on the point of view.<br>
The built in Matlab&copy; or Octave command `patch()` basically renders a set of polygons (=facets, patches). We can associate the FEM mesh elements with such a drawing primitive (here: triangle) and therefore fullfill all the previous stated constraints in order to implement a plot function. In the current implementation the meshdata (the PDE solution at the mesh nodes and the meshing triangles defined by the vertices) have to be written into a text file from within your FreeFem++ script. This file is then parsed and processed by `ffread2patch()` in order to be plot by the `patch()` command. Basically `ffread2patch()` is splitting and rearranging the vertice data because `patch()` expects its input data to be bundled patch wise. A detailed documentation of the `patch()` command can be found at [1](https://de.mathworks.com/help/matlab/ref/patch.html), [2](https://de.mathworks.com/help/matlab/visualize/introduction-to-patch-objects.html) and [3](https://de.mathworks.com/help/matlab/creating_plots/how-patch-data-relates-to-a-colormap.html).

![](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/3dsurf_10.png)

## Source code snippets

### Getting started

Two absolute beginner examples:

  * Run
    * FreeFem++ `demo1_getstarted.edp`
    * From within Matlab/Octave run `demo1_getstarted1.m`
    * From within Matlab/Octave run `demo1_getstarted2.m`

[Screenshot: minimum example](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/2dgetstarted1.png)  
[Screenshot: minimum example](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/2dgetstarted2.png)  

### 2d plot examples

The 2d plot examples focus on displaying functions of the type R<sup>2</sup> &rarr; R or 2d meshes:

  * Run
    * FreeFem++ `demo2_plot2d.edp`
    * From within Matlab/Octave run `demo2_plot2d.m`

[Screenshot: density plot](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/2ddensity.png)  
[Screenshot: surf plot](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/2dsurf.png)  
[Screenshot: 2d-mesh](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/2dmesh.png)  

### 3d plot examples

The 3d plot examples focus on displaying functions of the type R<sup>3</sup> &rarr; R (i.e. a 3d object boundary colored with a scalar value like a temperature) or 3d mesh surfaces:

  * Run
    * FreeFem++ `demo3_plot3d_cyl.edp`
    * FreeFem++ `demo3_plot3d_box.edp`
    * From within Matlab/Octave run `demo3_plot3dbd.m`

[Screenshot: surf plot](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/3dsurf_2.png)  
[Screenshot: surface of a 3d-mesh](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/3dmesh.png)

### Advanced: 3d slicing examples

To make the inside visible it is also possible to cut a 3d FreeFem++ simulation along a slicing plane. You have to write the mesh elements as well as the boundary information to use this feature:

  * Run
    * From within Matlab/Octave run `demo4_start_slicer_gui.m`
    * A minimum example: `demo4_slice3d.m`

[Screenshot: slice3d](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/3dsurf_slice7.png)  
[Screenshot: boundary](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/3dsurf_slice8.png)  
[Screenshot: crosssection](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/3dsurf_slice9.png)  
[Screenshot: GUI-Slicer](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/GUI_Slicer.png)  

### Advanced: 2d isovalues

For 2d problems it is sometimes helpful to have an isovalue - plot which will be created in the following example:<br>

  * Run
    * FreeFem++ `demo5_isovalues_2d.edp`
    * From within Matlab/Octave run `demo5_isovalues_2d.edp`

[Screenshot: isovalues](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/2dsurf_isovalues.png)  
[Screenshot: surf](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/2dsurf_surf.png)  

### Advanced: 2d vectorfields

For 2d problems it is sometimes helpful to have an vectorfield - plot as well. Such a plot  will be created in the following example:<br>

  * Run
    * FreeFem++ `demo6_vector_2d.edp`
    * From within Matlab/Octave run `demo6_vector_2d.m`

[Screenshot: vectorfields](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/2dvectorfield.png)  

## Implementation - FreeFem++

### Simple 2d problems

From within the FreeFem++ script write the mesh elements (triangles defined by its vertices) and the solution of the PDE at the nodes respectively:

```cpp
ofstream tridata ("export_tri.txt");
for (int i=0; i<Th.nt; i++){
  for (int j=0; j<3; j++){
    tridata << Th[i][j].x << ";"
            << Th[i][j].y << ";"
            << u[][Vh(i,j)] << "\n";
  }
}
```

### 3d boundary problems

If the domain boundary (surface) is to be displayed it is enough to write the boundary elements only:

```cpp
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

### Advanced: 3d problems (slicing)
If a crosssection is to be made it is necessary to write the mesh elements (tetrahedra) as well as the boundary data:

```cpp
ofstream tetdata ("export_tet.txt");
for (int i=0; i<Th3d.nt; i++){
  for (int j=0; j<4; j++){
    tetdata << Th3d[i][j].x << ";"
            << Th3d[i][j].y << ";"
            << Th3d[i][j].z << ";"
            << u[][Vh(i,j)] << "\n";
  }
}
```

### Advanced: 2d problems (vector fields)

In order to plot 2d vector fields just write multiple columns:

```cpp
ofstream file ("export_tri_ncols.txt");
for (int i = 0; i < Th.nt; i++){
  for (int j = 0; j < 3; j++){
    file << Th[i][j].x << ";"
         << Th[i][j].y << ";"
         << uh[][Vh(i,j)] << ";"
         << vh[][Vh(i,j)] << "\n";
  }
}
```

## Implementation - Matlab&copy;/Octave

Hint: Library functions can be found in the folder `ffmatlib`. Therefore use the `addpath(ffmatlib)` command if you are working from another directory.

### General purpose plots in 2d and 3d boundary

The Matlab&copy;/Octave function `ffread2patch()` reads and rearranges file content in order to be plot with the `patch()` command. It's arguments depend on the number of columns and on the separation character. `ffread2patch()` can process 2d mesh elements or 3d boundary data (must be triangle):

```cpp
[X,Y,C, ...] = ffread2patch('filename.txt','Delimiter',';','Format','auto');
```
The number of columns can be set via the format specifier explicitely:

```cpp
[X,Y,Z,C] = ffread2patch('filename.txt','Delimiter',';','Format','%f %f %f %f');
```
X, Y and C are matrices which can be fed to `patch()`.

Hint: You can split the reading and conversion process into two different entities:

```cpp
[tridata] = ffreadfile('File1','filename.txt','Delimiter',';','Format','%f %f %f');
[X,Y,C] = fftri2patch(tridata);
```

### Advanced: 3d slicing

In order to perform a slice from a 3d problem you have to load the boundary data as well as the mesh element data. This can be accomplished with the function call:

```cpp
[file1data,file2data] = ffreadfile ('File1','filename1.txt','File2','filename2.txt', ...
                                    'Delimiter',';','Format','%f %f %f %f')');
```

Next perform the slice with the data stored in the variables file1data, file2data: 

```cpp
[BX,BY,BZ,BC] = slicebd2patch (boundarydata,S1,S2,S3);
[SX,SY,SZ,SC] = slicetet2patch (meshdata,S1,S2,S3);
```

The input arguments `S1..S3` contain x,y,z coordinates of three points defining the slicing plane. The output data `(BX,BY,BZ,BC)` and `(SX,SY,SZ,SC)` or the superposition of both `([SX BX],[SY BY],[SZ BZ],[SC BC])` can be plot with the `patch()` command.

### Advanced: Slicer_GUI

There is a graphical user interface (working in both worlds, Matlab and Octave either) which can be used to create crossections of PDE problems more easily. The function call is

```cpp
slicer_gui(bddatafile,tetdatafile);
```
You may also have a look at `demo4_start_slicer_gui.m`.

### Advanced: 2d isovalues

To create isolevel-curveplots with the Matlab / Octave command `contour()` it is necessary to interpolate the data given on the mesh vertices on a rectangular grid. In a first step load the raw data:

```cpp
[tridata] = ffreadfile('File1','filename.txt','Delimiter',';','Format','%f %f %f');
```

Then interpolate the third column in `tridata` on a rectangular grid defined by X and Y. If grid points are outside from any triangles the function will return NaN's.

```cpp
C = fftri2grid(tridata,X,Y);
```

The result can be plot with the command:

```cpp
[c,h] = contour(X,Y,C,8);
```

### Advanced: 2d vectorfields

To create vector field plots with the Matlab / Octave command `quiver()` it is necessary to interpolate the data given on the mesh vertices on a rectangular grid. In a first step load the raw data:

```cpp
[tridata] = ffreadfile('File1','filename.txt','Delimiter',';','Format','%f %f %f');
```

Then interpolate the data columns in `tridata` on a rectangular grid defined by X and Y. If grid points are outside from any triangles the function will return NaN's.

```cpp
[U,V] = fftri2grid(tridata,X,Y);
```

The result can be plot with the command:

```cpp
quiver(X,Y,U,V);
```

## Files

  * `ffread2patch.m` Read FreeFem++ simulation results and convert vertex/triangle data to patch plot data.
  * `ffreadfile.m` Read one or two FreeFem++ simulation result files.
  * `fftri2patch.m` Convert FreeFem++ vertex/triangle data to patch plot data.
  * `slicetet2patch.m` Slice 3d mesh elements (tetrahedra) and convert to patch plot data.
  * `slicebd2patch.m` Slice 3d boundary (triangle) data and convert to patch plot data.
  * `demo4_start_slicer_gui.m` Starts the slicer graphical user interface (Slicer_GUI) to slice 3d data.
  * `slicer_gui.m` Slicer_GUI - Implementation.
  * `fftri2grid.m` Interpolates R<sup>2</sup> &rarr; R<sup>n</sup> FEM data on a rectangular mesh grid.

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

## Hardware acceleration

It should be emphasized that responsivity and smotheness of the plots will strongly depend on the degree of freedom of the PDE problem and on the capabilities of your graphics hardware. For larger problems (lots of thousand of vertices) you should prefer a dedicated graphics card instead on board graphics. Make use of hardware acceleration extensively. Some notes on trouble shooting and tweaking:<br><br>
 If `get(gcf,'RendererMode')` is set to auto Matlab/Octave will decide on its own which renderer is the best for the current graphic task.

  * `get(figure_handle,'Renderer')` returns the current figure() renderer
  * `set(figure_handle,'Renderer','OpenGL')` forces a figure() to switch to OpenGL
  * `set(figure_handle,'Renderer','painters')` forces a figure() to switch to vector graphics

Generally OpenGL can be considered to be faster than painters. To get an OpenGL info you can type `opengl info` within Matlab. Ensure the line `Software` shows `false` otherwise you are running in Software OpenGL. If Hardware-accelerated OpenGL is available on the system you may alter the mode manually with the `opengl software` and `opengl hardware` commands.

## The License

GPLv3+

Have fun ...

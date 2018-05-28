# How to plot FreeFem++ simulations in Matlab and Octave

Once you have successfully simulated a PDE problem using FreeFem++ you may want to have a look at your results from within Matlab&copy; or Octave. In this repository you will find some code snippets showing how to make this wish come true.

## Basic theory

The widely used Matlab&copy;/Octave commands `surf()` and `mesh()` plot measurement data or functions based on a rectangular grid. The command `patch()` plots polygons (=facets, patches) based on vertice coordinates which are allowed to be independent from a rectangular grid. We will associate those drawing primitives with FE-mesh elements (here: triangles) and hence enable plotting of irregular tesselation structures like FreeFem++ meshes.<br>To do this the meshdata (the PDE solution at the mesh nodes and the meshing triangles defined by the nodal coordinates) have to be written into a text file from within your FreeFem++ program. This file is then parsed and processed by `ffread2patch()` in order to be plot by the `patch()` command. Basically `ffread2patch()` is splitting and rearranging the continuous vertice data because `patch()` expects its drawing coordinates bundled patch wise. A detailed documentation of the `patch()` command can be found here: [1](https://de.mathworks.com/help/matlab/ref/patch.html), [2](https://de.mathworks.com/help/matlab/visualize/introduction-to-patch-objects.html) and [3](https://de.mathworks.com/help/matlab/creating_plots/how-patch-data-relates-to-a-colormap.html).

![](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/3dsurf_10.png)

## Getting started: A minimum example

Two simple 2d examples:

  * Run
    * FreeFem++ `demo1_getstarted.edp`
    * From within Matlab/Octave run `demo1_getstarted1.m`
    * From within Matlab/Octave run `demo1_getstarted2.m`

[Screenshot: minimum example](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/2dgetstarted1.png)  
[Screenshot: minimum example](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/2dgetstarted2.png)  

## Running the 2d plot examples

The 2d plot examples focus on displaying functions of the type R<sup>2</sup> &rarr; R or 2d meshes:

  * Run
    * FreeFem++ `demo2_plot2d.edp`
    * From within Matlab/Octave run `demo2_plot2d.m`

[Screenshot: density plot](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/2ddensity.png)  
[Screenshot: surf plot](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/2dsurf.png)  
[Screenshot: 2d-mesh](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/2dmesh.png)  

## Running the 3d plot examples

The 3d plot examples focus on displaying functions of the type R<sup>3</sup> &rarr; R (i.e. a 3d object boundary colored with a scalar value like a temperature) or 3d mesh surfaces:

  * Run
    * FreeFem++ `demo3_plot3d_cyl.edp`
    * FreeFem++ `demo3_plot3d_box.edp`
    * From within Matlab/Octave run `demo3_plot3dbd.m`

[Screenshot: surf plot](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/3dsurf_2.png)  
[Screenshot: surface of a 3d-mesh](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/3dmesh.png)

## Running the 3d slicing examples

To make the inside visible it is also possible to cut a 3d FreeFem++ simulation along a slicing plane. You have to write the mesh elements as well as the boundary information to use this feature.

  * Run
    * FreeFem++ `demo4_slice3d.edp`
    * From within Matlab/Octave run `slicer_gui.m`
    * A minimum example: `demo4_slice3d.m`

[Screenshot: slice3d](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/3dsurf_slice7.png)  
[Screenshot: boundary](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/3dsurf_slice8.png)  
[Screenshot: crosssection](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/3dsurf_slice9.png)  

## Implementation

**FreeFem++ 2d:** From within the FreeFem++ script write the mesh elements (triangles defined by its vertices) and the solution of the PDE at the nodes respectively:

```cpp
ofstream datamesh ("export_tri.txt");
for (int i=0; i<Th.nt; i++){
  for (int j=0; j<3; j++){
    datamesh << Th[i][j].x << ";"
             << Th[i][j].y << ";"
             << u[][Vh(i,j)] << "\n";
  }
}
```

**FreeFem++ 3d - boundaries:** If the domain boundary (surface) is to be displayed it is enough to write the boundary elements only:

```cpp
int idx;
int nbelement=Th3d.nbe;
ofstream bedata ("boundary_file.txt");
for (int k=0;k<nbelement;++k){
  for (int num=0;num<3;num++){
    idx=Th3d.be(k)[num];
    bedata << Th3d(idx).x << ";"
           << Th3d(idx).y << ";"
           << Th3d(idx).z << ";"
           <<  u(Th3d(idx).x,Th3d(idx).y,Th3d(idx).z) << "\n";
  }
}
```

**FreeFem++ 3d - slice:** If a crosssection is to be made it is necessary to write the mesh elements (tetrahedra) as well as the boundary data:

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

**Matlab&copy;/Octave - 2d and 3d boundary:**

The Matlab&copy;/Octave function `ffread2patch()` reads and rearranges the prior written file content in order to be plot using the `patch()` command. It's arguments depend on the number of columns and on the separation character. `ffread2patch()` can process both 2d mesh elements or 3d boundary data (triangle) and as many columns as you want:

```cpp
[X,Y,C, ...] = ffread2patch('filename.txt','Delimiter',';','Format','auto');
```
If you don't like to autodetect the number of columns you can give the format specifier explicitely:

```cpp
[X,Y,Z,C] = ffread2patch('filename.txt','Delimiter',';','Format','%f %f %f %f');
```
XX, YY and CC are matrices which can be fed to `patch()`.

Hint: You can split the reading and rearranging process into two different entities with the two commands:

```cpp
[M] = ffreadfile('File1','filename.txt','Delimiter',';','Format','%f %f %f');
[X,Y,C] = tri2patch(M);
```

**Matlab&copy;/Octave - 3d slice:**

The Matlab&copy;/Octave function `ffreadfile()` loads the boundary and / or mesh element data into a temporary variable, which is then processed with the slicing package `slicebd2patch()` (converts the boundary data), `slicetet2patch()` (creates crosssection data). The output of the two latter functions can be plot invoking the `patch()` command.

## Files

  * `ffread2patch.m` Read FreeFem++ simulation results and convert vertex/triangle data to patch plot data.
  * `ffreadfile.m` Read one or two FreeFem++ simulation result files.
  * `tri2patch.m` Convert FreeFem++ vertex/triangle data to patch plot data.
  * `slicetet2patch.m` Slice 3d mesh elements (tetrahedra) and convert to patch plot data.
  * `slicebd2patch.m` Slice 3d boundary (triangle) data and convert to patch plot data.

## Software

  * [FreeFem++][freefem]
  * [Octave][octave]

[freefem]:    http://www.freefem.org//
             "FreeFem++ solver for partial differential equations"
[octave]:     https://www.gnu.org/software/octave/
             "GNU Octave scientific programming language"

[matlab]:     https://de.mathworks.com/products/matlab.html
             "Matlab scientific programming language"

## Hardware acceleration

For larger degree of freedom problems speed may become a concern. If `get(gcf,'RendererMode')` is set to auto Matlab/Octave will decide on its own which renderer is the best for the current graphic task.

  * `get(figure_handle,'Renderer')` returns the current figure() renderer
  * `set(figure_handle,'Renderer','OpenGL')` forces a figure() to switch to OpenGL
  * `set(figure_handle,'Renderer','painters')` forces a figure() to switch to vector graphics

Older Matlab releases had a zbuffer-renderer (raster graphics) as well. Generally OpenGL can be considered to be faster than painters. To check if OpenGL is available on the system you can type `opengl info` within Matlab. Ensure the line `Software` shows `false` otherwise you are running in Software OpenGL. If Hardware-accelerated OpenGL is available on the system you may alter the mode manually with the `opengl software` and `opengl hardware` commands. 

## The License

GPLv3+


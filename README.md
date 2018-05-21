# How to plot FreeFem++ solutions in Matlab and Octave

This is a collection of minimal examples to show how FreeFem++ simulation results can be plotted in Matlab/Octave. It should be emphasized that it is not necessary to have the Matlab-PDEtools installed to run the examples.

## Basic theory

In contrast to the Matlab/Octave functions `surf()` and `mesh()` which do work on rectangular grids the function `patch()` basically plots polygons (=facets, patches) and hence enables plotting of irregular tesselation structures like FreeFem++ meshes. To do this the meshdata (triangles defined by vertex data) and the solution has to be written into a text file from within the FreeFem++ script. This file is then parsed and processed by `ff2patch()` in order to be plot with the `patch()` command. `ff2patch()` is doing nothing but splitting and rearranging the continuous vertice data into batches of three adjacent numbers because `patch()` expects its input bundled per triangle. Detailed documentation of the `patch()` command can be found here: [1](https://de.mathworks.com/help/matlab/ref/patch.html), [2](https://de.mathworks.com/help/matlab/visualize/introduction-to-patch-objects.html) and [3](https://de.mathworks.com/help/matlab/creating_plots/how-patch-data-relates-to-a-colormap.html).

![](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/3dmesh.png)

## Running the 2d plot examples

The 2d plot examples focus on displaying functions of the type R<sup>2</sup> &rarr; R or 2d meshes:

  * Seek into the folder `plot2d` and run
    * FreeFem++ ffgendata2ddisc.edp
    * From within Matlab/Octave run runplotdemo2d.m

[Screenshot: density plot](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/2ddensity.png)  
[Screenshot: surf plot](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/2dsurf.png)  
[Screenshot: 2d-mesh](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/2dmesh.png)  

## Running the 3d plot examples

The 3d plot examples focus on displaying functions of the type R<sup>3</sup> &rarr; R (i.e. a 3d object boundary colored with a scalar value like a temperature) or 3d mesh surfaces:

  * Seek into the folder `plot3d` and run
    * FreeFem++ ffgendata3dcyl.edp
    * FreeFem++ ffgendata3dbox.edp
    * From within Matlab/Octave run runplotdemo3d.m

[Screenshot: surf plot](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/3dsurf_2.png)  
[Screenshot: surface of a 3d-mesh](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/3dmesh.png)

The Matlab/Octave functions `ffslicebd3d.m` (cuts the boundary data) and `ffslicetet3d.m` (gives the crosssection data) cut a 3d FreeFem++ simulation along a plane and make its inside visble. You will need to write the complete mesh data as well as the boundary data from within FreeFem++ to use this function. It already does it's job but it is still experimental. That means i may change the calling syntax in future if i find it necessary. 

  * Seek into the folder `slice3d` and run
    * FreeFem++ ffgendata3d.edp
    * From within Matlab/Octave run  rundemoslice3d.m

[Screenshot: slice3d](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/3dsurf_slice3.png)  
[Screenshot: zoom](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/3dsurf_slice4.png)  

## Implementation

**2d:** Write the complete mesh data from within FreeFem++:

```cpp
ofstream datamesh ("tridata2d.txt");
for (int i=0; i<Th.nt; i++){
  for (int j=0; j<3; j++){
    datamesh << Th[i][j].x << ";"
             << Th[i][j].y << ";"
             << u[][Vh(i,j)] << "\n";
  }
}
```

**3d:** If the domain boundary is to be displayed it is sufficient to write the boundary elements only. If a full 3d sclice is to be made it is necessary to write the mesh data as well:

```cpp
int idx;
int nbelement=Th3d.nbe;
ofstream bedata ("nbtridata3d.txt");
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

The Matlab/Octave function `ff2patch()` rearranges the prior written file content in order to be plot using the `patch()` command. The arguments depend on the number of columns and on the separation character. `ff2patch()` does not care about content hence you can process complete net data as well as boundary data and as many columns as you want:

```cpp
[XX,YY,CC] = ff2patch('filename.txt','Delimiter',';','Format','auto');
```
If you don't like to autodetect the number of columns you can use the format specifier:

```cpp
[XX,YY,ZZ,CC] = ff2patch('filename.txt','Delimiter',';','Format','%f %f %f %f');
```
XX, YY and CC are matrices which can be fed to `patch()`.

## Files

  * `ff2patch.m` Library function which reads and rearranges the FreeFem++ output data in order to be plot within Matlab/Octave
  * `ffslice3d.m` Cuts 3d simulation data along a plane
  * `runplotdemo2d.m` Matlab/Octave file demonstrating 2d surf, 2d density and 2d meshplots
  * `runplotdemo3d.m` Matlab/Octave file demonstrating 3d surf including tesselation and 3d meshplots
  * `ffgendata2ddisc.edp` Creates some fantasy art data (diffusion of temperature in a 2d-sheet metal)
  * `ffgendata3dcyl.edp` Creates a cylindrical 3d mesh
  * `ffgendata3dbox.edp` Creates some fantasy art (poisson problem - spatial 3d temperature diffusion by heat conduction)

## Software and system requirements

  * [FreeFem++][freefem]
  * [Octave][octave]

[freefem]:    http://www.freefem.org//
             "FreeFem++ solver for partial differential equations"
[octave]:     https://www.gnu.org/software/octave/
             "GNU Octave scientific programming language"

## The License

GPLv3+


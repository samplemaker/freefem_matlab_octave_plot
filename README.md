# How to plot FreeFem++ solutions in Matlab and Octave

This is a collection of minimal examples to show how FreeFem++ simulation results can be plotted in Matlab/Octave. It should be emphasized that it is not necessary to have the Matlab-PDEtools installed to run the examples.

## Basic theory

In contrast to the Matlab/Octave functions `surf()` and `mesh()` which do work on cartesian meshes the function `patch()` basically plots polygons (=facets) and hence enables plotting of irregular tesselation structures like FreeFem++ meshes. To do this the meshdata (triangle/vertex data) and the solution has to be written prior into a text file from within the FreeFem++ script. This file is then parsed and processed by `ff2patch()` in order to be plot with the `patch()` command. `ff2patch()` is doing nothing else but splitting and rearranging the continuous vertice data into batches of three adjacent numbers because `patch()` expects its input bundled per triangle (=facet). A documentation of the `patch()` command can be found [here](https://de.mathworks.com/help/matlab/ref/patch.html).

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

## Implementation

**2d:** In FreeFem++ we have to write all data:

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

**3d:** To display the surface it is sufficient to write the boundary elements only. If a 3d sclice is to be done it is also necessary to write the complete net data.

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

The Matlab/Octave function `ff2patch()` rearranges the prior written file content in order to be plot using the `patch()` command. The arguments depend on the number of columns and on the separation character.  `ff2patch()` does not care about content hence you can process as many columns as you want:

```cpp
[XX,YY,CC] = ff2patch('filename.txt','Delimiter',';','Format','auto');
```
If you don't like to autodetect the number of columns use

```cpp
[XX,YY,CC] = ff2patch('filename.txt','Delimiter',';','Format','%f %f %f');
```

The Matlab/Octave function `ffslice3d()` is a proof of concept to cut 3d FreeFem++ simulation data along a plane. You will need to write boundary and tetrahedron data as well from within FreeFem++ to use this function.

[Screenshot: slice3d](https://raw.githubusercontent.com/samplemaker/freefem_matlab_octave_plot/public/screenshots/3dsurf_slice3.png)  

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


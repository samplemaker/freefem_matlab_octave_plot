# How to plot FreeFem++ solutions in Matlab and Octave

Examples demonstrating how to write FreeFem++ simulation results and how to create plots in Matlab/Octave from the results respectively. It should be emphasized that it is NOT necessary to have the Matlab-PDEtools installed to run these examples!

## Basic theory

In contrast to the Matlab/Octave functions `surf()` and `mesh()` which do work on cartesian meshes the built in function `patch()` basically plots triangles (polygons) and hence enables plotting of irregular tesselation structures like FreeFem++ meshes. To do this the meshdata (triangle/vertex data) and the solution has to be written prior into a text file from within the FreeFem++ script. This file is then parsed and processed by the `ff2patch()` function in order to be plot with the `patch()` command. `ff2patch()` is doing nothing else but splitting and rearranging the continuous vertice data stream into batches of three adjacent numbers each because `patch()` expects its input to be bundled triangle (=facet) wise. Documentation of `patch()` can be found here [patch objects](https://de.mathworks.com/help/matlab/ref/patch.html).

## Running the 2d plot examples

The 2d plot examples focus on displaying functions R<sup>2</sup> &rarr; R or 2d meshes:

  * Seek into folder **plot2d** and run
    * FreeFem++ ffgendata2ddisc.edp
    * From within Matlab/Octave run runplotdemo2d.m

[Screenshot: density plot](https://github.com/samplemaker/freefem_matlab_octave_plot/blob/master/screenshots/2ddensity.png)  
[Screenshot: surf plot](https://github.com/samplemaker/freefem_matlab_octave_plot/blob/master/screenshots/2dsurf.png)  
[Screenshot: 2d-mesh](https://github.com/samplemaker/freefem_matlab_octave_plot/blob/master/screenshots/2dmesh.png)  

## Running the 3d plot examples

The 3d plot examples focus on displaying functions R<sup>3</sup> &rarr; R (i.e. a 3d object boundary colored with a scalar value like temperature etc.) or 3d mesh surfaces:

  * Seek into folder **plot3d** and run
    * FreeFem++ ffgendata3dcyl.edp
    * FreeFem++ ffgendata3dbox.edp
    * From within Matlab/Octave run runplotdemo3d.m

[Screenshot: surf plot](https://github.com/samplemaker/freefem_matlab_octave_plot/blob/master/screenshots/3dsurf.png)  
[Screenshot: surface of a 3d-mesh](https://github.com/samplemaker/freefem_matlab_octave_plot/blob/master/screenshots/3dmesh.png)

## Files

  * `ff2patch.m` Library function which reads and rearranges the FreeFem++ output data in order to be plot within Matlab/Octave
  * `runplotdemo2d.m` Matlab/Octave file demonstrating 2d surf, 2d density and 2d meshplots
  * `runplotdemo3d.m` Matlab/Octave file demonstrating 3d surf including tesselation and 3d meshplots
  * `ffgendata2ddisc.edp` Creates some fantasy art data (diffusion of temperature in a 2d-sheet metal)
  * `ffgendata3dcyl.edp` Creates a cylindrical 3d mesh
  * `ffgendata3dbox.edp` Creates some fantasy art (poisson problem - spatial 3d temperature diffusion by heat conduction)

## Implementation

2d: In FreeFem++ we can simply write the complete mesh:

```javascript
ofstream datamesh ("tridata2ddisc.txt");
for (int i=0; i<Th.nt; i++){
  for (int j=0; j<3; j++){
    datamesh << Th[i][j].x << ";"
             << Th[i][j].y << ";"
             << u[][Vh(i,j)] << "\n";
  }
}
```

3d: In FreeFem++ boundary elements must be written in order to display the surface only

```javascript
int idx;
int nbelement=Th33d.nbe;
ofstream bedata ("tridata3dbox.txt");
for (int k=0;k<nbelement;++k){
  for (int num=0;num<3;num++){
    idx = Th33d.be(k)[num];
    bedata << Th33d(idx).x << ";"
           << Th33d(idx).y << ";"
           << Th33d(idx).z << ";"
           << u[][idx] << "\n";
  }
}
```

The Matlab/Octave function `ff2patch()` call depends on the number of columns and on the separation character. It does not care about the content, so you can process as many scalars or coordinates as you want:

```javascript
[XX,YY,CC] = ff2patch('filename.txt','Delimiter',';','Format','auto');
```
If you don't like to autodetect the number of columns use

```javascript
[XX,YY,CC] = ff2patch('filename.txt','Delimiter',';','Format','%f %f %f');
```

## Hacking

The folder `experimental` contains unfinished examples (meshes) and tests (3d slicing function)

## Software and system requirements

  * [FreeFem++][freefem]
  * [Octave][octave]

[freefem]:    http://www.freefem.org//
             "FreeFem++ solver for partial differential equations"
[octave]:     https://www.gnu.org/software/octave/
             "GNU Octave scientific programming language"

## The License

GPLv3+

 

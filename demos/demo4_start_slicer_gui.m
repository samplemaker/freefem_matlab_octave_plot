clear all;
addpath('ffmatlib');

%Boundary data file
bddatafile='temp_demo4_bddata3d_box.txt';
%File containing the mesh elements
tetdatafile='temp_demo4_tetdata3d_box.txt';

slicer_gui(bddatafile,tetdatafile);

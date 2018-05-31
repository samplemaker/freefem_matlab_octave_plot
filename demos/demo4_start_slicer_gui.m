clear all;
addpath('ffmatlib');

%boundary data file
bddatafile='temp_demo4_bddata3d_box.txt';
%file containing the mesh elements
tetdatafile='temp_demo4_tetdata3d_box.txt';

slicer_gui(bddatafile,tetdatafile);

clear all;
fprintf('Running FreeFem++ scripts ...\n');
tic;
fprintf('demo_pdeplot.edp\n');
system('FreeFem++ -ns -nw -v 0 demo_pdeplot.edp');
fprintf('demo1_getstarted.edp\n');
system('FreeFem++ -ns -nw -v 0 demo1_getstarted.edp');
fprintf('demo2_plot2d.edp\n');
system('FreeFem++ -ns -nw -v 0 demo2_plot2d.edp');
fprintf('demo3_plot3d_box.edp\n');
system('FreeFem++ -ns -nw -v 0 demo3_plot3d_box.edp');
fprintf('demo3_plot3d_cyl.edp\n');
system('FreeFem++ -ns -nw -v 0 demo3_plot3d_cyl.edp');
fprintf('demo4_slice3d.edp\n');
system('FreeFem++ -ns -nw -v 0 demo4_slice3d.edp');
fprintf('demo5_isovalues_2d.edp\n');
system('FreeFem++ -ns -nw -v 0 demo5_isovalues_2d.edp');
fprintf('demo6_vector_2d.edp\n');
system('FreeFem++ -ns -nw -v 0 demo6_vector_2d.edp');
fprintf('demo7_slice3d.edp\n');
system('FreeFem++ -ns -nw -v 0 demo7_slice3d.edp');
fprintf('demo8_slice3d_vectors.edp\n');
system('FreeFem++ -ns -nw -v 0 demo8_slice3d_vectors.edp');
toc;
fprintf('Starting demo ...\n');
demo_pdeplot
pause(5);
close all;
demo1_getstarted1
pause(5);
close all;
demo1_getstarted2
pause(5);
close all;
demo2_plot2d
pause(5);
close all;
demo3_plot3dbd
pause(5);
close all;
demo4_slice3d
pause(5);
close all;
demo4_start_slicer_gui
pause(5);
close all;
demo5_isovalues_2d
pause(5);
close all;
demo6_vector_2d
pause(5);
close all;
demo7_slice3d_2dgrid
pause(5);
close all;
demo8_slice3d_2dgrid_vectors
pause(5);
close all;

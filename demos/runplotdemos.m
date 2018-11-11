clear all;
fprintf('Running FreeFem++ scripts ...\n');
tic;
fprintf('capacitor_2d_p1.edp\n');
system('FreeFem++ -ns -nw -v 0 capacitor_2d_p1.edp');
fprintf('demo_pdeplot_2d_p1.edp\n');
system('FreeFem++ -ns -nw -v 0 demo_pdeplot_2d_p1.edp');
fprintf('convective_rolls.edp\n');
system('FreeFem++ -ns -nw -v 0 convective_rolls.edp');
fprintf('periodic_bc.edp\n');
system('FreeFem++ -ns -nw -v 0 periodic_bc.edp');
fprintf('complex_pde_v2.edp\n');
system('FreeFem++ -ns -nw -v 0 complex_pde_v2.edp');
toc;
fprintf('Starting demo ...\n');
capacitor_2d
pause(5);
close all;
demo_pdeplot
pause(5);
close all;
convective_rolls
pause(5);
close all;
testsubplots
pause(5);
close all;
periodic_bc
pause(5);
close all;
bdcoloring
pause(5);
close all;
complex_pde_v2
pause(5);
close all;
complex_pde_v3
pause(5);
close all;
fprintf('Running FreeFem++ scripts (3D) ...\n');
fprintf('capacitor_3d.edp\n');
system('FreeFem++ -ns -nw -v 0 capacitor_3d.edp');
capacitor_3d
pause(5);
close all;


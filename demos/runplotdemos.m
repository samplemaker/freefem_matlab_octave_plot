clear all;
fprintf('Running FreeFem++ scripts ...\n');
tic;
fprintf('capacitor_2d_p1.edp\n');
system('FreeFem++ -ns -nw -v 0 capacitor_2d_p1.edp');
fprintf('demo_pdeplot_2d_p1.edp\n');
system('FreeFem++ -ns -nw -v 0 demo_pdeplot_2d_p1.edp');
fprintf('convective_rolls.edp\n');
system('FreeFem++ -ns -nw -v 0 convective_rolls.edp');
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


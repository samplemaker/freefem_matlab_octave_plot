clear all;
fprintf('capacitor_2d.edp\n');
system('FreeFem++ -ns -nw -v 0 capacitor_2d.edp');
fprintf('interpolate_complex.edp\n');
system('FreeFem++ -ns -nw -v 0 interpolate_complex.edp');
fprintf('capacitor_3d.edp\n');
system('FreeFem++ -ns -nw -v 0 capacitor_3d.edp');
capacitor_2d
pause(5);
close all;
interpolate
pause(5);
close all;
capacitor_3d
pause(5);
close all;

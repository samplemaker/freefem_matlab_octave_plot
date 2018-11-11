clear all;
addpath('ffmatlib');


%Compute a conformal plot based on a PDE simulation result
[p,b,t]=ffreadmesh('complex_test.msh');
[u]=ffreaddata('complex_test.txt');
figure;
hold on;
%Z1: real part = const
%Z2: imag part = const
[Z1, Z2] = ffcplxmesh([0,0], [2*pi(),2*pi()], [3,3], [12,12]);
[w] = ffinterpolate(p,b,t,real(Z1),imag(Z1),u);
plot(real(w),imag(w),'b','LineWidth',2);
[w] = ffinterpolate(p,b,t,real(Z2),imag(Z2),u);
plot(real(w),imag(w),'r','LineWidth',2);
grid;
title('Conformal plot of a PDE solution');
daspect([1,1,1]);


%Project a polar plot onto a PDE solution
[p,b,t]=ffreadmesh('capacitorp1.msh');
[u]=ffreaddata('capacitor_potential_p1only.txt');
[Z1, Z2] = ffcplxmesh([0.5,0], [4,2*pi()], [10,10], [10,11]);
%Map to polar coordinates
strFunc='@(Z)(real(Z).*exp(1i*imag(Z)))';
f = str2func(strFunc);
U1 = f(Z1);
U2 = f(Z2);
figure;
hold on;
%constant radius
[w] = ffinterpolate(p,b,t,real(U1),imag(U1),u);
plot3(real(U1),imag(U1),w,'g','LineWidth',2);
%constant angle
[w] = ffinterpolate(p,b,t,real(U2),imag(U2),u);
plot3(real(U2),imag(U2),w,'g','LineWidth',2);
ffpdeplot(p,b,t, ...
          'XYData',u, ...
          'ZStyle','continuous', ...
          'Mesh','off', ...
          'CBTitle','U[V]', ...
          'Title','Curved Interpolation');
ylabel('y');
xlabel('x');
zlabel('u');
grid;
lighting gouraud;
view([-47,24]);
camlight('headlight');


%Interpolates on a curved line
N = 100;
s = linspace(0,2*pi(),N);
Z = 3.5*(cos(s)+1i*sin(s)).*sin(0.5*s);
w = ffinterpolate(p,b,t,real(Z),imag(Z),u);
figure('position', [0, 0, 800, 300])
subplot(1,2,1);
hold on;
plot3(real(Z),imag(Z),real(w),'g','LineWidth',2);
ffpdeplot(p,b,t, ...
          'XYData',u, ...
          'ZStyle','continuous', ...
          'Mesh','off', ...
          'ColorBar','off', ...
          'Title','Single curve Interpolation');
ylabel('y');
xlabel('x');
zlabel('u');
grid;
lighting gouraud;
view([-47,24]);
camlight('headlight');
subplot(1,2,2);
plot(s,real(W),'b');
grid;
xlabel('s');
ylabel('u');
title('Interpolation Values');
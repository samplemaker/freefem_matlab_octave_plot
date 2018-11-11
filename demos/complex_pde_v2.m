clear all;
addpath('ffmatlib');


%Compute a conformal plot based on a PDE simulation result
[p,b,t]=ffreadmesh('complex_test.msh');
[uC]=ffreaddata('complex_test.txt');
figure;
hold on;
%ZX: real part = const
%ZY: imag part = const
[ZX, ZY] = ffcplxmesh([0,0], [2*pi(),2*pi()], [3,3], [12,12]);
[WX, WY] = ffipocplx(p,b,t,uC,ZX,ZY);
plot(real(WX),imag(WX),'b','LineWidth',1);
plot(real(WY),imag(WY),'r','LineWidth',1);
grid;
title('Conformal plot of a PDE solution');
daspect([1,1,1]);


%Project a polar plot onto a PDE solution
[p,b,t]=ffreadmesh('capacitorp1.msh');
[uC]=ffreaddata('capacitor_potential_p1only.txt');
%ZX: real part = const
%ZY: imag part = const
[ZX, ZY] = ffcplxmesh([0.5,0], [4,2*pi()], [10,10], [10,11]);
%Map to polar coordinates
strFunc='@(Z)(real(Z).*exp(1i*imag(Z)))';
f = str2func(strFunc);
ZU = f(ZX);
ZV = f(ZY);
[WU, WV] = ffipocplx(p,b,t,uC,ZU,ZV);
figure;
hold on;
%constant radius
plot3(real(ZU),imag(ZU),real(WU),'g','LineWidth',1);
%constant angle
plot3(real(ZV),imag(ZV),real(WV),'g','LineWidth',1);
ffpdeplot(p,b,t, ...
          'XYData',uC, ...
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
W = ffipocplx(p,b,t,uC,Z);
figure('position', [0, 0, 800, 300])
subplot(1,2,1);
hold on;
plot3(real(Z),imag(Z),real(W),'g');
ffpdeplot(p,b,t, ...
          'XYData',uC, ...
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
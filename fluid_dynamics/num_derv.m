clear;
clc;

N=100;

x=linspace(0,2*pi,N);
dx=x(2)-x(1);
y=sin(x);
dydx=cos(x);
d2ydx2=-sin(x);

D=dervx(N)*0.5/dx;
D2=dervx2(N)/dx^2;
dydxn=D*y';
d2ydx2n=(D2)*y';

plot(x,[dydx;dydxn']);
plot(x,[d2ydx2;d2ydx2n']);
clear;
clc;

N=40;
k=1;
x=linspace(0,2*pi,N);
dx=x(2)-x(1);
dt=0.002;
u(:,1)=sin(k*x);
tf=10;

for t=1:tf/dt;
    for i=1:N-1
        u(i,t+1)=(-dt/dx)*(u(i+1,t)-u(i,t))*u(i,t)+u(i,t);
    end
    u(N,t+1)=u(1,t+1);
    
    plot(x,u(:,t+1));
    pause();
end
clear;
clc;

Lx=0.15;
Ly=0.30;
Nx=20;
Ny=40;
dt=1;
tf=60*60;
alpha=1E-5;

TN=20;
TS=90;
TW=90;
TE=90;
Ti=20;

x=linspace(0,Lx,Nx);
y=linspace(0,Ly,Ny);
dx=x(2)-x(1);
dy=y(2)-y(1);
[X,Y]=meshgrid(x,y);

for i=1:Nx
    for j=1:Ny
        T(i,j)=Ti;
    end
end

for i=1:Nx
   T(i,1)=TS;
   T(i,Ny)=TN;
end

for j=1:Ny
   T(1,j)=TE;
   T(Nx,j)=TW;
end

for t=1:tf/dt
    
[Tx,Ty]=gradient(T,dx,dy);
[Txx,Tyx]=gradient(Tx,dx,dy);
[Tyx,Tyy]=gradient(Ty,dx,dy);
    
T=alpha*dt*(Txx+Tyy)+T;

for i=1:Nx
   T(i,1)=TS;
   %T(i,Ny)=TN;
end

for j=1:Ny
   T(1,j)=TE;
   T(Nx,j)=TW;
end

contourf(X,Y,T');
title(['Temperatura [°C] - t = ' num2str(round((t)/dt/60,2),'%.2f') ' min']); 
axis([0 Lx 0 Ly]);
pbaspect([Lx Ly Ly]);
colorbar;
colormap('hot');
xlabel('x [m]');
ylabel('y [m]');
pause(0.1);

end
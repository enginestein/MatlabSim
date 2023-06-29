clear;
clc;

Lx=0.15;
Ly=0.30;
Nx=20;
Ny=40;
dt=10;
tf=60*60;
T0=20;
rho0=2400;
fck0=20;

x=linspace(0,Lx,Nx);
y=linspace(0,Ly,Ny);
dx=x(2)-x(1);
dy=y(2)-y(1);
[X,Y]=meshgrid(x,y);

for i=1:Nx
    for j=1:Ny
        T(i,j)=T0;
    end
end

v = VideoWriter('concreto1.mp4');
open(v);
for t=1:tf/dt
    
for i=1:Nx
    for j=1:Ny
        cp(i,j)=cp_c(T(i,j));
        k(i,j)=k_c(T(i,j),1);
        rho(i,j)=rho_c(rho0,T(i,j));
        fck(i,j)=fck0*1.004*exp(-1*((T(i,j)-104.3)/531.8)^2);
    end
end

[Tx,Ty]=gradient(T,dx,dy);
[Txx,Txy]=gradient(k.*Tx,dx,dy);
[Tyx,Tyy]=gradient(k.*Ty,dx,dy);
    
T=dt*(Txx+Tyy)./(rho.*cp)+T;

for i=1:Nx
   T(i,1)=345*log10(8*t*dt/60+1)+T0;
   %T(i,Ny)=TN;
end

for j=1:Ny
   T(1,j)=345*log10(8*t*dt/60+1)+T0;
   T(Nx,j)=345*log10(8*t*dt/60+1)+T0;
end

subplot(1,2,1);
contourf(100*X,100*Y,T');
axis([0 100*Lx 0 100*Ly]);
pbaspect([Lx Ly Ly]);
colorbar;
h = colorbar;
set(get(h,'title'),'string','T [°C]');
colormap('hot');
xlabel('x [cm]');
ylabel('y [cm]');

subplot(1,2,2);
contourf(100*X,100*Y,fck');
axis([0 100*Lx 0 100*Ly]);
pbaspect([Lx Ly Ly]);
colorbar;
h = colorbar;
set(get(h,'title'),'string','fck [MPa]');
colormap('jet');
xlabel('x [cm]');
ylabel('y [cm]');

suptitle(['Viga de concreto em situação de incêndio - t = ' num2str(round(t*dt/60,2),'%.2f') ' min']); 

%pause(0.01);
frame = getframe(gcf);
writeVideo(v,frame);

end
close(v);
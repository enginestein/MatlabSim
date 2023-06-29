% Solves 2D navier-stokes using the vorticity-streamfunction equation for
% the lid-driven cavity flow at Re = 100

clear;
clc;

N=32; % number of spatial points

x=linspace(0,1,N); % spatial mesh
y=linspace(0,1,N); % spatial mesh
[X,Y]=meshgrid(x,y);
dx=x(2)-x(1); % space discretization
dt=0.0005; % time step
tf=0.5; % final time
nu=1; % kinematic viscosity
Uinf=100; % velocity scale
Co=nu*dt/(dx*dx)
itmax=4000; % Max iteration for poisson solver

% Flags
f1=0; % flag for viscous term
f2=0; % Flag for recording video

u=zeros(N,N);
v=zeros(N,N);
w=zeros(N,N);
    
% v = VideoWriter('ns1.mp4');
% open(v);

for t=1:tf/dt;

    % Boundary Conditions
    
    for j=1:N
       u(N,j)=Uinf; 
       u(1,j)=0; 
       u(j,1)=0; 
       u(j,N)=0; 
       v(N,j)=0; 
       v(1,j)=0; 
       v(j,1)=0; 
       v(j,N)=0; 
    end

    [dudx,dudy] = gradient(u,dx);
    [dvdx,dvdy] = gradient(v,dx);

    w=dudy-dvdx;
    
    [dwdx,dwdy] = gradient(w,dx);
    [d2wdx2,dwdxdy] = gradient(dwdx,dx);
    [dwdydx,d2wdy2] = gradient(dwdy,dx);
            
    w2=dt*(-dwdx.*u-dwdy.*v+nu*(d2wdx2+d2wdy2))+w;
    w=w2;

    psi=poisson_solver(N,N,itmax,-w,dx,0.000001);
    
    [v,u] = gradient(psi,dx);
    v=-v;
    
  % Boundary Conditions
    
    for j=1:N
       u(N,j)=Uinf; 
       u(1,j)=0; 
       u(j,1)=0; 
       u(j,N)=0; 
       v(N,j)=0; 
       v(1,j)=0; 
       v(j,1)=0; 
       v(j,N)=0; 
    end

    [dudx,dudy] = gradient(u,dx);
    [dvdx,dvdy] = gradient(v,dx);

    w=dudy-dvdx;
    
    f=(dudx.^2+2*dudy.*dvdx+dvdy.^2);
    p=poisson_solver(N,N,itmax,f,dx,0.000001);
    
    contourf(X,Y,psi,20);
    %quiver(X,Y,u,v);
    set(gca,'FontSize',12);
    colormap('jet');
    title(['\psi - Re = ' num2str(Uinf) ' - t = ' num2str(t*dt)]); 
    axis([0,1,0,1]);
    xlabel('x');
    ylabel('y');
    axis('square');
    pause(0.05);
    
%     frame = getframe(gcf);
%     writeVideo(v,frame);
end

%close(v);
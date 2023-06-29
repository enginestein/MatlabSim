% Solves 2D navier-stokes using the vorticity-streamfunction equation for
% the lid-driven cavity flow at Re = 100

clear;
clc;

N=64; % number of spatial points

x=linspace(0,1,N); % spatial mesh
y=linspace(0,1,N); % spatial mesh
[X,Y]=meshgrid(x,y);
dx=x(2)-x(1); % space discretization
dt=0.001; % time step
tf=0.5; % final time
nu=1E-8; % kinematic viscosity
Uinf=1; % velocity scale
Co=nu*dt/(dx*dx)
itmax=4000; % Max iteration for poisson solver
tol=0.0000001;

% Flags
f1=0; % flag for viscous term
f2=0; % Flag for recording video

% Intial conditions

for i=1:N
    for j=1:N
        u(i,j)=sin((i-1)*dx);
        v(i,j)=cos((i-1*dx));
    end
end

[dudx,dudy] = gradient(u,dx);
[dvdx,dvdy] = gradient(v,dx);

w=dudy-dvdx;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V = VideoWriter('taylorgreen.mp4');
open(V);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t=1:tf/dt;

    % Boundary Conditions
    
    for j=1:N
       u(N,j)=u(1,j); 
       u(j,1)=u(j,N); 
       v(N,j)=v(1,j); 
       v(j,1)=v(j,N); 
    end

    [dudx,dudy] = gradient(u,dx);
    [dvdx,dvdy] = gradient(v,dx);

    w=dudy-dvdx;
    
    [dwdx,dwdy] = gradient(w,dx);
    [d2wdx2,dwdxdy] = gradient(dwdx,dx);
    [dwdydx,d2wdy2] = gradient(dwdy,dx);
            
    w2=dt*(-dwdx.*u-dwdy.*v+nu*(d2wdx2+d2wdy2))+w;
    w=w2;

    psi=poisson_solver(N,N,itmax,-w,dx,tol);
    
    [v,u] = gradient(psi,dx);
    v=-v;
    
  % Boundary Conditions
    
    for j=1:N
       u(N,j)=u(1,j); 
       u(j,1)=u(j,N); 
       v(N,j)=v(1,j); 
       v(j,1)=v(j,N); 
    end

    [dudx,dudy] = gradient(u,dx);
    [dvdx,dvdy] = gradient(v,dx);

    w=dudy-dvdx;
    
    f=(dudx.^2+2*dudy.*dvdx+dvdy.^2);
    p=poisson_solver(N,N,itmax,f,dx,tol);
    
    contourf(X,Y,psi,20);
    %quiver(X,Y,u,v);
    set(gca,'FontSize',12);
    colormap('jet');
    colorbar;
    title(['\psi - Re = ' num2str(Uinf/nu) ' - t = ' num2str(t*dt)]); 
    axis([0,1,0,1]);
    xlabel('x');
    ylabel('y');
    axis('square');
    %pause(0.05);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    frame = getframe(gcf);
    writeVideo(V,frame);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

close(V);
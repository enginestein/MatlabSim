% Solves 2D navier-stokes using the vorticity-streamfunction equation for
% the lid-driven cavity flow at Re = 100

clear;
clc;

Nx=20; % number of spatial points
Ny=20;
x=linspace(0,1,Nx); % spatial mesh
y=linspace(0,1,Ny); % spatial mesh
[X,Y]=meshgrid(x,y);
dx=x(2)-x(1); % space discretization
dy=y(2)-y(1);
dt=0.01; % time step
tf=0.5; % final time
Re=1; % Reynolds number

% Generate differentiation matrices
dX=dervx(Nx);
dX2=dervx2(Nx);
dY=dervx(Ny);
dY2=dervx2(Ny);

% Identity matrices
Ix=eye(Nx);
Iy=eye(Ny);

% Differential operators
% First derivatives
Dx=kron(Ix,dX);
Dy=kron(dY,Iy);

Dx=0.5*Dx./dx;
Dy=0.5*Dy./dy;

% Second derivatives
Dxx=kron(Ix,dX2);
Dyy=kron(dY2,Iy);

Dxx=Dxx./(dx^2);
Dyy=Dyy./(dy^2);

% Flags
f1=0; % flag for viscous term
f2=0; % Flag for recording video

u=zeros(Nx,Ny);
v=zeros(Nx,Ny);
w=zeros(Nx,Ny);
    
% v = VideoWriter('ns1.mp4');
% open(v);

for t=1:tf/dt;
    
    % Boundary Conditions
    
    for j=1:Ny
       u(Nx,j)=1; 
       u(1,j)=0; 
       v(Nx,j)=0; 
       v(1,j)=0; 
    end
    
    for i=1:Nx
       u(i,Ny)=0; 
       u(i,1)=0; 
       v(i,Ny)=0; 
       v(i,1)=0; 
    end
    
    u=reshape(u,Nx*Ny,1);
    v=reshape(v,Nx*Ny,1);
    
    dudx=Dx*u;
    dudy=Dy*u;
    dvdx=Dx*v;
    dvdy=Dy*v;

    w=dudy-dvdx;
    
    dwdx=Dx*w;
    dwdy=Dy*w;
    d2wdx2=Dxx*w;
    d2wdy2=Dyy*w;
            
    w2=dt*(-dwdx.*u-dwdy.*v+(1/Re)*(d2wdx2+d2wdy2))+w;
    w=w2;

    psi=(Dxx+Dyy)\w;
    
    u=Dy*psi;
    v=-Dx*psi;
    
  % Boundary Conditions
    
    u=reshape(u,Nx,Ny);
    v=reshape(v,Nx,Ny);
  
    for j=1:Ny
       u(Nx,j)=1; 
       u(1,j)=0; 
       v(Nx,j)=0; 
       v(1,j)=0; 
    end
    
    for i=1:Nx
       u(i,Ny)=0; 
       u(i,1)=0; 
       v(i,Ny)=0; 
       v(i,1)=0; 
    end
    
    u=reshape(u,Nx*Ny,1);
    v=reshape(v,Nx*Ny,1);

    dudx=Dx*u;
    dudy=Dy*u;
    dvdx=Dx*v;
    dvdy=Dy*v;

    w=dudy-dvdx;
    
    f=(dudx.^2+2*dudy.*dvdx+dvdy.^2);
    p=(Dxx+Dyy)\f;
    
    u=reshape(u,Nx,Ny);
    v=reshape(v,Nx,Ny);
    w=reshape(w,Nx,Ny);
    psi=reshape(psi,Nx,Ny);
    p=reshape(p,Nx,Ny);
    
    contourf(X,Y,psi,20);
    %quiver(X,Y,u,v);
    set(gca,'FontSize',12);
    colormap('jet');
    title(['\psi - Re = ' num2str(Re) ' - t = ' num2str(t*dt)]); 
    axis([0,1,0,1]);
    xlabel('x');
    ylabel('y');
    axis('square');
    pause(0.05);
    
    w=reshape(w,Nx*Ny,1);
    psi=reshape(psi,Nx*Ny,1);
    p=reshape(p,Nx*Ny,1);
    
%     frame = getframe(gcf);
%     writeVideo(v,frame);
end

%close(v);
% Solves 2D navier-stokes
clear;
clc;

Nx=32; % number of spatial points
Ny=32;
Lx=1; % domain size
Ly=1;
dt=0.01; % time step
tf=100; % final time
Re=1000; % Reynolds number
Ri=0; % Richardson number
Sc=1; % Schmidt number
itmax=4000; % Max iteration for poisson solver
tol=0.00001; % tolerance for poisson solver

x=linspace(0,Lx,Nx); % spatial mesh
y=linspace(0,Ly,Ny); % spatial mesh
[X,Y]=ndgrid(x,y);
dx=x(2)-x(1); % space discretization
dy=y(2)-y(1);
Co=(1/Re)*dt/(dx*dx)

% Flags
f1=0; % flag for viscous term
f2=0; % Flag for recording video

% Intial conditions

for i=1:Nx
    for j=1:Ny
        u(i,j)=0;
        v(i,j)=0;
        phi(i,j)=1;
    end
end

y1=linspace(-1,1,Ny);
for j=1:Ny
u(:,j)=tanh(10*y1(j));
end

[dudx,dudy] = gradient(u,dx,dy);
[dvdx,dvdy] = gradient(v,dx,dy);

w=dudy-dvdx;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V = VideoWriter('taylorgreen.mp4');
% open(V);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t=1:tf/dt;

    % Boundary Conditions
    
    for j=1:Ny
       u(Nx,j)=u(1,j); 
       v(Nx,j)=v(1,j); 
    end
    
    for i=1:Nx
       u(i,1)=0; 
       u(i,Ny)=0; 
       v(i,1)=0; 
       v(i,Ny)=0; 
    end

    [dudx,dudy] = gradient(u,dx);
    [dvdx,dvdy] = gradient(v,dx);

    w=dudy-dvdx;
    
    [dwdx,dwdy] = gradient(w,dx);
    [d2wdx2,dwdxdy] = gradient(dwdx,dx);
    [dwdydx,d2wdy2] = gradient(dwdy,dx);
    
    [dphidx,dphidy] = gradient(phi,dx);
    [d2phidx2,dphidxdy] = gradient(dphidx,dx);
    [dphidydx,d2phidy2] = gradient(dphidy,dx);
            
    w2=dt*(-dwdx.*u-dwdy.*v+(1/Re)*(d2wdx2+d2wdy2)+Ri.*dphidx)+w;
    w=w2;
    
    phi2=dt*(-dphidx.*u-dphidy.*v+(1/(Re*Sc))*(d2phidx2+d2phidy2))+phi;
    phi=phi2;

    psi=poisson_solver2(Nx,Ny,itmax,-w,dx,dy,tol);
    
    [v,u] = gradient(psi,dx);
    v=-v;
    
  % Boundary Conditions
    
    for j=1:Ny
       u(Nx,j)=u(1,j); 
       v(Nx,j)=v(1,j); 
    end
    
    for i=1:Nx
       u(i,1)=0; 
       u(i,Ny)=0; 
       v(i,1)=0; 
       v(i,Ny)=0; 
    end

    [dudx,dudy] = gradient(u,dx);
    [dvdx,dvdy] = gradient(v,dx);

    w=dudy-dvdx;
    
    f=(dudx.^2+2*dudy.*dvdx+dvdy.^2);
    p=poisson_solver2(Nx,Ny,itmax,f,dx,dy,tol);
    
    contourf(X,Y,u,20);
    %quiver(X,Y,u,v);
    set(gca,'FontSize',12);
    colormap('jet');
    colorbar;
    title(['\psi - Re = ' num2str(Re) ' - t = ' num2str(t*dt)]); 
    axis([0,Lx,0,Ly]);
    xlabel('x');
    ylabel('y');
    axis('square');
    pbaspect([Lx Ly 1]);
    pause(0.05);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     frame = getframe(gcf);
%     writeVideo(V,frame);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%close(V);
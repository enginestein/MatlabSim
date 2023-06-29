% Solves 2D navier-stokes
clear;
clc;

Nx=128; % number of spatial points
Ny=128;
Lx=10; % domain size
Ly=10;
cex=Lx/2;
cey=Ly/2;
dt=0.001; % time step
tf=100; % final time
Umag=1;
Re=100; % Reynolds number
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
        u(i,j)=1;
        v(i,j)=0;
    end
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
       u(1,j)=Umag; 
       u(Nx,j)=u(Nx-1,j); 
       v(1,j)=0; 
       v(Nx,j)=v(Nx-1,j); 
    end
    
    for i=1:Nx
       u(i,1)=0; 
       u(i,Ny)=0;
       v(i,1)=0; 
       v(i,Ny)=0;
    end
    
    for i=1:Nx
    for j=1:Ny
        x=dx*(i-1);
        y=dy*(j-1);
        
        if sqrt((x-cex)^2+(y-cey)^2)<=0.5
            u(i,j)=0;
            v(i,j)=0;
        end
    end
    end

    [dudx,dudy] = gradient(u,dx);
    [dvdx,dvdy] = gradient(v,dx);

    w=dudy-dvdx;
    
    [dwdx,dwdy] = gradient(w,dx);
    [d2wdx2,dwdxdy] = gradient(dwdx,dx);
    [dwdydx,d2wdy2] = gradient(dwdy,dx);
            
    w2=dt*(-dwdx.*u-dwdy.*v+(1/Re)*(d2wdx2+d2wdy2))+w;
    w=w2;
    
    psi=poisson_solver2(Nx,Ny,itmax,-w,dx,dy,tol);
    
    [v,u] = gradient(psi,dx);
    v=-v;
    
  % Boundary Conditions
    
    for j=1:Ny
       u(1,j)=u(2,j); 
       u(Nx,j)=u(Nx-1,j); 
       v(1,j)=v(2,j); 
       v(Nx,j)=v(Nx-1,j); 
    end
    
    for i=1:Nx
       u(i,1)=Umag; 
       u(i,Ny)=Umag;
       v(i,1)=0; 
       v(i,Ny)=v(i,Ny-1);
    end
    
        for i=1:Nx
    for j=1:Ny
        x=dx*(i-1);
        y=dy*(j-1);
        
        if sqrt((x-cex)^2+(y-cey)^2)<=0.5
            u(i,j)=0;
            v(i,j)=0;
        end
    end
    end

    [dudx,dudy] = gradient(u,dx);
    [dvdx,dvdy] = gradient(v,dx);

    w=dudy-dvdx;
    
    f=(dudx.^2+2*dudy.*dvdx+dvdy.^2);
    p=poisson_solver2(Nx,Ny,itmax,f,dx,dy,tol);
    
    umag=sqrt(u.^2+v.^2);
    contourf(X,Y,umag,20);
    %quiver(X,Y,u,v);
    set(gca,'FontSize',12);
    colormap('jet');
    colorbar;
    title(['U_{mag} - Re = ' num2str(Re) ' - t = ' num2str(t*dt)]); 
    axis([0,Lx,0,Ly]);
    xlabel('x');
    ylabel('y');
    axis('square');
    pbaspect([Lx Ly 1]);
    pause(0.5);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     frame = getframe(gcf);
%     writeVideo(V,frame);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%close(V);
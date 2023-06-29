clear;
clc;

N=11; % number of spatial points

x=linspace(0,1,N); % spatial mesh
y=linspace(0,1,N); % spatial mesh
[X,Y]=meshgrid(x,y);
dx=x(2)-x(1); % space discretization
eps=0.1; % magnitude of perturbation
dt=0.1; % time step
tf=10; % final time
Re=100; % Reynolds number
itmax=4000;

% Flags
f1=0; % flag for viscous term
f2=0; % Flag for recording video

u=zeros(N,N);
v=zeros(N,N);
p=zeros(N,N);

dudx=zeros(N,N);
dudy=zeros(N,N);
dvdx=zeros(N,N);
dvdy=zeros(N,N);

% Initial condition
for i=1:N
    for j=1:N
        u(i,j,1)=0;
        v(i,j,1)=0;
        p(i,j,1)=0;
    end
end

%Boundary conditions

for j=1:N
   u(N,j)=1; 
end
    
% v = VideoWriter('ns1.mp4');
% open(v);

for t=1:tf/dt;
    for i=2:N-1
        for j=2:N-1
            [dudx,dudy] = gradient(u,dx);
            [dvdx,dvdy] = gradient(v,dx);
            [dpdx,dpdy] = gradient(p,dx);
            
            [d2udx2,dudxdy] = gradient(dudx,dx);
            [dvdydx,d2vdy2] = gradient(dvdy,dx);
            [d2udydx,d2udy2] = gradient(dudy,dx);
            [d2vdx2,d2vdxdy] = gradient(dvdx,dx);
            
            u2(i,j)=dt*(-dudx(i,j)*u(i,j)-dudy(i,j)*v(i,j)-dpdx(i,j)+(1/Re)*(d2udx2(i,j)+d2udy2(i,j)))+u(i,j);
            v2(i,j)=dt*(-dvdx(i,j)*u(i,j)-dvdy(i,j)*v(i,j)-dpdy(i,j)+(1/Re)*(d2vdx2(i,j)+d2vdy2(i,j)))+v(i,j);
        end
    end
    
    %%% Solves poisson
    
    [dudx,dudy] = gradient(u2,dx);
    [dvdx,dvdy] = gradient(v2,dx);

    f=(dudx.^2+2*dudy.*dvdx+dvdy.^2);
    p=poisson_pressure(N,N,itmax,f,dx);
     
     %%%%%%%%%% Obtain new values
     
     for i=2:N-1
        for j=2:N-1
            [dudx,dudy] = gradient(u2,dx);
            [dvdx,dvdy] = gradient(v2,dx);
            [dpdx,dpdy] = gradient(p,dx);
            
            [d2udx2,dudxdy] = gradient(dudx,dx);
            [dvdydx,d2vdy2] = gradient(dvdy,dx);
            [d2udydx,d2udy2] = gradient(dudy,dx);
            [d2vdx2,d2vdxdy] = gradient(dvdx,dx);
            
            u3(i,j)=dt*(-dudx(i,j)*u(i,j)-dudy(i,j)*v(i,j)-dpdx(i,j)+(1/Re)*(d2udx2(i,j)+d2udy2(i,j)))+u2(i,j);
            v3(i,j)=dt*(-dvdx(i,j)*u(i,j)-dvdy(i,j)*v(i,j)-dpdy(i,j)+(1/Re)*(d2vdx2(i,j)+d2vdy2(i,j)))+v2(i,j);
        end
     end
     
     u=u3;
     v=v3;
     
    %Boundary conditions %%%%%%%%%%%%%%%%%%

    for i=1:N

       u(1,i)=0;
       u(i,1)=0;
       u(i,N)=0;
       
       v(N,i)=0;
       v(1,i)=0;
       v(i,1)=0;
       v(i,N)=0;
       u(N,i)=1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    contourf(X,Y,p,20);
    title(['t = ' num2str(t*dt)]); 
    axis([0,1,0,1]);
    axis('square');
    pause(0.1);
    
%     frame = getframe(gcf);
%     writeVideo(v,frame);
end

%close(v);
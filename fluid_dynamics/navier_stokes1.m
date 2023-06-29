clear;
clc;

N=11; % number of spatial points

x=linspace(0,1,N); % spatial mesh
y=linspace(0,1,N); % spatial mesh
[X,Y]=meshgrid(x,y);
dx=x(2)-x(1); % space discretization
eps=0.1; % magnitude of perturbation
dt=0.01; % time step
tf=10; % final time
Re=100; % Reynolds number

% Flags
f1=0; % flag for viscous term
f2=0; % Flag for recording video

u=zeros(N,N);
v=zeros(N,N);
p=zeros(N,N);

u1=zeros(N-1,N-1);
v1=zeros(N-1,N-1);
p1=zeros(N-1,N-1);

dudx=zeros(N,N);
dudy=zeros(N,N);
dvdx=zeros(N,N);
dvdy=zeros(N,N);

dudx1=zeros(N-1,N-1);
dudy1=zeros(N-1,N-1);
dvdx1=zeros(N-1,N-1);
dvdy1=zeros(N-1,N-1);


% Initial condition
for i=1:N
    for j=1:N
        u(i,j,1)=0;
        v(i,j,1)=0;
    end
end

%Boundary conditions

for j=1:N
   u(N,j)=1; 
end

% Separated grid
    
for i=1:N-1
    for j=1:N-1
        u1(i,j)=(u(i,j)+u(i,j+1)+u(i+1,j+1)+u(i+1,j))/4;
        v1(i,j)=(v(i,j)+v(i,j+1)+v(i+1,j+1)+v(i+1,j))/4;
        p1(i,j)=(p(i,j)+p(i,j+1)+p(i+1,j+1)+p(i+1,j))/4;
    end
end

% Solve poisson equation

[dudx1,dudy1] = gradient(u1,dx);
[dvdx1,dvdy1] = gradient(v1,dx);

for i=2:N-2
    for j=2:N-2
        p1(i,j)=(dx*dx*(dudx1(i,j).^2+2*dudy1(i,j)*dvdx1(i,j)+dvdy1(i,j).^2)+p1(i+1,j)+p1(i-1,j)+p1(i,j+1)+p1(i,j-1))/4;
    end
end

 % Staggered grid
     
     p(1,1)=p1(1,1);
     p(N,1)=p1(N-1,1);
     p(1,N)=p1(1,N-1);
     p(N,N)=p1(N-1,N-1);
     
     for i=2:N-1
         p(i,1)=0.5*(p1(i-1,1)+p1(i,1));
         p(i,N)=0.5*(p1(i-1,N-1)+p1(i,N-1));
         p(N,i)=0.5*(p1(N-1,i-1)+p1(N-1,i));
         p(1,i)=0.5*(p1(1,i-1)+p1(1,i));
     end
     
     for i=2:N-2
         for j=2:N-2
            p(i,j)=(p1(i+1,j)+p1(i,j+1)+p1(i,j)+p1(i,j))/4;
         end
     end
    
% v = VideoWriter('ns1.mp4');
% open(v);

for t=1:tf/dt;
    for i=1:N
        for j=1:N
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
            
    for i=1:N-1
       for j=1:N-1
          u1(i,j)=(u2(i,j)+u2(i,j+1)+u2(i+1,j+1)+v2(i+1,j))/4;
          v1(i,j)=(v2(i,j)+v2(i,j+1)+v2(i+1,j+1)+v2(i+1,j))/4;
          p1(i,j)=(p(i,j)+p(i,j+1)+p(i+1,j+1)+p(i+1,j))/4;
       end
    end
    
    [dudx1,dudy1] = gradient(u1,dx);
    [dvdx1,dvdy1] = gradient(v1,dx);
            
     for it=1:100
         for i=2:N-2
           for j=2:N-2
                p1(i,j)=(dx*dx*(dudx1(i,j).^2+2*dudy1(i,j)*dvdx1(i,j)+dvdy1(i,j).^2)+p1(i+1,j)+p1(i-1,j)+p1(i,j+1)+p1(i,j-1))/4;
            end
         end
     end
     
     % Staggered grid
     
     p2(1,1)=p1(1,1);
     p2(N,1)=p1(N-1,1);
     p2(1,N)=p1(1,N-1);
     p2(N,N)=p1(N-1,N-1);
     
     for i=2:N-1
         p2(i,1)=0.5*(p1(i-1,1)+p1(i,1));
         p2(i,N)=0.5*(p1(i-1,N-1)+p1(i,N-1));
         p2(N,i)=0.5*(p1(N-1,i-1)+p1(N-1,i));
         p2(1,i)=0.5*(p1(1,i-1)+p1(1,i));
     end
     
     for i=2:N-2
         for j=2:N-2
            p2(i,j)=(p1(i+1,j)+p1(i,j+1)+p1(i,j)+p1(i,j))/4;
         end
     end
     
     %%%%%%%%%%%%%%
     
     
     %%%%%%%%%% Obtain new values
     
     for i=1:N
        for j=1:N
            [dudx,dudy] = gradient(u2,dx);
            [dvdx,dvdy] = gradient(v2,dx);
            [dpdx,dpdy] = gradient(p2,dx);
            
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
     p=p2;
     
    %Boundary conditions %%%%%%%%%%%%%%%%%%

    for i=1:N
       u(N,i)=1;
       u(1,i)=0;
       u(i,1)=0;
       u(i,N)=0;
       
       v(N,i)=0;
       v(1,i)=0;
       v(i,1)=0;
       v(i,N)=0;
    end
    
    % Separated grid
    
for i=1:N-1
    for j=1:N-1
        u1(i,j)=(u(i,j)+u(i,j+1)+u(i+1,j+1)+u(i+1,j))/4;
        v1(i,j)=(v(i,j)+v(i,j+1)+v(i+1,j+1)+v(i+1,j))/4;
        p1(i,j)=(p(i,j)+p(i,j+1)+p(i+1,j+1)+p(i+1,j))/4;
    end
end

% Solve poisson equation

[dudx1,dudy1] = gradient(u1,dx);
[dvdx1,dvdy1] = gradient(v1,dx);

for i=2:N-2
    for j=2:N-2
        p1(i,j)=(dx*dx*(dudx1(i,j).^2+2*dudy1(i,j)*dvdx1(i,j)+dvdy1(i,j).^2)+p1(i+1,j)+p1(i-1,j)+p1(i,j+1)+p1(i,j-1))/4;
    end
end

 % Staggered grid
     
     p(1,1)=p1(1,1);
     p(N,1)=p1(N-1,1);
     p(1,N)=p1(1,N-1);
     p(N,N)=p1(N-1,N-1);
     
     for i=2:N-1
         p(i,1)=0.5*(p1(i-1,1)+p1(i,1));
         p(i,N)=0.5*(p1(i-1,N-1)+p1(i,N-1));
         p(N,i)=0.5*(p1(N-1,i-1)+p1(N-1,i));
         p(1,i)=0.5*(p1(1,i-1)+p1(1,i));
     end
     
     for i=2:N-2
         for j=2:N-2
            p(i,j)=(p1(i+1,j)+p1(i,j+1)+p1(i,j)+p1(i,j))/4;
         end
     end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    contour(X,Y,abs(dudy-dvdx),50);
    title(['t = ' num2str(t*dt)]); 
    axis([0,1,0,1]);
    axis('square');
    pause(0.1);
    
%     frame = getframe(gcf);
%     writeVideo(v,frame);
end

%close(v);
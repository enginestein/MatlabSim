clear;
clc;

dt=0.0005; % time step
itmax=1000; % number of iterations
N=10; % number of particles
Lx=1; % domain size in x
Ly=1; % domain size in y
Lz=1; % domain size in z
vmag=1; % velocity magnitude order

x0=rand(N,1)';
y0=rand(N,1)';
z0=rand(N,1)';

for n=1:N
    v0x(n)=vmag*rand(1)*(-1)^round(rand(1)*10,0);
    v0y(n)=vmag*rand(1)*(-1)^round(rand(1)*10,0);
    v0z(n)=vmag*rand(1)*(-1)^round(rand(1)*10,0);
end

count=0; %number of colisions

for i=0:1:itmax
    
    x=x0+v0x*dt;
    y=y0+v0y*dt;
    z=z0+v0z*dt;
    
    x0=x;
    y0=y;
    z0=z;
    
    % Check for colisions between particles
    for k=1:N-1
        if round(x(k),3)==round(x(k+1),3)&&round(y(k),3)==round(y(k+1),3)&&round(z(k),3)==round(z(k+1),3)
            v0x(k)=-(v0x(k)+v0x(k+1))*0.5*v0x(k)/abs(v0x(k));
            v0y(k)=-(v0y(k)+v0y(k+1))*0.5*v0y(k)/abs(v0y(k));
            v0z(k)=-(v0z(k)+v0z(k+1))*0.5*v0z(k)/abs(v0z(k));
            
            v0x(k+1)=-(v0x(k)+v0x(k+1))*0.5*v0x(k+1)/abs(v0x(k+1));
            v0y(k+1)=-(v0y(k)+v0y(k+1))*0.5*v0y(k+1)/abs(v0y(k+1));
            v0z(k+1)=-(v0z(k)+v0z(k+1))*0.5*v0z(k+1)/abs(v0z(k+1));
            count=count+1;
        end
    end
    
    col(i+1)=count;
    
    % Check if it particle is out of bounds
    for k=1:N
        if(x(k)<0||x(k)>Lx)
            v0x(k)=-v0x(k);
        end
        if(y(k)<0||y(k)>Ly)
            v0y(k)=-v0y(k); 
        end
        if(z(k)<0||z(k)>Ly)
            v0z(k)=-v0z(k); 
        end
    end
            
    % Draw each particle
    for k=1:N
        plot3(x(k),y(k),z(k),'-o','Color','b','MarkerSize',10,'MarkerFaceColor','b');
        title(mean(sqrt(v0x(:).^2+v0y(:).^2+v0z(:).^2)));
        axis('square');
        %filledCircle([x(k),y(k)],0.02,10,'k');
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        grid on
        hold on;
    end
    hold off;
    
    axis ([0 Lx 0 Ly 0 Lz]);
    pause(0.01);
end
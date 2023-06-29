clear;
clc;

dt=0.005; % time step
itmax=1000; % number of iterations
N=60; % number of particles
Lx=1; % domain size in x
Ly=1; % domain size in y
vmag=10; % velocity magnitude order
or=2; % collision order of sensitivty

x0=Lx*rand(N,1)';
y0=Ly*rand(N,1)';

for n=1:N
    v0x(n)=vmag*rand(1)*(-1)^round(rand(1)*10,0);
    v0y(n)=vmag*rand(1)*(-1)^round(rand(1)*10,0);
    color(n)=0;
end

count=0; %number of colisions
v = VideoWriter('test5.avi');
open(v);

for i=0:1:itmax
    
    x=x0+v0x*dt;
    y=y0+v0y*dt;
    
    x0=x;
    y0=y;
    
    % Check for colisions between particles
    for k=1:N-1
        if round(x(k),or)==round(x(k+1),or)&&round(y(k),or)==round(y(k+1),or)
            v0x(k)=-(v0x(k)+v0x(k+1))*0.5*v0x(k)/abs(v0x(k));
            v0y(k)=-(v0y(k)+v0y(k+1))*0.5*v0y(k)/abs(v0y(k));
            v0x(k+1)=-(v0x(k)+v0x(k+1))*0.5*v0x(k+1)/abs(v0x(k+1));
            v0y(k+1)=-(v0y(k)+v0y(k+1))*0.5*v0y(k+1)/abs(v0y(k+1));
            count=count+1;
            color(k)=1;
            color(k+1)=1;
        end
    end
    
    nump(i+1)=0;
    for k=1:N
        if (x(k)>Lx/4&&x(k)<3*Lx/4&&y(k)>Ly/4&&y(k)<3*Ly/4)
            nump(i+1)=nump(i+1)+1;
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
    end
            
    % Draw each particle
    for k=1:N
        title(['Number of particles: ' num2str(nump(i+1))]);
        if color(k)==1
            plot(x(k),y(k),'-o','Color','r','MarkerSize',10,'MarkerFaceColor','r');
        else
            plot(x(k),y(k),'-o','Color','b','MarkerSize',10,'MarkerFaceColor','b');            
        end
        axis('square');
        %filledCircle([x(k),y(k)],0.02,10,'k');
        hold on;
    end
    hold off;
    axis ([Lx/4 3*Lx/4 Ly/4 3*Ly/4]);
    frame = getframe(gcf);
    writeVideo(v,frame);
    %pause(0.01);

end
t=linspace(0,(itmax),itmax+1);
plot(t,nump);

close(v);
clear;
clc;

dt=0.005; % time step
itmax=500; % number of time iterations
randit=20; % number of stocastic iterations
N=100; % number of particles
Lx=1; % domain size in x
Ly=1; % domain size in y
vmag=1; % velocity magnitude order
or=2; % collision order of sensitivty

count=0; %number of colisions

for r=1:randit
    for n=2:N

    x0=rand(n,1)';
    y0=rand(n,1)';

    for j=1:n
        v0x(j)=vmag*rand(1)*(-1)^round(rand(1)*10,0);
        v0y(j)=vmag*rand(1)*(-1)^round(rand(1)*10,0);
    end

        for i=0:1:itmax

            x=x0+v0x*dt;
            y=y0+v0y*dt;

            x0=x;
            y0=y;

            % Check for colisions between particles
        for k=1:n-1
            if round(x(k),or)==round(x(k+1),or)&&round(y(k),or)==round(y(k+1),or)
                v0x(k)=-(v0x(k)+v0x(k+1))*0.5*v0x(k)/abs(v0x(k));
                v0y(k)=-(v0y(k)+v0y(k+1))*0.5*v0y(k)/abs(v0y(k));
                v0x(k+1)=-(v0x(k)+v0x(k+1))*0.5*v0x(k+1)/abs(v0x(k+1));
                v0y(k+1)=-(v0y(k)+v0y(k+1))*0.5*v0y(k+1)/abs(v0y(k+1));
                count=count+1;
            end
        end

            col(i+1,n)=count;

            % Check if it particle is out of bounds
            for k=1:n
                if(x(k)<0||x(k)>Lx)
                    v0x(k)=-v0x(k);
                end
                if(y(k)<0||y(k)>Ly)
                    v0y(k)=-v0y(k); 
                end
            end

            % Draw each particle
    %         for k=1:N
    %             filledCircle([x(k),y(k)],0.02,10,'k');
    %             hold on;
    %         end
    %         hold off;
    % 
    %         axis ([0 Lx 0 Ly]);
    %         pause();
        end
        count=0;
    end
    clear x0;
    clear y0;
    clear v0x;
    clear v0y;
    maxcol(r,:)=max(col);
end
plot(maxcol');
plot(mean(maxcol))
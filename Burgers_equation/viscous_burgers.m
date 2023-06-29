clear;
clc;

N=100; % number of spatial points
k=2; % wavenumber
x=linspace(0,2*pi,N); % spatial mesh
dx=x(2)-x(1); % space discretization
eps=0; % magnitude of perturbation
dt=0.02; % time step
tf=10; % final time
Re=30; % Reynolds number
f1=1; % flag for viscous term

% Initial condition
for i=1:N
    u(i,1)=sin(k*x(i))+eps*rand(1);
end

% Evaluate blow-up time
dudx=k*cos(k*x);
tb=-1/min(dudx);

for t=1:tf/dt;
    for i=2:N-1
        u(i,t+1)=(-0.5*dt/dx)*(u(i+1,t)-u(i-1,t))*u(i,t)+u(i,t)+f1*(dt/Re)*(u(i+1,t)-2*u(i,t)+u(i-1,t))/dx^2;
    end
    
    plot(x,u(:,t+1),'LineWidth',2);
    title(['Viscous Burgers equation - ','t = ' num2str(t*dt) ', t_b = ' num2str(tb) ', Re = ' num2str(Re)]);
    axis('square');
    xlabel('x');
    ylabel('u(x)');
    axis([0 2*pi -2 2]);
    pause();
end
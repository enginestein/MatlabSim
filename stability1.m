clear;
clc;

l=1;
k=1;
Rc=1;
dt=0.01;
U0=1;

u(1)=0;
for R=0.9:0.01:10

    for i=1:100
        a=k*(R-Rc);
        u(i+1)=(a-l*(u(i)-U0)^2)*dt+u(i);
    end
    plot(u);
    hold on;
    pause();
end
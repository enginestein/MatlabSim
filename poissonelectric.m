%------------------------------------------------------------------------
%---Poisson Equation Solver using Successive Overrelaxation Method-------
%------------------------------------------------------------------------
clc
clear all;
close all;
%-------------------------------------------------------------------------
%--Define dimension of the 2-d box and -----------------------------------
%-------------------------------------------------------------------------
loop=1;
Nx=50;
Ny=50;
M=4000; % maximum iteration value
%--------------------------------------------------------------------------
%--Initialize V-matrix and Rho-matrix--------------------------------------
%---V(i,j)=potential inside the 2-D box------------------------------------
%---rho(i,j)=charge density inside the 2-D box-----------------------------
%---i is along X-axis and j is along Y-axis--------------------------------
%--------------------------------------------------------------------------
V(1:Nx,1:Ny)=0.0;
rho(1:Nx,1:Ny)=0.0;
rho(25,25)=25; % charge at the center of the box
%--------------------------------------------------------------------------
% Boundary conditions
V(1,:)=50; 
V(Nx,:)=0; 
V(:,1)=0; 
V(:,Ny)=0; 
%-------------------------------------------------------------------------
w=cos(pi/Nx)+cos(pi/Ny); % Converging Term
Ncount=0;
loop=1;
while loop==1;
   Rmin=0; 
for i=2:Nx-1
for j=2:Ny-1
   Residue=w.*(0.25.*(V(i-1,j)+V(i+1,j)+V(i,j-1)+V(i,j+1)+rho(i,j))-V(i,j));
   Rmin=Rmin+abs(Residue);
   V(i,j)=V(i,j)+Residue;
end
end
Rmin=Rmin/(Nx*Ny); % Average Residue per grid point
if(Rmin>=0.00001)
    Ncount=Ncount+1;
    if(Ncount>M)
        loop=0;
        disp(['solution doesnt converge in ',num2str(M),' iter'])
    end
else
    loop=0;
    disp(['solution converges in ',num2str(Ncount) ,' iteration'])
end
end
%------------------------------------------------------------------------
% Plot the result
%------------------------------------------------------------------------
 X=1:Nx;
 Y=1:Ny;
 contour(X,Y,V,0:1:50,'linewidth',2)
 h=gca; 
get(h,'FontSize') 
set(h,'FontSize',12)
colorbar('location','eastoutside','fontsize',12);
axis([1 Nx 1 Ny])
xlabel('X','fontSize',12);
ylabel('Y','fontSize',12);
title('Electric Potential Distribution, V(X,Y)','fontsize',12);
fh = figure(1);
set(fh, 'color', 'white'); 

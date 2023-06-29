function p=poisson_pressure(Nx,Ny,M,f,dx,tol)

%------------------------------------------------------------------------
%---Poisson Equation Solver using Successive Overrelaxation Method-------
%------------------------------------------------------------------------
%-------------------------------------------------------------------------
%--Define dimension of the 2-d box and -----------------------------------
%-------------------------------------------------------------------------
loop=1;
%--------------------------------------------------------------------------
%--Initialize V-matrix and Rho-matrix--------------------------------------
%---V(i,j)=potential inside the 2-D box------------------------------------
%---rho(i,j)=charge density inside the 2-D box-----------------------------
%---i is along X-axis and j is along Y-axis--------------------------------
%--------------------------------------------------------------------------
p(1:Nx,1:Ny)=0.0;
%--------------------------------------------------------------------------
% Boundary conditions
p(1,:)=0; 
p(Nx,:)=0; 
p(:,1)=0; 
p(:,Ny)=0; 
%-------------------------------------------------------------------------
Ncount=0;
loop=1;
while loop==1;
   Rmin=0; 
for i=2:Nx-1
for j=2:Ny-1
   Residue=(0.25.*(p(i-1,j)+p(i+1,j)+p(i,j-1)+p(i,j+1)+f(i,j)*dx.^2)-p(i,j));
   Rmin=Rmin+abs(Residue);
   p(i,j)=p(i,j)+Residue;
end
end
Rmin=Rmin/(Nx*Ny); % Average Residue per grid point
if(Rmin>=tol)
    Ncount=Ncount+1;
    if(Ncount>M)
        loop=0;
    end
else
    loop=0;
end
end
end
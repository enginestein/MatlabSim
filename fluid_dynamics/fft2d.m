function differentiated_data=fft2d(data_spacedomain,Nx,Ny,dx,dy)

x = 0 : dx : (Nx-1)*dx; % x-distance
y = 0 : dy : (Ny-1)*dy; % y-distance
Nyq_kx = 1/(2*dx); % Nyquist of data in first dimension
Nyq_ky = 1/(2*dy); % Nyquist of data in second dimension
dkx = 1/(Nx*dx);   % x-Wavenumber increment
dky = 1/(Ny*dy);   % y-Wavenumber increment
kx = -Nyq_kx : dkx : Nyq_kx-dkx; % x-wavenumber
ky = -Nyq_ky : dky : Nyq_ky-dky; % y-wavenumber
data_wavenumberdomain = zeros(size(data_spacedomain)); % initialize data
% Compute 2D Discrete Fourier Transform
for i1 = 1 : Nx
for j1 = 1 : Ny
for i2 = 1 : Nx
for j2 = 1 : Ny
data_wavenumberdomain(j1,i1) = data_wavenumberdomain(j1,i1) + ...
(2i*pi*ky(j1))*(2i*pi*kx(i1))* ...
(data_spacedomain(j2,i2)*exp(-1i*(2*pi)*(kx(i1)*x(i2)+ky(j1)*y(j2))));
end
end
end
end
differentiated_data = ifft2(ifftshift(data_wavenumberdomain));
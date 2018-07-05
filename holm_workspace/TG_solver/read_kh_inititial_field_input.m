% Purpose: read the kh_initial_field.in that I had been using for both my
% JFM papers. I got the input file from Ali. 
% Note: 
% It seems that I cannot re-produce the same eigenfunctions for density.

fname='kh_initial_field.in';
A = zeros(350,4);
fid=fopen(fname);
fgetl(fid)
for i=1:350;
    A(i,:) = str2num(fgetl(fid));
end
z=A(:,1);
psi = A(:,2);
dpsi= A(:,3);
rho = A(:,4);

nx = 100;
ij = sqrt(-1);
k_fgm = 0.441;
Lx = 2*pi/k_fgm;
x=linspace(-Lx/2,Lx/2,nx);
rho_fgm=real(rho*exp(ij*k_fgm*x));
[XX,ZZ]=meshgrid(x,z);
figure;
contourf(XX,ZZ,rho_fgm)
axis equal;

rho_tot = 0.01*rho_fgm + (1-tanh(z))*ones(1,nx);
contourf(XX,ZZ,rho_tot)
axis equal;
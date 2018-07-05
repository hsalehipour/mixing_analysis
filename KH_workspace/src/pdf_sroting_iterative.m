clear all;
nz= 500;
zmax=15;    zmin=-15;
z=linspace(zmin,zmax,nz);
rho=1-tanh(z);


% Sort using bin pdf approach
nbin_max=50;
rmax=max(rho);
rmin=min(rho);
dr=rmax-rmin;

   
      
%=====STEP 1
zb = zeros(1,nbin_max+1);
rhob = zeros(1,nbin_max+1);
for i=1:nz
    dzi = z(2)-z(1);
    rho_index = 2 + nbin_max*(rmax-rho(i))/dr;  % Bin the density
    j = floor(rho_index);
    j = min(j,nbin_max);
    zb(j) = zb(j) + dzi;
end


k=1;
for i=1:nbin_max+1      % Running sum of z - eliminate voids
if (zb(i)>0) 
    k=k+1;
    zb(k)=zb(k-1)+zb(i);
end
end
nbin = k;
rhob(1)=rmax;
rhob(2:nbin+1) = rmax-([2:nbin+1]*dr)/(nbin+1);

ratio = abs(zb(2:nbin)-zb(1:nbin-1))/(zmax-zmin);

set_ratio = 0.1;
find(ratio>set_ratio)
zb1 = zeros(1,nbin_max+1);
dr = rhob(2)-rhob(1);
for i=1:nz
    dzi = z(2)-z(1);
    rho_index = 2 + nbin_max*(rmax-rho(i))/dr;  % Bin the density
    j = floor(rho_index);
    j = min(j,nbin_max);
    zb1(j) = zb1(j) + dzi;
end

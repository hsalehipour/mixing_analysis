clear all;
nz= 500;
z=linspace(-15,15,nz);
rho=1-tanh(z);


% Sort using bin pdf approach
nbin_max=50;
rmax=max(rho);
rmin=min(rho);
dr=rmax-rmin;



%=====STEP 1
zb = zeros(1,nbin_max);
rhob = zeros(1,nbin_max);
for i=1:nz
    dzi = z(2)-z(1);
    rho_index = 2 + nbin_max*(rmax-rho(i))/dr;  % Bin the density
    j = floor(rho_index);
    j = min(j,nbin_max);
    zb(j) = zb(j) + dzi;
end

zbb=zeros(1,nbin_max+1);
for i=2:nbin_max+1      % Running sum of z - eliminate voids
    zbb(i)=zbb(i-1)+zb(i-1);
end

rhob(1)=rmax;
rhob(2:nbin_max) = rmax-([2:nbin_max]*dr)/nbin_max;



%=====STEP 2
r2 = 1.7;%rhob(2);
r1 = 0.3;%rhob(nbins-1);
nbin_max = 3*nbin_max-2;
zb = zeros(1,nbin_max);
rhob = zeros(1,nbin_max);

for i=1:nz
    
    %Map \rho(theta) ===> \theta(rho) using Chebyshev transf
    rr = rho(i);
    if     (rr>=r2 && rr <= rmax)
        theta = acos( (rr-r2)/(rmax-r2) );
    elseif (rr>r1 && rr < r2)
        theta = pi - asin( 2.0*(rr-r1)/(r2-r1)-1 );
    else
        theta = 2*pi -acos( (r1-rr)/(r1-rmin) );
    end
    
    rho_index = 2 + nbin_max*(1.0-theta/2.0/pi);  % Bin the density
    j = floor(rho_index);
    j = min(j,nbin_max);
    zb(j) = zb(j) + dzi;
end
         
k=1;
for i=1:nbin_max      % Running sum of z - eliminate voids
if (zb(i)>0) 
    k=k+1;
    zb(k)=zb(k-1)+zb(i);
end
end
nbins = k;



% zone 1
N = floor(nbins/3);
theta=linspace(0,pi/2,N);
rho1=rmax*cos(theta)+r2*(1-cos(theta));


% zone 2
N = nbins-2*floor(nbins/3);
theta=linspace(pi/2,3*pi/2,N);
rho2 = 0.5*(r1*(1-sin(theta))+r2*(1+sin(theta)));


% zone 3
N = floor(nbins/3);
theta=linspace(3*pi/2,2*pi,N);
rho3=r1*(1-cos(theta)) + rmin*cos(theta);
rhob = [rho1 rho2 rho3];



plot(z,rho,'r'); hold on;
plot(zb(1:nbins)-15,rhob(1:nbins),'o')

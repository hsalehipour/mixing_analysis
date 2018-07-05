clear all; clc; close all;

Re  = 48000;
Pr  = 1;
J   = 1.2;%0.01:0.02:0.32;
% Ri0 = 0:0.05:0.75;
% R = Ri0/J;
R   = 5;% sqrt(Pr);



k = 2.11;%2:0.01:2.2;
%k = 1.155:0.001:1.165;         % streamwise wavenumber
l = 0;                  % spanwise   wavenumber

% control parameters
iBC1 = [0 1];   % frictionless, insulated lower wall
iBCN = iBC1;
imode = 0;      % 0=all modes, 1=FGM 

% vertical profiles
nz = 500;
z = linspace(-5,5,nz)';
U = tanh(z);
rho = 1-tanh(R*z);
% a = 0.5;
% rho = 1-tanh(R*(z-a));


% initialization
nu      = 1/Re;
kappa   = nu/Pr;

if imode==0
   tic
   [sig,w,r] = my_SSF(z,U,rho,k,l,nu,kappa,iBC1,iBCN,imode,J);
   toc
else    
sig = zeros(length(k),length(J));
for ir=1:length(J);
    for ik=1:length(k);
        tic
        sig(ik,ir) = fast_SSF(z,U,rho,k(ik),l,nu,kappa,iBC1,iBCN,imode,J(ir));
        toc
    end
end
end


% superposition of both left-propagating mode and right-propagating modes
% (See Smyth et al 1988) for more info.
% fld = r(:,1)+r(:,2); 
% % 
% % % plotting commands
% figure;
% nx = 100;
% ij = sqrt(-1);
% k_fgm = k;
% Lx = 2*pi/k_fgm;
% x=linspace(-Lx/2,Lx/2,nx);
% fld_fgm=real(fld*exp(ij*k_fgm*x));
% [XX,ZZ]=meshgrid(x,z);
% contourf(XX,ZZ,fld_fgm)
% axis equal;
% 
% fld = cos(2*pi*XX/Lx).*sech(ZZ).*tanh(ZZ);
% figure;
% contourf(XX,ZZ,fld)
% axis equal;

% rho_fgm = fld_fgm;
% rho_tot = 0.01*rho_fgm + rho*ones(1,nx);
% contourf(XX,ZZ,rho_tot)
% axis equal;

% 
% [KK, RR] = meshgrid(k,R);
% contourf(KK,RR,sig');

% command for outputing the initial_field.in.
% output_initial_eigv(fname, k_fgm,z,w(:,1)+w(:,2),r(:,1)+r(:,2));

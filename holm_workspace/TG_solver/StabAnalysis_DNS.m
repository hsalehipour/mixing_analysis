clear all; clc; close all;

% load DNS profile data which should include Re, Pr, J and k=2*pi/Lx;
load('DNS_profile.mat');
Re  = stab.Re;
Pr  = stab.Pr;
J   = stab.Rib;
k   = stab.k;
% k   = 0.9:0.05:2;%2:0.1:1.1;%stab.k;
% k = 1.4;

% input parameters
nz = 500;
z = linspace(-5,5,nz)';

% read the profiles
U   = interp1(stab.z,stab.u  ,z,'pchip');
rho = interp1(stab.z,stab.rho,z,'pchip');


% control parameters
iBC1 = [0 1];   % frictionless, insulated lower wall
iBCN = iBC1;
imode = 0;      % 0=all modes, 1=FGM 
l = 0;          % spanwise   wavenumber

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

clear all; clc; close all;

Re  = 410;
Pr  = 700;
J   = 0.2;
R   = 5;
a   = -0.25;      % offset


k = 0.4:0.1:1;%1:25;
%k = 1.155:0.001:1.165;         % streamwise wavenumber
l = 0;                  % spanwise   wavenumber
% k=1.051;
% k=1.41;
% For Re=410, Pr=700, J=0.2, R=5, a=-0.3:  k=0.668 (FGM), k=0.757 (2nd FGM)

% control parameters
iBC1 = [0 1];   % frictionless (rigid is 1), insulated walls
iBCN = iBC1;
imode = 0;      % 0=all modes, 1=FGM 

% vertical profiles
nz = 200;
z = linspace(-5,5,nz)';
U = -tanh(z);   %-tanh(z);   % -sin(pi/2*z);
rho = 1-tanh(R*(z-a));



% initialization
nu      = 1/Re;
kappa   = nu/Pr;

if imode==0
    sig = zeros(length(k),2);
    for ik=1:length(k);
        tic
        [sig_all,w,r] = my_SSF(z,U,rho,k(ik),l,nu,kappa,iBC1,iBCN,imode,J);
        sig(ik,:) = sig_all(1:2);
        toc
    end  
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
% fld = r(:,1);% +r(:,2); 
% % 
% % % plotting commands
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



% 2) Analysis
ntime = length(time);
nz = length(z);

%% 2-1) RHS vs LHS of sigma2D
sig_LHS = ddx(KE,time);
sig_RHS = sigma.*KE*2;
err = abs(sig_RHS-sig_LHS)./sig_LHS*100;



%% 2-2) RHS vs LHS of sigma3D 
sig3_LHS = ddx(KE3d,time);
sig3_RHS = sigma3d.*KE3d*2;
err3 = abs(sig3_RHS-sig3_LHS)./sig3_LHS*100;


%% 2-3) Resolution comparison vs. Kolmogroph scale
Lk = (nu^3./eps_klmg).^0.25;
Lb = Lk/sqrt(Pr);
dx = Lx/nelx/norder;
dy = Ly/nely/norder;
dz = Lz/nelz/norder; 
rsln = max([dx./Lb, dy./Lb, dz./Lb],[],2);


%% 2-4) Calculate N2, Shear and mean values within shear layer
gradrho = ddx(rhobar,z);
shear   = ddx(ubar,z);
N2      = -g/rho0*gradrho;


% find shear layer for future averaging (not needed since in we do <f>/<g>)    
% Selected tol = 0.001*(initial density gradient across shear layer)
tol=1e-3;   
shl_zmax  = zeros(ntime,1);
shl_zmin  = zeros(ntime,1);
for i=1:ntime
    indz = abs(gradrho(:,i))>tol;
%     indz = abs(shear(:,i))>tol;
%     indz = z>=-1 & z<=1;
    shl_zmax(i) = max(z(indz),[],1)';
    shl_zmin(i) = min(z(indz),[],1)';
%     shl_ind(i,:)  = (z<=shl_zmax(i) & z>=shl_zmin(i));
end
shl_thk = shl_zmax - shl_zmin;

N2avg        = mean(N2,1)';
shear_avg    = mean(shear,1)';
gradrho_avg  = -N2avg*rho0/g;

% Calculate avg Rich. number, Ri = <N2>/<S2>
Ri_avg = N2avg./mean(shear.^2)';

% Calculate a bulk measure of Rich. number, Ri = <N2>/<S>^2
% Note: this definition is Lz dependent, so averaging should be done over
% the shear layer only.
Ri_bul = N2avg./(mean(shear,1).^2)'.*shl_thk/Lz;

%% 2-5) calculate Buoyancy freq
% NOTE 1: Re_b = dissp/nu/N2
% NOTE 2: in D/N2avg/nu, since D is averaged over the whole domain, N2avg
%         should also be averaged similarly (i.e. mean(N2,1) 
%         not trapz()/shl_thk )

% Find Reb
D = sig_Dv.*KE*2;
Reb = abs(D)./N2avg/nu;

% Reb_profile = epsbar./N2/nu;
% for i=1:ntime
%     shl_ind  = z<=shl_zmax(i) & z>=shl_zmin(i);
%     Reb(i) = mean(Reb_profile(shl_ind,i));
% end



%% 2-6) t2d and t3d
[~, indt2d] = max(KE2d);
[~, indt3d] = max(KE3d);
t2d = time(indt2d);
t3d = time(indt3d);

% find the start the turbulent phase
% check to ensure KE3d > KE2d for all t>tfo (added July 27,2014)
% MUST HAVE to avoid 2d-like turbulence into the turb-phase-averaged data
for ii=1:ntime-1;
    if(KE3d(ii)<KE2d(ii)); 
        indtfo=ii+1;     
    end;
end
tfo = time(indtfo);
% BUT FOR REB plots, the definition below gives higher Reb in low-Ri cases 
% BUT gives wrong values of mixing regime parameter
% ttemp = time(KE3d>=KE2d);
% tfo = ttemp(1);
% [~,indtfo] = min(abs(time-tfo));

% ensure t3d>tfo
if (t3d<tfo);
    [~, ind] = max(KE3d(time>=tfo));
    indt3d = ind  + indtfo;
    t3d = time(indt3d);
end

% find the end of turbulent phase or ...
% time associated with the start of relaminarization phase
cr = Reb;
tol = 1;
[~,ind] = min(abs(abs(cr(time>t3d))-tol));
indtrl = ind + indt3d;
trl = time(indtrl);
indtfe = indtrl;
tfe = trl;

%% 2-7) Potential Energy calculation
zz = z*ones(1,ntime);
RPE = Ri/R*trapz(z,rhob.*zz,1)'/Lz;
PE  = Ri/R*trapz(z,rhobar.*zz,1)'/Lz;
Dp = Ri*2/Lz/R/Re/Pr;
Dpt = Dp*time;
APE = PE-RPE;




%% 2-8) Time derivative calculations 
%ttime=linspace(time(1),time(end),ntime*100);
%PE2 = spline(time,PE,ttime);
%PEb2 = spline(time,PEb,ttime);
ddt_PE  = ddx(PE,  time);
ddt_RPE = ddx(RPE, time);




%% 2-9) Mixing calculation
M = ddt_RPE - Dp;
D = sig_Dv.*KE*2;
H = sig_Hb.*KE*2;
S = H - M;

% Mixing efficiency 
eff = M./(M-D);
%eff_c = trapz(time,M)./trapz(time,M-D);
eff_c = zeros(1,ntime);
eff_i = zeros(1,ntime);     % initial growth 
eff_f = zeros(1,ntime);     % fully turbulent phase
eff_d = zeros(1,ntime);     % decaying turbulence

for i=2:ntime
     eff_c(i)  = trapz(time(1:i),M(1:i))/...
                 trapz(time(1:i),M(1:i)-D(1:i));
end
for i=2:indtfo
    ind = 1;
    eff_i(i) = trapz(time(ind:i),M(ind:i))/...
               trapz(time(ind:i),M(ind:i)-D(ind:i));
end
for i=indt3d+1:indtfe
    ind = indt3d;    
    eff_f(i) = trapz(time(ind:i),M(ind:i))/...
               trapz(time(ind:i),M(ind:i)-D(ind:i));
end
for i=indtfe+1:ntime
    ind = indtfe;
    eff_d(i) = trapz(time(ind:i),M(ind:i))/...
               trapz(time(ind:i),M(ind:i)-D(ind:i));
end

ind = eff_i~=0;         eff_i_avg(ic) = mean(eff_i(ind));
ind = eff_f~=0;         eff_f_avg(ic) = mean(eff_f(ind));
ind = eff_d~=0;         eff_d_avg(ic) = mean(eff_d(ind));
ind = eff_c~=0;         eff_c_avg(ic) = mean(eff_c(ind));

ind = time>=t3d & time<=tfe;
eff_avg (ic) = mean(eff(ind));
eff_err(ic) = std(eff(ind));

eff_f_avg (ic) = mean(eff_f(ind));
eff_f_err(ic) = std(eff_f(ind));


%% 2-10) Calculate approxmiate mixing
fac = 2*KE3d;
H3d  = fac.*sig3_Hb;
Sh2d = fac.*sig3_Sh;
Shbar= fac.*sig3_Rs;
D3d  = fac.*sig3_Dv;
An   = fac.*sig3_An;

% viscous dissipation
Dbar = -trapz(z,shear.^2./Re)'./Lz;
D2d = D - D3d-Dbar;

% Rates in different reservoirs
dKE2d  = ddx(KE2d,time);
dKEbar = ddx(KEbar,time);


% New bouyancy fluxes
H2d = H - H3d;
XK = dKEbar + Shbar - Dbar;
XP = H2d - S;

indh = H3d>=0;
UP = XP*0;
UP(indh) = H2d(indh)-S(indh);
UP(~indh) = H(~indh)-S(~indh);
XP = UP;

% find beta using nl-least square fit: XP = beta*An
% modelfun = @(b,x)(b*x);
% x1 = An;    x=x1*0;
% y1 = XP;    y=y1*0;
% for i=2:ntime
%      x(i)  = trapz(time(1:i),x1(1:i));
%      y(i)  = trapz(time(1:i),y1(1:i));
% end
% beta0 = 1;
% beta(ic) = nlinfit(x,y,modelfun,beta0);
% beta = beta';

% find beta using nl-least square fit: XP = beta*An
modelfun = @(b,x)(b(1)*x(:,1)+b(2)*x(:,2));
x = zeros(ntime,2);     y = zeros(ntime,1);
x1 = An;
x2 = D2d;
y1 = XP;    
for i=2:ntime
     x(i,1)  = trapz(time(1:i),x1(1:i));
     x(i,2)  = trapz(time(1:i),x2(1:i));
     y(i)    = trapz(time(1:i),y1(1:i));
end
% x = [An D2d];
% y = XP;
beta0 = [0 0];
beta(:,ic) = nlinfit(x,y,modelfun,beta0);


%% 2-11 mixing recipe
ind = time>=t3d & time<=tfe;

Ri_turb(ic) = mean(Ri_avg(ind));

intH3d(ic)  = trapz(time(ind),H3d(ind));
intAn(ic)   = trapz(time(ind),An(ind));
intXP(ic)   = trapz(time(ind),XP(ind));
intM(ic)    = trapz(time(ind),M(ind));
intD(ic)    =-trapz(time(ind),D(ind));
intS(ic)    = trapz(time(ind),S(ind));
intDp(ic)   = trapz(time(ind),time(ind)*0+Dp);
intXK(ic)   = trapz(time(ind),XK(ind));
intShbar(ic) = trapz(time(ind),Shbar(ind));
intH2d(ic)   = trapz(time(ind),H2d(ind));
intSh2d(ic)  = trapz(time(ind),Sh2d(ind));
intD3d(ic)   = -trapz(time(ind),D3d(ind));
intD2d(ic)   = -trapz(time(ind),D2d(ind));

eff_cum = intM./(intM+intD);
gamma = eff./(1-eff);

%% 2-12 Cross-correlation
ind = time>=t3d & time<=tfe;
[ccXPAn  ,lags] = xcorr(XP,  smooth(An,8)  ,'coeff');        % < M, Hb_3d >
% [ccXPAn  ,lags] = xcorr(XP(ind),  An(ind)  ,'coeff');        % < M, Hb_3d >
% [ccXPAn  ,lags] = xcorr(XP,  An  ,'coeff');        % < M, Hb_3d >



%% 2-13) Kappa
% (1) Diascalar kappa : the most precise kappa (Winters & D'Asaro 1996)
% (2) Osborn-Cox formula
% (3) Conventional definition using turbulent buoyancy flux
% (4) Osborn formula based on mixing efficiency


% (1) kappa = kappa_0* <|grad rho|^2>S* / (ddz_rho*)^2
% where (*) refers to the adiabatically restratified surface
% NEW (July 29, 2014)
% In principle M = g phi_d , 
% where              phi_d(z) = <|grad rho|^2>S*/(ddz_rho*)
% But we noticed that my calculation of <|grad rho|^2>S* is not the most 
% accurate, and thus the above equality is not well satisfied for lower Ri
% So I simply use phi_d = M/g;
% NOTE: kappa0 = 1/(Re*Pr);
% at t=0, M=Dp (or as in Winters 1995, phi_d = phi_i)
gradrhob = ddx(rhob,z);
phi_d = -M./Ri;
kappa_star = phi_d./mean(gradrhob)'/kappa0;
% kappa_star = mean(chi_ref)./mean(gradrhob.^2);


% (2) kappa = kappa_0* <|grad rho'|^2> / (ddz_rhobar)^2
% gradrho  = ddx(rhobar,z);
% kappa_cox  = mean(chi3d)'  ./mean(gradrho.^2)';
% Note: chi_ref in the following formula, 
% kappa_coxfull  = mean(chi_pr)'./mean(gradrho.^2)';


% (3) kappa = <rho'w'>/ddz_rhobar;
% NOTE: kappa0=nu/Pr (where nu = 1/Re) is absolutely important to be used
% here in order to non-dimensionalize our derived kappa and kappa_3d based
% on buoyancy flu formula just to be consistent with our own general
% non-dimensionalization.
% IN GENERAL : kappa_nodim = kappa_dim (m^2/s)/kappa_ref
rhow    = 2*R/Ri*sig_Hb.*KE;
rhow3d  = 2*R/Ri*sig3_Hb.*KE3d;
kappa    = -rhow./gradrho_avg/kappa0;
kappa_3d = -rhow3d./gradrho_avg/kappa0;

% (4) kappa_osborn = 0.2*eps/N2
kappa_osb = 0.2*Pr*Reb;

%% 2-14) Turbulent PR

% irreversible turbulent production (needs a lot of elaboration, see JFM#3)
PP = M + abs(D);

% Calculate kappa_m (non-dimensional) and turbulent Pr
% kappa_m = -uw./shear_avg/nu;
kappa_m = PP./mean(shear.^2)'/nu;

% Pr_t = (kappa_m/nu0)*nu / (kappa_star*kappa0)* kappa 
% Pr_t = kappa_m/kappa_star*Pr
Pr_t = kappa_m./kappa_star*Pr;
ind = time>=t3d & time<=trl;
Prt_avg(ic) = mean(Pr_t(ind));
 % Kundu pg 587, flux Ri, Ri_g = Pr_t.*Rf

 
 
%% 2-15) Taylor microscale Reynolds number
% Based on Smyth 2001 (Eq. 29)
q = sqrt(2*(KE2d+KE3d));
% q = sqrt(2*KE3d);
lam = sqrt(5*nu*q.^2./abs(D));
Re_taylor = q.*lam./nu;
ind = time>=t3d & time<=trl;
Re_taylor0(ic) = Re_taylor(indt3d);


%% 2-15) Length scales
% NOTE: both nominator and denominator must have been averaged along the
%       same height (z) as long as its value is zero outside the shear
%       layer
% STILL SOME ISSUES with definition of N
% Note: To circumvent the problem one may follow Smyth & Moum (2000) in
%       defining integral length scales and define Nb but that
%       overestimates N.
%       So I can re-calculate the averaging of all parameters overt the
%       patch thickness and not the entire domain. Here the stratified
%       layer thickness is defined based on|drho/dz|>=1e-3; 
%       This is somewhat similar to the alternative approach of 
%       Smyth & Moum (2000, pg9)
% Note: 
%       Reb = (Lo/Lk)^4/3;
%       Ri  = (Lc/Lo)^4/3;
%       Lsc = (rho_r*U0^2)/ (4*g*rho0) = d/(4*Rib)
%       Rib = Ri_0/R;

% Define Turbulent Dissipation (i.e. D or D3d or D2d+D3d?)
turb_disp = abs(D);

% Integral length scales for shear and density layer thickness
Iu   = trapz(z,1-ubar.^2,1);
Irho = trapz(z,1-(rhobar-1).^2,1);
Lsc  = 1/4/g;

% Define a bulkd Richardson number based on the integral length scales
Ri_I = (g*2./Irho)./(2./Iu).^2;

% averaging factor to be pre-multiplied by domain-averaged values;
% avg_fac = shl_thk./Lz;
% 
% % 0) Re-scale the Thorpe length scale to represent the rms value within the
% % turbulent patch (i.e. not averaged over the entire computational domain).
% LT3d_patch = LT3d./sqrt(avg_fac);
% 
% % 1) Kolmogrov length scale, recalculated;
% Lk       = (nu^3./turb_disp).^0.25;
% Lk_patch = Lk.*avg_fac.^0.25;
% 
% % 2) Batchelor length scale
% Lb       = Lk/sqrt(Pr);
% Lb_patch = Lk_patch/sqrt(Pr);
% 
% % 3) Ozmidov length scale
% Lo       = (turb_disp./N2avg.^1.5).^0.5;
% Lo_patch = Lo.*avg_fac.^0.25;
% % N3avg = mean(N2.^1.5,1)';
% % Nb = sqrt(g*rho0./Irho');
% 
% % 4) Corssin length scale 
% % Note: in order for Ri_bul = (Lc/Lo)^4/3, Lc should be defined using
% % mean(shear)^3;
% Lc       = (turb_disp./shear_avg.^3).^0.5;
% Lc_patch = Lc.*avg_fac;
% % shear3_avg    = mean(shear.^3,1)';
% % Lc = (turb_disp./shear3_avg).^0.5;
% % Lc = (turb_disp./shear_avg.^3).^0.5;
% 
% % 5) Ellison length scale (Ivey Imberger 1991)
% % Added in Oct 2015: Le = <(rho')^2>_{xyz}^(0.5)/|drhobar/dz|;
% % Note:  
% % since the numerator of Le is first averaged and then taken to power 0.5,
% % the result would be Lz dependent unless calculated over the depth of the
% % turbulent patch (parameterized here by shl_thk). Hence the factor appears
% % to resolve this issue.
% % Le       = rhop_rms./abs(gradrho_avg);
% % Le_patch = Le.*avg_fac.^0.5;
% 
% % 6) Energy-containing length scales (see pg 9-10, Smyth_etal_2001)
% q2 = 2/3*KE3d;
% Len = q2.^(1.5)./turb_disp;






%% 2-16) calculate different non-dimensional numbers
% % NOTE 1: Rf =~sig3_Hb./sig3_Rs because sig3_Rs = <uw dU/dz> =~ uw*dU/dz
% 
% Ri_profile = N2./shear.^2;
% % Rf   = (-g*rhow3d)./(uw.*shear_avg);       
% % above expression is wrong b/c double division by Lz 
% %uwshear = -2*sig3_Rs.*KE3d;
% uwshear = -2*KE3d.*(sig3_Rs+sig3_Sh);
% Rf   = (+g*rhow3d)./(-uwshear);      

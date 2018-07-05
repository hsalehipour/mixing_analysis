%% Read Bills data
% Load data for the 2 cruises
load('../data_Bill/events_TIWE');       A_TIWE=A;
load('../data_Bill/events_FLX91');      A_FLX91=A;

% We'll combine the TIWE and FLX91 datasets, since I assume the difference
% doesn't matter here.
Lt=[A_TIWE.Lt; A_FLX91.Lt]; % Thorpe scale
Gd=[A_TIWE.Gd; A_FLX91.Gd]; % Gamma_d
E=Gd./(1+Gd);

eps=[A_TIWE.eps; A_FLX91.eps]; % epsilon
chi=[A_TIWE.chi; A_FLX91.chi]; % chi

% Adjust N^2 as in Smyth et al (2001)
r=[A_TIWE.Lt; A_FLX91.Lt]./[A_TIWE.Le; A_FLX91.Le]; 
% To override this adjustment, uncomment the next line
%r=ones(size(r));

% Find Lo and buoy Reynolds
Lo=[A_TIWE.Lo; A_FLX91.Lo].*(r.^.75);% adjusted Ozmidov scale 
Rb=[A_TIWE.Rb; A_FLX91.Rb].*r;% adjusted buoyancy Reynolds number 
R_OT=Lo./Lt;


% FLX91 old patches    
ind_old = A_FLX91.Lt<=A_FLX91.Lo;
FLX91o_Rb=A_FLX91.Rb(ind_old);
FLX91o_Gd=A_FLX91.Gd(ind_old); 
FLX91o_Jd=A_FLX91.Gd(ind_old).*A_FLX91.eps(ind_old); 
% FLX91o_Kap = nu*FLX91o_Gd.*FLX91o_Rb;

% FLX91 young patches    
FLX91y_Rb=A_FLX91.Rb(~ind_old);
FLX91y_Gd=A_FLX91.Gd(~ind_old); 
FLX91y_Jd=A_FLX91.Gd(~ind_old).*A_FLX91.eps(~ind_old); 
% FLX91y_Kap = nu*FLX91y_Gd.*FLX91y_Rb;

% TIWE old patches
ind_old = A_TIWE.Lt<=A_TIWE.Lo;
TIWEo_Rb=A_TIWE.Rb(ind_old); 
TIWEo_Gd=A_TIWE.Gd(ind_old); 
TIWEo_Jd=A_TIWE.Gd(ind_old).*A_TIWE.eps(ind_old); 
% TIWEo_Kap = nu*TIWEo_Gd.*TIWEo_Rb;

% TIWE young patches
TIWEy_Rb=A_TIWE.Rb(~ind_old); 
TIWEy_Gd=A_TIWE.Gd(~ind_old); 
TIWEy_Jd=A_TIWE.Gd(~ind_old).*A_TIWE.eps(~ind_old);
% TIWEy_Kap = nu*TIWEy_Gd.*TIWEy_Rb;

% nbins is arbitrary. The number of patches in each bin appear effectively
% in the standard deviation
bin_min = 1;
bin_max = 1e5; 
nbins   = 50;  
bins=log10(bin_min):(log10(bin_max)-log10(bin_min))/nbins:log10(bin_max);
bins=10.^bins;  

% Bin Gamma data : old patches
[FLX91o_bined_Gam, FLX91o_npatch, FLX91o_sdev] = ...
    bindata(FLX91o_Rb,FLX91o_Gd,bins);
[TIWEo_bined_Gam, TIWEo_npatch, TIWEo_sdev] = ...
    bindata(TIWEo_Rb,TIWEo_Gd,bins);
FLX91o_err = FLX91o_sdev./sqrt(FLX91o_npatch);
TIWEo_err  =  TIWEo_sdev./sqrt( TIWEo_npatch);

% Bin Gamma data : young patches
[FLX91y_bined_Gam, FLX91y_npatch, FLX91y_sdev] = ...
    bindata(FLX91y_Rb,FLX91y_Gd,bins);
[TIWEy_bined_Gam, TIWEy_npatch, TIWEy_sdev] = ...
    bindata(TIWEy_Rb,TIWEy_Gd,bins);
FLX91y_err = FLX91y_sdev./sqrt(FLX91y_npatch);
TIWEy_err  =  TIWEy_sdev./sqrt( TIWEy_npatch);

%% Plot scattered data of Bill's data
figure(1);
oldcolor   = [1 0.5 0.5];
youngcolor = [0.5 0.5 1];

loglog(Rb(R_OT>=1), Gd(R_OT>=1), '.', 'Color',oldcolor  , 'MarkerSize',18); 
hold all;
loglog(Rb(R_OT<1) , Gd(R_OT<1) , '.', 'Color',youngcolor, 'MarkerSize',18);
xlabel('$Re_b$'  , 'interpreter','latex')
ylabel('$\Gamma$', 'interpreter','latex')


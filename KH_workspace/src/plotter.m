% % plotter (A.M. modified by H.S.)
% Date: 9 Dec 2015
% Objective: This function plots Gamma vs. Reb for KH DNS data of HS in
% bins and overlays on its top some observed measurments
%
%
%
%% Startup
% Figure 1  : will plot Gamma-Reb for all DNS (binned), and scattered Bill.
%             Both DNS and observational data are classified based on "young" 
%             vs. "mature" distinction.
% Figrue 2,3  will plot histograms of Reb and Gamma of Bill's data
% Figure 4    will plot Gamma-Reb for all "mature" DNS (binned), Bill's
%             "mature" data, BB curve, and the two parameterization curve.

close all; clear all;
addpath(genpath('/usr/local/MATLAB/matlab_extra'));
% set(0, 'DefaultAxesFontName', 'Nimbus Roman No9 L')

errlim = 5000;

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

% % Bin Kappa data
% [FLX91_bined_Kap, FLX91_npatch, FLX91_sdev] = ...
%     bindata(FLX91_Rb,FLX91_Kap,bins);
% [TIWE_bined_Kap, TIWE_npatch, TIWE_sdev] = ...
%     bindata(TIWE_Rb,TIWE_Kap,bins);
% FLX91_err = FLX91_sdev./sqrt(FLX91_npatch);
% TIWE_err  =  TIWE_sdev./sqrt( TIWE_npatch);

%% Read DNS data of of Salehipour & Peltier(2015) and Salehipour(GRL,2016)
load('../data/DNS_KHI');

bin_min = 1;
bin_max = 1e5; 
nbins   = 125;  
dnsbins=log10(bin_min):(log10(bin_max)-log10(bin_min))/nbins:log10(bin_max);
dnsbins=10.^dnsbins;  

% DNS mature (old) patches
DNSo_Reb = DNS.old.Reb;
DNSo_Gam = DNS.old.Gam;
[DNSo_Gam_bined, DNSo_npatch, DNSo_Gam_sdev] =...
    bindata(DNSo_Reb,DNSo_Gam,dnsbins);
DNSo_err  =  DNSo_Gam_sdev./sqrt(DNSo_npatch);

% DNS young patches
DNSy_Reb = DNS.young.Reb;
DNSy_Gam = DNS.young.Gam;
[DNSy_Gam_bined, DNSy_npatch, DNSy_Gam_sdev] =...
    bindata(DNSy_Reb,DNSy_Gam,dnsbins);
DNSy_err  =  DNSy_Gam_sdev./sqrt(DNSy_npatch);

%% Plot scattered data of Bill's data
figure(1);
oldcolor   = [1 0.5 0.5];
youngcolor = [0.5 0.5 1];

loglog(Rb(R_OT>=1), Gd(R_OT>=1), '.', 'Color',oldcolor  , 'MarkerSize',18); 
hold all;
loglog(Rb(R_OT<1) , Gd(R_OT<1) , '.', 'Color',youngcolor, 'MarkerSize',18);
xlabel('$Re_b$'  , 'interpreter','latex')
ylabel('$\Gamma$', 'interpreter','latex')


%% Plot 1d histograms Bill's data 
oldcolor   = [1 0.5 0.5];
youngcolor = [0.5 0.5 1];

figure (2); % histograms of Gamma
nbin = 75;
hndl1 = plotHist(Gd(R_OT>=1),nbin,oldcolor,oldcolor,1);
hndl2 = plotHist(Gd(R_OT<1),nbin,youngcolor,youngcolor,1);

figure (3); % histograms of Reb
hndl3 = plotHist(Rb(R_OT>=1),nbin,oldcolor,oldcolor,1);
hndl4 = plotHist(Rb(R_OT<1),nbin,youngcolor,youngcolor,1);


%% plot DNS data of the KH-ansatz (younge and old) bined based on Reb 

figure(1);

% DNS "mature" patches: Fill the shaded region of error bars
fcolor = [1 0 0];
ecolor = [1 0.5 0.5];
indfill = DNSo_err>=0 & DNSo_err<=errlim; 
h_dns_old = plot_filled_binned(dnsbins(indfill),DNSo_Gam_bined(indfill),...
    DNSo_err(indfill),ecolor,fcolor);

% DNS "mature" patches: Fill the shaded region of error bars
fcolor = [0 0 1];
ecolor = [0.5 0.5 1];
indfill = DNSy_err>=0 & DNSy_err<=errlim; 
h_dns_young = plot_filled_binned(dnsbins(indfill),DNSy_Gam_bined(indfill),...
    DNSy_err(indfill),ecolor,fcolor);
 
uistack(h_dns_old,'top');


%% Plot a 2d histogram of all Bill's data
% gam_bins = 10.^(linspace(-2,1,100));
% hist3([[TIWE_Rb; FLX91_Rb], [TIWE_Gd; FLX91_Gd]],'Edges',{bins gam_bins},'FaceAlpha',.65);
% 
% % [NN, CC] = hist3([[TIWE_Rb; FLX91_Rb], [TIWE_Gd; FLX91_Gd]],{bins gam_bins});
% % imagesc(CC{1},CC{2},NN'); 
% % axis xy;
% set(get(gca,'child'),'FaceColor','interp',...
%     'LineStyle', 'None',...
%     'CDataMode','auto');
% view(2)
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');
% colormap(brewermap(256,'Purples'))
% uistack(get(gca,'child'),'bottom');


%% parameterization
param = @(x,x0,y0)  2*y0*(x/x0).^(0.5)./(1+x/x0);
% (x0,y0) : the coordinates of eff-peak in the x-y space (x=Reb, y=Gamma);
% Gamma_param=2*Gammastar*(Reb_smooth/Rebstar).^.5./(1+Reb_smooth/Rebstar);







%% plot Bill data : FLX91
% fcolor = 'm';
% ecolor = [0.75 0 0.75];
fcolor = [0.68 0.92 1];
ecolor = [0.04 0.52 0.78];

% scatter plot
% loglog(FLX91_Rb, FLX91_Gd, 'r.', 'MarkerSize',18);

figure;
%=== old patches
% Fill the shaded region of error bars
% indfill = FLX91_err>=0 & FLX91_err<=errlim; 
indfill = FLX91o_npatch>=5; 
h_FLX91o = plot_filled_binned(...
    bins(indfill),...
    FLX91o_bined_Gam(indfill),...
    FLX91o_err(indfill),ecolor,fcolor);
         

%=== young patches
% Fill the shaded region of error bars
% indfill = FLX91_err>=0 & FLX91_err<=errlim; 
indfill = FLX91y_npatch>=5; 
h_FLX91y = plot_filled_binned(...
    bins(indfill),...
    FLX91y_bined_Gam(indfill),...
    FLX91y_err(indfill),0.5*ecolor,0.5*fcolor);                 



%% plot Bill data : TIWE
fcolor = [0.5 1 0.5];
ecolor = [0.17 0.51 0.34];

% scatter plot
% loglog(TIWE_Rb, TIWE_Gd, 'k.', 'MarkerSize',18);

%=== old patches
% Fill the shaded region of error bars
% indfill = TIWE_err>=0 & TIWE_err<=errlim;    
indfill = TIWEo_npatch>=5; 
h_TIWEo = plot_filled_binned(...
    bins(indfill),...
    TIWEo_bined_Gam(indfill),...
    TIWEo_err(indfill),ecolor,fcolor);                 


%=== young patches
% Fill the shaded region of error bars
% indfill = TIWE_err>=0 & TIWE_err<=errlim;    
indfill = TIWEy_npatch>=5; 
h_TIWEy = plot_filled_binned(...
    bins(indfill),...
    TIWEy_bined_Gam(indfill),...
    TIWEy_err(indfill),0.5*ecolor,0.5*fcolor); 



%% Plot Bouffard data
fcolor = [0.8 0.8 0.8];
ecolor = [0.5 0.5 1];

% load('../data_Bouffard/Gamma_Reb.mat')
load('../data_Bouffard/Gamma_Reb_STD.mat')
loglog(BB_Reb,BB_Gamma,'-','Color',ecolor,'LineWidth',3); 

% % options for shaded region of errors
indfill = BB_Gamma>0; 
opts={'EdgeColor', ecolor,...       
      'FaceColor', fcolor};
[y1handle, y2handle, h_BB] = fill_between(BB_Reb(indfill), ...
    BB_GammaUP(indfill), BB_GammaDOWN(indfill), [], opts{:});  
set(y1handle,'Visible','off');
set(y2handle,'Visible','off');                   





%% plot the parameterizaqtion curve
Reb_smooth=[.1:.01:1 1:1:1e2 1e2:10:1e3 1e3:1e2:5e3 5e3:1e3:2e5];
Gamma_param_lb = param(Reb_smooth,100,0.2);
Gamma_param_ub = param(Reb_smooth,300,0.5);

h_parlb=loglog(Reb_smooth,Gamma_param_lb,'--k','linewidth',4);
h_parub=loglog(Reb_smooth,Gamma_param_ub,'--k','linewidth',4);
h_osb  = plot(Reb_smooth,Reb_smooth*0+0.2,'k:','linewidth',3);


%%
axis([3 1e5 0. 0.75]);
box on;
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'linear');
xlabel('$Re_b$'  , 'interpreter','latex')
ylabel('$\Gamma$', 'interpreter','latex')

legend([h_dns, h_FLX91o, h_TIWE, h_BB, h_osb, h_parub, h_parlb],...
    {'DNS [Salehipour \& Peltier (2015)]',...
     'Obs. FLX91 [Smyth et al. (2001)]',...
     'Obs. TIWE~ [Smyth et al. (2001)]',...
     'Obs. [Bouffard \& Boegman (2013)]', ...
     'Osborn (1980)',...
     'Upper bound param. $Re_b^{\star}=300$, $\Gamma~~^{\star}=0.5$',...
     'Lower bound param. $Re_b^{\star}=100$, $\Gamma~~^{\star}=0.2$',...
},...
    'interpreter','latex');

%exportfig(gcf, 'fig/tst.eps', 'FontMode', 'fixed','FontSize', 30,'width',20,'height',12);
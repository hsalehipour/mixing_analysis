%% Read Data
clear all;      
%close all;
fadrs = '/home/hesam/kh-workspace/Argo_data_CWhalen/av_epsilon_Oct_2014/';
%fname = [fadrs, 'av_world_K_250to500m_3_1p5deg_2014.mat'];
%fname = [fadrs, 'av_world_K_500to1000m_3_1p5deg_2014.mat'];
fname = [fadrs, 'av_world_K_1000to2000m_3_1p5deg_2014.mat'];
load(fname);
nu = 1e-6;

% Mapping toolbox settings
DataRef = [1/1.5 -90 180];           % Reference vector for mapping toolbox
projection = 'robinson';


%% Find the buoyancy Reynolds number from the Argo float data
lon=x;       lat=y;
K_osb = c;
Reb = K_osb/0.2/nu;

% Plotting input and colorbar setting
DataMin = 0;        DataMax = 4;
tick_range = [0 1 2 3 4];

% figure;
% [a,b,h] = geoplot(log10(Reb),DataRef, DataMin, DataMax, projection);
% set(h, 'YTick', (tick_range-b)/a);
% set(h, 'YTickLabel', tick_range);
% title(h, '$log_{10} Re_b$','interpreter', 'latex');

% figure;
% h = imgplot(lon,lat,log10(Reb), DataMin, DataMax);
% set(h, 'YTick', tick_range);
% title(h, '$log Re_b$','interpreter', 'latex');

% pcolor(x,y,log10(Reb)); shading flat; colorbar;
% ylim([-90 90]);
% axis equal;


%% Mixing efficiency
% fit1: a power law to eff-vs-Reb                   for Reb<=100
% fit2: a line to the  log10(eff) vs. log10(Reb)    for Reb>100
% (made sure the two fits are connected at Reb=100)
%   y1=a*x.^b+c;
%   y2=2*x.^(-0.5);
%=========Fit1
% General model Power2:
%      f(x) = a*x^b+c
% Coefficients (with 95% confidence bounds):
%        a =      -0.427  (-0.4585, -0.3954)
%        b =     -0.8887  (-0.9716, -0.8058)
%        c =      0.2039  (0.1987, 0.2091)
% 
% Goodness of fit:
%   SSE: 2.157
%   R-square: 0.7214
%   Adjusted R-square: 0.7211
%   RMSE: 0.03591
%=========Fit2
% Linear model Poly1:
%      f(x) = p1*x + p2
% Coefficients (with 95% confidence bounds):
%        p1 =     -0.5113  (-0.532, -0.4906)
%        p2 =      0.2958  (0.2429, 0.3487)
% 
% Goodness of fit:
%   SSE: 16.02
%   R-square: 0.7093
%   Adjusted R-square: 0.709
%   RMSE: 0.1291


% c1 =-0.427;      c2 =-0.8887;     c3 = 0.2071;
c1 =-0.4;      c2 =-0.9;     c3 = 0.2;      % adjusted

eff=Reb./0;
for i=1:length(y)
for j=1:length(x)
    if(Reb(i,j)<100)
        eff(i,j) = c1*Reb(i,j).^c2+c3;
    else
        eff(i,j) = 2*Reb(i,j).^(-0.5);
    end
end
end  


% Plotting input and colorbar setting
DataMin = 0;        DataMax = 0.25;

% tick_range = [0 0.1 0.2 0.3];
% figure;
% [a,b,h] = geoplot(eff./(1-eff), DataRef, DataMin, DataMax, projection);
% set(h, 'YTick', (tick_range-b)/a);
% set(h, 'YTickLabel', tick_range);
% title(h, '$\mathcal{E}/(1-\mathcal{E})$','interpreter', 'latex');

% figure;
% h = imgplot(lon,lat,eff./(1-eff), DataMin, DataMax);
% %title(h, '$\mathcal{E}$/(1-$\mathcal{E})$','interpreter', 'latex');
% ylabel(h, '$\mathcal{E}$/(1-$\mathcal{E})$','interpreter', 'latex','Rotation',0);


% figure;
% pcolor(x,y,eff./(1-eff)); shading flat; colorbar;
% ylim([-90 90]);
% axis equal;




%% Turbulent Pr
% General model: (fit for 0<Prt<20)
%      f(x) = a*x^b+c*x^d
% Coefficients (with 95% confidence bounds):
%        a =       167.4  (140.9, 193.8)
%        b =      -2.895  (-3.061, -2.728)
%        c =       3.597  (3.234, 3.96)
%        d =      -0.209  (-0.2336, -0.1844)
% 
% Goodness of fit:
%   SSE: 2502
%   R-square: 0.8598
%   Adjusted R-square: 0.8596
%   RMSE: 0.9899



%p1 = 1.0;    p2=10;   q1 =-2;
%Prt = (p1*Reb+p2)./(Reb + q1);

% Fit 1: a*x^b+c*x^d
a = 167.4;    b = -2.9;      c=3.6;        d=-0.2;
Prt = a*Reb.^b+c*Reb.^d;
Prt(Reb>10^3) = 1;

% Plotting input and colorbar setting
DataMin = 0;        DataMax = 5;
tick_range = 0:5; 

% figure;
% [a,b,h] = geoplot(Prt, DataRef, DataMin, DataMax, projection);
% set(h, 'YTick', (tick_range-b)/a);
% set(h, 'YTickLabel', tick_range);
% title(h, '$Pr_t$','interpreter', 'latex');

% figure;
% h = imgplot(lon,lat,Prt, DataMin, DataMax);
% set(h, 'YTick', tick_range);
% % title(h, '$Pr_t$','interpreter', 'latex');
% ylabel(h, '$Pr_t$','interpreter', 'latex','Rotation',0);

% % figure;
% % pcolor(x,y,Prt); shading flat; colorbar;
% % ylim([-90 90]);
% % axis equal;



%% Kappa based on New Eff
K_new = eff./(1-eff).*Reb*nu;

% % Plotting input and colorbar setting
% % DataMin = 0;        DataMax = 10;
% % tick_range = 0:5; 
% DataMin = -6;        DataMax = -3.5;
% tick_range = -6:-4;
% figure;
% %[a,b,h] = geoplot(K_new./K_osb, DataRef, DataMin, DataMax, projection);
% [a,b,h] = geoplot(log10(K_osb), DataRef, DataMin, DataMax, projection);
% set(h, 'YTick', (tick_range-b)/a);
% set(h, 'YTickLabel', tick_range);
% title(h, '$K^*_{\rho}/K_{\rho}$','interpreter', 'latex');

% figure;
% pcolor(x,y,log10(K_osb)); shading flat; colorbar;
% ylim([-90 90]);
% axis equal;


%% Comparison with Bryan and Lewis (1979) profile and POP1 profiles
% calculate the global averages
Kosb_avg = mean(K_osb(~isnan(K_osb)));
Knew_avg = mean(K_new(~isnan(K_new)));

z =-5000:0;

% Bryan and Lewis (1979)
kappa_BL   = 0.8 + 1.05/pi*atan((abs(z)-2500)/222.2);

% Large (1994) kpp_ccsm3:
bckgrnd_vdc1    = 0.524;
bckgrnd_vdc2    = 0.313;
bckgrnd_vdc_dpth= 1000.0e02;
bckgrnd_vdc_linv= 4.5e-05;
Prandtl         = 10.0;
rich_mix        = 50.0;
kappa_POP1 = bckgrnd_vdc1 + ...
             bckgrnd_vdc2*atan((abs(z)*100-bckgrnd_vdc_dpth)*bckgrnd_vdc_linv);

% Read netcdf files containing other vertical profiles of kappa
fadrs = '/home/hesam/kh-workspace/guido_kappa_profiles/';
fname1 = [fadrs, 'kappa_obs.nc'];
fname2 = [fadrs, 'kappa_ccsm4_globalavg.nc'];

% file 1
ncid = netcdf.open(fname1,'NC_NOWRITE');
depth1   = netcdf.getVar(ncid,0,'double')';     % depth levels 
depth2   = netcdf.getVar(ncid,1,'double')';     % depth levels
kappa_ls = netcdf.getVar(ncid,2,'double')';     % Lampkin & Speer (2007)
kappa_g  = netcdf.getVar(ncid,3,'double')';     % Ganachaud (2003
kappa_dl = netcdf.getVar(ncid,4,'double')';     % Decloedt & Luther (2012)
kappa_nf = netcdf.getVar(ncid,5,'double')';     % Nikurashin & Ferrari (2013) 
netcdf.close(ncid);


% file 2
ncid = netcdf.open(fname2,'NC_NOWRITE');
depth3     = netcdf.getVar(ncid,0,'double')';     % depth levels 
kappa_cesm = netcdf.getVar(ncid,1,'double')';     % ccsm4 global avg model
netcdf.close(ncid);

fac = 1e-4;    % conversion factor from cm^2/s to m^2/s
figure;
plot(kappa_BL  *fac,z,'--','LineWidth',2); hold all;
plot(kappa_POP1*fac,z,'-','LineWidth',2);
plot(kappa_cesm*fac,depth3,'-','LineWidth',2);
plot(kappa_g ,depth1,'--','LineWidth',4);
plot(kappa_ls(11:end)  ,depth1(11:end),'--','LineWidth',4);
plot(kappa_dl,depth1,'-.' ,'LineWidth',4);
%plot(kappa_nf,depth2,'-.' ,'LineWidth',4);

figure(1); hold all;          
plot([Kosb_avg Kosb_avg], -Zlim1,'g-o','LineWidth',2,'MarkerFace','g')
plot([Knew_avg Knew_avg], -Zlim1,'r-o','LineWidth',2,'MarkerFace','r')

xlim([5e-6 1e-3]);
ylim([-5000 0]);
set(gca, 'XScale', 'log');
grid on; box on;
xlabel('$K_{\rho} (m^2/s)$', 'interpreter','latex')
ylabel('$depth (m)$', 'interpreter','latex')
hndl = legend('BL79', ...
              'NCAR CCSM3',...
              'NCAR CESM1',...
              'G03',...
              'LS07',...
              'DL12',...
              'Argo-Osborn',...
              'Argo-DNS');
              %'NF13',...



% set(hndl, 'interpreter','latex'); 
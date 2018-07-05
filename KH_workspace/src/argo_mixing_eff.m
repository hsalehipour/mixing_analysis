%% Read Data
%clear all;     clc; 
%close all;
fadrs = '/home/hesam/kh-workspace/Argo_data_CWhalen/av_epsilon_Oct_2014/';

fname = [fadrs, 'av_world_K_250to500m_3_1p5deg_2014.mat'];    ic=1;
fname = [fadrs, 'av_world_K_500to1000m_3_1p5deg_2014.mat'];     ic=2;
fname = [fadrs, 'av_world_K_1000to2000m_3_1p5deg_2014.mat'];  ic=3;
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
% =======================
% A multiparameter Parameterization of mixing efficiency
% =======================
Ri= 0.25;
eff=Reb./0;
for i=1:length(y)
for j=1:length(x)
    eff(i,j) = mixeff_param(Reb(i,j),Ri);    
end
end  



% =======================
% Ri = 0:0.05:1;
% nri = length(Ri);
% Knew_avg = zeros(nri,1);
% 
% for kk=1:length(Ri);
% eff=Reb./0;
% for i=1:length(y)
% for j=1:length(x)
%     eff(i,j) = mixeff_param(Reb(i,j),Ri(kk));    
% end
% end  
% K_new = eff./(1-eff).*Reb*nu;
% Knew_avg(kk) = mean(K_new(~isnan(K_new)));
% end


% Plotting input and colorbar setting
% DataMin = 0;        DataMax = 0.5;

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


%% Calculate the global averages
if (ic==1) 
    Kosb_avg=[]; 
    Knew_avg=[]; 
    zlevel1 = [];
    zlevel2 = [];
end;
Kosb_avg = [Kosb_avg; mean_weighted(lat,K_osb)];
Knew_avg = [Knew_avg; mean_weighted(lat,K_new)];
zlevel1 = [zlevel1; -Zlim1'];
zlevel2 = [zlevel2; -mean(Zlim1)];

% Kosb_avg = mean(K_osb(~isnan(K_osb)));
% Knew_avg = mean(K_new(~isnan(K_new)));

Kappaosb_avg = repmat(Kosb_avg,[1 2])';
Kappanew_avg = repmat(Knew_avg,[1 2])';
if (ic==3);
    figure(1); hold all;
    plot(Kappaosb_avg(:), zlevel1,'ko-','LineWidth',2,'MarkerFace','r')
    plot(Kappanew_avg(:), zlevel1,'go--','LineWidth',2,'MarkerFace','g')
end
% plot(Kosb_avg, zlevel2,'go','LineWidth',2,'MarkerFace','g')
% plot(Knew_avg, zlevel2,'ro','LineWidth',2,'MarkerFace','g')
% 
% hndl = legend('Argo-Osborn',...
%               'Argo-DNS, $Ri=0.25$',...
%               'Argo-DNS, $Ri=0.5$',...
%               'Argo-DNS, $Ri=0.75$',...
%               'Argo-DNS, $Ri=0.95$');          
% set(hndl, 'interpreter','latex'); 
% gtext(['$Ri$ =', num2str(Ri)],'interpreter','latex');


%% Comparison with Bryan and Lewis (1979) profile and POP1 profiles

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
fname3 = [fadrs, 'weighted_argo_kappa_ccsm4_globalavg.nc'];
% 
% % file 1
ncid = netcdf.open(fname1,'NC_NOWRITE');
depth1   = netcdf.getVar(ncid,0,'double')';     % depth levels 
depth2   = netcdf.getVar(ncid,1,'double')';     % depth levels
kappa_ls = netcdf.getVar(ncid,2,'double')';     % Lampkin & Speer (2007)
kappa_g  = netcdf.getVar(ncid,3,'double')';     % Ganachaud (2003
kappa_dl = netcdf.getVar(ncid,4,'double')';     % Decloedt & Luther (2012)
kappa_nf = netcdf.getVar(ncid,5,'double')';     % Nikurashin & Ferrari (2013) 
netcdf.close(ncid);
% 
% 
% file 2
ncid = netcdf.open(fname2,'NC_NOWRITE');
depth3     = netcdf.getVar(ncid,0,'double')';     % depth levels 
kappa_cesm = netcdf.getVar(ncid,1,'double')';     % ccsm4 global avg model
netcdf.close(ncid);

% % file 3
% ncid = netcdf.open(fname3,'NC_NOWRITE');
% kappa_argo_cesm = netcdf.getVar(ncid,1,'double')';     % argo-loc ccsm4 global avg model
% netcdf.close(ncid);



% ===plot vertical profiles
if(ic==3);
fac = 1e-4;    % conversion factor from cm^2/s to m^2/s
figure(1);  hold all;
plot(kappa_BL  *fac,z,'--','LineWidth',2); 
plot(kappa_POP1*fac,z,'-','LineWidth',2);
plot(kappa_cesm*fac,depth3,'-','LineWidth',2);
%plot(kappa_argo_cesm*fac,depth3,'-','LineWidth',2);
plot(kappa_g ,depth1,'--','LineWidth',4);
plot(kappa_ls(11:end)  ,depth1(11:end),'--','LineWidth',4);
plot(kappa_dl,depth1,'-.' ,'LineWidth',4);
% %plot(kappa_nf,depth2,'-.' ,'LineWidth',4);
% 
xlim([5e-6 1e-3]);
ylim([-5000 0]);
set(gca, 'XScale', 'log');
grid on; box on;
xlabel('$K_{\rho} (m^2/s)$', 'interpreter','latex')
ylabel('$depth (m)$', 'interpreter','latex')
hndl = legend('Argo-Osborn',...
              'Argo-DNS',...
              'BL79', ...
              'NCAR CCSM3',...
              'NCAR CESM1',...
              'G03',...
              'LS07',...
              'DL10');
              %'NF13',...
end;
% set(hndl, 'interpreter','latex'); 



%% Write Argo location file in .nc 
% Purpose: to compute CESM global average only on the Argo locations 
% argo_loc = Reb;
% argo_loc(~isnan(Reb))=1;
% argo_loc(isnan(Reb))=0;
% nlat = length(lat);
% nlon = length(lon);
% 
% ncid = netcdf.create('argo_latlon.nc','NC_WRITE');
% x_dimid = netcdf.defDim(ncid,'nlat',nlat);
% y_dimid = netcdf.defDim(ncid,'nlon',nlon);
% dimid = [x_dimid, y_dimid];
% lat_varID = netcdf.defVar(ncid,'lat','double',x_dimid);
% lon_varID = netcdf.defVar(ncid,'lon','double',y_dimid);
% loc_varID = netcdf.defVar(ncid,'location','double',dimid);
% netcdf.endDef(ncid);
% % 
% netcdf.putVar(ncid,lat_varID,lat);
% netcdf.putVar(ncid,lon_varID,lon);
% netcdf.putVar(ncid,loc_varID,argo_loc);
% netcdf.close(ncid);

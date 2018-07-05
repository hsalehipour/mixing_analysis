

%% Plot Reynolds number from the Argo float data

% Plotting input and colorbar setting
% DataMin = 0;        DataMax = 4;
% tick_range = [0 1 2 3 4];
% figure;
% h = imgplot(lon,lat,log10(Reb(:,:,1)), DataMin, DataMax);
% % set(h, 'YTick', tick_range);
% title(h, '$log Re_b$','interpreter', 'latex');

% figure;
% [a,b,h] = geoplot(log10(Reb),DataRef, DataMin, DataMax, projection);
% set(h, 'YTick', (tick_range-b)/a);
% set(h, 'YTickLabel', tick_range);
% title(h, '$log_{10} Re_b$','interpreter', 'latex');

% pcolor(x,y,log10(Reb)); shading flat; colorbar;
% ylim([-90 90]);
% axis equal;



%% Mixing efficiency
% =======================


% Plotting input and colorbar setting
% DataMin = 0;        DataMax = 0.5;
% figure;
% h = imgplot(lon,lat,eff(:,:,1)./(1-eff(:,:,1)), DataMin, DataMax);
% title(h, '\Gamma=E/(1-E)');%,'interpreter', 'latex');

% tick_range = [0 0.1 0.2 0.3];
% figure;
% [a,b,h] = geoplot(eff./(1-eff), DataRef, DataMin, DataMax, projection);
% set(h, 'YTick', (tick_range-b)/a);
% set(h, 'YTickLabel', tick_range);
% title(h, '$\mathcal{E}/(1-\mathcal{E})$','interpreter', 'latex');
%ylabel(h, '$\mathcal{E}$/(1-$\mathcal{E})$','interpreter', 'latex','Rotation',0);

% figure;
% pcolor(x,y,eff./(1-eff)); shading flat; colorbar;
% ylim([-90 90]);
% axis equal;




%% Turbulent Pr

% Plotting input and colorbar setting
% DataMin = 0;        DataMax = 3;
% figure;
% h = imgplot(lon,lat,Prt(:,:,1), DataMin, DataMax);
% title(h, '$Pr_t$','interpreter', 'latex');

%tick_range = 0:5;
%ylabel(h, '$Pr_t$','interpreter', 'latex','Rotation',0);
% set(h, 'YTick', tick_range);


% figure;
% [a,b,h] = geoplot(Prt, DataRef, DataMin, DataMax, projection);
% set(h, 'YTick', (tick_range-b)/a);
% set(h, 'YTickLabel', tick_range);
% title(h, '$Pr_t$','interpreter', 'latex');


% % figure;
% % pcolor(x,y,Prt); shading flat; colorbar;
% % ylim([-90 90]);
% % axis equal;



%% Kappa based on New Eff

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
% aa = repmat(Kappa_osb_avg,[2 1]);       % averaged num. and denom. indiv
% bb = repmat(Kappa_new_avg,[2 1]);       

%aa = repmat(K_osb_avg,[2 1]);           % averaged the whole ratio
%bb = repmat(K_new_avg,[2 1]);

% figure(1); hold all;
% plot(aa(:), zlev,'b','LineWidth',3)
% plot(bb(:), zlev,'r','LineWidth',3)
% gtext(['$Ri$ =', num2str(Ri)],'interpreter','latex');


% plot(Kosb_avg, zlevel2,'go','LineWidth',2,'MarkerFace','g')
% plot(Knew_avg, zlevel2,'ro','LineWidth',2,'MarkerFace','g')
% 
% hndl = legend('Argo-Osborn',...
%               'Argo-DNS, $Ri=0.25$',...
%               'Argo-DNS, $Ri=0.5$',...
%               'Argo-DNS, $Ri=0.75$',...
%               'Argo-DNS, $Ri=0.95$');          
% set(hndl, 'interpreter','latex'); 

%% Plot the histograms of data
nbin = 200;

% Gamma histograms
% figure(1); hold on;
% h_Gamma = argo_PlotHist(Gamma(:),nbin,'k','k',0);
% xlabel('$\mathcal{E}/(1-\mathcal{E})$', 'interpreter','latex')
% 
% 
% % Reb and Kappa histograms
% figure; 
% h_Reb = argo_PlotHist(Reb(:),nbin,'b','b',1);
% xlabel('$Re_b$', 'interpreter','latex')
% 
% figure; hold on;
% h_Kosb = argo_PlotHist(K_osb(:),nbin,'g','g',1);
% h_Knew = argo_PlotHist(K_new(:),nbin,'k','k',1);
% xlabel('$K_{\rho}$ $(m^2/s)$', 'interpreter','latex')
% legend([h_Kosb h_Knew],{'Osborn (1980)~','gen. Osborn'},'interpreter','latex');
% xlim([1e-8 1e-2]); ylim([0 1500]); box on;
% % 
% figure; hold on;
% h_Kmnew = argo_PlotHist(Km_new(:),nbin,'k','k',1);
% xlabel('$K_{m}$ $(m^2/s)$', 'interpreter','latex')
% %ylabel('$count$', 'interpreter','latex')
% xlim([1e-8 1e-2]); ylim([0 1500]); box on;



% 
% %% green
% yy = K_new(:);      % Ri=0.25
% fcolor = 'g';
% ecolor = 'g';
% 
% %% red
% yy = K_new(:);      % Ri=0.44
% ecolor = 'r';
% fcolor = 'r';

%% Comparison with Bryan and Lewis (1979) profile and POP1 profiles

% % ===plot vertical profiles
% if(ic==ncase);
% fac = 1e-4;    % conversion factor from cm^2/s to m^2/s
% figure(1);  hold all;
% plot(kappa_BL  *fac,zz,'--','LineWidth',2); 
% plot(kappa_POP1*fac,zz,'--','LineWidth',2);
% plot(kappa_cesm*fac,depth3,'--','LineWidth',2);
% %plot(kappa_argo_cesm*fac,depth3,'--','LineWidth',2);
% plot(kappa_g ,depth1,'--','LineWidth',4);
% plot(kappa_ls(11:end)  ,depth1(11:end),'--','LineWidth',4);
% plot(kappa_dl,depth1,'--' ,'LineWidth',4);
% % %plot(kappa_nf,depth2,'--' ,'LineWidth',4);
% % 
% xlim([5e-6 1e-3]);
% ylim([-5000 0]);
% set(gca, 'XScale', 'log');
% grid on; box on;
% xlabel('$K_{\rho} (m^2/s)$', 'interpreter','latex')
% ylabel('$depth (m)$', 'interpreter','latex')
% hndl = legend('Argo-Osborn',...
%               'Argo-DNS',...
%               'BL79', ...
%               'NCAR CCSM3',...
%               'NCAR CESM1',...
%               'G03',...
%               'LS07',...
%               'DL10');
%               %'NF13',...
% end;
% set(hndl, 'interpreter','latex'); 
% colorCurve;

%exportfig(gcf, 'fig/Gamma_Ri044_250to450_new.pdf', 'FontMode', 'fixed','FontSize', 28,'width',18,'height',16);

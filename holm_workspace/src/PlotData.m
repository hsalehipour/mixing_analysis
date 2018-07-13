% 2) Ploting

%% ploting the evolution of KE
% xminplot = 0;       xmaxplot=400;
% yminplot = 0.89;     ymaxplot= 1;
% figure(1);     
% plot(time, KE/KE(1));               hold on;
% ylabel('$\mathcal{K}/\mathcal{K}_0$', 'interpreter','latex');
% xlabel('t');
% 
% figure(2);     
% semilogy(time,KE2d/KE(1));          hold on;
% ylabel('$\mathcal{K}_{2d}/\mathcal{K}_0$', 'interpreter','latex');
% xlabel('t');
% % 
% % 
% figure(3);     
% semilogy(time,KE3d/KE(1));          hold on;
% semilogy(t2d ,(KE3d(indt2d)/KE(1)),'*')
% % semilogy(tfo ,(KE3d(indtfo)/KE(1)),'o','MarkerSize',4)
% % semilogy(tfe ,(KE3d(indtfe)/KE(1)),'s','MarkerSize',4)
% ylabel('$\mathcal{K}_{3d}/\mathcal{K}_0$', 'interpreter','latex');
% xlabel('t');

% figure(4);     plot(time,KEbar/KE(1));       hold on;
% ylabel('$\overline{\mathcal{K}}/\mathcal{K}_0$', 'interpreter','latex');



%% ploting the error between sigma_2d/3d LHS and RHS
% figure;     hold on;
% plot(time,sig_LHS,'ro','MarkerSize',2)
% plot(time,sig_RHS,'b-');
% ylabel('${d\over dt}\mathcal{K}$', 'interpreter','latex');
% xlabel('t');
% hndl = legend(   '${d\over dt}\mathcal{K}$',...
%                  '-$\mathcal{H} + \mathcal{D}$');
% set(hndl, 'interpreter','latex');
% grid on; box on;
% % 
% % sigma3d
% figure;     hold on;
% plot(time,sig3_LHS,'ro','MarkerSize',2)
% plot(time,sig3_RHS,'b');
% ylabel('${d\over dt}\mathcal{K}_{3d}$', 'interpreter','latex');
% xlabel('t');
% hndl = legend(   '${d\over dt}\mathcal{K}_{3d}$',...
%                  '-$\mathcal{H}_{3d} + \mathcal{D}_{3d} + \overline{\mathcal{S}h} + \mathcal{S}h_{2d} + \mathcal{A}$');
% set(hndl, 'interpreter','latex');
% grid on; box on;
% % 
% % ploting the error between LHS and RHS of dP/dt = 2K*Hb+Dp
% figure;     hold on;
% plot(time,ddt_PE,'ro','MarkerSize',2)
% plot(time,2*KE.*sig_Hb+Dp);
% ylabel('${d\mathcal{P} \over dt}$', 'interpreter','latex');
% xlabel('t');
% hndl = legend(   '${d\over dt}\mathcal{P}$',...
%                  '$\mathcal{H} + \mathcal{D}_P$');
% set(hndl, 'interpreter','latex');
% grid on; box on;

% figure(10);     hold on;
% plot(time,sig_LHS,'bo','MarkerSize',2);
% plot(time,ddt_PE,'ro','MarkerSize',2)
% plot(time,sig3_LHS,'go','MarkerSize',2)
% 
% plot(time,sig_RHS,'b-');
% plot(time,2*KE.*sig_Hb+Dp,'r-');
% plot(time,sig3_RHS,'g-');
% 
% xlabel('t');
% hndl = legend(   '${d\over dt}\mathcal{K}$',...
%                  '${d\over dt}\mathcal{P}$',...
%                  '${d\over dt}\mathcal{K}_{3d}$',...
%                  '-$\mathcal{H} + \mathcal{D}$',...
%                  ' $\mathcal{H} + \mathcal{D}_P$',...
%                  '-$\mathcal{H}_{3d} + \mathcal{D}_{3d} + \overline{\mathcal{S}h} + \mathcal{S}h_{2d} + \mathcal{A}$');
% set(hndl, 'interpreter','latex');
% grid on; box on;


%% ploting the components sigma3d
% figure; hold on;
% % fac = KE3d./KE;
% fac = 2*KE3d;
% tstart = 0;
% yminplot = -4e-4;       ymaxplot =  4e-4;
% xminplot = 50;          xmaxplot=300;
% plot(time-tstart,sigma3d.*KE3d*2,'b--','LineWidth',1); 
% plot(time-tstart,-H3d,'g');
% plot(time-tstart,D3d,'r');
% plot(time-tstart,Shbar,'m');
% plot(time-tstart,Sh2d,'b');
% plot(time-tstart,An,'k');
% hndl = legend(   '${d\over dt}\mathcal{K}_{3d}$',...
%                  '-$\mathcal{H}_{3d}$', ...
%                  '$\mathcal{D}_{3d}$', ...
%                  '$\overline{\mathcal{S}h}$', ...
%                  '$\mathcal{S}h_{2d}$', ...
%                  '$\mathcal{A}$');
% set(hndl, 'interpreter','latex');
% xlabel('t');
% colorCurve;
% grid on; box on;

%% ploting the components sigma3d normalized by KE3d
% figure; hold on;
% fac = 1;
% yminplot = -0.1;  ymaxplot = 0.3;
% xminplot = 0;     xmaxplot = 350;
% tstart = 0;
% plot(time-tstart,fac.*(-sig3_Hb),'g');
% plot(time-tstart,fac.*sig3_Rs,'m');
% plot(time-tstart,fac.*sig3_Sh,'b');
% plot(time-tstart,fac.*sig3_An,'r');
% hndl = legend(   '-$\hat\mathcal{H}_{3d}$', ...
%                  '$\hat\mathcal{S}h$', ...
%                  '$\hat\mathcal{S}h_{2d}$', ...
%                  '$\hat\mathcal{A}$');
% set(hndl, 'interpreter','latex','location','northeast');
% xlabel('t');
% box('on')
% grid on;

%% ploting the components of sigma
% figure; hold on;
% tstart = 0;
% yminplot = -5e-4;   ymaxplot =  8e-4;
% xminplot = 0;       xmaxplot = 350;
% plot(time-tstart,D,'r');     hold on;
% plot(time-tstart,H,'g');
% plot(time-tstart,M,'b');
% %plot(time-tstart,S,'k');
% hndl = legend('$\mathcal{D}$',...
%                '$\mathcal{H}$',...
%                '$\mathcal{M}$');
% set(hndl, 'interpreter','latex');
% xlabel('t');
% box('on')
% grid on;

%% ploting the resolution vs. Kolmogroph scale
% Lz = 10;
% figure(9);     plot(time,rsln);     hold on;


%% plot N2 profile at the end
% figure(10);    plot(N2(:,end),z);   hold on;


%% ploting the potential energy budget
% figure;
% xminplot = 0;   xmaxplot = 350;
% yminplot = 0;   ymaxplot = 0.02;
% plot(time,PE-PE(1)-Dpt,time,RPE-RPE(1)-Dpt,time,APE-APE(1));         
% hold on;
% hndl = legend('$\mathcal{P}$ - $\mathcal{P}^0$ - $\mathcal{D}_p$',...
%               '$\mathcal{P}_B$ - $\mathcal{P}^0_{B}$ - $\mathcal{D}_p$',...
%               '$\mathcal{P}_A$ - $\mathcal{P}^0_{A}$');
% set(hndl, 'interpreter','latex');
% xlabel('t')
% 
% figure(26);     hold on;
% plot(time,eff,time,eff_c,time,eff3_c)


%% plot efficiency
% figure
% indt = time>=25;
% plot(time(indt),smooth(eff(indt)));     hold on;
% indt = time>=tfo;
% %plot(time(indt),smooth(Rf(indt)),'r')
% plot(time(indt),smooth(Rf(indt)./(1-Rf(indt))),'r')
% 
% yminplot = 0;       ymaxplot = 1;
% xminplot = 0;       xmaxplot = 350;
% xlabel('t')
% hndl = legend('$\mathcal{E}$',...
%               '$R_f$');
% set(hndl, 'interpreter','latex');

% % 
%  figure(5); hold on;
% plot(time,smooth(eff));
%  plot(time,smooth(eff)./(1-smooth(eff)));
% yminplot = 0;        ymaxplot = 1;
% xminplot = 25;       xmaxplot = 400;
%  xlabel('t');
% ylabel('$\mathcal{E}_i$','interpreter','latex');
% grid on; box on;
% hndl = legend('$Pr=1$','$Pr=2$','$Pr=4$','$Pr=8$','$Pr=16$');
% set(hndl, 'interpreter','latex');

% Calculating the standard deviation from mean during the fully turbulent
% phase
% Ri_range = 0.12:0.02:0.2;
% Pr_range = [1 2 4 8 16];
% xx = Pr_range;
% figure;
% errorbar(xx,eff_avg,eff_err,'k-o')
% % figure
% % errorbar(xx,eff_f_avg,eff_f_err,'k-o')
% xlabel('$Pr$', 'interpreter','latex');
% ylabel('$\mathcal{E}_i^f$','interpreter','latex');
% grid on; box on;
% hold on;
% plot(xx,xx*0+0.2,'k--');

% Plot cumulative efficiency for [t3d tfe]
% figure;
% plot(xx,intM./(intM+intD),'k-o')
% xlabel('$Pr$', 'interpreter','latex');
% ylabel('$\mathcal{E}_c$','interpreter','latex');
% grid on; box on;
% hold on;
% plot(xx,xx*0+0.2,'k--');

%% plot density profiles
% figure(33)
%plot(N2(:,end),z);     hold on;
% plot(ubar(:,end),z,'r');     hold on;
% plot(rhobar(:,end)-1,z,'g');     hold on;
% plot(rhob(:,end).*z,z,'g');     hold on;


 
%% Buoyancy Reynolds 
% figure(44);     hold on;
% xminplot = 0;       xmaxplot = 300;
% % indt = time>=xminplot & time<=xmaxplot;
% plot(time,Reb);              hold on;
% xlabel('t');
% ylabel('$Re_b$','interpreter','latex');
% grid on; box on;
% hndl = legend('$Pr=1$','$Pr=2$','$Pr=4$','$Pr=8$','$Pr=16$');
% set(hndl, 'interpreter','latex');

%% Taylor microscale Reynolds number
% figure(55);     hold on;
% xminplot = 0;       xmaxplot = 300;
% plot(time,Re_taylor);              hold on;
% xlabel('t');
% ylabel('$Re_{\lambda}$','interpreter','latex');
% grid on; box on;
% hndl = legend('$Pr=1$','$Pr=2$','$Pr=4$','$Pr=8$','$Pr=16$');
% set(hndl, 'interpreter','latex');

%% kappa_m, kappa_s, Turbulent PR
% figure;
% xminplot = t3d;     xmaxplot = trl;
% indt = time>=xminplot & time<=xmaxplot;
% semilogy(time(indt)-xminplot,kappa_3d(indt),'r');       hold on;
% [haxes,hline1,hline2] = plotyy(...
%                    time(indt)-xminplot,kappa_m(indt),...
%                    time(indt)-xminplot,Pr_t(indt),...
%                    'semilogy','semilogy');
% set(haxes(1),...
%     'xlim',[0 xmaxplot-xminplot],...
%     'ylim',[1e-5 1e-1],...
%     'box', 'off',...
%     'Ycolor','k');   
% ylabel(haxes(1),'$K_{\rho} ,K_m$','interpreter','latex');
% set(haxes(2),...
%     'xlim',[0 xmaxplot-xminplot],...
%     'ylim',[0 20],...
%     'Yscale', 'linear',...
%     'Ycolor','k');
% ylabel(haxes(2),'$Pr_t$','interpreter','latex');
% xlabel('$t$-$t_{3d}$','interpreter','latex');
% gtext('$K_m$', 'interpreter','latex','FontSize',16);
% gtext('$K_{\rho}$', 'interpreter','latex','FontSize',16);
% gtext('$Pr_t$', 'interpreter','latex','FontSize',16);



% figure;
% Pr_range = [1 2 4 8 16];
% plot(Pr_range,Prt_avg,'bo');


%% Pr-specific plots
% plot(time(indt)-t3d,smooth(Pr_t(indt)))
% plot(ssfit(time(indt)-t3d,Pr_t(indt),1e-2))
% 
% figure(25);     hold on;
% ind = time>=t3d & time<=trl;
% plot(Pr,mean(Pr_t(ind)),'bo');
% 
% figure(26);     hold on;
% plot(Pr,mean(kappa(ind)),'bo');
% plot(Pr,mean(kappa_3d(ind)),'b*');
% plot(Pr,mean(kappa_m(ind)),'bs');

% 


% Movie of ref density variation
% figure
% for i=1:ntime
%     subplot(121); plot(z,rhob(:,i),'r');    ylim([0 2])
%     title(['t=',num2str(round(time(i)))]);
% 
%     subplot(122); plot(z,rhobar(:,i));      ylim([0 2])
%     title(['t=',num2str(round(time(i)))]);
%     pause(0.1)
% end


% figure (35);    hold on;
% r = 0.001;
% cf1 = ssfit (time,sig3_Sh.*KE3d./KE,r);
% cf2 = ssfit (time,sig3_An.*KE3d./KE,r);
% %plot(time,feval(cf1,time))
% plot(time,feval(cf2,time))


% figure(63)
% L_en = (2*KE).^1.5./mean(epsbar,1)';
% L_o  = (mean(epsbar,1)'./N2avg.^1.5).^0.5;
% plot(time, L_en./L_o);      hold on;


%% Mixing factors
% xminplot = 0;         xmaxplot = 300;
% dKE = ddx(KE,time);
% dKE2d = ddx(KE2d,time);
% dKE3d = ddx(KE3d,time);

% (1)
% figure
% ymaxplot = 1e-1;     yminplot = 1e-6;
% semilogy(time,KE2d,time,KE3d);
% hndl = legend(  '$\mathcal{K}_{2d}$', ...
%                 '$\mathcal{K}_{3d}$');
% set(hndl, 'interpreter','latex');
% xlabel('t');
% box('on'); grid on;

% % (2)
% figure
% ymaxplot = 8e-4;     yminplot = -8e-4;
% plot(time,dKE,time,dKE2d,time,dKE3d);
% hndl = legend(  '${d\over dt} \mathcal{K}$', ...
%                 '${d\over dt} \mathcal{K}_{2d}$', ...
%                 '${d\over dt} \mathcal{K}_{3d}$');
% set(hndl, 'interpreter','latex');
% xlabel('t');
% box('on'); grid on;

% (3) H2d H3d
% ymaxplot =  8e-4;     yminplot = -4e-4;
% figure
% hold on;
% plot(time,H2d,'b');
% plot(time,H3d,'r')
% plot(time,H,'g--');
% hndl = legend(  '$\mathcal{H}_{2d}$', ...
%                 '$\mathcal{H}_{3d}$', ...
%                 '$\mathcal{H}$');
% set(hndl, 'interpreter','latex');
% xlabel('t');
% box('on')
% grid on;


%% Plot all evolution Eqns
% xminplot = 0;         xmaxplot = 300;
% % 
% % (1) KEbar evolution eqn
% ymaxplot =  6e-4;     yminplot = -8e-4;
% figure
% hold on;
% plot(time,dKEbar,'--', ...
%      time,-Shbar,...
%      time, XK,...
%      time, Dbar);
% hndl = legend(  ' ${d \over dt} \overline{\mathcal{K}}$', ...
%                 '-$\overline{\mathcal{S}h}$', ...
%                 ' $\mathcal{U}_{K}$',...
%                 ' $\overline{\mathcal{D}}$');
% set(hndl, 'interpreter','latex');
% xlabel('t');
% box('on')
% grid on;
% colorCurve;
% % 
% % (2) KE2d evolution eqn
% ymaxplot =  8e-4;     yminplot = -6e-4;
% figure
% hold on;
% plot(time,dKE2d, '--',...
%      time,-H2d,...
%      time,-Sh2d,...
%      time,-An,...
%      time,D2d,...
%      time,-XK);
% hndl = legend(  '${d \over dt} \mathcal{K}_{2d}$', ...
%                 '-$\mathcal{H}_{2d}$', ...
%                 '-$\mathcal{S}h_{2d}$',...
%                 '-$\mathcal{A}$',...
%                 ' $\mathcal{D}_{2d}$',...
%                 '-$\mathcal{U}_K$');
% set(hndl, 'interpreter','latex');
% xlabel('t');
% box('on')
% grid on;
% colorCurve
% 
% 
% (3) APE evolution eqn
% ymaxplot =  6e-4;     yminplot = -8e-4;
% figure
% hold on;
% plot(time,S, '--',...
%      time,H2d,...
%      time,-XP)
% hndl = legend(  ' ${d \over dt} \mathcal{P}_{A}$', ...
%                 ' $\mathcal{H}_{2d}$', ...             
%                 '-$\mathcal{U}_P$');
% set(hndl, 'interpreter','latex');
% xlabel('t');
% box('on')
% grid on;
% colorCurve


% (4) BPE evolution eqn
% ymaxplot =  4e-4;     yminplot = -3e-4;
% figure
% hold on;
% plot(time,M, '--',...
%      time,H3d,...
%      time,XP,...
%      time,Dp+time*0)
% hndl = legend(  ' ${d\over dt} \mathcal{P}_B$', ...
%                 ' $\mathcal{H}_{3d}$', ...
%                 ' $\mathcal{U}_P$',...
%                 ' $\mathcal{D}_{P}$');
% set(hndl, 'interpreter','latex');
% xlabel('t');
% box('on')
% grid on;
% colorCurve



%% Mixing: approximate formulation: M = H3d + beta*An
% xminplot = 0;         xmaxplot = 300;
% ymaxplot =  7e-4;     yminplot = -2e-4;
% figure
% hold on;
% plot(time,M, ...
%      time,H3d+beta(ic).*An,'--');
% hndl = legend(  '$\mathcal{M}$', ...
%                 '$\mathcal{M}^{\prime} = \mathcal{H}_{3d} +  \beta\mathcal{A}$');
% set(hndl, 'interpreter','latex');
% xlabel('t');
% box('on')
% grid on;
% colorCurve;
% 

%% Mixing: Cross Correlation plots
% figure; hold on;
% plot(lags, ccXPAn,'k')
% hndl = legend('$\langle \mathcal{U}_P$, $\mathcal{A} \rangle$');
% set(hndl, 'interpreter','latex');
% xlabel('\tau');          
% grid on; box on;
% ylim([-0.1 1]);
% % xlim([-length(lags)/2, length(lags)/2])
% xlim([-150, 150])


% figure;
% xminplot = 0;         xmaxplot = 300;
% ymaxplot = 1e-4;      yminplot = -0.1e-4;
% % ymaxplot = 6e-4;      yminplot = -1e-4;
% plot(time,XP,time,An)
% hndl = legend(  '$\mathcal{U}_P$', ...
%                 '$\mathcal{A}$');
% set(hndl, 'interpreter','latex');
% xlabel('t');
% box('on'); grid on;
% colorCurve;



%% Entrainement
% Plot the density profile
% figure(10)
% plot(time,RPE-RPE(1))
% ylabel('$\mathcal{P}_B$', 'interpreter','latex');
% xlabel('t');
% hold on;
% % 
% figure(11); hold on;
% [~, ind] = min(abs(time-300));
% plot(rhobar(:,ind),z);
% ylabel('z');
% hndl = xlabel('$\overline{\rho}$');
% set(hndl, 'interpreter','latex');
% hndl = legend('$Pr=1$','$Pr=2$','$Pr=4$','$Pr=8$','$Pr=16$');
% set(hndl, 'interpreter','latex');


%% Mixing ingredients
% Pr_range = [1 2 4 8 16];      xx = Pr_range;
% Ri_range = 0.10:0.02:0.2;       xx = Ri_range;
% Ri_range = [0.01 0.02 0.04 0.08 0.16];       xx = Ri_range;

% Re_range = [4000 6000 8000 12000];  xx = Re_range;

% Reb_range(ic) = mean(Reb(ind));
% % 
% % Plot the plain ingredients
% figure; hold on;
% plot(xx, intM);
% plot(xx, intXP);
% plot(xx, intH3d);
% plot(xx, intAn);
% grid on; box on;
% hndl = legend('$\langle \mathcal{M}     \rangle_f$',...
%               '$\langle \mathcal{U}_P   \rangle_f$',...
%               '$\langle \mathcal{H}_{3d}\rangle_f$',...
%               '$\langle \mathcal{A}\rangle_f$');
% set(hndl, 'interpreter','latex');        
% % xlabel('$Pr$','interpreter','latex');
% xlabel('$Ri_0$','interpreter','latex');
% % xlabel('$Re$','interpreter','latex');
% 
% colorCurve;
% gtext('$\langle \mathcal{M}     \rangle_f$','interpreter','latex');
% gtext('$\langle \mathcal{U}_P   \rangle_f$','interpreter','latex');
% gtext('$\langle \mathcal{H}_{3d}\rangle_f$','interpreter','latex');
% gtext('$\langle \mathcal{A}     \rangle_f$','interpreter','latex');


% 
% 
% Plot the ingredients normalized by D
% figure; hold on;
% plot(xx, intM./intD);
% plot(xx, intXP./intD);
% plot(xx, intH3d./intD);
% plot(xx, intAn./intD);
% colorCurve;
% plot(xx, intM./(intM+intD));
% grid on; box on;
% hndl = legend('$\langle \mathcal{M}     \rangle_f/ \langle\mathcal{D} \rangle_f$',...
%               '$\langle \mathcal{U}_P   \rangle_f/ \langle\mathcal{D} \rangle_f$',...
%               '$\langle \mathcal{H}_{3d}\rangle_f/ \langle\mathcal{D} \rangle_f$',...
%               '$\langle \mathcal{A}\rangle_f/ \langle\mathcal{D}  \rangle_f$',...
%               '$\mathcal{E}_c$');
% set(hndl, 'interpreter','latex');        
% % % xlabel('$Pr$','interpreter','latex');
% xlabel('$Ri_0$','interpreter','latex');
% xlabel('$Re$','interpreter','latex');
% gtext('$\langle \mathcal{M}     \rangle_f/ \langle\mathcal{D} \rangle_f$','interpreter','latex');
% gtext('$\langle \mathcal{U}_P   \rangle_f/ \langle\mathcal{D} \rangle_f$','interpreter','latex');
% gtext('$\langle \mathcal{H}_{3d}\rangle_f/ \langle\mathcal{D} \rangle_f$','interpreter','latex');
% gtext('$\langle \mathcal{A}\rangle_f/ \langle\mathcal{D}  \rangle_f$','interpreter','latex');
% gtext('$\mathcal{E}_c$','interpreter','latex');

% % plot beta
% figure(53)
% plot(Pr_range, beta);
% xlabel('$Pr$','interpreter','latex');
% ylabel('$\beta$','interpreter','latex');
% grid on; box on;

%% Regimes of mixing

% figure(14);
% PI_m = intXP./intH3d;
% p = polyfit(Ri_turb,PI_m,2);
% plot(Ri_turb,PI_m,'k^');     hold all;
% xx = linspace(min(Ri_turb), max(Ri_turb), 100);
% plot(xx,polyval(p,xx),'k');
% plot([0.3 0.45], [1 1], 'k--');
% ylabel('$\Pi_m = \langle\mathcal{U}_P\rangle_f / \langle\mathcal{H}_{3d}\rangle_f$','interpreter','latex')
% xlabel('$\langle Ri \rangle_z $', 'interpreter','latex')



%% New plots: Kappa JFM paper

%{
% if Pr==1;   MarkerChar='o';      end;
% if Pr==2;   MarkerChar='x';      end;
% if Pr==4;   MarkerChar='+';      end;
% if Pr==8;   MarkerChar='^';      end;
% if Pr==16;  MarkerChar='*';      end;
% MarkerSize      = -50./log10(Ri_avg(ind));
% MarkerFColor    = [1 1 1];

% figure (1);
% hold all;
% scatter(kappa_star(ind),kappa_osb(ind),...
%      MarkerSize, MarkerChar, ...
%     'MarkerEdgeColor',[0 0 0],...
%     'MarkerFaceColor',[1 1 1]);
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');
% xlabel('$\widetilde{K}^*_{\rho}$', 'interpreter','latex')
% ylabel('$\widetilde{K}^{osb}_{\rho}$', 'interpreter','latex')
% xminplot = 0.1;      xmaxplot = 1e3;
% yminplot = 0.1;      ymaxplot = 1e3;
% axis equal; box on;
% % loglog([0.1 2e3], [0.1 2e3],'k-')

%
% figure (2);
% hold all;
% scatter(kappa_star(ind),kappa_cox(ind),...
%      MarkerSize, MarkerChar, ...
%     'MarkerEdgeColor',[0 0 0],...
%     'MarkerFaceColor',[1 1 1]);
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');
% xlabel('$\widetilde{K}^*_{\rho}$', 'interpreter','latex')
% ylabel('$\widetilde{K}^{cox}_{\rho}$', 'interpreter','latex')
% xminplot = 0.1;      xmaxplot = 1e3;
% yminplot = 0.1;      ymaxplot = 1e3;
% axis equal; box on;
% % loglog([0.1 2e3], [0.1 2e3],'k-')


% figure (3);
% hold all;
% %scatter(kappa_star(ind),kappa_3d(ind),...
% scatter(kappa_star(ind),kappa(ind),...
%      MarkerSize, MarkerChar, ...
%     'MarkerEdgeColor',[0 0 0],...
%     'MarkerFaceColor',[1 1 1]);
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');
% xlabel('$\widetilde{K}^*_{\rho}$', 'interpreter','latex')
% ylabel('$\widetilde{K}^{ml}_{\rho''}$', 'interpreter','latex')
% % ylabel('$\widetilde{K}^{ml}_{\rho_{3d}}$', 'interpreter','latex')
% xminplot = 0.1;      xmaxplot = 1e3;
% yminplot = 0.1;      ymaxplot = 1e3;
% axis equal; box on;
% % loglog([0.1 2e3], [0.1 2e3],'k-')

%====================================
% K_rho, K_m and eff plots
%====================================
% figure (4);
% % scatter(Reb(ind),kappa_star(ind),... 
% scatter(Reb(ind),kappa_star(ind)/Pr,...
%       MarkerSize, MarkerChar, ...
%      'MarkerEdgeColor',[0 0 0],...
%      'MarkerFaceColor',[1 1 1]);
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');
% hold all;
% xlabel('$Re_b$', 'interpreter','latex')
% ylabel('$\widetilde{K}_{\rho}/ Pr$', 'interpreter','latex')
% % ylabel('$\widetilde{K}_{\rho}= K_{\rho}/\kappa$', 'interpreter','latex')
% xminplot = 1;        xmaxplot = 1e4;
% yminplot = 0.1;      ymaxplot = 1e3;
% loglog([xminplot xmaxplot], 0.2*[xminplot xmaxplot],'k--')
% loglog([10 xmaxplot], 2*[10 xmaxplot].^0.5,'k--')
% p1 = 10^(2/3)*Pr^(-0.5);
% p2 = (3*log(sqrt(Pr)))^2;
% loglog([p1 p2],  0.2*[p1 p2].^(1.5),'k--')
% loglog([p2 100], 0.2*[p2 100],'r--')
% loglog([100 xmaxplot], 2*[100 xmaxplot],'b--')


% figure (2);
% gamma = smooth(eff(ind),'rlowess')./(1-smooth(eff(ind),'rlowess'));
% % scatter(Reb(ind),gamma,...
% scatter(Reb(ind),smooth(eff(ind),'rlowess'),...
%        MarkerSize, MarkerChar, ...
%       'MarkerEdgeColor',[0 0 0],...
%       'MarkerFaceColor',MarkerFColor);
% set(gca, 'XScale', 'log');
% hold all;
% xlabel('$Re_b$', 'interpreter','latex')
% ylabel('$\mathcal{E}$', 'interpreter','latex')
% % ylabel('$\Gamma \hspace{5pt} = \mathcal{E}$/(1-$\mathcal{E}$)', 'interpreter','latex')
% box on;
% xminplot = 1;            xmaxplot = 1e4;
% yminplot = 1e-2;         ymaxplot = 0.5;
% loglog([10 xmaxplot], 1.5*[10 xmaxplot].^(-0.5),'k--')
% loglog([1 1e3], [0.2 0.2],'k--')
% 
% % % Note: to address 2d-residu in response to Mash et al issue of
% Gamma=1/5
% MarkerSize = 35;
% if Pr>1;   
%     MarkerChar='^';
%     MarkerFColor = [1 0 0];
% end;
% if Ri<0.12;   
%     MarkerChar='s';
%     MarkerFColor = [0 1 0];
% end;
% % [~, indtt3d] = max(KE3d);
% % tt3d = time(indtt3d);
% % tind = time >= tt3d & time<t3d;
% % gamma = smooth(eff(tind),'rlowess')./(1-smooth(eff(tind),'rlowess'));
% % MarkerSize = -50./log10(Ri_avg(tind));
% % scatter(Reb(tind),gamma,...
% %        MarkerSize, MarkerChar, ...
% %       'MarkerEdgeColor',[0 0 0],...
% %       'MarkerFaceColor',[1 0 0]);

% figure (3);
% scatter(Ri_avg(ind),smooth(eff(ind),'rlowess'),...
%       MarkerSize, MarkerChar, ...
%      'MarkerEdgeColor',[0 0 0],...
%      'MarkerFaceColor',[1 1 1]);
% set(gca, 'XScale', 'log');
% hold all;
% xlabel('$Ri_b = N^2/\langle \left( d\overline{u}/ dz\right)^2\rangle$', 'interpreter','latex')
% ylabel('$\mathcal{E}$', 'interpreter','latex')
% box on;
% xminplot = 0.01;            xmaxplot = 0.6;
% yminplot = 1e-2;         ymaxplot = 0.5;
% loglog([xminplot xmaxplot], 0.6*[xminplot xmaxplot].^(0.6),'k--')
% loglog([xminplot xmaxplot], [0.2 0.2],'k--')


% figure (3);
% scatter(Reb(ind),kappa_m(ind),...
%       MarkerSize, MarkerChar, ...
%      'MarkerEdgeColor',[0 0 0],...
%      'MarkerFaceColor',[1 1 1]);
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');
% hold all;
% xlabel('$Re_b$', 'interpreter','latex')
% ylabel('$\widetilde{K}_m$', 'interpreter','latex')
% box on;
% % xminplot = 1;         xmaxplot = 1e4;
% % yminplot = 0.1;         ymaxplot = 1e3;
% % loglog([xminplot xmaxplot], [xminplot xmaxplot],'k--')
% 
% 
% figure (4);
% scatter(Reb(ind),kappa_m(ind)./Ri_avg(ind),...
%       MarkerSize, MarkerChar, ...
%      'MarkerEdgeColor',[0 0 0],...
%      'MarkerFaceColor',[1 1 1]);
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');
% hold all;
% xlabel('$Re_b$', 'interpreter','latex')
% ylabel('$\widetilde{K}_m/Ri_b$', 'interpreter','latex')
% box on;




%====================================
% Pr_t plots
%====================================
% figure (4);
% scatter(Reb(ind),Ri_avg(ind)./smooth(eff(ind),'rlowess'),...
%      MarkerSize, MarkerChar, ...
%     'MarkerFaceColor', [1 1 1],...
%     'MarkerEdgeColor', [0 0 0]);
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');
% hold all;
% xminplot = 1;         xmaxplot = 1e4;
% yminplot = 0.1;       ymaxplot = 100;
% xlabel('$Re_b$', 'interpreter','latex')
% ylabel('$Pr_t$', 'interpreter','latex')
% box on;

% figure (5);
% scatter(Ri_avg(ind),Ri_avg(ind)./smooth(eff(ind),'rlowess'),...
%     MarkerSize, MarkerChar, ...
%     'MarkerFaceColor', [1 1 1],...
%     'MarkerEdgeColor', [0 0 0]);
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');
% hold all;
% xminplot = 0.01;         xmaxplot = 1;
% yminplot = 0.5;         ymaxplot = 100;
% xlabel('$Ri_b = N^2/\langle \left( d\overline{u}/ dz\right)^2\rangle$', 'interpreter','latex')
% ylabel('$Pr_t$', 'interpreter','latex')
% box on;


% Appendix Figures
% figure (1);
% scatter(time(ind),Reb(ind),...
%       MarkerSize, MarkerChar, ...
%      'MarkerEdgeColor',[0 0 0],...
%      'MarkerFaceColor',[1 1 1]);
% set(gca, 'XScale', 'log');
% hold all;
% xlabel('$Re_b$', 'interpreter','latex')
% ylabel('$\mathcal{E}$', 'interpreter','latex')
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');
% hold all;
% xlabel('$time$', 'interpreter','latex')
% ylabel('$Re_b$', 'interpreter','latex')
% xminplot = 1e2;         xmaxplot = 2e4;
% yminplot = 1;           ymaxplot = 1e4;
% box on;

% figure (2);
% scatter(Ri_avg(ind),Reb(ind),...
%       MarkerSize, MarkerChar, ...
%      'MarkerEdgeColor',[0 0 0],...
%      'MarkerFaceColor',[1 1 1]);
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');
% hold all;
% % xlabel('$Ri_b = \langle\overline{N^2}\rangle/\langle\overline{S^2}\rangle$', 'interpreter','latex')
% xlabel('$Ri_b = N^2/\langle \left( d\overline{u}/ dz\right)^2\rangle$', 'interpreter','latex')
% ylabel('$Re_b$', 'interpreter','latex')
% xminplot = 0;         xmaxplot = 0.7;
% yminplot = 1;         ymaxplot = 1e4;


% Plotting relationship function between Rib=N2/S2 and N2/<S>^2
% figure (3);           
% scatter(Ri_bul(ind),Ri_avg(ind),...
%       MarkerSize, MarkerChar, ...
%      'MarkerEdgeColor',[0 0 0],...
%      'MarkerFaceColor',[1 1 1]);
% hold all;
% xlabel('$Ri_b^{\dagger}$', 'interpreter','latex')
% ylabel('$Ri_b$', 'interpreter','latex')
% %xlabel('$Ri_b^{\ddagger} = g/\rho_0\langle d\overline{\rho}/ dz \rangle_d/ \langle d\overline{u}/ dz \rangle_d^2  $', 'interpreter','latex')
% %ylabel('$Ri_b^{\dagger } = g/\rho_0\langle d\overline{\rho}/ dz \rangle  / \langle \left( d\overline{u}/dz\right)^2\rangle$', 'interpreter','latex')
% box on; grid on;
% xx = linspace(0, 3.5, 100);
% plot(xx, 0.9*xx./(xx+0.8),'r-');
%----
% Fit info
% General model Rat11:
%      f(x) = (p1*x + p2) / (x + q1)
% Coefficients (with 95% confidence bounds):
%        p1 =      0.8933  (0.8667, 0.9198)
%        p2 =    -0.01021  (-0.0125, -0.00792)
%        q1 =       0.777  (0.7312, 0.8229)
% 
% Goodness of fit:
%   SSE: 1.294
%   R-square: 0.9744
%   Adjusted R-square: 0.9744
%   RMSE: 0.02216

 
  
% figure (3);
% % scatter(Re_taylor(ind)./Ri_avg(ind),kappa_star(ind)/Pr,...
% % scatter((10/3).^2./Ri_avg(ind).^(1.6),Reb(ind),...
% scatter(Ri_avg(ind),0.2*Reb(ind).^((-0.2+Ri_avg(ind))),...
%       MarkerSize, MarkerChar, ...
%      'MarkerEdgeColor',[0 0 0],...
%      'MarkerFaceColor',[1 1 1]);
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');
% hold all;
% xlabel('$Re_{\lambda}/Ri_b$', 'interpreter','latex')
% ylabel('$Re_b$', 'interpreter','latex')
% xminplot = 1;         xmaxplot = 1e4;
% yminplot = 1;         ymaxplot = 1e4;
% grid on;
%}

%% JFM: Holmboe vs KH mixing

% Figure 1: Isosurfaces produced by VisIt

% Figure 2: includes four panels:
if R==1;   % KH
    LineThickness=1;
    LineColor = 'k';
    FaceColor = [0.5 0.5 0.5];
end   
if R~=1;   % Holmboe
    LineThickness=4;
    LineColor = 'k';
    FaceColor = [1 1 1];
end      
% MarkerSize      = 55;%6e3./(time-t2d);
% % 
%exportfig(gcf, 'fig/test.eps', 'FontMode', 'fixed','FontSize', 16,'color','cmyk','width',16,'height',14);
%exportfig(gcf, 'fig/test.eps', 'FontMode', 'fixed','FontSize', 16,'color','cmyk','width',6.5,'height',5.5);
     

%=======
% figure(1); hold all;
% plot(time-t2d,KE2d,'-' ,'Color',LineColor,'LineWidth',LineThickness)
% plot(time-t2d,KE3d,'--','Color',LineColor,'LineWidth',LineThickness)
% xlabel('$t$-$t_{2d}$', 'interpreter','latex')
% box on;
% yminplot = 1e-6;        ymaxplot =  2e-2;
% xminplot = -150;        xmaxplot=300;
% xlim([xminplot xmaxplot]);
% ylim([yminplot ymaxplot]);
% ylabel('$\hspace{5pt}$','interpreter','latex');
% set(gca, 'YScale', 'log');
% hndl = legend('$H : \mathcal{K}_{2d}$',...
%               '$H : \mathcal{K}_{3d}$',...
%               '$K: \mathcal{K}_{2d}$',...
%               '$K: \mathcal{K}_{3d}$');
% set(hndl, 'interpreter','latex');      
% 
%========
% % exportfig(gcf, 'fig/BPE_APE.eps', 'FontMode', 'fixed','FontSize', 18,'width',4*8/3,'height',5);
% figure(2); hold all; 
% xminplot = -150;        xmaxplot=300;
% yminplot = 1e-5;        ymaxplot=2e-2;
% plot(time-t2d,APE-APE(1),'k-','LineWidth',LineThickness)
% plot(time-t2d,RPE-RPE(1)-Dpt,'k--','LineWidth',LineThickness)
% xlabel('$t$-$t_{2d}$', 'interpreter','latex')
% hndl = legend('$H: \Delta \mathcal{P}_A$',...
%               '$H: \Delta \mathcal{P}_B$ - $D_pt$' ,...
%               '$K: \Delta \mathcal{P}_A$',...
%               '$K: \Delta \mathcal{P}_B$ - $D_pt$');              
% set(hndl, 'interpreter','latex');  
% set(gca,'Yscale','log');
% xlim([xminplot xmaxplot]);
% ylim([yminplot ymaxplot]);
% box on;

% %=======
% figure(2); hold all;
% yminplot =  1e-6;        ymaxplot = 1e-3;
% xminplot = -150;         xmaxplot = 300;
% plot(time-t2d,  M ,'-'   ,'Color',LineColor,'LineWidth',LineThickness);
% plot(time-t2d, -D ,'--'  ,'Color',LineColor,'LineWidth',LineThickness);
% set(gca,...
%     'xlim',[xminplot xmaxplot],...
%     'ylim',[yminplot ymaxplot],...
%     'yscale','log')
% xlabel('$t$-$t_{2d}$','interpreter','latex');
% ylabel('$\mathcal{M}, \varepsilon$','interpreter','latex');
% hndl = legend('$H : \mathcal{M}$',...
%               '$H : \varepsilon$',...
%               '$K: \mathcal{M}$',...
%               '$K: \varepsilon$');
% set(hndl, 'interpreter','latex');  
% box on; grid on;
% 
% %=======
% figure(3); hold all;
% xminplot = -150;     xmaxplot = 300;
% yminplot = 1;        ymaxplot= 2e2;
% plot(time-t2d,  Reb ,'-'   ,'Color',LineColor,'LineWidth',LineThickness);
% set(gca,...
%     'xlim',[xminplot xmaxplot],...
%     'ylim',[yminplot ymaxplot],...
%     'yscale','log');
% xlabel('$t$-$t_{2d}$','interpreter','latex');
% ylabel('$Re^*_b$','interpreter','latex');
% hndl = legend('$H$','$K$');
% set(hndl, 'interpreter','latex');  
% box on; grid on;

%
%=======
% SAMPLE for double axis (for future reference)
% figure(2); hold all;
% yminplot =  1e-6;        ymaxplot = 2e-3;
% xminplot = -150;         xmaxplot = 300;
% y2minplot = 1e-1;        y2maxplot= 2e2;
% 
% plot(time-t2d, M ,'-'   ,'Color',LineColor,'LineWidth',LineThickness);
% [haxes,hline1,hline2] = plotyy( time-t2d,-D , time-t2d,Reb,'semilogy','semilogy');
% set(hline1, 'LineStyle','--','Color',LineColor,'LineWidth',LineThickness);
% set(hline2, 'LineStyle',':','Color',LineColor,'LineWidth',LineThickness+2);
% set(haxes(1),...
%     'xlim',[xminplot xmaxplot],...
%     'ylim',[yminplot ymaxplot],...
%     'yscale','log',...
%     'box', 'off',...
%     'ycolor','k');   
% set(haxes(2),...
%     'xlim',[xminplot xmaxplot],...
%     'ylim',[y2minplot y2maxplot],...
%     'yscale','log',...
%     'box', 'off',...
%     'ycolor','k');   
% 
% xlabel('$t$-$t_{2d}$','interpreter','latex');
% if ic==ncase
%     ylabel(haxes(1),'$\mathcal{M}, \varepsilon$','interpreter','latex');
%     ylabel(haxes(2),'$Re^*_b$','interpreter','latex');    
% end
% gtext('$\mathcal{M}$', 'interpreter','latex','FontSize',16);
% gtext('$\varepsilon$', 'interpreter','latex','FontSize',16);
% gtext('$Re^*_b$', 'interpreter','latex','FontSize',16);

% 
% 
% %=======
% figure(4); hold all;     
% yminplot = 0;            ymaxplot = 0.6;
% xminplot = 1e0;          xmaxplot = 2e2;
% tind = Reb>=3;
% scatter(Reb(tind),eff(tind),...
%     MarkerSize,'o',...
%     'MarkerFaceColor', FaceColor,...
%     'MarkerEdgeColor','k');
% grid on; box on;
% xlim([xminplot xmaxplot]);
% ylim([yminplot ymaxplot]);
% set(gca, 'XScale', 'log');
% % set(gca, 'YScale', 'log');
% xlabel('$Re^*_b$', 'interpreter','latex')
% ylabel('$\eta$', 'interpreter','latex')
% hndl = legend('$H$','$K$');
% set(hndl, 'interpreter','latex');
% % plot([xminplot xmaxplot], 1/6+0*[xminplot xmaxplot],'k--');
% [~, myind] = min(abs(time-t2d-100));
% scatter(Reb(indt2d),eff(indt2d),...
%     3*MarkerSize,'+',...
%     'LineWidth', 2,...
%     'MarkerEdgeColor','r');
% scatter(Reb(indt3d),eff(indt3d),...
%     3*MarkerSize,'*',...
%     'LineWidth', 2,...
%     'MarkerEdgeColor','r');
% scatter(Reb(myind),eff(myind),...
%     3*MarkerSize,'x',...
%     'LineWidth', 2,...
%     'MarkerEdgeColor','r');
% %=======
% figure(5); hold all;      
% yminplot = 1e0;          ymaxplot = 1e3;
% xminplot = 1e0;          xmaxplot = 2e2;
% tind = time>=t3d & Reb>=3;
% scatter(Reb(tind),kappa_star(tind),...
%      MarkerSize,'o',...
%     'MarkerFaceColor', FaceColor,...
%     'MarkerEdgeColor','k');
% %plot(Reb(ind),kappa_cox(ind) );
% grid on; box on;
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');
% xlim([xminplot xmaxplot]);
% ylim([yminplot ymaxplot]);
% xlabel('$Re^*_b$', 'interpreter','latex')
% ylabel('$K^*_{\rho}/\kappa$', 'interpreter','latex')
% % plot([xminplot xmaxplot], 0.2*Pr*[xminplot xmaxplot],'k--');
% scatter(Reb(indt3d),kappa_star(indt3d),...
%     3*MarkerSize,'*',...
%     'LineWidth', 2,...
%     'MarkerEdgeColor','r');
% scatter(Reb(myind),kappa_star(myind),...
%     3*MarkerSize,'x',...
%     'LineWidth', 2,...
%     'MarkerEdgeColor','r');


% % Figures 3 & 4 : include the following z-t contour plots
% %============== z-t plot =============
% exportfig(gcf, 'fig/H_zt_cox.eps', 'FontMode', 'fixed','FontSize', 12,'width',24,'height',8);
% exportfig(gcf, 'fig/KH_zt_chi.eps', 'FontMode', 'fixed','FontSize', 17,'width',36,'height',12);
 
%== Snippet to change the width of the colorbar
% axpos = get(gca,'Position');
% cpos = get(hndl,'Position');
% cpos(3) = 1.5*cpos(3);
% set(hndl,'Position',cpos);
% set(gca,'Position',axpos);
%==

% % addpath(genpath('/home/hesam/holm_workspace/colormaps/'))
% ncnt        = 1024;        % contour levels
% ncnt_chi    = 1024;        % for plots with chi
% ncmp        = 1024;        % colormap levels
% indz   = z>=-5 & z<=5;
% tmax = 300;
% [TTIME, TZZ]=meshgrid(time-t2d,z(indz));
% 
% % for some reasons chi_pr has some grid-dependant noise, BUT I checked by
% % increasing nz in .usr and these did not improve.
% chi_smoothed = chi_pr*0;        
% for i=1:ntime
%     %chi_smoothed(:,i) = smooth(z,chi_pr(:,i),'sgolay',4);
%     chi_smoothed(:,i) = smooth(z,chi_pr(:,i),17);
% end
% indeps0= log10(epsbar(indz,:))<-4.5;
% indchi0= log10(2*kappa0*chi_smoothed(indz,:))<-6;
% indzero= abs(gradrhob(indz,:))<1e-3;

% % N2
% figure;
% contourf(TTIME,TZZ,N2(indz,:),ncnt,'LineStyle','none')
% hndl = colorbar;%       caxis([0 1])
% colormap(cubehelix(ncmp));
% title(hndl, '$\overline{N^2}$','interpreter', 'latex');
% xlabel('$t$-$t_{2d}$', 'interpreter','latex')
% ylabel('$z$', 'interpreter','latex')
% xlim([min(time-t2d) tmax])
% 
% espbar
% figure;
% contourf(TTIME,TZZ,log10(epsbar(indz,:)),ncnt,'LineStyle','none')
% hndl = colorbar;       caxis([-4.5 -2.5])
% colormap(cubehelix(ncmp));
% title(hndl, '$log (\overline{\varepsilon}) $','interpreter', 'latex');
% xlabel('$t$-$t_{2d}$', 'interpreter','latex')
% ylabel('$z$', 'interpreter','latex')
% xlim([min(time-t2d) tmax])
% 
% 
% % chi = | grad \rho'|^2 where \rho' = rho_2d + rho_3d
% figure;
% contourf(TTIME,TZZ,log10(2*kappa0*chi_smoothed(indz,:)),ncnt_chi,'LineStyle','none')
% hndl = colorbar;       caxis([-6 -2])
% colormap(cubehelix(ncmp))
% title(hndl, '$log (\overline{\chi}) $','interpreter', 'latex');
% xlabel('$t$-$t_{2d}$', 'interpreter','latex')
% ylabel('$z$', 'interpreter','latex')
% xlim([min(time-t2d) tmax])
% 
% % Re*_bar = epsbar/nu/N2
% figure;
% Reb_bar = epsbar(indz,:)./N2(indz,:)/nu;
% Reb_bar(indeps0) = 1e-10  ;
% Reb_bar(indzero) = 1e-10  ;
% contourf(TTIME,TZZ,log10(Reb_bar),ncnt,'LineStyle','none')
% hndl = colorbar;       caxis([1 3])
% colormap(cubehelix(ncmp))
% title(hndl, '$log (\overline{Re_b}) $','interpreter', 'latex');
% xlabel('$t$-$t_{2d}$', 'interpreter','latex')
% ylabel('$z$', 'interpreter','latex')
% xlim([min(time-t2d) tmax])
% 
% % Cox number = chi/grad(rhobar)^2
% figure;
% kappa_profile_cox = chi_smoothed(indz,:)./gradrho(indz,:).^2;
% kappa_profile_cox(indchi0) = 1e-10  ;
% kappa_profile_cox(indzero) = 1e-10  ;
% contourf(TTIME,TZZ,log10(kappa_profile_cox),ncnt_chi,'LineStyle','none')
% hndl = colorbar;       caxis([0 3])
% colormap(cubehelix(ncmp))
% title(hndl, '$log (K^{cox}_{\rho}/\kappa) $','interpreter', 'latex');
% xlabel('$t$-$t_{2d}$', 'interpreter','latex')
% ylabel('$z$', 'interpreter','latex')
% xlim([min(time-t2d) tmax])


% Buoyancy flux: Hbar= (rho w)
% figure;
% indH0= log10(Hbar(indz,:))<-8;  % small values
% indHn= Hbar(indz,:)<0;          % negative values
% kappa_profile_ml = Hbar(indz,:)./N2(indz,:)/kappa0;
% kappa_profile_ml(indH0) = 1e-10;
% kappa_profile_ml(indHn) = NaN;
% contourf(TTIME,TZZ,log10(kappa_profile_ml),ncnt,'LineStyle','none')
% hndl = colorbar;       caxis([0 3])
% colormap(parula(ncmp));
% title(hndl, '$\overline{\rho w}$','interpreter', 'latex');
% xlabel('$t$-$t_{2d}$', 'interpreter','latex')
% ylabel('$z$', 'interpreter','latex')
% xlim([min(time-t2d) tmax])
% % gamma_profile = Hbar(indz,:)./epsbar(indz,:);        % gamma = R_f/(1-R_f)
% % Mchi = -kappa0*g*chi_ref./gradrhob;

%=======
% added  section on Reynolds number effect on HWI (revisied paper)
%=======
% Added for illustrating the sigmoidal dependence of sigma3d on Re
% myind = sigma3d>0 & KE3d>=1e-8*KE(1) & KE3d<=1e-4*KE(1);
% aa(ic) = mean(sigma3d(myind));
% myfac = 1e-4./max(exp(2*aa(ic)*(time(myind))));
% figure(3); hold all;
% plot(time(myind), myfac*exp(2*aa(ic)*(time(myind))),'k--')

% figure(4); 
% plot([1 3 5 10 20 40 60]*100, [0, aa])
% xlabel('$Re$', 'interpreter','latex')
% ylabel('$\sigma_{3d}$', 'interpreter','latex')

%=======
% plots of rho, U profiles 
%=======
% figure(1); hold all;
% plot(rhobar(:,indt2d),z,'k-' ,'LineWidth',LineThickness);
% plot(ubar(:,indt2d)+1,z,'k--','LineWidth',LineThickness);
% hndl = legend('$H: \overline{\rho}(z)$',...
%               '$H: \overline{U}(z)+1$',...
%               '$K: \overline{\rho}(z)$',...
%               '$K: \overline{U}(z)+1$');
% set(hndl, 'interpreter','latex');  
% ylim([-5 5]);   xlim([-0.1 2.1]);
% if ic==1; gtext('$t = t_{2d}$','interpreter','latex'); end;
% xlabel('$\overline{\rho},\overline{U}+1$', 'interpreter','latex')
% ylabel('$z$', 'interpreter','latex')
% box on;
% 
% 
% figure(2);  hold all;
% plot(rhobar(:,indt3d),z,'k-' ,'LineWidth',LineThickness)
% plot(ubar(:,indt3d)+1,z,'k--','LineWidth',LineThickness);
% hndl = legend('$H: \overline{\rho}(z)$',...
%               '$H: \overline{U}(z)+1$',...
%               '$K: \overline{\rho}(z)$',...
%               '$K: \overline{U}(z)+1$');
% set(hndl, 'interpreter','latex');  
% ylim([-5 5]);   xlim([-0.1 2.1]);
% if ic==1; gtext('$t = t_{3d}$','interpreter','latex'); end
% xlabel('$\overline{\rho},\overline{U}+1$', 'interpreter','latex')
% ylabel('$z$', 'interpreter','latex')
% box on;
% 
% figure(3);  hold all;
% plot(rhobar(:,indt2d+50),z,'k-' ,'LineWidth',LineThickness)
% plot(ubar(:,indt2d+50)+1,z,'k--','LineWidth',LineThickness);
% hndl = legend('$H: \overline{\rho}(z)$',...
%               '$H: \overline{U}(z)+1$',...
%               '$K: \overline{\rho}(z)$',...
%               '$K: \overline{U}(z)+1$');
% set(hndl, 'interpreter','latex');  
% ylim([-5 5]);   xlim([-0.1 2.1]);
% if ic==1; gtext('$t =t_{2d}+100$','interpreter','latex'); end
% xlabel('$\overline{\rho},\overline{U}+1$', 'interpreter','latex')
% ylabel('$z$', 'interpreter','latex')
% box on;
% exportfig(gcf, 'fig/tst.eps', 'FontMode', 'fixed','FontSize', 18,'width',8,'height',8);

%=======
% plots of R, Rib 
%=======
% figure(4); hold all;
% plot(time,Iu./Irho,'k-' ,'LineWidth',LineThickness);
% hndl = legend('$H$','$K$');
% set(hndl, 'interpreter','latex');  
% xlabel('$time$', 'interpreter','latex')
% ylabel('$R=\ell_u/\ell_{\rho}$', 'interpreter','latex')
% box on;
% xlim([0 500]);
% scatter(t2d,Iu(indt2d)./Irho(indt2d),...
%     3*MarkerSize,'+',...
%     'LineWidth', 2,...
%     'MarkerEdgeColor','r');
% scatter(t3d,Iu(indt3d)./Irho(indt3d),...
%     3*MarkerSize,'*',...
%     'LineWidth', 2,...
%     'MarkerEdgeColor','r');
% scatter(t2d+100,Iu(indt2d+50)./Irho(indt2d+50),...
%     3*MarkerSize,'x',...
%     'LineWidth', 2,...
%     'MarkerEdgeColor','r');

% figure(5);  hold all;
% plot(time,Ri_I,'k-' ,'LineWidth',LineThickness)
% hndl = legend('$H$','$K$');
% set(hndl, 'interpreter','latex'); 
% xlabel('$time$', 'interpreter','latex')
% ylabel('$Ri_g(0)$', 'interpreter','latex')
% box on;
% xlim([0 500]);
% scatter(t2d,Ri_I(indt2d),...
%     3*MarkerSize,'+',...
%     'LineWidth', 2,...
%     'MarkerEdgeColor','r');
% scatter(t3d,Ri_I(indt3d),...
%     3*MarkerSize,'*',...
%     'LineWidth', 2,...
%     'MarkerEdgeColor','r');
% scatter(t2d+100,Ri_I(indt2d+50),...
%     3*MarkerSize,'x',...
%     'LineWidth', 2,...
%     'MarkerEdgeColor','r');
% 
% figure(6);  hold all;
% Rib = Ri_I./(Iu./Irho);
% plot(time,Rib,'k-' ,'LineWidth',LineThickness)
% hndl = legend('$H$','$K$');
% set(hndl, 'interpreter','latex'); 
% xlabel('$time$', 'interpreter','latex')
% ylabel('$Ri_b$', 'interpreter','latex')
% box on;
% xlim([0 500]);
% scatter(t2d,Rib(indt2d),...
%     3*MarkerSize,'+',...
%     'LineWidth', 2,...
%     'MarkerEdgeColor','r');
% scatter(t3d,Rib(indt3d),...
%     3*MarkerSize,'*',...
%     'LineWidth', 2,...
%     'MarkerEdgeColor','r');
% scatter(t2d+100,Rib(indt2d+50),...
%     3*MarkerSize,'x',...
%     'LineWidth', 2,...
%     'MarkerEdgeColor','r');
% 
% 
% figure(7);  hold all;
% Re_I = Iu./2/nu;
% plot(time,Re_I,'k-' ,'LineWidth',LineThickness)
% hndl = legend('$H$','$K$');
% set(hndl, 'interpreter','latex'); 
% xlabel('$time$', 'interpreter','latex')
% ylabel('$Re$', 'interpreter','latex')
% box on;
% xlim([0 500]);
% scatter(t2d,Re_I(indt2d),...
%     3*MarkerSize,'+',...
%     'LineWidth', 2,...
%     'MarkerEdgeColor','r');
% scatter(t3d,Re_I(indt3d),...
%     3*MarkerSize,'*',...
%     'LineWidth', 2,...
%     'MarkerEdgeColor','r');
% scatter(t2d+100,Re_I(indt2d+50),...
%     3*MarkerSize,'x',...
%     'LineWidth', 2,...
%     'MarkerEdgeColor','r');

%=======
% plots instantaneous and cumulative Mixing for low-Re and high-Re HWI
%=======
% figure(1); hold all; 
% xminplot = -150;        xmaxplot=250;
% M_cum = RPE-RPE(1)-Dpt;
% plot(time-t2d,M_cum./Dpt,'k-','LineWidth',LineThickness)
% plot(t3d-t2d,M_cum(indt3d)/Dpt(indt3d),'r*','MarkerSize',12)
% xlabel('$t$-$t_{2d}$', 'interpreter','latex')
% ylabel('$\mathcal{M}_c/D_p$', 'interpreter','latex')
% xlim([xminplot xmaxplot]);
% box on;
% 
% 
% figure(2); hold all; 
% xminplot = -150;        xmaxplot=250;
% plot(time-t2d,M./Dp,'k-','LineWidth',LineThickness)
% plot(t3d-t2d,M(indt3d)./Dp,'r*','MarkerSize',12)
% xlabel('$t$-$t_{2d}$', 'interpreter','latex')
% ylabel('$\mathcal{M}/D_p=K^*_{\rho}/\kappa $', 'interpreter','latex')
% xlim([xminplot xmaxplot]);
% box on;
% hndl = legend('$Re=\hspace{5.1pt}500:H05$-$hr~$',...
%                 '$Re=4000:H40~$',...
%                 '$Re=6000:H60~$');
% set(hndl, 'interpreter','latex');

%============ Removed plots
% figure(1); hold all;
% plot(time-t2d,Reb,'k-','LineWidth',LineThickness)
% plot(t3d-t2d,Reb(indt3d),'ko','MarkerFaceColor',[1 1 1])
% xlabel('$t-t_{2d}$', 'interpreter','latex')
% ylabel('$Re^*_b$', 'interpreter','latex')
% grid on; box on;
% xminplot = -150;        xmaxplot=300;
% xlim([xminplot xmaxplot]);
% % hndl = legend('$H_L:Re=~500$','$H~:Re=4000$');
% % set(hndl, 'interpreter','latex');  
% 
% M_acc = zeros(1,ntime);
% N2_acc = zeros(1,ntime);
% for i=2:ntime
%      M_acc(i) = trapz(time(1:i),M(1:i));
%      N2_acc(i)= trapz(time(1:i),N2avg(1:i));
% end
% for i=indt3d+1:ntime
%      M_acc(i) = trapz(time(indt3d:i),M(indt3d:i));
%      N2_acc(i)= trapz(time(indt3d:i),N2avg(indt3d:i));
% end
% kappa_acc  = M_acc./N2_acc;
% plot(time-t2d,M_acc,'k-','LineWidth',LineThickness)
% plot(t3d-t2d,M_acc(indt3d),'ko','MarkerFaceColor',[1 1 1])
% grid on; box on;
% ylabel('$\int_0^t \mathcal{M} \hspace{2pt} dt$', 'interpreter','latex')
% hndl = legend('$H$','$K$');
% set(hndl, 'interpreter','latex');  
% box on; grid on;
% hndl = legend('$H_0:Re=~500$','$H~:Re=4000$');
% set(hndl, 'interpreter','latex'); 




%% Spectral analysis of kinetic energy
% The following code has now been slightly modified for the purpose of SOC
% paper

%{
for i=1:5;
A = loadtxt([fadrs,'s.000000', num2str(i)','.dat'], 9, 2);
nkx = A(2,1);
nky = A(3,1);
ktime = A(4,1);
kh  = A(5,:)';
KEx = A(6,:)';
KEy = A(7,:)';
KEz = A(8,:)';
PEx = A(9,:)';      % not really PE, it is spectra of (rho')^2

     
delkx = 2*pi/Lx;
delky = 2*pi/Ly;

kx = 0:delkx:nkx*delkx;
ky = 0:delky:nky*delky;

[KX, KY] = meshgrid(kx,ky);
KH = reshape(kh,[nky+1 nkx+1]);
KEX = reshape(KEx,[nky+1 nkx+1]);
KEY = reshape(KEy,[nky+1 nkx+1]);
KEZ = reshape(KEz,[nky+1 nkx+1]);
PEX = reshape(PEx,[nky+1 nkx+1]);

% Define some wavenumber and time indecies
[~, indtime]= min(abs(time-ktime));
[~, indlkx]  = min(abs(kx-1./Lk_patch(indtime)));
[~, indlox]  = min(abs(kx-1./Lo_patch(indtime)));
[~, indlenx] = min(abs(kx-1./Len_patch(indtime)));

[~, indlky]  = min(abs(ky-1./Lk_patch(indtime)));
[~, indloy]  = min(abs(ky-1./Lo_patch(indtime)));
[~, indleny] = min(abs(ky-1./Len_patch(indtime)));

indknqx = min(indlkx+round(0.8*(nkx-indlkx)), nkx+1);
indknqy = min(indlky+round(0.8*(nky-indlky)), nky+1);

sxKEx = sum(KEX(:,1:indknqx),1)*delky;      % x-spectra of KEx
sxKEy = sum(KEY(:,1:indknqx),1)*delky;      % x-spectra of KEy
sxKEz = sum(KEZ(:,1:indknqx),1)*delky;      % x-spectra of KEz
sxPE  = sum(PEX(:,1:indknqx),1)*delky;      % x-spectra of KEz

syKEx = sum(KEX(1:indknqy,:),2)*delkx;      % y-spectra of KEx
syKEy = sum(KEY(1:indknqy,:),2)*delkx;      % y-spectra of KEy
fac = 1;
% delkh = KX./KH*delkx + KY./KH*delky;
% A(9,:)  = KX(:);
% A(10,:) = KY(:);
% A(11,:) = delkh(:);
% sortedA = sortrows(A',5);


% % Kx spectra
% figure;
% fac = 1;%Lo_patch(indtime);
% % scaling = 1/(sum(sxKEx)*delkx);
% scaling = 1;
% loglog(fac*kx(1:indknqx), sxKEx.*scaling,'k-','LineWidth',LineThickness);
% hold all;
% loglog(fac*kx(indlkx),sxKEx(indlkx)*scaling,'ks','LineWidth',1,'MarkerFaceColor',[1 1 1],'MarkerSize',12)
% loglog(fac*kx(indlox),sxKEx(indlox)*scaling,'ko','LineWidth',1,'MarkerFaceColor',[1 1 1],'MarkerSize',12)

% Kx compensated spectra
figure(1);
% scaling = 1/(sum(sxKEx)*delkx)*kx.^(5/3); %abs(D(indtime)).^(-2/3).*kx.^(5/3);
scaling = kx.^(5/3);%*0+1; %abs(D(indtime)).^(-2/3).*kx.^(5/3);
loglog(fac*kx(1:indknqx), sxKEx.*scaling(1:indknqx),'k-','LineWidth',LineThickness);   
hold all;
loglog(fac*kx(indlkx),sxKEx(indlkx)*scaling(indlkx),'ks','LineWidth',1,'MarkerFaceColor',[1 1 1],'MarkerSize',12)
loglog(fac*kx(indlox),sxKEx(indlox)*scaling(indlox),'ko','LineWidth',1,'MarkerFaceColor',[1 1 1],'MarkerSize',12)
loglog(fac*kx(indlenx),sxKEx(indlenx)*scaling(indlenx),'k*','LineWidth',1,'MarkerFaceColor',[1 1 1],'MarkerSize',12)


% figure;
% fac = 1;%Lo_patch(indtime);
% scaling = 1;
% loglog(fac*kx(1:indknqx), sxKEz.*scaling,'k-','LineWidth',LineThickness);
% hold all;
% loglog(fac*kx(indlkx),sxKEz(indlkx)*scaling,'ks','LineWidth',1,'MarkerFaceColor',[1 1 1],'MarkerSize',12)
% loglog(fac*kx(indlox),sxKEz(indlox)*scaling,'ko','LineWidth',1,'MarkerFaceColor',[1 1 1],'MarkerSize',12)

% Kz compensated spectra
figure(2);
scaling = kx;%*0+1;
loglog(fac*kx(1:indknqx), sxKEz.*scaling(1:indknqx),'k-','LineWidth',LineThickness);   
hold all;
loglog(fac*kx(indlkx),sxKEz(indlkx)*scaling(indlkx),'ks','LineWidth',1,'MarkerFaceColor',[1 1 1],'MarkerSize',12)
loglog(fac*kx(indlox),sxKEz(indlox)*scaling(indlox),'ko','LineWidth',1,'MarkerFaceColor',[1 1 1],'MarkerSize',12)
loglog(fac*kx(indlenx),sxKEz(indlenx)*scaling(indlenx),'k*','LineWidth',1,'MarkerFaceColor',[1 1 1],'MarkerSize',12)

end
%}

% hndl = legend('$H03:Re=~300$',...
%               '$H05:Re=~500$',...
%               '$H10:Re=1000~$',...
%               '$H20:Re=2000~$',...
%               '$H40:Re=4000~$',...
%               '$H60:Re=6000~$');
% set(hndl, 'interpreter','latex'); 

% 
% hndl = legend('$Re=\hspace{5.1pt}500:H05$-$hr$',...
%                 '$Re=4000:H40~$',...
%                 '$Re=6000:H60~$',...
%                 '$Re=4000:K$');
% set(hndl, 'interpreter','latex');
% 
% hndl = legend('$Re=~500:H05~$',...
%               '$Re=~500:H05$-$hr$');
% set(hndl, 'interpreter','latex');
% xlabel('$k_y$', 'interpreter','latex')
% xlabel('$k_x\ell_O$', 'interpreter','latex')
% ylabel('$\widehat{\mathcal{K}''_x}(k_x) \varepsilon^{{-}2/3}k_x^{5/3}$', 'interpreter','latex')
% ylabel('$\widehat{\mathcal{K}''_x}(k_x)k_x^{5/3}$', 'interpreter','latex')
% ylabel('$\widehat{\mathcal{K}''_x}(k_x)$', 'interpreter','latex')
% ylabel('$\widehat{\mathcal{K}''_x}(k_x)k_x^{5/3}/\mathcal{K}''_x$', 'interpreter','latex')
% ylabel('$\widehat{\mathcal{K}''_y}(k_y)/\mathcal{K}''_y$','interpreter','latex')
% xlim([0.1 100])
% exportfig(gcf, 'fig/tst.eps', 'FontMode', 'fixed','FontSize', 20,'width',8,'height',6.5);


%% SOC paper: Self-regulation Holmboe (part of that discussed ISSF/ICTAM)
%
ncnt        = 256;        % contour levels
ncmp        = 1024;        % colormap levels
nbin        = 512;         % number of bins
tmax = max(time)-t2d;
indz   = z>=-5 & z<=5;
indeps0= log10(epsbar(indz,:)/Dp/Pr)<0;%-6;
[TTIME, TZZ]=meshgrid(time-t2d,z(indz));
Ri_profile2 = Ri_profile(indz,:);
Ri_profile2(abs(Ri_profile2)>5) = NaN;
Ri_profile2(indeps0) = NaN  ;
if ic==1; profile_data_all = []; end;
aa = Ri_profile2(:,time>=t3d);
profile_data_all =[profile_data_all; aa(:)];

% hndl = legend('$~~R3$-$J016$','$~~R5$-$J016$','$R10$-$J016$','$R25$-$J016$',...
%     '$~~R5$-$J008$','$R10$-$J008$',...
%     '$~~R5$-$J032$','$R10$-$J032$');
% set(hndl, 'interpreter','latex');  


%{
% Holmboe burns, KH flares!
ncmp   = 1024;        % colormap levels
ncnt   = 256;
indz   = z>=-5 & z<=5;
tmax   = 200;
[TTIME, TZZ]=meshgrid(time-t2d,z(indz));
figure;
fld = abs(eps_prt(indz,:));%./(ones(length(find(indz)),1)*sum(eps_prt(indz,:)));
contourf(TTIME,TZZ,log10(fld),ncnt,'LineStyle','none')
hndl = colorbar;       caxis([-4.5 -2.5])
colormap(cubehelix(ncmp));
title(hndl, '$log(\varepsilon''_z)$','interpreter', 'latex');
xlabel('$t$-$t_{2d}$', 'interpreter','latex')
ylabel('$z$', 'interpreter','latex')
xlim([0 tmax])
%}

% figure(1); hold all;
% plot(time,Iu./Irho,'k-' ,'LineWidth',LineThickness);
% xlabel('$time$', 'interpreter','latex')
% ylabel('$R=\ell_u/\ell_{\rho}$', 'interpreter','latex')
% box on;
% xlim([0 350]);
% 
% figure(2);  hold all;
% plot(time,Irho,'k-' ,'LineWidth',LineThickness)
% xlabel('$time$', 'interpreter','latex')
% ylabel('$\ell_{\rho}$', 'interpreter','latex')
% box on;
% xlim([0 350]);
% 
% figure(3);  hold all;
% Ri0 = interp1(z,Ri_profile,0);
% plot(time-t3d,Ri0,'k-' ,'LineWidth',LineThickness);
% xlabel('$t$-$t_{3d}$', 'interpreter','latex')
% ylabel('$Ri_g(0)$', 'interpreter','latex')
% box on;
% xlim([-10 200]);
% 
% 
% figure;
% contourf(TTIME,TZZ,Ri_profile2,ncnt,'LineStyle','none')
% hndl = colorbar;       caxis([-1 1])
% colormap(brewermap([],'Spectral'))
% hold all;
% plot(time-t2d,[-Irho/2; Irho/2],'m'  ,'LineWidth',2);
% plot(time-t2d,[-Iu/2; Iu/2]    ,'--k','LineWidth',2);
% title(hndl, '$Ri_g$','interpreter', 'latex');
% xlabel('$t$-$t_{2d}$', 'interpreter','latex')
% ylabel('$z$', 'interpreter','latex')
% xlim([min(time-t2d) tmax])
% 
% 
% figure(51);  hold all;
% binmin = -1;     binmax = max(Ri_profile2(:));
% plotHist(Ri_profile2(:,time>=t3d),nbin,'r','w',binmin,binmax,0);
% xlabel('$Ri_g$', 'interpreter','latex')
% ylabel('PDF' , 'interpreter','latex')
% xlim([-0.5 1]);
% box on;
% xpatch = [0.2 0.25 0.25 0.2];
% ypatch = [0 0 8 8];
% % patch(xpatch,ypatch,[1 0.4 0.4],'LineStyle','None')
% colorCurve;


% if ic==ncase
% figure(61);  hold all;
% binmin = -1;     binmax = max(profile_data_all(:));
% plotHist(profile_data_all,nbin,'r','w',binmin,binmax,0);
% xlabel('$Ri_g$', 'interpreter','latex')
% ylabel('PDF' , 'interpreter','latex')
% xlim([-0.5 1]);
% box on;
% end

% figure(6);  hold all;
% RR = -(Shbar-XK-M)./(D3d+D2d);
% R_eq(ic) = mean(RR(ind));
% plot(Ri,R_eq(ic),'ko','MarkerFace',[0.75 0.75 0.75],'MarkerSize',8);
% xlabel('$Ri_g$', 'interpreter','latex')
% ylabel('$\mathcal{F}$' , 'interpreter','latex')
% box on;


%% SOC: Plot 2d histogram
%{
nbin_fld = 150;

fld = Ri_profile;
bin_fld   = linspace(0,0.75,nbin_fld);

%fld = M_sur./epsbar;


% fld = eps_prt;
% bin_fld   = linspace(-1,2,nbin_fld);
% bin_fld  = 10.^(linspace(-6,-2,nbin_fld));

time_eddy = 4;     % eddy turnover time
binmax = max(bin_fld(:));
binmin = min(bin_fld(:));
indz   = z>=-5 & z<=5;
indtime = time>=t2d;% & time<=trl;
indeps0= log10(epsbar(indz,:)/Dp/Pr)<0;
fld_profile = fld(indz,:);
fld_profile(fld_profile>binmax) = NaN;
fld_profile(fld_profile<binmin) = NaN;
fld_profile(indeps0) = NaN  ;
[TTIME, TZZ]=meshgrid(time(indtime)-min(time(indtime)),z(indz));
% TZZ = TZZ./(ones(length(z(indz)),1)*Irho(ind)/2);     % normalize by Irho
 
nbin_z   = round(nz/20);
nbin_time= round((max(TTIME(:))-min(TTIME(:)))/time_eddy);
bin_time = linspace(0,max(time(indtime))-min(time(indtime)),nbin_time);
bin_zz   = linspace(-5,5,nbin_z);
 

% plot 2D histogram
h1 = plot_histogram2d(TTIME,fld_profile(:,indtime),[],bin_time,bin_fld);
hold(h1,'on');  
plot(h1,[min(TTIME(:)), max(TTIME(:))], [0.25 0.25],'k--','LineWidth',2);


% colormap(h2,brewermap([],'*BuPu'))
% caxis(h1,[0 50])
% caxis(h2,[0 2])

% plot Ri_g values at the interface, upper and lower ends of Irho and Iu
fld0 = interp1(z,fld,0);
fld_topr = fld0*0;    fld_botr=fld0*0;
fld_topu = fld0*0;    fld_botu=fld0*0;
for i=1:ntime; fld_topr(i) = interp1(z,fld(:,i), Irho(i)/2); end
for i=1:ntime; fld_botr(i) = interp1(z,fld(:,i),-Irho(i)/2); end
for i=1:ntime; fld_topu(i) = interp1(z,fld(:,i), Iu(i)/2); end
for i=1:ntime; fld_botu(i) = interp1(z,fld(:,i),-Iu(i)/2); end

if (ic-1)*(ic-3)*(ic-5)==0
plot(h1,time(indtime)-min(time(indtime)),fld0(indtime),'.','MarkerSize',8,'MarkerEdgeColor',[0 0 0]);
plot(h1,time(indtime)-min(time(indtime)),fld_topr(indtime),'.','MarkerSize',14,'MarkerEdgeColor',[0 0 0]);
plot(h1,time(indtime)-min(time(indtime)),fld_topu(indtime),'o','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0.5 0.5 0.5]);
end
drawvline(t3d-t2d ,binmax,binmin)
drawvline(tfe-t2d ,binmax,binmin)
ylim([binmin binmax]);

% add labels 
hndl = colorbar;       caxis([0 50])
title(hndl, '$count$','interpreter', 'latex');
xlabel('$t$-$t_{2d}$', 'interpreter','latex')
ylabel('$Ri_g$', 'interpreter','latex')
% set(gca,'YScale','log');
% colormap(h1, brewermap([],'RdPu'))
colormap(h1,brewermap([],'PuBuGn'))
%}

%% SOC paper: Terminal velocity 
%%{
% figure(77); hold all;
% plot(time,(g*Irho).^0.5);
% xlabel('$time$', 'interpreter','latex')
% ylabel('$\mathcal{V} = \sqrt{g''\ell_\rho}$', 'interpreter','latex')
%plot(mean(Ri_avg(ind)), mean(eff(ind)),'o')
%plot(max(Ri_avg(ind)), max(eff(ind)),'s')
%%}

%% SOC paper: quasi-equilibrium
% Sh_prod = Shbar-XK;
% D_tot   = abs(D3d+D2d);
% figure(66);  hold all;
% %RR = -(Shbar-XK-M)./(D3d+D2d);
% RR = smooth(Sh_prod-M)./smooth(D_tot);
% plot(time-t3d,RR,'k','LineWidth',2);
% plot(trl-t3d, RR(indtrl), 'ko','MarkerSize',10,'MarkerFaceColor',[1 1 1]);
% % plot(time(ind)-t3d,RR(ind),'k','LineWidth',2);
% % plot(time(time>=trl)-t3d,RR(time>=trl),'k:','LineWidth',2);
% xlabel('$t$-$t_{3d}$', 'interpreter','latex')
% ylabel('$\mathcal{F}$' , 'interpreter','latex')
% xlim([0 time(end)-t3d]);
% box on;


%% SOC: 1/f analysis

%{
figure(22);
xx = linspace(-Lx/2,Lx/2,nxp);
nfreq=2^floor(log2(nxp/2));
f = 1/Lx*(0:nfreq-1);
% time_range = indt3d:indtrl;
time_range = [indt3d, indt3d+20, indtrl, indtrl+20, indtrl+40];
nspectra = length(time_range);
PSD = zeros(nspectra,nfreq);

figure; hold all;
for i=1:nspectra
    
    indtime = time_range(i);
    % interpolate into a evenly-spaced mesh
    yy_top = interp1(xsoc(:,indtime),vort3d_top(:,indtime),xx,'pchip');
    yy_bot = interp1(xsoc(:,indtime),vort3d_bot(:,indtime),xx,'pchip');
    
    % spectral analysis using fft
    YY_top = fft(sqrt(yy_top));
    YY_bot = fft(sqrt(yy_bot));    
    YY_top = YY_top(1:nfreq);
    YY_bot = YY_bot(1:nfreq);
    
    % compute power spectral density at top and bottom flanks
    PSD_top = YY_top.*conj(YY_top)/nxp;%     PSD_top = PSD_top/sum(PSD_top);
    PSD_bot = YY_bot.*conj(YY_bot)/nxp;%     PSD_bot = PSD_bot/sum(PSD_bot);
    PSD(i,:) = 0.5*(PSD_top+PSD_bot); 
    
    f_en = 1./Len_patch(indtime)/10;
    [~, ii] = min(abs(f-f_en));
    loglog(f(1:ii), PSD(i,1:ii),'k');   hold on
    set(gca, 'YScale', 'log','XScale', 'log');

%     hold off;
%     loglog(f, PSD/k); 
%     hold on;
%     plot(f, 1e1./f ,'r:');
% %     ylim([1e-8 1e0]);
%     title(['time = ', num2str(time(i)), ', t_{2d}=', num2str(t2d), ', t_{3d}=', num2str(t3d)])
%     pause(0.5); 
end
% spatial frequency of energy-containing eddies

% loglog(f(ii+1:end), PSD(ii+1:end)/k,'k:');
% loglog(f_ellison, PSD(ii)/k,'ko','MarkerFaceColor',[1 1 1]);
% loglog(f, 1e1./f ,'r:');
%}

%{ 
Not used anymore
[TTIME, TXX]=meshgrid(time,xsoc(:,1));
contourf(TXX,TTIME,vort2d_bot,128,'LineStyle','none')

wavelet decomposition
N = 16;                               % Decomposition Level 
wname = 'coif5';                     % Near symmetric wavelet
[c, s] = wavedec2(fld,N,wname);    % Multilevel 2-D wavelet decomposition.
opt = 'gbl';                    % Global threshold
%thr = 2;                       % Threshold
Z = 0.5*mean(fld);
thr = sqrt(4/3*Z*(log2(numel(fld))));
% thr = var(vorty_out(:))*sqrt(2*log(numel(vorty_out)));
sorh = 'h';                     % Hard thresholding
keepapp = 1;                    % Approximation coefficients cannot be thresholded
[fld_cmp,cden,sden,perf0,perfl2] = wdencmp(opt,c,s,wname,N,thr,sorh,keepapp);
perf0
perfl2

figure;
binrange = 10.^(linspace(-2,3,200));
% fld = [vort2d_top(:,ind);vort2d_bot(:,ind)];
hndl=histogram(fld,binrange); set(gca, 'YScale', 'log','XScale', 'log');
hndl.DisplayStyle = 'stairs';
hndl.Orientation = 'horizontal';
%}


%% SOC paper: mixing efficiency plots

%{
NOT USED
figure(4);  hold all;
plot(time,PE-PE(1)-Dpt,'k-' ,'LineWidth',LineThickness);  
xlabel('$time$', 'interpreter','latex')
ylabel('$M_c = \int_0^t M \hspace{2pt} dt$', 'interpreter','latex');
xlim([0 350]);
%}


% Ri_I_soc_avg(ic) = mean(Ri_I(ind));
% Ri_I_soc_std(ic) = std (Ri_I(ind));
% 
% Rib_soc_avg(ic) = mean(Ri_bul(ind));
% Rib_soc_std(ic) = std (Ri_bul(ind));
% if ic==8
%     figure(1); hold all;
%     errorbar(Ri_I_soc_avg,eff_avg,eff_err,'ko','MarkerFace',[0.5 0.25 0.5],'MarkerSize',10);
%     plot([0 2], [0.15 0.15], 'k--','LineWidth',2);
%     xlim([0 2]); ylim([0 0.4]);
%     xlabel('$Ri_{I}$', 'interpreter','latex')
%     ylabel('$E$', 'interpreter','latex')
% end


% figure(1); hold all;
% plot(Reb,eff,'o')
% figure(2); hold all;
% plot(Ri_I,eff,'o')
% figure(5); hold all;
% indtime = time>=t3d;% & time<=trl;
% Ri0 = interp1(z,Ri_profile,0);
% eff_sm = smooth(eff(indtime),'rlowess');
% myRi = Ri_avg;    % Ri_avg
% [effmx, indmx] = max(eff_sm);
% Rimx = myRi(indt3d+indmx-1);
% plot(myRi(indtime),eff_sm,'ko','MarkerFace',[1 1 1],'MarkerSize',8,'MarkerEdgeColor',[0 0 0]);
% plot(Rimx,effmx,'ko','MarkerFace',[0.5 0.5 0.5],'MarkerSize',8);
% xlabel('$Ri_b$', 'interpreter','latex')
% ylabel('$E$', 'interpreter','latex')
% 
% % % 
% figure(6); hold all;
% plot(Reb(indtime),eff_sm,'ko','MarkerFace',[1 1 1],'MarkerSize',8,'MarkerEdgeColor',[0 0 0]);
% plot(Reb(indt3d+indmx-1),effmx,'ko','MarkerFace',[0.5 0.5 0.5],'MarkerSize',8);
% xlabel('$Re_b$', 'interpreter','latex')
% ylabel('$E$', 'interpreter','latex')

%% SOC: definition of N2 and S2 (connection with Prt)
%
indz   = z>=-5 & z<=5;
indeps0= log10(epsbar(indz,:)/Dp/Pr)<0;%-6;
N2_profile = N2(indz,:);
N2_profile(indeps0) = NaN;
S2_profile = shear(indz,:).^2;
S2_profile(indeps0) = NaN;

N2_ts = max(0,mode(N2_profile)');   % time series
S2_ts = mode(S2_profile)';

%% SOC: cummulative time-dependent mixing efficiency
% stratified turbulence energy
%{
EST = KE3d+KE2d+APE;
D_cum = D3d*0;
Sh_prod_cum= Shbar*0;
Sh_prod = Shbar-XK;
BF = H;     % total buoyancy flux
BF_cum = BF*0;
S2avg = mean(shear.^2)';
N2_cum = N2_ts*0;   % time series
S2_cum = S2_ts*0;

for i=2:ntime;
    D_cum(i) = trapz(time(1:i),-D2d(1:i)-D3d(1:i));
    Sh_prod_cum(i) = trapz(time(1:i),Sh_prod(1:i));
    N2_cum(i) = trapz(time(1:i),N2_ts(1:i));
    S2_cum(i) = trapz(time(1:i),S2_ts(1:i));
    BF_cum(i) = trapz(time(1:i),BF(1:i));
end
M_cum = RPE-RPE(1)-Dpt;
Ri_cum = N2_cum./S2_cum;
eff_cum = M_cum./(M_cum + D_cum);
Prt_cum = Ri_cum./eff_cum;

% figure(110); hold all;
% plot(time-t2d,M_cum./D_cum);
% xlabel('$t-t_{2d}$', 'interpreter','latex')
% ylabel('$\Gamma_c$', 'interpreter','latex')
% 
% figure(111); hold all;
% plot(time-t2d,M_cum./(M_cum + D_cum));
% xlabel('$t-t_{2d}$', 'interpreter','latex')
% ylabel('$E_c$', 'interpreter','latex')
% 
% figure(112); hold all;
% plot(time-t2d,D_cum./Sh_prod_cum);
% xlabel('$t-t_{2d}$', 'interpreter','latex')
% ylabel('$\frac{\mathcal{D}_c}{\mathcal{P}}$', 'interpreter','latex')
% 
figure(113); hold all;
plot(time-t2d,Ri_cum);
xlabel('$t-t_{2d}$', 'interpreter','latex')
ylabel('$Ri^c_b$', 'interpreter','latex')
% 
figure(114); hold all;
plot(time-t2d,Prt_cum);
xlabel('$t-t_{2d}$', 'interpreter','latex')
ylabel('$Pr^c_t$', 'interpreter','latex')

% Assess E_ST energy balance equation 
% figure;
% plot(time, ddx(EST,time), 'k', time, Sh_prod - M + D2d+D3d, 'go');
% xlabel('$time$', 'interpreter','latex')
% hndl = legend(   '${d\over dt}E_{ST}$','$\mathbb{P} - \mathcal{M} - \epsilon''$');
% set(hndl, 'interpreter','latex');
% drawvline(tfe,2e-4,-2e-4)
% drawvline(t3d,2e-4,-2e-4)
% drawvline(t2d,2e-4,-2e-4)
% ylim([-2e-4, 2e-4])
% xlim([0 400])
% 
% plot components of E_st
% figure;
% plot(time, Sh_prod, time, M, time, abs(D2d+D3d));
% xlabel('$time$', 'interpreter','latex')
% hndl = legend( '$\mathbb{P}$', '$\mathcal{M}$', '$\epsilon''$');
% set(hndl, 'interpreter','latex');
% drawvline(tfe,2e-4,-2e-4)
% drawvline(t3d,2e-4,-2e-4)
% drawvline(t2d,2e-4,-2e-4)
% ylim([-2e-4, 2e-4])
% xlim([0 400])
%}

%% Check linear stability of the end state
% [~, indtime] = min(abs(time-t2d-100));  %(t2d+100);
% kappa_rho = kappa_star*kappa0+kappa0;
% myPrt = nu*kappa_m./kappa_rho;
% myRet = 1./nu./kappa_m;
% 
% stab.z   = zz(:,indtime);
% stab.rho = rhobar(:,indtime);
% stab.u   = ubar(:,indtime);
% stab.Re  = myRet(indtime);  %Iu(end)./2./nu;
% stab.Pr  = myPrt(indtime);  %Pr;
% stab.Rib = Ri_I(indtime)./(Iu(indtime)./Irho(indtime));
% stab.k   = 2*pi/Lx;
% save('../TG_solver/DNS_profile.mat','stab');

%% Save Lu and Lrho to in .dat files 
% Purpose: to be used by visit python driver to extract soc related
% quantities
% fid = fopen([fadrs,'flanks.dat'], 'w');   
% fprintf(fid, '%17.9e %17.9e %17.9e \r\n', [time';Irho;Iu]);
% fclose(fid);

%% SOC: plot of energy-containing length scales
%{
figure(201); hold all;
plot(time(ind)-t3d,Len_patch(ind),'k-','LineWidth',2)
plot(time(time>=trl)-t3d,Len_patch(time>=trl),'k:','LineWidth',2);
xlabel('$t-t_{3d}$', 'interpreter','latex')
ylabel('$\ell_{en}$', 'interpreter','latex')
xlim([0, 180]);
grid on; box on;
%}

%% SOC: cross spectrum of rho' and w'
% This analysis is added on July 10, 2018 and for the purpose of my SOC
% paper. Notice that the necessary data for this analysis must have been
% already prepared through the python code in pySrc which drives visits and
% calculates rho' and w' on various horizontal lines. 
% 
nfiles = length(dir([fadrs '/cc.*.dat']));
for i=2:nfiles;
    A = loadtxt([fadrs,'cc.000', num2str(i-1)','.dat'], 13, 1);
    
     % rho' at different z-levels
    rho_fluct_Irho_bot  = A(2,:)';   
    rho_fluct_Iu_bot    = A(3,:)';
    rho_fluct_interf    = A(4,:)';
    rho_fluct_Irho_top  = A(5,:)';   
    rho_fluct_Iu_top    = A(6,:)';
    
    % w' at different z-levels
    w_fluct_Irho_bot  = A(7,:)';   
    w_fluct_Iu_bot    = A(8,:)';
    w_fluct_interf    = A(9,:)';
    w_fluct_Irho_top  = A(10,:)';   
    w_fluct_Iu_top    = A(11,:)';
    cctime = A(12,1);
    xmesh = A(13,:)';
    nxp = size(A,2);

    % uniform mesho for interpolation
    xx = linspace(-Lx/2,Lx/2,nxp);
      

    window_type = hanning(nxp/3); % type of filtering window
    noverlap    = [];           % noverlap, defaul = length(window_type)/2
    nfft        = [];           % number of freq (default = max(256, 2^p), p=nextpow2(nxp/length(window_type))
    Fs          = nxp/(Lx/2/pi);% sampling frequency
   
    
    % 1) cross correlation at top flank
    figure(1);
    yy1 = interp1(xmesh,rho_fluct_Iu_top,xx,'pchip');
    yy2 = interp1(xmesh,w_fluct_Iu_top,xx,'pchip');
    [Pxy,F] = cpsd(yy1, yy2, window_type , noverlap, nfft,Fs);
    loglog(F, abs(Pxy))
    xlim([1, max(F)]);  ylim([1e-11, 1e-2]);    grid on;
    
    % add various wavenumbers associated with different wavelength to the
    % spectra.
    [~, indtime]= min(abs(time-cctime));
    [~, indlkx]  = min(abs(F-1./Lk_patch(indtime)));
    [~, indlox]  = min(abs(F-1./Lo_patch(indtime)));
    [~, indlenx] = min(abs(F-1./Len_patch(indtime)));
    hold all;
    loglog(F(indlkx),abs(Pxy(indlkx)),'ks')
    loglog(F(indlox),abs(Pxy(indlox)),'ko')
    loglog(F(indlenx),abs(Pxy(indlenx)),'k*')

    % 2) bottom flank
    figure(2)
    yy1 = interp1(xmesh,rho_fluct_Iu_bot,xx,'pchip');
    yy2 = interp1(xmesh,w_fluct_Iu_bot,xx,'pchip');
    [Pxy,F] = cpsd(yy1, yy2, window_type, noverlap, nfft, Fs);
    loglog(F, abs(Pxy))
    xlim([1, max(F)]);  ylim([1e-11, 1e-2]);     grid on;
    hold all;
    % add other wavenumbers
    loglog(F(indlkx),abs(Pxy(indlkx)),'ks')
    loglog(F(indlox),abs(Pxy(indlox)),'ko')
    loglog(F(indlenx),abs(Pxy(indlenx)),'k*')
    
    figure(3)
    yy1 = interp1(xmesh,rho_fluct_interf,xx,'pchip');
    yy2 = interp1(xmesh,w_fluct_interf,xx,'pchip');
    [Pxy,F] = cpsd(yy1, yy2, window_type, noverlap, nfft, Fs);
    loglog(F, abs(Pxy))
    xlim([1, max(F)]);  ylim([1e-11, 1e-2]);     grid on;
    hold all;
    % add other wavenumbers
    loglog(F(indlkx),abs(Pxy(indlkx)),'ks')
    loglog(F(indlox),abs(Pxy(indlox)),'ko')
    loglog(F(indlenx),abs(Pxy(indlenx)),'k*')
    


end

%% ML Prediction data
% Idea #1: given chi,eps and N2 --> predict eff!

% NOTE:
% The sorting procedure is not so accurate for laminar flows just at the 
% beginning and at the end of the simulation and thus leads to scattered 
% oscilations in efficiency. We need to smooth out those episodes by:
% (1) eff-->0 when apply a sigmoid function to the initial episode
% (2) eff=smooth(eff) when Reb_pert < 7


Reb_pert = abs(D2d+D3d)/nu./N2avg;
ind_end = Reb_pert<7 & time>=trl;
time_end = time(ind_end);
my_eff = eff./(1+exp(-(time-t2d/2)));
my_eff = max(my_eff,0);
my_eff(ind_end)=smooth(my_eff(ind_end),round(length(time_end)/10));

%% ML: Osbron-Cox method: generating prediction HWI data 
%{
if ic==1
    fid1 = fopen('prediction_data_features_chi_eps_N2.dat', 'w');   
    num_training_data = 0;
end
nz_profile = 512;
zft = linspace(min(z), max(z), nz_profile);
XX = 2*kappa0*chi_pr;    XX = max(0., XX);      % X  , scalar dissip.
DD = epsbar;             DD = max(0., DD);      % eps, turbul dissip.
for i=1:ntime;
    features = zeros(nz_profile, 4);
    ft1 = interp1(z, XX(:,i), zft);
    ft2 = interp1(z, DD(:,i), zft);
    ft3 = interp1(z, N2(:,i), zft);
    features(:,1:3) = [ft1', ft2', ft3']; 
    features(:,4) = my_eff(i);
    fprintf(fid1, '%17.9e, \t %17.9e, \t %17.9e, \t %17.9e \r\n', features');
end


num_training_data = num_training_data + ntime;
if (ic==ncase)
    fid2 = fopen('training_data_info.dat', 'w');  
    fprintf(fid2, '%6i, \t %6i \r\n', nz_profile, num_training_data);
    fclose(fid1);
    fclose(fid2);
end
%}


%% ML: Osbron-Cox method: generating NORMALIZED prediction HWI data 
% NOTE: chi is not normalied here (I need to understand how to do that!)
%{
if ic==1
    fid1 = fopen('prediction_data_normalized_features_chi_eps_N2.dat', 'w');   
    num_training_data = 0;
end
nz_profile = 512;
zft = linspace(min(z), max(z), nz_profile);
XX = 2*kappa0*chi_pr;    XX = max(0., XX);      % X  , scalar dissip.
DD = epsbar;             DD = max(0., DD);      % eps, turbul dissip.
for i=1:ntime;
    features = zeros(nz_profile, 4);
    ft1 = interp1(z, XX(:,i), zft);
    ft2 = interp1(z, DD(:,i)/Dp, zft);
    ft3 = interp1(z, N2(:,i)/N2avg(i), zft);
    features(:,1:3) = [ft1', ft2', ft3']; 
    features(:,4) = my_eff(i);
    fprintf(fid1, '%17.9e, \t %17.9e, \t %17.9e, \t %17.9e \r\n', features');
end


num_training_data = num_training_data + ntime;
if (ic==ncase)
    fid2 = fopen('prediction_data_info.dat', 'w');  
    fprintf(fid2, '%6i, \t %6i \r\n', nz_profile, num_training_data);
    fclose(fid1);
    fclose(fid2);
end
%}

%% Comparison b/w ML predictions and Cox-Osborn formula

% some usefull func handles
r2_score = @(label,prediction) ...
    1.0-sum((label-prediction).^2)/(sum((label-mean(label)).^2) + 1e-18); 
mean_square_err = @(label,prediction) mean((label-prediction).^2);

chi = 2*kappa0*chi_pr;      chi = max(0., chi);
B_cox = 0.5*mean(chi).*mean(N2)./mean(gradrho.^2);
dissp = mean(epsbar);
eff_cox = B_cox./(dissp + B_cox);
eff_true = my_eff';

if ic==1;    
    eff_true_all = [];   
    eff_cox_all  =[];
end
eff_true_all = [eff_true_all eff_true];
eff_cox_all = [eff_cox_all eff_cox];

% figure;
% plot(time,eff_true,time,eff_cox);



% Ri_cum(ic) = mean(N2avg)/mean(mean(shear.^2));
% delrho = Irho*0.;
% delU   = Iu*0.;
% for i=1:ntime; 
%     delrho(i) = interp1(z,rhobar(:,i), Irho(i)/2) - interp1(z,rhobar(:,i), -Irho(i)/2);
%     delU  (i) = interp1(z,  ubar(:,i), Irho(i)/2) - interp1(z,  ubar(:,i), -Irho(i)/2);
% end
% % negative sign is for delrho
% Ri_rho = -(g.*delrho./delU.^2).*Irho;
% figure(66); hold all;
% plot(time,Ri_rho);


% aa=xsoc(:,1);
% bb=vort3d_top(:,time_range(1));
% for i=time_range; 
%     aa=[aa;xsoc(2:end,i)+Lx*(i-1)]; 
%     bb=[bb;vort3d_top(2:end,i)];
% end
% plot(aa,bb)


%% Draw vertical lines for t2d, t3d
% drawvline(t2d ,ymaxplot,yminplot)
% drawvline(tfo ,ymaxplot,yminplot)
% drawvline(t3d ,ymaxplot,yminplot)
% drawvline(tfe ,ymaxplot,yminplot)
% % % 
% % % shadedpatch(tfo,tfe,yminplot,ymaxplot,'r')
% ylim([yminplot,ymaxplot])
% xlim([xminplot,xmaxplot])
% % 
% if xminplot<t2d && xmaxplot>t2d;        
%     gtext('$t^s_{2d}$', 'interpreter','latex'); 
% end
% if xminplot<tfo && xmaxplot>tfo;
%     gtext('$t^o_{f}$' , 'interpreter','latex'); 
% end
% if xminplot<t3d && xmaxplot>t3d;
%     gtext('$t^s_{3d}$', 'interpreter','latex'); 
% end
% if xminplot<tfe && xmaxplot>tfe;       
%     gtext('$t^e_{f}$' , 'interpreter','latex'); 
% end
% gtext(['$Pr$ =', num2str(Pr)],'interpreter','latex');
% gtext(['$Ri_0$ =', num2str(Ri)],'interpreter','latex');
% 

% addpath(genpath('/home/hesam/holm_workspace/colormaps/'))
% addpath(genpath('/home/hesam/software/MATLAB2010a/figuremaker/'))
% exportfig(gcf, 'fig2.eps', 'FontMode', 'fixed', 'FontSize', 10, 'color', 'cmyk', 'height', 4, 'width',5);
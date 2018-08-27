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
% 
% 
% figure(3);     
% semilogy(time,KE3d/KE(1));          hold on;
% semilogy(t3d ,(KE3d(indt3d)/KE(1)),'*')
% % semilogy(tfo ,(KE3d(indtfo)/KE(1)),'o','MarkerSize',4)
% % semilogy(tfe ,(KE3d(indtfe)/KE(1)),'s','MarkerSize',4)
% ylabel('$\mathcal{K}_{3d}/\mathcal{K}_0$', 'interpreter','latex');
% xlabel('t');

% figure(4);     plot(time,KEbar/KE(1));       hold on;
% ylabel('$\overline{\mathcal{K}}/\mathcal{K}_0$', 'interpreter','latex');



%% ploting the error between sigma_2d/3d LHS and RHS
% figure(1);     hold on;
% plot(time,sig_LHS,'ro','MarkerSize',2)
% plot(time,sig_RHS,'b-');
% ylabel('${d\over dt}\mathcal{K}$', 'interpreter','latex');
% xlabel('t');
% hndl = legend(   '${d\over dt}\mathcal{K}$',...
%                  '-$\mathcal{H} + \mathcal{D}$');
% set(hndl, 'interpreter','latex');
% grid on; box on;
% 
% % sigma3d
% figure(1);     hold on;
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
% figure(1);     hold on;
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
% yminplot = 0.;     ymaxplot = 0.35;
% xminplot = 50;     xmaxplot = 100;
% tstart = 0;
% plot(time-tstart,fac.*smooth(-sig3_Hb),'g');
% %plot(time-tstart,fac.*smooth(sig3_Rs),'m');
% plot(time-tstart,fac.*smooth(sig3_Sh),'b');
% plot(time-tstart,fac.*smooth(sig3_An),'r');
% hndl = legend(   '-$\hat\mathcal{H}_{3d}$', ...
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

%% Revised GRL parameterization and plots

% ========== Construction steps of the mixing efficiency param.============
% Find the right flank
% [effmx, indmx] = max(eff(ind));
% Rebmx = Reb(indt3d+indmx);
% ind_rf = time>=t3d & Reb>=Rebmx;
% ind_lf = time>=t3d & Reb<=Rebmx;


% % Fit RIGHT flank: Only DNS with `a' right flank
% % ======================
% %   eff/eff_peak = (1+2p)(Reb/Reb_peak)^p/(1+2p(Reb/Reb_peak)^(p+0.5))
% %   p = 1;
% % ======================
% % Step 1: Construct eff_peak and Reb_peak based on fit only to DNS data
% %         with right flank
% % Step 2: Take eff(Reb,Ri) in the limit of Reb->infty which gives:
% %         Psi = 3/2*eff_peak*sqrt(Reb_peak);
% % Step 3: Find Ri formulas for eff_peak and Psi_peak
% % Step 4: Find the best "p" value for other DNS data with the left flank
% %         using the formula for eff_peak and Reb_peak where:
% %         Reb_peak = 4/9*(Psi./eff_peak).^2;
% % FIT : eff/eff_peak = (1+2p)(Reb/Reb_peak)^p/(1+2p(Reb/Reb_peak)^(p+0.5))
% xdata = Reb(ind_rf);
% ydata = eff(ind_rf);
% fo = fitoptions('Method','NonlinearLeastSquares','start',100);
% ft = fittype('y0*(1+2*p)*(x/x0)^p/(1+2*p*(x/x0)^(p+0.5))',...
%     'problem',{'p','y0'},'options',fo);
% [curve1,gof] = fit(xdata,ydata,ft,'problem',{1, effmx});
% 
% % Plot left flank data against their parameterized versions:
% figure(1); hold all;
% plot(curve1,xdata,ydata,'ko');
% set(gca,'Xscale','log');

% if(ic==1); fit1_Ri=[];  fit_Rebp=[];  fit_effp=[];   end;
% fit1_Ri = [fit1_Ri  ; Ri_avg(ind_rf)];
% fit_Rebp =[fit_Rebp ; Ri_avg(ind_rf)*0+curve1.x0];
% fit_effp = [fit_effp; Ri_avg(ind_rf)*0+curve1.y0];
% fit_psi = 3/2*fit_effp.*sqrt(fit_Rebp);

% if(ic==1); fit2_Ri=[];  fit_effp=[]; end;
% fit_effp = [fit_effp ; effmx];
% fit2_Ri = [fit2_Ri ; Ri_avg(indt3d+indmx)];


% fit_Rebp(ic) = curve1.x0;
% fit_effp(ic) = curve1.y0;
% fit_Ri(ic)  = Ri_avg(indt3d+indmx);
% fit_p(ic)  = curve1.p;
% 


% 
% 
% Plot left flank data against their parameterized versions:
% xdata = Reb(ind);
% ydata = eff(ind);
% Rimx = Ri_avg(indt3d+indmx);
% Psi  = mixeff_param_Psi(Rimx);
% eff_peak = mixeff_param_effpeak(Rimx);
% Reb_peak = 4/9*(Psi./eff_peak).^2;
% figure(2); hold all;
% param = @(x,p)  (1+2*p)*x.^p./(1+2*p*x.^(p+0.5));
% plot(xdata,ydata,'o');
% plot(xdata,eff_peak*param(xdata/Reb_peak,0.55))


%% Old GRL parameterization

% ========== Construction steps of the mixing efficiency param.============
% % % Find the right flank
% [effmx, indmx] = max(eff(ind));
% Rebmx = Reb(indt3d+indmx);
% ind_rf = time>=t3d & Reb>=Rebmx;
% ind_lf = time>=t3d & Reb<=Rebmx;
% % tMarkerSize = -50./log10(Ri_avg(ind_rf));
% % 
% % Fit right flank: eff = c*Reb^(-0.5); Following Monismith et al
% xdata = Reb(ind_rf);
% ydata = eff(ind_rf);
% fo = fitoptions('Method','NonlinearLeastSquares');             
% ft = fittype('a*x^n','problem','n','options',fo);
% if (length(xdata) >=2)
%     [curve1,gof] = fit(xdata,ydata,ft,'problem',-0.5);
% else
%     curve1.a = NaN;
% end

% (Not required) Fit left flank: eff = a*Reb^b+c;
% xdata = Reb(ind_lf);
% ydata = eff(ind_lf);
% ft = fittype('power2');
% if (length(xdata) >=3); [curve2,gof] = fit(xdata,ydata,ft);     end;

% (If required) Plots the fits
% xx = linspace(min(Reb(ind)),max(Reb(ind)),length(Reb(ind))');
% figure;
% loglog(Reb(ind_rf),smooth(eff(ind_rf),'rlowess'),'ro');      hold all;
% loglog(Reb(ind_lf),smooth(eff(ind_lf),'rlowess'),'go');
% plot(curve1,'r');
% plot(curve2,'g');


% % Find the Ri dependence of Psi, eff=Psi(Ri)Reb^(-0.5)
% % NOTE: only data with 'a' right flank.
% if(ic==1); fit1_Ri=[]; fit_psi=[]; end;
% % fit_x = [fit_x ; mean(Ri_bul(ind_rf))];
% % fit_x = [fit_x ; mean(Ri_avg(ind_rf))];
% % fit_y = [fit_y ; curve1.a];
% fit1_Ri = [fit1_Ri ; Ri_avg(ind_rf)];
% fit_psi = [fit_psi ; Ri_avg(ind_rf)*0+curve1.a];

% % % Find possible relation b/w eff_peak and Ri
% % NOTE: all data
% if(ic==1); fit2_Ri=[]; fit_effp=[]; end;
% %fit_x = [fit_x ; Ri_bul(indmx)];
% fit2_Ri = [fit2_Ri ; Ri_avg(indt3d+indmx)];
% % fit2_Ri = [fit2_Ri ; Ri_I(indt3d+indmx)];
% fit_effp = [fit_effp ; effmx];
% % % [curve,gof] = fit(fit_x,fit_y,'rat23');


% Find the Reb-dependence of left flank: Theta(Reb)=eff/eff_peak
% NOTE: only using Pr=1 data as the scatter increases for higher Pr.
% Rimx = Ri_avg(indt3d+indmx);
% eff_peak = mixeff_param_effpeak(Rimx);
% Psi = mixeff_param_Psi(Rimx);
% Reb_peak = (Psi./eff_peak).^2;
% if(ic==1); fit_x=[]; fit_y=[]; end;
% fit_x = [fit_x ; Reb(ind_lf)./Reb_peak];
% fit_y = [fit_y ; eff(ind_lf)./eff_peak];


%% GRL plots
% ======== Figures of the multiparameter param. of mixing efficiency=======
% 
if Pr==1;   MarkerChar='o';      end;
if Pr~=1;   MarkerChar='s';      end;

% addpath(genpath('/home/hesam/kh-workspace/src/extra/'))
% addpath(genpath('/home/hesam/holm_workspace/colormaps/'))
% colormap(brewermap(ncmp,'*RdPu'))
% exportfig(gcf, 'fig/kappa_m_param.eps', 'FontMode', 'fixed','FontSize', 16,'width',8,'height',6);


% ncnt        = 256;        % contour levels
% ncmp        = 256;        % colormap levels
% ri=0.0:5e-4:1;
% reb=[2:0.01:20 20.1:0.1:300 350:50:1e4];
% [eff_par, Reb_par,Ri_par] = mixeff_param(reb,ri);
% Prt_par = Ri_par./eff_par;          
% Prt_par(eff_par<1e-4)=NaN;  Prt_par(Prt_par>100)=NaN;
% Psi = mixeff_param_Psi(Ri_par);
% eff_peak = mixeff_param_effpeak(Ri_par);
% Reb_peak = 4/9*(Psi./eff_peak).^2;
% 
% Contour plots of eff(Ri,Reb)
% figure;
% contourf(Reb_par,Ri_par,eff_par,ncnt,'LineStyle','none');
% % contour(Reb_par,Ri_par,eff_par,ncnt);
% set(gca, 'XScale', 'log');
% hold on;
% plot(Reb_peak(:,1),Ri_par(:,1),'--w','LineWidth',2);
% xlabel('$Re_b$', 'interpreter','latex')
% ylabel('$Ri$', 'interpreter','latex')
% colormap(paruly(ncmp));
% hndl = colorbar;
% title(hndl, '$\mathcal{E}$','interpreter', 'latex');
% % set(gca,'view',[48,8]);
% % c = smooth(eff(ind),'rlowess');
% % scatter(Reb(ind),Ri_bul(ind),72,...
% %         c,'filled', 'o');
% % caxis([0 1/3]);


% % Contour plots of Pr_t(Ri,Reb)
% figure;
% contourf(Reb_par,Ri_par,Prt_par,ncnt,'LineStyle','none');
% set(gca, 'XScale', 'log');
% xlabel('$Re_b$', 'interpreter','latex')
% ylabel('$Ri$', 'interpreter','latex')
% hndl = colorbar;
% title(hndl, '$Pr_t$','interpreter', 'latex');
% colormap(brewermap(ncmp,'*PuBuGn'))
% hold on;
% plot(Reb_peak(:,1),Ri_par(:,1),'--k','LineWidth',2);
% caxis([1 10])

% % K_rho: plot DNS data and compare them with their parameterized equivalents.
% figure(3);
% hold all;
% eff_DNS     = smooth(eff(ind),'rlowess');
% [eff_DNS_par, Reb_DNS_par,Ri_DNS_par] = mixeff_param(Reb(ind),Ri_avg(ind));
% Gamma_DNS     = eff_DNS./(1-eff_DNS);
% Gamma_DNS_par = eff_DNS_par./(1-eff_DNS_par);
% kappa_DNS     = Pr*Gamma_DNS.*Reb(ind);
% kappa_DNS_par = Pr*Gamma_DNS_par.*Reb_DNS_par;
% % if (ic==1); kappa_rmse=0;  nn=0;   end;
% % kappa_rmse = kappa_rmse + sum((kappa_DNS-kappa_DNS_par).^2);
% % nn = nn + length(kappa_DNS);
% % if (ic==ncase); kappa_rmse = sqrt(kappa_rmse/nn);    end;    
% % scatter(Gamma_DNS,diag(Gamma_DNS_par), [], MarkerChar, ...
% scatter(kappa_DNS,diag(kappa_DNS_par), [], MarkerChar, ...
%      'MarkerEdgeColor',[0 0 0],...
%      'MarkerFaceColor',[1 1 1]);
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');
% xlabel('$K^{dns}_{\rho}/\kappa$', 'interpreter','latex')
% ylabel('$K^{par}_{\rho}/\kappa$', 'interpreter','latex')
% xminplot = 1;      xmaxplot = 1e3;
% yminplot = 1;      ymaxplot = 1e3;
% % axis equal; box on;
% % loglog([0.1 2e3], [0.1 2e3],'k-')
% % 
% % 
% % % K_m: plot DNS data and compare them with their parameterized equivalents.
% figure(4);
% hold all;
% kappa_m_DNS     = Ri_avg(ind).*Reb(ind)./(1-eff_DNS);
% kappa_m_DNS_par = Ri_DNS_par.*Reb_DNS_par./(1-eff_DNS_par);
% 
% scatter(kappa_m_DNS,diag(kappa_m_DNS_par), [], MarkerChar, ...
%      'MarkerEdgeColor',[0 0 0],...
%      'MarkerFaceColor',[1 1 1]);
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');
% xlabel('$K^{dns}_{m}/\nu$', 'interpreter','latex')
% ylabel('$K^{par}_{m}/\nu$', 'interpreter','latex')
% xminplot = 1;      xmaxplot = 1e3;
% yminplot = 1;      ymaxplot = 1e3;
% % axis equal; box on;
% % loglog([0.1 2e3], [0.1 2e3],'k-')


% % Psi and Eff_peakc figure
% figure(5);
% [haxes,hline1,hline2] = plotyy(fit2_Ri,fit_effp,fit1_Ri,fit_psi);
% hold(haxes(1),'on');
% hold(haxes(2),'on');
% 
% ri=0.0:0.0001:1;
% Psi = mixeff_param_Psi(ri);
% eff_peak = mixeff_param_effpeak(ri);
% plot(haxes(2),ri,Psi,'k--','LineWidth',1);
% plot(haxes(1),ri,eff_peak,'k--','LineWidth',1);
% set(hline2,'LineStyle','None', 'Marker', 'o','MarkerSize', 8, ...
%     'MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerFaceColor',[0.5 0.5 0.5]);
% set(hline1,'LineStyle','None', 'Marker','p','MarkerSize', 8, ...
%     'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor',[0 0 0]);
% set(haxes(2),...
%     'xscale','log',...
%     'xlim',[5e-3 1e0],...
%     'ylim',[0 16],...
%     'YTick',0:2:16,...
%     'Ycolor','k');
% set(haxes(1),...
%     'xscale','log',...
%     'xlim',[5e-3 1e0],...
%     'ylim',[0 0.4],...
%     'YTick',0:0.05:0.4,...
%     'box', 'on',...
%     'Ycolor','k');   
% xlabel(haxes(1),'$Ri$','interpreter','latex');
% ylabel(haxes(1),'$E^{\star}$','interpreter','latex');
% ylabel(haxes(2),'$\Psi$','interpreter','latex');
% 
% % creat the uncertainty patch
% x = [0.46 1 1 0.46];
% y = [0 0 0.4 0.4];
% ptch_hndl = patch(x,y,[1 0.4 0.4],'Parent', haxes(1), 'LineStyle','None');
% uistack(ptch_hndl,'bottom');


%============= Old Reb-based Pr_t plot
% Prt parametrization
% reb=[2:0.01:20 20.1:0.1:300 350:50:1e4];
% a = 167.4;    b = -2.9;      c=3.6;        d=-0.2;
% Prt_par = a*reb.^b+c*reb.^d;
% Prt_par(reb>10^3) = 1;
% 
% figure (5);
% scatter(Reb(ind),Ri_avg(ind)./smooth(eff(ind),'rlowess'),...
%      48, MarkerChar, ...
%     'MarkerFaceColor', [1 1 1],...
%     'MarkerEdgeColor', [0 0 0]);
% set(gca, 'XScale', 'log');
% hold all;
% xminplot = 1;         xmaxplot = 1e4;
% yminplot = 0;         ymaxplot = 10;
% xlabel('$Re_b$', 'interpreter','latex')
% ylabel('$Pr_t = K_m/K_{\rho}$', 'interpreter','latex')
% box on;
% if (ic==ncase); plot(reb,Prt_par,'--r','LineWidth',2); end;

% ============== Old plots of Psi and E_star (individual)
% figure(5); hold all;
% plot(fit_x,fit_y,'p','MarkerSize', 8, ...
%     'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor',[0 0 0]);
% ri=0.0:0.001:0.5;
% eff_peak = mixeff_param_effpeak(ri);
% plot(ri,eff_peak,'k--');
% xlim([0 0.5]);    ylim([0 0.4]);
% box on; grid on;
% xlabel('$Ri$', 'interpreter','latex')
% ylabel('$\mathcal{E}^{\star}$', 'interpreter','latex')


% figure(6); hold all;
% x = [0.44 0.5 0.5 0.44];
% y = [0 0 15 15];
% patch(x,y,[1 0.4 0.4],'LineStyle','None')
% plot(fit_x,fit_y,'o','MarkerSize', 8, ...
%     'MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerFaceColor',[0.5 0.5 0.5]);
% ri=0.0:0.001:0.5;
% Psi = mixeff_param_Psi(ri);
% plot(ri,Psi,'k--','LineWidth',1);
% xlim([0 0.5]);    ylim([0 14]);
% box on; grid on;
% xlabel('$Ri$', 'interpreter','latex')
% ylabel('$\Psi = \mathcal{E}Re_b^{0.5}$', 'interpreter','latex')


% ========= Supplimentary plots
% figure (7);
% scatter(Reb(ind),Ri_avg(ind),48,eff(ind),'filled');
% hold all;
% set(gca, 'XScale', 'log');
% xlabel('$Re_b$', 'interpreter','latex')
% ylabel('$Ri$', 'interpreter','latex')
% colormap(paruly);
% hndl = colorbar;
% title(hndl, '$E$','interpreter', 'latex');
% caxis([0 1/3]);
% box on;
% % ri=0.0:5e-4:1;
% % Psi = mixeff_param_Psi(ri);
% % eff_peak = mixeff_param_effpeak(ri);
% % Reb_peak = 4/9*(Psi./eff_peak).^2;
% % plot(Reb_peak,ri,'--k');

% 
% figure (1); hold all;
% MarkerSize      = -50./log10(Ri_avg(ind));
% scatter3(Ri_bul(ind),Reb(ind),smooth(eff(ind),'rlowess'),...
%     MarkerSize, MarkerChar, ...
%     'MarkerEdgeColor',[0 0 0],...
%     'MarkerFaceColor',[1 1 1]);
% hold all;
% set(gca, 'YScale', 'log');
% set(gca, 'ZScale', 'log');
% xlabel('$Ri_b$', 'interpreter','latex')
% ylabel('$Re_b$', 'interpreter','latex')
% zlabel('$\mathcal{E}$', 'interpreter','latex')
% box on;

%% Paper with Ali
if(ic==1); DNSdata_Reb=[]; DNSdata_Gam=[]; DNSdata_M=[]; end;
Gam = eff./(1-eff);

% DNSdata_Reb = [DNSdata_Reb ; Reb(ind)];
% DNSdata_Gam = [DNSdata_Gam ; Gam(ind)];
% DNSdata_M = [DNSdata_M ; M(ind)];
% DNS.old.M = DNSdata_M;
% DNS.old.Gam = DNSdata_Gam;
% DNS.old.Reb = DNSdata_Reb;

% myind = eff>=0 & time<t3d;
% DNSdata_Reb = [DNSdata_Reb ; Reb(myind)];
% DNSdata_Gam = [DNSdata_Gam ; Gam(myind)];
% DNSdata_M = [DNSdata_M ; M(myind)];
% DNS.young.M = DNSdata_M;
% DNS.young.Gam = DNSdata_Gam;
% DNS.young.Reb = DNSdata_Reb;

% myind = eff>=0;
% DNSdata_Reb = [DNSdata_Reb ; Reb(myind)];
% DNSdata_Gam = [DNSdata_Gam ; Gam(myind)];
% DNSdata_M = [DNSdata_M ; M(myind)];
% DNS.all.M = DNSdata_M;
% DNS.all.Gam = DNSdata_Gam;
% DNS.all.Reb = DNSdata_Reb;
% plotter;
% 


%% SOC: for mixing efficiency vs Ri plot
% %  using all KH DNS except Pr=16

%{ 
NOt USED ANY MORE
% figure(5); hold all;
% [effmx, indmx] = max(eff(ind));
% Rimx = Ri_I(indt3d+indmx-1);
% plot(Rimx,effmx,'bp','MarkerFace',[0.68 0.92 1],'MarkerSize',10);
%}

% Ri_I_soc_avg(ic) = mean(Ri_I(ind));
% Ri_I_soc_std(ic) = std (Ri_I(ind));
% 
% if ic==27
%     figure(1); hold all;
%     errorbar(Ri_I_soc_avg,eff_avg,eff_err,'ko', 'MarkerFace',[0.68 0.92 1],'MarkerSize',10);
%     plot([0 2], [0.15 0.15], 'k--','LineWidth',2);
%     xlim([0 2]); ylim([0 0.4]);
%     xlabel('$Ri_{I}$', 'interpreter','latex')
%     ylabel('$E$', 'interpreter','latex')
% end
% 
% % Terminal velocity and quasi-equilibrium
% figure(77); hold all;
% plot(time,(g*Irho).^0.5);
% xlabel('$time$', 'interpreter','latex')
% ylabel('$\mathcal{V} = \sqrt{g''\ell_\rho}$', 'interpreter','latex')
% %plot(mean(Ri_avg(ind)), mean(eff(ind)),'o')
% %plot(max(Ri_avg(ind)), max(eff(ind)),'s')


%% SOC: Plot 2d histogram

%{
M_sur = chi_pr./abs(gradrhob)*Ri/R/Re/Pr - Dp*Lz;
% fld = Ri_profile;
% bin_fld   = linspace(0,0.75,nbin_fld);
fld = M_sur./epsbar;
%bin_fld   = linspace(-1,2,nbin_fld);
nbin_fld = 150;
bin_fld  = 10.^(linspace(-4,0,nbin_fld));

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
% hold(h2,'on');  plot(h2,[min(TTIME(:)), max(TTIME(:))], [0.25 0.25],'k--');
% colormap(h1,brewermap([],'PuBuGn'))
hndl = colorbar;       caxis([0 50])
title(hndl, '$count$','interpreter', 'latex');
xlabel('$t$-$t_{2d}$', 'interpreter','latex')
ylabel('$Ri_g$', 'interpreter','latex')

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

% if (ic-1)*(ic-3)*(ic-5)==0
plot(h1,time(indtime)-min(time(indtime)),fld0(indtime),'d','MarkerSize',8,'MarkerEdgeColor',[1 1 1]);
plot(h1,time(indtime)-min(time(indtime)),fld_topr(indtime),'s','MarkerSize',8,'MarkerEdgeColor',[1 1 1]);
plot(h1,time(indtime)-min(time(indtime)),fld_topu(indtime),'o','MarkerSize',8,'MarkerEdgeColor',[1 1 1]);
% end
drawvline(t3d-t2d ,0.75,0)
drawvline(tfe-t2d ,0.75,0)
ylim([binmin binmax]);
%}

%% SOC: Self-regulation Holmboe (see PlotData.m in holm_workspace)

% tmax = max(time)-t2d;
% indz   = z>=-5 & z<=5;
% [TTIME, TZZ]=meshgrid(time-t2d,z(indz));
% Ri_profile2 = N2(indz,:)./shear(indz,:).^2;
% Ri_profile2(abs(Ri_profile2)>5) = NaN;
% indeps0= log10(epsbar(indz,:))<-6;
% Ri_profile2(indeps0) = NaN  ;
% 
% figure(5);  hold all;
% binmin = -1;     binmax = max(Ri_profile2(:));
% plotHist(Ri_profile2(:,ind),256,'r','w',binmin,binmax,0);
% xlabel('$Ri_g$', 'interpreter','latex')
% ylabel('PDF' , 'interpreter','latex')
% xlim([-0.5 1]);
% box on;

% figure(6);  hold all;
% RR = -(Shbar-XK-M)./(D3d+D2d);
% R_eq(ic) = mean(RR(ind));
% plot(Ri,R_eq(ic),'kp','MarkerFace',[0.75 0.75 0.75],'MarkerSize',8);
% xlabel('$Ri_g(0)$', 'interpreter','latex')
% ylabel('\mathcal{R}' , 'interpreter','latex')
% box on;

%% SOC: cummulative time-dependent mixing efficiency
% D_cum = D3d*0;
% Sh_cum= Shbar*0;
% Sh_prod = Shbar-XK;
% for i=2:ntime;
%     D_cum(i) = trapz(time(1:i),-D2d(1:i)-D3d(1:i));
%     Sh_cum(i) = trapz(time(1:i),Sh_prod(1:i));
% end
% M_cum = RPE-RPE(1)-Dpt;
% figure(110); hold all;
% plot(time-t2d,M_cum./D_cum);
% 
% figure(112); hold all;
% % stratified turbulence energy
% EST = KE3d+KE2d+APE;
% % plot(time-t2d,(Sh_cum-M_cum-D_cum));%,time-t2d,,'ro')
% %plot(time-t2d,(EST-EST(1))./D_cum, time-t2d,(Sh_cum-M_cum-D_cum)./D_cum,'ro')
% plot(time-t2d,Sh_cum./D_cum);


%% ML: Osbron-Cox method: generating training data 
%{
if ic==1
    fid1 = fopen('training_data_features_chi_eps_N2.dat', 'w');   
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

%% ML: Osbron-Cox method: generating NORMALIZED training data 
% NOTE: chi is not normalied here (I need to understand how to do that!)
%{
if ic==1
    fid1 = fopen('training_data_normalized_features_chi_eps_N2.dat', 'w');   
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
    fid2 = fopen('training_data_info.dat', 'w');  
    fprintf(fid2, '%6i, \t %6i \r\n', nz_profile, num_training_data);
    fclose(fid1);
    fclose(fid2);
end
%}

%% ML: Osborn-Cox method: Increasing the size of training data by adding random gaussian noise
% NOT used yet!
%{
num_rand_profile = 10;
amp_noise = 0.01;
if ic==1
    fid1 = fopen('training_data_extra_features_chi_eps_N2.dat', 'w');   
    num_training_data = 0;
end
nz_profile = 512;
zft = linspace(min(z), max(z), nz_profile);
XX = 2*kappa0*chi_pr;    XX = max(0., XX);      % X  , scalar dissip.
DD = epsbar;             DD = max(0., DD);      % eps, turbul dissip.
for i=1:ntime;
    features = zeros(nz_profile, 4);
    features_extra = zeros(nz_profile*num_rand_profile, 4);
    ft1 = interp1(z, XX(:,i), zft);
    ft2 = interp1(z, DD(:,i), zft);
    ft3 = interp1(z, N2(:,i), zft);
    features(:,1:3) = [ft1', ft2', ft3']; 
    features(:,4) = my_eff(i);
    for nr=1:num_rand_profile
        ft1_extra = add_rand_lognormal_noise(ft1);
        ft2_extra = add_rand_lognormal_noise(ft2);
        ft3_extra = add_rand_lognormal_noise(ft3);
        ft4_extra = zeros(1,nz_profile) + my_eff(i);
        istart = 1+(nr-1)*nz_profile;
        iend   = istart + nz_profile - 1;
        features_extra(istart:iend,:) = [ft1_extra', ft2_extra', ft3_extra', ft4_extra']; 
    end
    fprintf(fid1, '%17.9e, \t %17.9e, \t %17.9e, \t %17.9e \r\n', [features', features_extra']);
end

num_training_data = num_training_data + ntime;
if (ic==ncase)
    fid2 = fopen('training_data_info.dat', 'w');  
    fprintf(fid2, '%6i, \t %6i \r\n', nz_profile, num_training_data);
    fclose(fid1);
    fclose(fid2);
end
%}


%% ML: Osbron-Cox method: alternative data file for nn_mixing (may not end up being useful)
%{ 
if ic==1
    fid1 = fopen('training_data_features_N2_chi_mixing_eps.dat', 'w');   
    num_training_data = 0;
end
nz_profile = 512;
zft = linspace(min(z), max(z), nz_profile);
for i=1:ntime;
    features = zeros(nz_profile, 4);
    ft1 = interp1(z, N2(:,i), zft);
    ft2 = interp1(z, 2*kappa0*chi_pr(:,i), zft);
    features(:,1:2) = [ft1', ft2']; 
    features(:,3) = mixing(i);
    features(:,4) = dissipation(i);
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

%% dl_Osborn routines

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
my_eff(ind_end)=smooth(my_eff(ind_end),round(length(time_end)));



%% dl_Osborn: Comparison b/w ML predictions and Cox-Osborn formula
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

%% dl_Osborn: generate the training data
epsbar_turb = epsbar - shear.^2./Re;

dl_fname = 'KHI_training_data.dat';
% dl_fname = 'KHI_testing_data.dat';
dl_osborn_write_train(dl_fname, my_eff, epsbar,N2, kappa0,z,shl_zmax,shl_zmin)

%% test
% [effmx, indmx] = max(eff(ind));
% Rimx = Ri_I(indt3d+indmx);
% Fr_t = (Reb.^2.*Iu'/(2*g*Re^2*Lx^4)).^(1/6);
% Re_t = (2*Re^2*Reb*g*Lx^4./Iu').^(1/3);
% figure(1); hold all;
% scatter(Re_t(ind),1./Fr_t(ind),...
%      48, MarkerChar, ...
%     'MarkerFaceColor', [1 1 1],...
%     'MarkerEdgeColor', [0 0 0]);
% plot(Re_t(indt3d+indmx),1/Fr_t(indt3d+indmx),'rp','MarkerFace',[1 0.5 0.5],'MarkerSize',12);
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');

% figure(1); hold all;
% plot(Re,Ri,'ko','MarkerSize',12);
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');

% figure(12); hold all;
% plot(-D(ind),M(ind));


% indz   = z>=-20 & z<=20;
% Ri_profile2 = N2(indz,ind)./shear(indz,ind).^2;
% indeps0= log10(epsbar(indz,ind))<-6;
% Ri_profile2(indeps0) = NaN  ;
% Ri_profile2(abs(Ri_profile2)>1) = NaN;
% figure(22); hold all;
% plotHist(Ri_profile2(:),200,'r','w',0);
% % xlabel('$Ri_g$', 'interpreter','latex')
% % ylabel('probability' , 'interpreter','latex')
% xlim([-0.5 1]);
% grid on; box on;
% % colorCurve;


% ind = time>=t3d & time<=trl;
% ind2=indt3d-1;
% time_eddy = 10;     % eddy turnover time
% while (ind2<=ntime-time_eddy);
%     ind1 = ind2+1;
%     ind2 = ind1+time_eddy;
%     Ri_profile2 = N2(indz,ind1:ind2)./shear(indz,ind1:ind2).^2;
%     indeps0= log10(epsbar(indz,ind1:ind2))<-6;
%     Ri_profile2(indeps0) = NaN  ;
%     Ri_profile2(abs(Ri_profile2)>1) = NaN;
%     figure(22); hold all;
%     plotHist(Ri_profile2(:),100,'r','w',0);
% %     plot(mean(RR(ind1:ind2)),nanmean(Ri_profile2(:)),'o');
%     colorCurve;
%     xlim([-0.5 1]);%    ylim([0 0.5]);
%     grid on; box on;
%     pause(0.5);
% end

% figure;
% tmax = max(time)-t2d;
% [TTIME, TZZ]=meshgrid(time-t2d,z(indz));
% Ri_profile2 = N2(indz,:)./shear(indz,:).^2;
% Ri_profile2(abs(Ri_profile2)>2) = NaN;
% indeps0= log10(epsbar(indz,:))<-6;
% Ri_profile2(indeps0) = NaN  ;
% contourf(TTIME,TZZ,Ri_profile2,128,'LineStyle','none')
% hndl = colorbar;       caxis([-1 1])
% % colormap(brewermap([],'RdBu'))
% colormap(brewermap([],'Spectral'))
% hold all;
% plot(time-t2d,[-Irho/2; Irho/2],'m'  ,'LineWidth',2);
% plot(time-t2d,[-Iu/2; Iu/2]    ,'--k','LineWidth',2);
% title(hndl, '$\overline{Ri}_g$','interpreter', 'latex');
% xlabel('$t$-$t_{2d}$', 'interpreter','latex')
% ylabel('$z$', 'interpreter','latex')
% xlim([min(time-t2d) tmax])
% 
% figure;
% contourf(TTIME,TZZ,N2(indz,:),128,'LineStyle','none')
% hndl = colorbar;%       caxis([0 1])
% colormap(cubehelix(128));
% title(hndl, '$\overline{N^2}$','interpreter', 'latex');
% xlabel('$t$-$t_{2d}$', 'interpreter','latex')
% ylabel('$z$', 'interpreter','latex')
% xlim([min(time-t2d) tmax])


% ri=0.0:5e-4:1;
% nu0=50e-4;
% kappa_KPP=nu0*(1-(ri./0.7).^2).^3;
% kappa_KPP(ri>0.7)=0;
% figure(1); hold all;
% myRi = N2avg./(mean(shear,1).^2)'.*Iu'/Lz;
% myRib = g*Iu/2;
% myRig0 = myRib.*Iu./Irho;   % same as Ri_I

% scatter(myRi(ind), kappa_star(ind)/Pr*1e-6,48,log10(Reb(ind)),'filled','o')
% plot(Ri_I(ind),Ri_avg(ind),'o')
% Ri_bul = N2avg./(mean(shear,1).^2)'.*avg_fac;
% plot(Ri_bul, Ri_avg,'o')
% scatter(Reb(ind),smooth(eff(ind),'rlowess'),72,Ri_avg(ind),'filled', 'o');
% set(gca,'Xscale','log');

% figure(2); hold all;
% % plot(time(time>=t2d)-t3d,Lo(time>=t2d)./LT3d(time>=t2d),'o');
% indyoung = time>=t2d & time<t3d;
% ReT3d = (LT3d./Lk).^(4/3);
% plot(Reb(indyoung),ReT3d(indyoung),'o');

% figure(1); hold all;
% Fr_t = (Lo_patch./LT3d_patch).^(2/3);
% Re_t = (LT3d_patch./Lk_patch).^(4/3);
% scatter(Reb(ind),Fr_t(ind),48,eff(ind),'filled');
%


% scatter(Reb(ind),1./Fr_t(ind).^2,48,eff(ind),'filled');
% loglog(Reb(~ind),eff(~ind)./(1-eff(~ind)),'k.');    hold on;
% loglog(Reb(ind) ,eff(ind)./(1-eff(ind)),'ko','MarkerSize',6,'MarkerFaceColor','w');

% nbins is arbitrary: number of patches in each bin need to be reported

%figure(2); hold all;
%plot(Reb,Lo./LT3d)

% figure(1); hold all;
% plot(Reb(ind),eff(ind))
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');
% 
% figure(2); hold all;
% plot(Reb(ind),M(ind))
% %set(gca, 'XScale', 'log');
% %set(gca, 'YScale', 'log');
% 
% figure(3); hold all;
% plot(Reb(ind),Lo(ind)./LT3d(ind))
%set(gca, 'XScale', 'log');
%set(gca, 'YScale', 'log');

% figure(3); hold all;
% scatter(LT3d(ind),Lo(ind),48,eff(ind),'filled')

% figure(2); hold all;
% plot([1 1e4],linspace(1e-8,1e-5,10)'*[1 1e4])

% plot(Reb(ind_lf), eff(ind_lf)./eff_peak, 'o');
% if(ic==1); fit_x=[]; fit_y=[]; end;
% fit_x = [fit_x ; Reb(ind)];
% fit_y = [fit_y ; smooth(eff(ind),'rlowess')];

% 
% fit_x = [fit_x ; Ri_bul(ind_rf)];
% % fit_x = [fit_x ; log(Reb(ind))];
% % fit_y = [fit_y ; smooth(eff(ind_rf),'rlowess')];
% fit_y = [fit_y ; smooth(eff(ind_rf),'rlowess')./Reb(ind_rf).^(-0.5)];

% % plot(Reb(ind),Pr_t(ind));
% myind = time<=tfe;
% plot(time(myind)-t3d,eff(myind));
% plot(time,KE2d,time,KE3d);

% ylim([0 1]);
% grid on;
% Ri_range = 0.10:0.02:0.2;       xx = Ri_range;
% plot(xx,intH3d./(intH3d+intD),xx,eff_cum)
% figure(1); plot(Ri_p(ind),eff(ind)); hold all;

% figure(1);
% loglog(Reb(ind),eff(ind)); hold all;
% plot(mean(Ri_p),eff_avg(ic),'ro'); hold all;
% plot(Ri_p(ind),XP(ind)./H3d(ind),'o'); hold all;

% figure(2);
% loglog(Lo(ind)./LT3d(ind), eff(ind),'ro');      hold all;
%semilogx(Reb(ind),Len(ind)./Lo(ind),'o');    hold all;
%semilogx(time(ind)-t3d,Len(ind)./Lo(ind));          hold all;
% plot(Reb(ind), Len(ind)./LT3d(ind),'*');       hold all;
% loglog(Reb(ind), kappa_star(ind)./kappa_cox(ind),'o');         hold all;
% loglog(kappa_star(ind), kappa_osb(ind),'o');         hold all;

% % Ri_t(ic) = mean(Ri_p(ind));
% % plot(Ri_t(ic), eff_cum(ic),'o')
% loglog(Reb(ind), tmp(ind),'r*'); hold all;

% semilogx(Reb(ind),eff(ind));        hold all;
% semilogx(Reb(ind),H3d(ind)./(H3d(ind)-D(ind))); hold all


% shear2_avg    = mean(shear.^2,1)';
% Re_shr = abs(D)./shear2_avg/nu;
% figure(1); semilogy(Reb(ind), Re_shr(ind));  hold all;

% ind = time>=t3d & time<=tfe;
% eps_tmp = 15/2*nu*mean(shear.^2,1);
% inteps_tmp = trapz(time(ind),eps_tmp(ind));
% figure;
% plot(time,-D3d,time,eps_tmp);
% ppcrt(ic) = mean(D3d(ind))/mean(eps_tmp(ind));
% figure(1);
% semilogy(time,kappa_3d  ,'b');        hold all;
% semilogy(time,kappa_osb ,'k');
% semilogy(time,kappa_cox ,'g');
% semilogy(time,kappa_star,'r-*');          hold all;
% hndl = legend('Mixing Length', ...
%               'Osborn 1980',...
%               'Osborn-Cox 1977',...
%               'Winters-D''Asaro 1996');
% set(hndl, 'interpreter','latex');    
% 

% figure(10);
% plot(time,XP,time,An);        hold all;

% figure(2);
% loglog(Reb(ind),kappa_3d(ind)  ,'b');        hold all;
% loglog(Reb(ind),kappa_osb(ind) ,'k');      hold all;
% loglog(Reb(ind),kappa_cox(ind) ,'g');       hold all;
% % kappa_tmp = 0.5*eff./(1-eff)*Pr.*Reb;
% loglog(Reb(ind),kappa_star(ind),'r-*');          hold all;
% hndl = legend('Mixing Length', ...
%               'Osborn 1980',...
%               'Osborn-Cox 1977',...
%               'Winters-D''Asaro 1996');
% set(hndl, 'interpreter','latex');  

% figure(3);
% semilogx(Reb(ind),gamma(ind),'*');          hold all;
% 
% figure(4);
% plot(Ri_avg(ind),gamma(ind),'*');          hold all;

% figure(5);
% loglog(Reb(ind),-H3d(ind)./D(ind),'*');          hold all;
% 
% figure(6);
% loglog(Reb(ind),-XP(ind)./D(ind),'*');          hold all;

% figure(7);
% loglog(Reb(ind),-D(ind),'*');          hold all;

% figure(8);
% loglog(Reb(ind),M(ind),'*');          hold all;
% 
% figure(9);
% if Ri<=0.02;
% %     plot(Ri_avg(ind),M(ind)./(Reb(ind).^0.7),'*');
%       nnn = 0.7;
%       fPr = Pr;
%       fRiPr = gamma(ind)./Reb(ind).^(nnn-1).*fPr;
%       plot(mean(Ri_avg(ind)),mean(fRiPr./fPr),'*');          
%     hold all    
% else
% %     plot(Ri_avg(ind),M(ind)./Reb(ind),'*');
%       nnn = 1.2;
%       fPr = Pr;
%       fRiPr = gamma(ind)./Reb(ind).^(nnn-1).*fPr;   
%       plot(mean(Ri_avg(ind)),mean(fRiPr./fPr),'*');
%     hold all    
% end





% 
% ROT =Lo./LT3d;
% figure(15);
% loglog(ROT(ind),eff(ind),'r*');     hold all;
% loglog(ROT(ind),0.16*ROT(ind).^(-0.63),'b-');

% plot(time,kappa_tmp,time,kappa_star)
% phi_d = -Ri*kappa_star.*mean(gradrhob)*kappa0;
% plot(time,phi_d,time,M)
% tmp(ic) = mean(kappa_tmp(ind)./Reb(ind)');
% loglog(Reb(ind),kappa_osb(ind) ,'k*');
% loglog(Reb(ind),eff_cum/0.2*kappa_osb(ind) ,'b*');

% hndl = legend('Flux formula', ...
%               'Osborn 1980',...
%               'Osborn-Cox',...
%               'Diascalar');
% set(hndl, 'interpreter','latex'); 

% loglog(Reb,0.2*Reb.^2*Pr,'-k')

% IMPORTANT COULD BE A PAPER on PARAMETERIZATION OF XP and EFF
% p = polyfit(Reb(ind),XP(ind),1);
% tintXP = XP*0;
% tintH3d= XP*0;
% tintD  = XP*0;
% tintM  = XP*0;
% tintAn = XP*0;
% 
% % % for i=2:ntime
% % %      tintXP(i)   = trapz(time(1:i),XP(1:i));
% % %      tintH3d(i)  = trapz(time(1:i),H3d(1:i));
% % % end
% for i=indt3d+1:indtfe
%     ind = indt3d;    
%     tintXP(i)  = trapz(time(ind:i),XP (ind:i));
%     tintH3d(i) = trapz(time(ind:i),H3d(ind:i));
%     tintD(i) = trapz(time(ind:i),D(ind:i));
%     tintM(i) = trapz(time(ind:i),M(ind:i));
%     tintAn(i)= trapz(time(ind:i),An (ind:i));
% end
% ind = time>t3d & time<=tfe;
% figure(12); hold all; 
% % plot(Reb(ind),tintXP(ind)./abs(tintD(ind)));grid on;
% % plot(Reb(ind),tintXP(ind));grid on;
% plot(Reb(ind), XP(ind));grid on;
% 
% figure(14); hold all; 
% plot(Reb(ind),tintH3d(ind));grid on;
% 
% figure(16); hold all; 
% plot(Reb(ind),tintAn(ind));grid on;

% plot(Reb(ind),polyval(p,Reb(ind)));grid on;

% figure(1); hold all;
% % plot(time,APE-APE(1));
% plot(time,XP);
% plot(time(indt3d),XP(indt3d),'*');
% 
% figure(2); hold all;
% plot(Ri,APE(indt3d)-APE(1),'o');
% temp = intDpAPE;%integrate(time,S);
% plot(Ri,temp(indt3d)-temp(indtfe),'o');
% plot(Ri,temp(indt3d),'*');
% figure(2); hold all;
% plot(Ri, APE(indt3d)-APE(indtfe),'o');
% plot(Ri, KE3d(indt3d)-KE3d(indtfe),'*');

% plot(Ri,intH2d(ic),'*k');
% plot(time-t2d,H2d);
% plot(t3d-t2d,H2d(indt3d));

% plot(Ri,intShbar(ic),'sk'); 
% plot(Ri,-intD3d(ic),'ok'); 


% plot(Ri,intShbar(indt3d)-intShbar(indtfe),'s'); 

% plot(Ri,KE2d(indt3d)-KE2d(indtfe),'*'); 
% plot(Ri,t3d-t2d,'o'); 
% 
% plot(t3d,KE(indt3d)/KE(1),'ko')

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


% addpath(genpath('/home/hesam/software/MATLAB2010a/figuremaker/'))
% exportfig(gcf, 'fig2.eps', 'FontMode', 'fixed', 'FontSize', 10, 'color', 'cmyk', 'height', 4, 'width',5);
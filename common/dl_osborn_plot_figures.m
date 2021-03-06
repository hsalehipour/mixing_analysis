
% Appendix figure: showing RMSE and R2 for all validation set based in
% increasing number of convNet layers.
test_rsme = [0.083206415, 0.07324314, 0.066356234, 0.06448664, 0.060592752, 0.05907017, 0.06259513];
test_r2  = [0.69783092, 0.76586290533, 0.8078237064, 0.81850034169, 0.83975750326, 0.8477095069477277, 0.828991600800043];
np_trainable=[2.0989e6, 1.0831e6, 5.9168e5, 3.6237e5, 2.6413e5, 2.3142e5, 2.3149e5];

nlayer = 1:7;
[haxes,hline1,hline2] = plotyy(nlayer, test_rsme.^2, nlayer, test_r2);
set(hline1, 'LineStyle','-', ...
    'Color',[0 0 0],...
    'MarkerFaceColor',[0.8 0.8 0.8],...
    'MarkerSize',6,...
    'Marker','square');
set(hline2, 'LineStyle','-',...
    'Color',[0 0 0],...
    'MarkerFaceColor',[0.8 0.8 0.8],...
    'MarkerSize',6,...
    'Marker','o');
set(haxes(1),...
    'xlim',[0 8],...
    'ycolor','k');   
set(haxes(2),...
    'xlim',[0 8],...
    'ycolor','k');   

ylabel(haxes(1),'$MSE$','interpreter','latex');
ylabel(haxes(2),'$r^2$','interpreter','latex');
xlabel('Number of ConvNet layers');
hndl = legend('MSE: Validation set', '$r^2$, Validation set');
set(hndl, 'interpreter','latex');

%{
% Main figures of the paper: comparing CNN predictions and Cox-Osborn
% methods with original DNS data
A = loadtxt([fadrs, 'dl_output.dat'], 2, 1); 
eff_nn = A(2,:)';

% Parameterization of Salehipour etal 2016 GRL
eff_par = mixeff_param(Reb,Ri_avg);

figure;
hndl = plot(time,eff_true, time,eff_nn, time, diag(eff_par), time, eff_cox);
set(hndl(1),'LineWidth',2,'Color',[0 0 0]);
set(hndl(2),...
    'LineStyle', 'None', ...
    'Color', [0 0 0], ...
    'MarkerFaceColor',[1 0.6 0.8],...
    'Marker','o',...
    'MarkerEdgeColor',[1 0. 0.]);
set(hndl(3),'LineWidth',2,'Color',[0.65 0.65 0.65],'LineStyle', ':');
set(hndl(4),'LineWidth',2,'Color',[0.65 0.65 0.65],'LineStyle', '--');



hndl = legend('$\mathcal{E}_{dns}$', '$E_{cnn}$', '$E_{par}$', '$R_{f,cox}$');
set(hndl, 'interpreter','latex');
%}

%{
% adding a bar graph on mean suqre error of different estimates
if ~exist('dl_err','var'); dl_err = []; end
L2_error = [mean_square_err(eff_true,eff_nn'), ...
            mean_square_err(eff_true,diag(eff_par)'), ...
            mean_square_err(eff_true,eff_cox) ];
dl_err = [dl_err; L2_error];
%}
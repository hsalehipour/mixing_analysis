% Read and analyze some data recieved from Lozovatsky and Fernando (2013)
% Is the Turbulent Prandtl number equal to unity at high Reynolds numbers?!

clear all; clc; close all;

fname='/home/hesam/kh-workspace/data_Lozovatsky_Fernando_2013/';
load([fname 'LS13.mat']);


Rf  = gamma0./(1+gamma0);
Prt = Ri./Rf;
ind=Ri<=0.5;

Psi = mixeff_param_Psi(Ri);
eff_peak = mixeff_param_effpeak(Ri);
Reb_peak = (Psi./eff_peak).^2;

figure
% line(Ri,Prt,'Color','r');
% ax1_pos = get(gca,'Position');
% ax2 = axes('Position',ax1_pos,...
%     'XAxisLocation','top',...
%     'YAxisLocation','right',...
%     'Color','none');
% line(flipud(Reb),Prt,'Parent',ax2,'Color','k')

scatter(Reb(ind)/300,Ri(ind),25,Rf(ind),'filled')
set(gca,'XScale','log');
hold all;
plot(Reb_peak(ind),Ri(ind))
clc; clear all; close all;
n=100;
rho = linspace(0.5,-0.5,n);
z = linspace(-1,1,2*n);
rho1 = rho;     rho2 = rho;

% rho = linspace(0,-0.5,n);
% z   = linspace(-0.5,0.5,2*n);
% rho1 = rho;     rho2 = rho+0.5;

%% PE
%initial
rhoi = [rho1 rho2];
PE = rhoi.*z;
figure(1);  plot(rhoi,z);   hold all;
figure(2);  plot(PE,z);    hold all;
PEi = trapz(z,PE)/(max(z)-min(z));

% final
rhof = [rho1 rho2];
ind = z>=-0.5 & z<=0.5;
rhof(ind) = 0;
PE = rhof.*z;
figure(1);      plot(rhof,z);
figure(2);      plot(PE,z);
PEf = trapz(z,PE)/(max(z)-min(z));

dPE = PEf - PEi;

%% BPE
rhosi = sort(rhoi,'descend');
PE = rhosi.*z;
figure(1);  plot(rhosi,z);   hold all;
figure(2);  plot(PE,z);    hold all;
BPEi = trapz(z,PE)/(max(z)-min(z));

% final
rhosf = sort(rhof,'descend');
PE = rhosf.*z;
figure(1);  plot(rhosf,z);   hold all;
figure(2);  plot(PE,z);    hold all;
BPEf = trapz(z,PE)/(max(z)-min(z));

dBPE = BPEf - BPEi;

eta = dBPE./abs((dPE-dBPE))
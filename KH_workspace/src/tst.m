% A rebuttal of Colm assertion about +-h vertical averaging for epsbar
% close all; clc;
ddz = 0;
dz = z(2)-z(1);
% eps_avg_time = epsbar*0;
% for i=1:length(z);
%     ddz = i*dz;
%     indz = z<=ddz & z>=-ddz;
%     LLz = max(z(indz))-min(z(indz));
%     eps_avg_time(i,:) = trapz(z(indz),epsbar(indz,:),1)./LLz;
%     plot(time,eps_avg_time(i,:)); hold all;
%     % pause;
% end

eps_avg_time = zeros(1,ntime);
N2_avg_time  = zeros(1,ntime);
N2_star      = -g/rho0*gradrhob;

for i=1:ntime;
    Lscale = Iu;
%     Lscale = shl_thk;

    indz = z<=Lscale(i)/2 & z>=-Lscale(i)/2;
    LLz = Lscale(i);
%     eps_avg_time(i) = trapz(z(indz),epsbar(indz,i),1)./LLz;
    eps_avg_time(i) = sum(dz*epsbar(indz,i),1)./LLz;

    
%    Lscale = Irho;
%     Lscale = shl_thk;
    indz = z<=Lscale(i)/2 & z>=-Lscale(i)/2;
    LLz = Lscale(i);
    N2_avg_time(i)  = sum(dz*N2_star(indz,i),1)./LLz;

    % pause;
end
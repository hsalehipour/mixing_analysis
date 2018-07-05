%==========================================================================
%                  Reading the output files for KH turbulence
%                     for spectral decomposition (PSD)
% Written by :              Hesam Salehipour
%                          January 14th, 2014
%==========================================================================
% This is a rather automated code to analyse the data;
% Multiple notes to which one should pay attention:
%   1) Make sure the data.in file contains the correct path to your desired
%      cases 
%   2) In each path, there should exist (i)   
%                                       (ii)  FILE.LIST
%   3) FILE.LIST has to include number of h*.dat files and their absolute
%      paths
%
%   4) All the neccessary analysis should be performed inside AnalyseData.m
%   5) All the ploting commands should be invoked inside PlotData.m
%==========================================================================


clear all;  close all;  clc;
%close all;  clc;


% REMEMBER that fh.0000001.dat is associated with t=60 and the rest are
% spaced with t=2 interval
time0 = 60;

%clear all;  close all;  clc;
%close all;  clc;

% Read number of cases to be compared (default = 1)
fname = 'data2.in';
fidm = fopen(fname);     
ncase=fscanf(fidm, '%f');



for ic=1:ncase
    % Reads the address for the current case
    fadrs = fgetl(fidm);
    
    % Read history files
    fid = fopen([fadrs, 'FILE.LIST']);
    h0=1;       hf=fscanf(fid,'%g');
    nhis = hf-h0+1;
    %     PSD_eyelid= zeros(129,nhis);
    %     PSD_core  = zeros(129,nhis);
    %     PSD_sp    = zeros(129,nhis);
    %     PSD_braid = zeros(129,nhis);
    %     PSD_tot   = zeros(129,nhis);
    
    h0 = 8;
    time = time0 + (h0-1)*2;
    %fname = sprintf([fadrs,'h.%07d.dat'],i);
    for i=1:h0
        fname = fgetl(fid);
    end
    A = loadtxt(fname,6,2);
    N  = A(1,:)';
    PSD_eyelid = A(2,:)';
    PSD_core   = A(3,:)';
    PSD_sp     = A(4,:)';
    PSD_braid  = A(5,:)';
    PSD_tot    = A(6,:)';
    sprintf('Reading file %s done!', fname)
    
    if (ic==1)
        PSD_eyelidPr1 = PSD_tot;
        PSD_corePr1   = PSD_tot;
        PSD_spPr1     = PSD_tot;
        PSD_braidPr1  = PSD_tot;
        PSD_totPr1    = PSD_tot;
    end
    ind = 1:25;
%     nPSD_eyelid = PSD_eyelid(ind)/sum(PSD_eyelid(ind));
%     nPSD_core   = PSD_core  (ind)/sum(PSD_core  (ind));
%     nPSD_sp     = PSD_sp    (ind)/sum(PSD_sp    (ind));
%     nPSD_braid  = PSD_braid (ind)/sum(PSD_braid (ind));
%     nPSD_tot    = PSD_tot   (ind)/sum(PSD_tot   (ind));
    
    nPSD_eyelid = PSD_eyelid(ind)/sum(PSD_eyelidPr1(ind));
    nPSD_core   = PSD_core  (ind)/sum(PSD_corePr1  (ind));
    nPSD_sp     = PSD_sp    (ind)/sum(PSD_spPr1    (ind));
    nPSD_braid  = PSD_braid (ind)/sum(PSD_braidPr1 (ind));
    nPSD_tot    = PSD_tot   (ind)/sum(PSD_totPr1   (ind));
%     figure(1);     
%     semilogy(N(ind), PSD_eyelid(ind));      hold on;
%     semilogy(N(ind), PSD_core(ind),'--')
%    
%     figure(2);     
%     semilogy(N(ind), PSD_braid(ind));        hold on;
%     semilogy(N(ind), PSD_sp(ind),'--');
    yminplot = 1e-8;    ymaxplot = 1e1;
    xminplot = 0;       xmaxplot = 24;

    figure(10);     
    semilogy(N(ind), N(ind).*nPSD_eyelid); hold on;
    ylabel('$\hat{P}(\gamma)$', 'interpreter','latex');
    xlabel('$\gamma$', 'interpreter','latex');
    xlim([xminplot xmaxplot]);       ylim([yminplot ymaxplot]);
    grid on; box on;


    figure(20);
    semilogy(N(ind), N(ind).*nPSD_core);   hold on;
    ylabel('$\hat{P}(\gamma)$', 'interpreter','latex');
    xlabel('$\gamma$', 'interpreter','latex')
    xlim([xminplot xmaxplot]);       ylim([yminplot ymaxplot]);
    grid on; box on;
      
    figure(30);     
    semilogy(N(ind), N(ind).*nPSD_braid);     hold on;
    ylabel('$\hat{P}(\gamma)$', 'interpreter','latex');
    xlabel('$\gamma$', 'interpreter','latex');
    xlim([xminplot xmaxplot]);       ylim([yminplot ymaxplot]);
    grid on; box on;

    
    figure(40);     
    semilogy(N(ind), N(ind).*nPSD_sp);  hold on;
    ylabel('$\hat{P}(\gamma)$', 'interpreter','latex');
    xlabel('$\gamma$', 'interpreter','latex');  
    xlim([xminplot xmaxplot]);       ylim([yminplot ymaxplot]);
    grid on; box on;

    
    figure(50);     
    semilogy(N(ind), N(ind).*nPSD_tot);    hold on;
    ylabel('$\hat{P}(\gamma)$', 'interpreter','latex');
    xlabel('$\gamma$', 'interpreter','latex');
    xlim([xminplot xmaxplot]);       ylim([yminplot ymaxplot]);
    grid on; box on;


%     figure(1);     ylim([1e-5 1e0]);   xlim([0,ind(end)]);
%     figure(2);     ylim([1e-5 1e0]);   xlim([0,ind(end)]);
%     figure(3);     ylim([1e-5 1e0]);   xlim([0,ind(end)]);
%     figure(4);     ylim([1e-5 1e0]);   xlim([0,ind(end)]);
%     figure(5);     ylim([1e-5 1e0]);   xlim([0,ind(end)]);
    
%     figure(6);     
%     semilogy(N(ind),nPSD_eyelid);    
%     hold on;
%     semilogy(N(ind),nPSD_sp,'r');    
%     semilogy(N(ind),nPSD_braid,'g');    
%     semilogy(N(ind),nPSD_tot,'y');    

    
    %[maxPSD(i), indx] = max(PSD(:,i));
    %iwn(i) = N(indx);

fid = fclose(fid);
    
end
fclose(fidm); 









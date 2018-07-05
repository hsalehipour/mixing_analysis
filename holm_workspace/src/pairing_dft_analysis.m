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


% REMEMBER that fh.0000001.dat is associated with t=2 and the rest are
% spaced with t=2 interval
time0 = 2;
nwn = 6;        
% number of wavenumbers investigated, 
% wave number k=1 for pairing mode
% wave number k=2 for KH mode
% others useless

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
    PSD  = zeros(nwn,nhis);
    time = zeros(nhis,1);
    %fname = sprintf([fadrs,'h.%07d.dat'],i);
    for i=h0:hf
        time(i) = time0 + (i-1)*2;
        fname = fgetl(fid);
        A = loadtxt(fname,2,2);
        N  = A(1,:)';
        PSD(:,i)= A(2,:)';
        sprintf('Reading file %s done!', fname)
    end

    figure(1);  hold all; % pairing
    plot(time,PSD(2,:)./PSD(3,:));
    
    figure(2);  hold all; % KH wave
    plot(time,PSD(3,:));

fid = fclose(fid);
    
end
fclose(fidm); 









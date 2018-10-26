%==========================================================================
%                  Reading the output files for KH turbulence
% Written by :              Hesam Salehipour
%                           September 8th, 2013
%==========================================================================
% This is a rather automated code to analyse the data;
% Multiple notes to which one should pay attention:
%   1) Make sure the data.in file contains the correct path to your desired
%      cases 
%   2) In each path, there should exist (i)   siminfo.dat file
%                                       (ii)  FILE.LIST
%                                       (iii) stat2d.dat
%                                       (iv)  stat3d.dat
%   3) FILE.LIST has to include number of h*.dat files and their absolute
%      paths
%
%   4) All the neccessary analysis should be performed inside AnalyseData.m
%   5) All the ploting commands should be invoked inside PlotData.m
%==========================================================================


% clear all;  close all;  clc;
% close all;  clc;

% Read number of cases to be compared (default = 1)
fname = 'data.in';
fidm = fopen(fname);     
ncase=fscanf(fidm, '%f');


for ic=1:ncase
    % Reads the address for the current case
    fadrs = fgetl(fidm);
    
    % 1) Startup
    StartUp;
    
    % 1) Reading data
    if (Pr==16)
        ReadData_old;      
    else
        ReadData;
    end
    
    
    % 2) Analysis
    if (Pr==16)
        AnalyseData_old;   
    else
        AnalyseData;
    end

    % 3) Plot data
    PlotData;

end
fclose(fidm); 









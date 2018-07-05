%==========================================================================
%                    Reading and Analysis the Argo data
% Written by :              Hesam Salehipour
%                           August 21st, 2015
%==========================================================================

clear all;  close all;  clc;
% clear all;  clc;

% Assume the value of gradient Richardson number
Ri= 0.25;

% Read number of cases to be compared (default = 1)
fname = 'argo_data.in';
fidm   = fopen(fname);

while 1
    
    % Read different datasets
    fadrs0 = fgetl(fidm);
    if isempty(fadrs0); break; end;
    ncase  = fscanf(fidm, '%f');
    
    for ic=1:ncase
        % Reads the address for the current case
        fadrs = [fadrs0,fgetl(fidm)];
        
        % 1) Startup & Reading data
        argo_StartUp;
    end 
end
fclose(fidm);               

% 2) Analysis
argo_AnalyseData

% 3) Plot data
argo_PlotData;
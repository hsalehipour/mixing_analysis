%==========================================================================
%                  Reading the output files of 2D DNS from NEK 
%                to be used in the sec. stability analysis code
% Written by :              Hesam Salehipour
%                          April 4th, 2014
%==========================================================================
% This is a rather automated code to start up the fortran secondary
% stability analysis code.
%
% Multiple notes to which one should pay attention:
% 
% BEFORE USING THIS CODE:
%   0) -- Read the 2D DNS results in VisIt
%      -- Export the fields coordinates, temp, x-velocity and y-velocity as
%         Xmdv format
%      -- Copy the generated .okc files from login.scinet into your desired
%         directory
%   1) Make sure the data_stab.in file contains the correct path to your
%      desire cases 
%   2) In each path, there should exist "FILE.LIST" (generated by
%      okcgenflist)
%   3) FILE.LIST has to include number of .okc files and their absolute
%      paths
%   4) Ensure the fortran code is compiled inside the "src" folder.
%   5) All the neccessary analysis should be performed inside secAnalyseData.m
%   6) All the ploting commands should be invoked inside secPlotData.m
%==========================================================================


clear all;  close all;  clc;
%close all;  clc;


%================ IMPORTANT INPUT==========================================
% H2:  the height of middle box (always taken to be 10 of H=30 where mesh is
%      more refined)
% RSC: ratio of shear to density layer
time = 72;

% PR=1
% Pr = 1;         Ri = 0.12;      Re = 6000;
% nelx = 64;      nelz = 66;      nelz2 = 44;

% % PR=2
% Pr = 2;         Ri = 0.12;      Re = 6000;
% nelx = 64;      nelz = 66;      nelz2 = 44;

% PR=4
% Pr = 4;         Ri = 0.12;      Re = 6000;
% nelx = 90;      nelz = 88;      nelz2 = 62;

% PR=8
% Pr = 8;         Ri = 0.12;      Re = 6000;
% nelx = 127;     nelz = 116;     nelz2 = 88;

% % % PR=16
Pr = 16;        Ri = 0.12;      Re = 6000;
nelx = 180;     nelz = 156;     nelz2 = 124;

% desired spanwise wave number range
dwn1 = 2;
dwn2 = 20;
ndwn = 10;

RSC = 1;        % ratio of shear to density layer
H2 = 10;     
N = 10;             % order
%==========================================================================



% Read number of cases to be compared (default = 1)
fname = 'data_stab.in';
fidm = fopen(fname);     
ncase=fscanf(fidm, '%f');


delz  = H2/(N*nelz2);

for ic=1:ncase
    % Reads the address for the current case
    fadrs = fgetl(fidm);
    
    % Read data files
    fid = fopen([fadrs, 'data_in/FILE.LIST']);
    h0=1;       hf=fscanf(fid,'%g');
    nhis = hf-h0+1;
    
    x_input   = [];       z_input  = [];
    rho_input = [];       ux_input = [];    uz_input = [];
   
    for i=h0:hf
        fname = fgetl(fid);
        A = loadtxt(fname, 6, 13);
        x_input  = [x_input; A(1,:)'];
        z_input  = [z_input; A(2,:)'];
        rho_input = [rho_input; A(4,:)'];
        ux_input  = [ux_input ; A(5,:)'];
        uz_input  = [uz_input ; A(6,:)'];       
        sprintf('Reading file %s done!', fname)
    end
    fclose(fid);

    data = [x_input, z_input, rho_input, ux_input uz_input];
    data_sorted = sortrows(data,[2 1]);
    
    x  = reshape(data_sorted(:,1),nelx*N,nelz*N);
    z  = reshape(data_sorted(:,2),nelx*N,nelz*N);
    rho= reshape(data_sorted(:,3),nelx*N,nelz*N);
    ux = reshape(data_sorted(:,4),nelx*N,nelz*N);
    uz = reshape(data_sorted(:,5),nelx*N,nelz*N);
    
    
    % get rid of multiple coincident points
    itr = 0;
    for i=N:N:N*nelz-1
        icnt = i-itr;
        arr1 = [x(:,icnt), z(:,icnt), rho(:,icnt), ...
                ux(:,icnt), uz(:,icnt)];
        arr2 = [x(:,icnt+1), z(:,icnt+1), rho(:,icnt+1), ...
                ux(:,icnt+1), uz(:,icnt+1)];
        % x-grid
        vec_part1 = reshape(arr1(:,1),2,nelx*N/2);
        vec_part2 = reshape(arr2(:,1),2,nelx*N/2);
        x(:,icnt+1) = [vec_part1(1,:)';vec_part2(1,:)'];
        
        % z-grid
        vec_part1 = reshape(arr1(:,2),2,nelx*N/2);
        vec_part2 = reshape(arr2(:,2),2,nelx*N/2);
        z(:,icnt+1) = [vec_part1(1,:)';vec_part2(1,:)'];
        
        % rho
        vec_part1 = reshape(arr1(:,3),2,nelx*N/2);
        vec_part2 = reshape(arr2(:,3),2,nelx*N/2);
        rho(:,icnt+1) = [vec_part1(1,:)';vec_part2(1,:)'];
        
        % ux
        vec_part1 = reshape(arr1(:,4),2,nelx*N/2);
        vec_part2 = reshape(arr2(:,4),2,nelx*N/2);
        ux(:,icnt+1) = [vec_part1(1,:)';vec_part2(1,:)'];
        
        % uy
        vec_part1 = reshape(arr1(:,5),2,nelx*N/2);
        vec_part2 = reshape(arr2(:,5),2,nelx*N/2);
        uz(:,icnt+1) = [vec_part1(1,:)';vec_part2(1,:)'];
        
        
        x(:,icnt) = [];
        z(:,icnt) = [];
        rho(:,icnt) = [];
        ux (:,icnt) = [];
        uz (:,icnt) = [];
        
        itr = itr+1;
    end
        
    itr = 0;
    for i=N:N:N*nelx-1
        icnt = i-itr;        
        x(icnt,:) = [];
        z(icnt,:) = [];
        rho(icnt,:) = [];
        ux (icnt,:) = [];
        uz (icnt,:) = [];        
        itr = itr+1;
    end
    
    
   % Interpolate to uniform grid.
   % Parameters used in sec stability analysis  code
    nxpnts=nelx*N;             % the factor 1/2 is arbitrary (HS)
    nzpnts=nelz2*N;
    xmax=max(max(x));     xmin=min(min(x));
    zmax = 5;             zmin=-5;
    xp = linspace(xmin,xmax,nxpnts);
    zp = linspace(zmin,zmax,nzpnts);
    [XX, ZZ] = meshgrid(xp,zp);
    RHO = interp2(x',z',rho',XX,ZZ);
    UX  = interp2(x',z',ux',XX,ZZ);
    UZ  = interp2(x',z',uz',XX,ZZ);
    
    
    % compute shear and stretch
    % shear   = ddz_ux + ddx_uz
    % stretch = ddx_ux - ddz_uz
    shear   = ddx(UX,zp') + ddx(UZ,xp');
    stretch = ddx(UX,xp') - ddx(UZ,zp');
    
    % draw to make sure
    mesh(XX,ZZ,RHO); view(0,90); axis equal;
    
         
    % write siminfo.dat
    fid = fopen([fadrs,'siminfo.dat'], 'w');    
    fprintf(fid, '%17.9e', time);
    fprintf(fid, '%17.9e', Re);
    fprintf(fid, '%17.9e', Ri);
    fprintf(fid, '%17.9e', Pr);
    fprintf(fid, '%17.9e', RSC);
    fprintf(fid, '%17.9e', xmin);
    fprintf(fid, '%17.9e', xmax);
    fprintf(fid, '%17.9e', zmin);
    fprintf(fid, '%17.9e', zmax);
    fprintf(fid, '%5i', nxpnts);
    fprintf(fid, '%5i', nzpnts);
    fprintf(fid, '%17.9e', dwn1);
    fprintf(fid, '%17.9e', dwn2);
    fprintf(fid, '%5i'   , ndwn);
    fclose(fid);

    % write the "data.in" file
    fid = fopen([fadrs,'data.in'], 'w');  
    fprintf(fid, 'XX    ZZ   rho  ux   uz   shear   stretch \n');
    for ij=1:nxpnts*nzpnts
    fprintf(fid, '%17.9e%17.9e%17.9e%17.9e%17.9e%17.9e%17.9e\n'...
                ,  XX(ij), ZZ(ij), ...
                   RHO(ij), UX(ij), UZ(ij), ...
                   shear(ij), stretch(ij));
    end
    fclose(fid);
    
    % Run the fortran compiled "sec_stab" code.
%     !./sec_stab
end
fclose(fidm); 









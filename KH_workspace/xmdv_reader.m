%==========================================================================
%                  Reading the sliced output files of 3D DNS from NEK 
%                to be used for further analysis in MATLAB
% Written by :              Hesam Salehipour
%                          Nov 13th, 2017
%==========================================================================



function [scfld_out, x_out,z_out] = xmdv_reader(scfld_in, x_in, z_in, nelx, nelz)
% nfld: number of scalar fields in .okc file besides the xyz grid.
%       eg: if exporting density only, nfld=1, exporting density and vel
%       vectorm nfld = 4

% I have always used N=10 as the polynomila order  (Hardcoded here)
N = 10;

% sort based on x,z corrdinates
data = [x_in, z_in, scfld_in];
data_sorted = sortrows(data,[2 1]);

x  = reshape(data_sorted(:,1),nelx*N,nelz*N);
z  = reshape(data_sorted(:,2),nelx*N,nelz*N);
fld= reshape(data_sorted(:,3),nelx*N,nelz*N);


% get rid of multiple coincident points
itr = 0;
for i=N:N:N*nelz-1
    icnt = i-itr;
    arr1 = [x(:,icnt)  , z(:,icnt)  , fld(:,icnt)];
    arr2 = [x(:,icnt+1), z(:,icnt+1), fld(:,icnt+1)];
    
    % x-grid
    vec_part1 = reshape(arr1(:,1),2,nelx*N/2);
    vec_part2 = reshape(arr2(:,1),2,nelx*N/2);
    x(:,icnt+1) = [vec_part1(1,:)';vec_part2(1,:)'];
    
    % z-grid
    vec_part1 = reshape(arr1(:,2),2,nelx*N/2);
    vec_part2 = reshape(arr2(:,2),2,nelx*N/2);
    z(:,icnt+1) = [vec_part1(1,:)';vec_part2(1,:)'];
    
    % vort_y
    vec_part1 = reshape(arr1(:,3),2,nelx*N/2);
    vec_part2 = reshape(arr2(:,3),2,nelx*N/2);
    fld(:,icnt+1) = [vec_part1(1,:)';vec_part2(1,:)'];
    
    % set repeated values to void
    x(:,icnt) = [];
    z(:,icnt) = [];
    fld(:,icnt) = [];
    
    itr = itr+1;
end

itr = 0;
for i=N:N:N*nelx-1
    icnt = i-itr;
    x(icnt,:) = [];
    z(icnt,:) = [];
    fld(icnt,:) = [];
    itr = itr+1;
end


% Interpolate to uniform grid of power 2
nxpnts=2^ceil(log2(min(nelx,nelz)*N));
nzpnts=nxpnts;
xmax=max(max(x));     xmin=min(min(x));
zmax = 5;             zmin=-5;
xp = linspace(xmin,xmax,nxpnts);
zp = linspace(zmin,zmax,nzpnts);
[x_out, z_out] = meshgrid(xp,zp);
scfld_out = interp2(x',z',fld',x_out,z_out);

end



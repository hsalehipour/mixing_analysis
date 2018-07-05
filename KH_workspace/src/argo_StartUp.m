%% 1) analysis startup


% Problem Setup
fname = [fadrs0,'siminfo.dat'];
fid = fopen(fname);     S=fscanf(fid, '%s');    fclose(fid);


%% Read data
load(fadrs);
lon=x;       lat=y; 
ndepth = ncase;

% setup initial values
nu = 1e-6;

nlon = length(lon);     
nlat = length(lat);
if ~exist('N2','var')
    N2    = zeros(nlat,nlon,ncase);
    disp  = N2;  
    K_osb = N2;
    zlev  = [];    
end

tf1 = strcmp(S,'bf');        if (tf1); N2   (:,:,ic)=c;   end;
tf2 = strcmp(S,'eps');       if (tf2); disp (:,:,ic)=c;   end;
tf3 = strcmp(S,'kappa');     if (tf3); K_osb(:,:,ic)=c;   end;

if (tf2);
    zlev = [zlev; -Zlim1'];
end

%% Find the buoyancy Reynolds number from the Argo float data

if (tf3) 
    Reb = K_osb/0.2/nu;
else
    Reb = disp./N2/nu;
    K_osb = 0.2*nu*Reb;
end

%% Mixing efficiency
% =======================
% A multiparameter Parameterization of mixing efficiency
% =======================
eff=Reb./0;
for k=1:ndepth
for i=1:nlat
for j=1:nlon
    eff(i,j,k) = mixeff_param(Reb(i,j,k),Ri); 
%      effp = 0.2; Rep=100;
%     eff(i,j,k) = 2*effp*(Reb(i,j,k)/Rep).^0.5./(1+Reb(i,j,k)/Rep);
%      if Reb(i,j,k)<Rep
%          eff(i,j,k)=effp*(Reb(i,j,k)/Rep).^0.5;
%      else
%          eff(i,j,k)=effp*(Reb(i,j,k)/Rep).^(-0.5);
%      end
end
end
end


%% Turbulent Pr
%
Prt = Ri./eff;          
Prt(eff<1e-4) = NaN;  
Prt(Prt>100)  = NaN;



%% Kappa based on New Eff
Gamma = eff./(1-eff);
K_new = Gamma.*Reb*nu;


%% Kappa_m based on New Eff
Km_new = Ri./(1-eff).*Reb*nu;


%% Calculate the global averages

Kosb_avg=zeros(1,ndepth); 
Knew_avg=zeros(1,ndepth);
%
Kappa_osb_avg=zeros(1,ndepth);
Kappa_new_avg=zeros(1,ndepth);

for k=1:ndepth
    Kosb_avg(k) = mean_weighted(lat,K_osb(:,:,k));
    Knew_avg(k) = mean_weighted(lat,K_new(:,:,k));
    Kappa_new_avg(k) = mean_weighted(lat,Gamma(:,:,k).*disp(:,:,k))./...
                       mean_weighted(lat,N2(:,:,k));
    Kappa_osb_avg(k) = mean_weighted(lat,          0.2*disp(:,:,k))./...
                       mean_weighted(lat,N2(:,:,k));
end


%% Comparison with Bryan and Lewis (1979) profile and POP1 profiles

zz =-5000:0;

% Bryan and Lewis (1979)
kappa_BL   = 0.8 + 1.05/pi*atan((abs(zz)-2500)/222.2);

% Large (1994) kpp_ccsm3:
bckgrnd_vdc1    = 0.524;
bckgrnd_vdc2    = 0.313;
bckgrnd_vdc_dpth= 1000.0e02;
bckgrnd_vdc_linv= 4.5e-05;
Prandtl         = 10.0;
rich_mix        = 50.0;
kappa_POP1 = bckgrnd_vdc1 + ...
             bckgrnd_vdc2*atan((abs(zz)*100-bckgrnd_vdc_dpth)*bckgrnd_vdc_linv);

% Read netcdf files containing other vertical profiles of kappa
fadrs = '/home/hesam/kh-workspace/guido_kappa_profiles/';
fname1 = [fadrs, 'kappa_obs.nc'];
fname2 = [fadrs, 'kappa_ccsm4_globalavg.nc'];
fname3 = [fadrs, 'weighted_argo_kappa_ccsm4_globalavg.nc'];
% 
% % file 1
ncid = netcdf.open(fname1,'NC_NOWRITE');
depth1   = netcdf.getVar(ncid,0,'double')';     % depth levels 
depth2   = netcdf.getVar(ncid,1,'double')';     % depth levels
kappa_ls = netcdf.getVar(ncid,2,'double')';     % Lampkin & Speer (2007)
kappa_g  = netcdf.getVar(ncid,3,'double')';     % Ganachaud (2003
kappa_dl = netcdf.getVar(ncid,4,'double')';     % Decloedt & Luther (2012)
kappa_nf = netcdf.getVar(ncid,5,'double')';     % Nikurashin & Ferrari (2013) 
netcdf.close(ncid);
% 
% 
% file 2
ncid = netcdf.open(fname2,'NC_NOWRITE');
depth3     = netcdf.getVar(ncid,0,'double')';     % depth levels 
kappa_cesm = netcdf.getVar(ncid,1,'double')';     % ccsm4 global avg model
netcdf.close(ncid);

% % file 3
% ncid = netcdf.open(fname3,'NC_NOWRITE');
% kappa_argo_cesm = netcdf.getVar(ncid,1,'double')';     % argo-loc ccsm4 global avg model
% netcdf.close(ncid);



%% Write Argo location file in .nc 
% Purpose: to compute CESM global average only on the Argo locations 
% argo_loc = Reb;
% argo_loc(~isnan(Reb))=1;
% argo_loc(isnan(Reb))=0;
% nlat = length(lat);
% nlon = length(lon);
% 
% ncid = netcdf.create('argo_latlon.nc','NC_WRITE');
% x_dimid = netcdf.defDim(ncid,'nlat',nlat);
% y_dimid = netcdf.defDim(ncid,'nlon',nlon);
% dimid = [x_dimid, y_dimid];
% lat_varID = netcdf.defVar(ncid,'lat','double',x_dimid);
% lon_varID = netcdf.defVar(ncid,'lon','double',y_dimid);
% loc_varID = netcdf.defVar(ncid,'location','double',dimid);
% netcdf.endDef(ncid);
% % 
% netcdf.putVar(ncid,lat_varID,lat);
% netcdf.putVar(ncid,lon_varID,lon);
% netcdf.putVar(ncid,loc_varID,argo_loc);
% netcdf.close(ncid);

% 1) Reading data

%% 1-1) Read stat2d file
fname = [fadrs,'stat2d.dat'];
A = loadtxt(fname,8,1);
time   =  A(1,:)';
KE     =  A(2,:)';
KEbar  =  A(3,:)';
KE2d   =  A(4,:)';
sig_Hb = -A(5,:)';      % based on my "reasonable" definition of buoy flux
sig_Dv =  A(6,:)';
%sigma  =  A(7,:)';     % not used because sigma = sig_Hb + sig_Dv
sigma  = -sig_Hb + sig_Dv;
eps_klmg= A(8,:)'*nu;   
% because it is not devided by Re in ppnek.usr



%% 1-2 Read stat3d file
fname = [fadrs,'stat3d.dat'];
% A = loadtxt(fname,9,1);
A = loadtxt(fname,10,1);
KE3d     =  A(2,:)';
sig3_Rs  =  A(3,:)';
sig3_Sh  =  A(4,:)';
sig3_An  =  A(5,:)';
sig3_Hb  = -A(6,:)';
sig3_Dv  =  A(7,:)';
uw       =  A(8,:)';
LT3d     =  A(9,:)';
rhop_rms =  A(10,:)';       % rms of rho'=rho2d+rho3d= <rho'^2>^(0.5)

%sigma3d = sig3_Rs + sig3_Sh + sig3_An + sig3_Hb + sig3_Dv;
sigma3d = sig3_Rs + sig3_Sh + sig3_An - sig3_Hb + sig3_Dv;




%% 1-3 Read history files
fid = fopen([fadrs, 'FILE.LIST']);
h0=1;       hf=fscanf(fid,'%g');
nhis = hf-h0+1;
rhobar  = zeros(nzp,nhis);
rhob    = zeros(nzp,nhis);
ubar    = zeros(nzp,nhis);
epsbar  = zeros(nzp,nhis);
chi_pr  = zeros(nzp,nhis);      % chi calculated using rho'=rho_2d + rho_3d 
chi_ref = zeros(nzp,nhis);      % chi calculated following rho_base
Hbar    = zeros(nzp,nhis);      % Hbar = <rho' w'>_{xy}= buoyflux(z)

for i=h0:hf
    %fname = sprintf([fadrs,'h.%07d.dat'],i);
    fname = fgetl(fid);
    A = loadtxt(fname, 10, 2);
    z = A(3,:)';
    rhobar (:,i) = A(4,:)';
    rhob   (:,i) = A(5,:)';
    ubar   (:,i) = A(6,:)';
    epsbar (:,i) = A(7,:)'*nu;
    chi_pr (:,i) = A(8,:)';
    chi_ref(:,i) = A(9,:)';
    Hbar   (:,i) = g*A(10,:)';      % it's not multiplied by g in .usr

    sprintf('Reading file %s done!', fname)
end
fid = fclose(fid);

%% 1-4 Read soc.*.dat files (if present)
%{
fadrs2 = [fadrs, 'soc\'];
ifsoc = exist(fadrs2,'dir')==7;

% Read the soc-related extracted quantitis at z=0 and z=+-Irho/2 
% these data are caclulated by python-visit routines.
if ifsoc
fid = fopen([fadrs2, 'FILE.LIST']);
h0=1;       hf=fscanf(fid,'%g');
A = loadtxt([fadrs2,'soc.0000.dat'], 8, 1);
nxp = length(A(1,:));
nhis = hf-h0+1;
xsoc = zeros(nxp,nhis);
zsoc = zeros(nxp,nhis);
vort2d_top = zeros(nxp,nhis);
vort2d_bot = zeros(nxp,nhis);
vort3d_top = zeros(nxp,nhis);
vort3d_bot = zeros(nxp,nhis);
vort2d_z0  = zeros(nxp,nhis);
vort3d_z0  = zeros(nxp,nhis);

for i=h0:hf
    fname = fgetl(fid);
    A = loadtxt(fname, 8, 1);
    xsoc(:,i) = A(1,:)';
    zsoc(:,i) = A(2,:)';
    vort2d_top(:,i) = A(3,:)';
    vort2d_bot(:,i) = A(4,:)';
    vort3d_top(:,i) = A(5,:)';
    vort3d_bot(:,i) = A(6,:)';
    vort2d_z0 (:,i) = A(7,:)';
    vort3d_z0 (:,i) = A(8,:)';

    sprintf('Reading file %s done!', fname)
end
fid = fclose(fid);
end
%}

%% 1-4 Read .okc files (if present) for x-z slices
% fadrs2 = [fadrs, 'soc\'];
% ifokc = exist(fadrs2,'dir')==7;
% 
% % calculate shear and density length scales to reduce 2D slices with dim
% %(x,z,t) into (x,t) time series.
% Iu   = trapz(z,1-ubar.^2,1);
% Irho = trapz(z,1-(rhobar-1).^2,1);
% 
% if ifokc    
%     fid = fopen([fadrs2, 'FILE.LIST']);
%     h0=1;       hf=fscanf(fid,'%g');
%     nhis = hf-h0+1;    
%     for i=h0:hf        
%         % read the number of outputed fields in .okc
%         fname = fgetl(fid);
%         fid2 = fopen(fname); 
%         C=textscan(fid2,'%d%d%d'); 
%         nfld = C{1};
%         fclose(fid2);
%         
%         % Read .okc data              
%         if nfld>4
%             nskip = 2*nfld+1;
%             A = loadtxt(fname, nfld, nskip);
%             x_in  = A(1,:)';
%             z_in  = A(2,:)';
%         else
%             nskip = 2*nfl+1;
%             A = loadtxt(fname, nfld, nskip);
%         end
% 
%         % read the scalar fields
%         vort3d_x_in  = A(nfld-3,:)';
%         vort3d_y_in  = A(nfld-2,:)';
%         vort3d_z_in  = A(nfld-1,:)';
%         vort2d_in    = A(nfld,:)';
%         
%         % get rid of repeated values in .okc files
%         [vort3d_x, x_okc, z_okc] = xmdv_reader(vort3d_x_in, x_in, z_in, nelx, nelz);
%         vort3d_y = xmdv_reader(vort3d_y_in, x_in, z_in, nelx, nelz);
%         vort3d_z = xmdv_reader(vort3d_z_in, x_in, z_in, nelx, nelz);
%         vort2d_z = xmdv_reader(vort2d_in, x_in, z_in, nelx, nelz);
%     end
% end

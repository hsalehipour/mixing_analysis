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
A = dlmread(fname, '', 1,0)';
KE3d     =  A(2,:)';
sig3_Rs  =  A(3,:)';
sig3_Sh  =  A(4,:)';
sig3_An  =  A(5,:)';
sig3_Hb  = -A(6,:)';
sig3_Dv  =  A(7,:)';
uw       =  A(8,:)';
LT3d     =  A(9,:)';
if size(A,2)>9
    rhop_rms =  A(10,:)';       % rms of rho'=rho2d+rho3d= <rho'^2>^(0.5)
end
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
chi3d   = zeros(nzp,nhis);      % chi calculated using rho_3d
chi_pr  = zeros(nzp,nhis);      % chi calculated using rho'=rho_2d + rho_3d
chi_ref  = zeros(nzp,nhis);     % or chi_base (not sure this is accurate)
Hbar    = zeros(nzp,nhis);      % Hbar = <rho' w'>_{xy}= buoyflux(z)

% Note: 26-Oct-2018: unfortunately I have not computed chi_pr consistently 
% for all KHI cases and now it is costly to rerun all post-processing.
% for the purpose of DL-paper I have had to re-run post-processing just for
% the 3 of the validation test cases. This enables comparison with
% Osborn-Cox formulation. So I will leave chi_pr here to be zero so that
% anything computed based on this becomes meaningless and don't cause the
% issue that I did not catch in the initial submission of my Rapids DLpaper

for i=h0:hf
    %fname = sprintf([fadrs,'h.%07d.dat'],i);
    fname = fgetl(fid);
    % A = loadtxt(fname, 9, 2);
    A = dlmread(fname, '', 2,0)';
    z = A(3,:)';
    rhobar (:,i) = A(4,:)';
    rhob   (:,i) = A(5,:)';
    ubar   (:,i) = A(6,:)';
    epsbar (:,i) = A(7,:)'*nu;
    if size(A,1)>9
        chi_pr (:,i) = A(8,:)';
        chi_ref(:,i) = A(9,:)';
        Hbar   (:,i) = g*A(10,:)';      % it's not multiplied by g in .usr
    else
        chi3d  (:,i) = A(8,:)';
        chi_ref (:,i) = A(9,:)';
    end

    sprintf('Reading file %s done!', fname)
end
fid = fclose(fid);


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
A = loadtxt(fname,8,1);
KE3d     =  A(2,:)';
sig3_Rs  =  A(3,:)';
sig3_Sh  =  A(4,:)';
sig3_An  =  A(5,:)';
sig3_Hb  = -A(6,:)';
sig3_Dv  =  A(7,:)';
uw       =  A(8,:)';
%LT3d     = A(8,:)';
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
for i=h0:hf
    %fname = sprintf([fadrs,'h.%07d.dat'],i);
    fname = fgetl(fid);
    A = loadtxt(fname, 7, 2);
    z = A(3,:)';
    rhobar (:,i) = A(4,:)';
    rhob   (:,i) = A(5,:)';
    ubar   (:,i) = A(6,:)';
    epsbar (:,i) = A(7,:)'*nu;
    sprintf('Reading file %s done!', fname)
end
fid = fclose(fid);


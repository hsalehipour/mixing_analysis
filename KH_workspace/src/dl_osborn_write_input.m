function dl_osborn_write_input(fname, my_eff, epsbar,N2, kappa0,z,Iu,nrand_extra)
% Purpose: write the input data for dl_osborn code.
    
    % Random number generator which generates N rand numbers within [a, b]
    rand_gen = @(a,b,N) a + (b-a).*rand(N,1);
    
    % check to see if file exists to be appended 
    if exist(fname,'file') == 2
        fid1 = fopen(fname, 'a');
    else
        fid1 = fopen(fname, 'w');
        fprintf(fid1, '%s, \t %s, \t %s \r\n', 'eps/Dp[Iu]','N2/N2avg[Iu]', 'eff_DNS');
    end
    
    ntime = size(N2,2);
    nz_profile = 512;
    DD = epsbar;             
    DD = max(0., DD);      % eps, turbul dissip.
    
    % fix the vertical resolution for all cases (Lz=30 often)
    Lz = max(z)-min(z);
    dz = Lz/(nz_profile-1);
    
    for i=1:ntime;
        Lbot = rand_gen(min(z),-Iu(i),nrand_extra);
        Ltop = rand_gen(Iu(i),max(z) ,nrand_extra);
        for nn=1:nrand_extra;
            % select only a portion of the profile
            indz = z<=Ltop(nn) & z>=Lbot(nn);
            
            % number of data points in the selected range
            nzr = round((max(z(indz))-min(z(indz)))/dz)+1;
            N2tot = trapz(z(indz),N2(indz,i),1);
            Dp    = kappa0*N2tot;
            features = zeros(nz_profile, 3);
            
            % interpolate to the selected z-range
            zft = linspace(min(z(indz)),max(z(indz)),nzr);
            ft1 = interp1(z(indz), DD(indz,i)/Dp, zft);
            ft2 = interp1(z(indz), N2(indz,i)/N2tot, zft);
            features(1:nzr,1:2) = [ft1', ft2'];
            features(:,3) = my_eff(i);
            fprintf(fid1, '%17.9e, \t %17.9e, \t %17.9e \r\n', features');
        end
    end
    fclose(fid1);
end
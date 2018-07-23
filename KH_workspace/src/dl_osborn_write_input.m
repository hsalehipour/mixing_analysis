function dl_osborn_write_input(fname, my_eff, epsbar,N2, kappa0,z,Iu)
% Purpose: write the input data for dl_osborn code.

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
    
    for i=1:ntime;       
        features = zeros(nz_profile, 3);
        
        % select only a portion of the profile
        indz = z<=Iu(i)/2. & z>=-Iu(i)/2.;
        N2avg = mean(N2(indz,i));
        Dp    = kappa0*N2avg;
        
        % interpolate to the selected z-range
        zft = linspace(min(z(indz)), max(z(indz)), nz_profile);
        ft1 = interp1(z(indz), DD(indz,i)/Dp, zft);
        ft2 = interp1(z(indz), N2(indz,i)/N2avg, zft);
        features(:,1:2) = [ft1', ft2']; 
        features(:,3) = my_eff(i);
        fprintf(fid1, '%17.9e, \t %17.9e, \t %17.9e \r\n', features');
    end
    fclose(fid1);
end
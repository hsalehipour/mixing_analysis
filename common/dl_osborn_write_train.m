function dl_osborn_write_train(fname, my_eff, M, epsbar,N2, kappa0,z,zmax_shl,zmin_shl)
% Purpose: write the input data for dl_osborn code.
        
    % check to see if file exists to be appended 
    if exist(fname,'file') == 2
        fid1 = fopen(fname, 'a');
    else
        fid1 = fopen(fname, 'w');
        fprintf(fid1, '%s, \t %s, \t %s, \t %s, \t %s \r\n', 'z/L_shl', 'eps/Dp','N2/N2tot', 'eff_DNS', 'M_DNS/Dp');
    end
    
    ntime = size(N2,2);
    nz_profile = 512;
    DD = epsbar;             
    DD = max(0., DD);     % eps, turbul dissip.
    MM = max(0., M);      % M, irreversible mixing
   
    for i=1:ntime;
        % remove the deadregions in the profile bsaed on L_shrl
        indz = z<=zmax_shl(i) & z>=zmin_shl(i);
        L_shl = max(z(indz)) - min(z(indz));
        
        % number of data points in the selected range
        N2tot = trapz(z(indz),N2(indz,i),1);
        Dp    = kappa0*N2tot;
        features = zeros(nz_profile, 5);
        
        % interpolate to the selected z-range
        zft = linspace(min(z(indz)),max(z(indz)),nz_profile);
        ft1 = zft./L_shl;
        ft2 = interp1(z(indz), DD(indz,i)/Dp, zft);
        ft3 = interp1(z(indz), N2(indz,i)/N2tot, zft);
        features(:, 1)  = ft1;
        features(:,2:3) = [ft2', ft3'];
        features(:,4) = my_eff(i);
        features(:,5) = MM(i)/Dp;
        fprintf(fid1, '%17.9e, \t %17.9e, \t %17.9e, \t %17.9e, \t %17.9e \r\n', features');
    end
    fclose(fid1);
end
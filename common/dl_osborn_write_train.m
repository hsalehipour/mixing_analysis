function dl_osborn_write_train(fname, my_eff, epsbar,N2, kappa0,z)
% Purpose: write the input data for dl_osborn code.
        
    % check to see if file exists to be appended 
    if exist(fname,'file') == 2
        fid1 = fopen(fname, 'a');
    else
        fid1 = fopen(fname, 'w');
        fprintf(fid1, '%s, \t %s, \t %s, \t %s \r\n', 'z/Lo_patch', 'eps/Dp','N2/N2tot', 'eff_DNS');
    end
    
    ntime = size(N2,2);
    nz_profile = 512;
    DD = epsbar;             
    DD = max(0., DD);      % eps, turbul dissip.
    Lz = max(z)-min(z);
   
    for i=1:ntime;                
        % number of data points in the selected range
        N2tot = trapz(z,N2(:,i),1);
        Dp    = kappa0*N2tot;
        features = zeros(nz_profile, 3);
        
        % interpolate to the selected z-range
        zft = linspace(min(z),max(z),nz_profile);
        ft1 = zft./Lz;
        ft2 = interp1(z, DD(:,i)/Dp, zft);
        ft3 = interp1(z, N2(:,i)/N2tot, zft);
        features(:, 1)  = ft1;
        features(:,2:3) = [ft2', ft3'];
        features(:,4) = my_eff(i);
        fprintf(fid1, '%17.9e, \t %17.9e, \t %17.9e, \t %17.9e \r\n', features');
    end
    fclose(fid1);
end
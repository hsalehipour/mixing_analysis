function dl_osborn_write_test(fname, my_eff, epsbar,N2, kappa0,z,Lo_patch, Iu, nrand_extra)
% Purpose: write the input data for dl_osborn code.
    
    % Random number generator which generates N rand numbers within [a, b]
    rand_gen = @(a,b,N) a + (b-a).*rand(N,1);
    
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
    
    for i=1:ntime;
        Lbot = rand_gen(min(z),-Iu(i),nrand_extra);
        Ltop = rand_gen(Iu(i),max(z) ,nrand_extra);
        for nn=1:nrand_extra;
            % select only a portion of the profile
            indz = z<=Ltop(nn) & z>=Lbot(nn);
            
            % number of data points in the selected range
            N2tot = trapz(z(indz),N2(indz,i),1);
            Dp    = kappa0*N2tot;
            features = zeros(nz_profile, 4);
            
            % interpolate to the selected z-range
            zft = linspace(min(z(indz)),max(z(indz)),nz_profile);
            ft1 = zft./Lo_patch(i);
            ft2 = interp1(z(indz), DD(indz,i)/Dp, zft);
            ft3 = interp1(z(indz), N2(indz,i)/N2tot, zft);
            features(:,1)   = ft1;
            features(:,2:3) = [ft2', ft3'];
            features(:,4)   = my_eff(i);
            fprintf(fid1, '%17.9e, \t %17.9e, \t %17.9e, \t %17.9e \r\n', features');
        end
    end
    fclose(fid1);
end
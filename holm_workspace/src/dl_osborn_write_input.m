function dl_osborn_write_input(my_eff, epsbar,chi_pr,N2,kappa0,z)
    N2avg = mean(N2,1)';
    Dp    = kappa0*N2avg;
    ntime = size(N2,2);

    fid1 = fopen('dl_input.dat', 'w');   
    nz_profile = 512;
    zft = linspace(min(z), max(z), nz_profile);
    XX = 2*kappa0*chi_pr;    
    DD = epsbar;
    XX = max(0., XX);      % X  , scalar dissip.
    DD = max(0., DD);      % eps, turbul dissip.
    for i=1:ntime;
        features = zeros(nz_profile, 4);
        ft1 = interp1(z, XX(:,i), zft);
        ft2 = interp1(z, DD(:,i)/Dp(i), zft);
        ft3 = interp1(z, N2(:,i)/N2avg(i), zft);
        features(:,1:3) = [ft1', ft2', ft3']; 
        features(:,4) = my_eff(i);
        fprintf(fid1, '%17.9e, \t %17.9e, \t %17.9e, \t %17.9e \r\n', features');
    end
end
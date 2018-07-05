function f_noised = add_rand_lognormal_noise(f, amp)
    switch nargin
        case 1
            amp=0.01;
    end            
    m = mean(f);
    v = var(f);
    mu = log((m^2)/sqrt(v+m^2));
    sigma = sqrt(log(v/(m^2)+1));
    f_noised = f + amp*lognrnd(mu, sigma,size(f));
end    
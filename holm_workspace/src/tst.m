x_inp = fit_x;
y_inp = fit_y;
nbins = 100;
nx = length(fit_x);

x_range = linspace(min(x_inp),max(x_inp),nbins+1);
x_out = zeros(nbins,1);
y_out = zeros(nbins,1);
cbin  = zeros(nbins,1);

for ibin=1:nbins
    for i=1:nx
        if (x_inp(i) >= x_range(ibin) && x_inp(i) < x_range(ibin+1))
            cbin(ibin) = cbin(ibin)+1;
            x_out(ibin) = x_out(ibin)+x_inp(i);
            y_out(ibin) = y_out(ibin)+y_inp(i);
        end
    end
end
    
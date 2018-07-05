% Purpose: Plot the histogram of Argo data in log scale with shaded areas.
function hndl = argo_PlotHist(yy,nbin,ecolor,fcolor,iflog)

if (iflog)
    binrange = 10.^(linspace(round(min(log10(yy))),round(max(log10(yy))),nbin));
else
    binrange = linspace(min(yy),max(yy),nbin);
end

[n,xout] = hist(yy(:),binrange,nbin);

if (iflog); set(gca, 'XScale', 'log');  end;
hold all;
opts={'EdgeColor', ecolor,...
    'FaceColor', fcolor};
plot(xout,n,'Color',ecolor,'LineWidth',2)
[y1handle, y2handle, hndl] = fill_between(xout,n,n*0,[],opts{:});
set(y1handle,'Visible','off');
set(y2handle,'Visible','off');

end



  
% Purpose: Plot the histogram of data in log scale with shaded areas.
function hndl = plotHist(yy,nbin,ecolor,fcolor,binmin,binmax,iflog)

if (iflog)
    binrange = 10.^(linspace(binmin,binmax,nbin));
else
    binrange = linspace(binmin,binmax,nbin);
end

[n,xout] = hist(yy(:),binrange,nbin);

if (iflog); set(gca, 'XScale', 'log');  end;
hold all;
opts={'EdgeColor', ecolor,...
    'FaceColor', fcolor};
plot(xout,n/trapz(xout,n),'Color',ecolor,'LineWidth',2)
% [y1handle, y2handle, hndl] = fill_between(xout,n,n*0,[],opts{:});
% set(y1handle,'Visible','off');
% set(y2handle,'Visible','off');

end



  
function [h1, h2] = plot_histogram2d(Xi,Yi,Zi,xbins,ybins)

nxbins = length(xbins);
nybins = length(ybins);

% map X/Y values to bin indices
ind_Xi = round( interp1(xbins, 1:nxbins, Xi(:), 'linear', 'extrap') );
ind_Yi = round( interp1(ybins, 1:nybins, Yi(:), 'linear', 'extrap') );
 
% limit indices to the range [1,numBins]
ind_Xi = max( min(ind_Xi,nxbins), 1);
ind_Yi = max( min(ind_Yi,nybins), 1);
 
% count number of elements in each bin
ncount= accumarray([ind_Yi(:) ind_Xi(:)], 1    , [nybins nxbins]);
if (~isempty(Zi))
zbins = accumarray([ind_Yi(:) ind_Xi(:)], Zi(:), [nybins nxbins],@mode);
end

% plot 2D histogram of intensity
figure;
ncount(ncount==0)=nan;
pcolor(xbins, ybins, ncount);   shading flat;
h1 = gca;

% plot 2D histograom of a 3rd parameter
if (~isempty(Zi))
figure;
zbins(zbins==0)=nan;
pcolor(xbins, ybins, zbins);   shading flat;
h2 = gca;
end

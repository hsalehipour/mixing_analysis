function [cf] = ssfit (x,y,r)
% Purpose : creat a smoothed spline fit with smoothing parameter "r"

% --- Create fit 
fopt= fitoptions('method','SmoothingSpline','SmoothingParam',r);
ifok = isfinite(x) & isfinite(y);
if ~all( ifok )
    warning( 'GenerateMFile:IgnoringNansAndInfs',...
        'Ignoring NaNs and Infs in data.' );
end
ftype = fittype('smoothingspline');

% Fit this model using new data
cf = fit(x(ifok),y(ifok),ftype,fopt);

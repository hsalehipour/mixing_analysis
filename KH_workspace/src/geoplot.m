% Purpose : Use the mapping toolbox to plot a series of data overlaid on
% top of the topography in grayscale given. The data is given in "Z" with a
% Reference vector "ZRef". The specified "projection" is also used.
% "ZMin" and "ZMax" determine the lower and upper bounds of the colorbar
% associated with the data.
% Written by: Hesam Salehipour, Feb 2015.
%
% Useful sources:
% 1- Understanding Raster Geodata in MATLAB help
% 2-
% http://www.mathworks.com/matlabcentral/answers/101346-how-do-i-use-multip
% le-colormaps-in-a-single-figure


function [a_coef, b_coef, hndl] = geoplot(Z,ZRef, ZMin, ZMax, projection)

load topo
bg = -topo;
% bg(bg<0)=0;

% Set up the axis properties
% projection = 'robinson';
axesm(projection)
gridm; framem;
axis off;

% Plot the background topographic image
h(1) = geoshow(bg,topolegend,'DisplayType','texturemap');
hold on;

% Plot the data
h(2) = geoshow(Z, ZRef, ...
    'DisplayType', 'surface',...
    'CData',Z);
set(h(2),'ZData',10 + 0*Z)
hold off;

% Initially, both CDatas are equal to Z.
m = 64;  % 256-elements is each colormap
colormap([gray(m); jet(m)])

% CData for background
bgmin = min(bg(:));
bgmax = max(bg(:));
C1 = min(m,round((m-1)*(bg-bgmin)/(bgmax-bgmin))+1);

% CData for data
% cmin = min(Z(:));
% cmax = max(Z(:));
cmin = ZMin;            cmax=ZMax;
C2 = round((m-1)*(Z-cmin)/(cmax-cmin))+m+1;
C2(C2<m+1) = m+1;       C2(C2>2*m) = 2*m;

% Update the CDatas for each object.
set(h(1),'CData',C1);
set(h(2),'CData',C2);

% Change the CLim property of axes so that it spans the
% CDatas of both objects.
caxis([min(C1(:)) max(C2(:))])

% Set the colorbar properties
hndl = colorbar;
set(hndl, 'ylim', [min(C2(:)) max(C2(:))])

% prefactor coefficients for linearly interpolating the colorbar into its actual value.
a_coef = (cmax-cmin)/(m-1);      
b_coef = cmin - (cmax-cmin)*(m+1)/(m-1);

% mask the land surfaces in black 
% geoshow('landareas.shp', 'FaceColor', 'black');

% Set the axis properties.
axesm('robinson',...
    'MeridianLabel','on',...
    'MLabelLocation',60,...
    'MLabelParallel','south',...
    'ParallelLabel', 'on',...
    'PLabelLocation', 30,...
    'PLabelMeridian', 'west');

end
%exportfig(gcf, 'fig/tst.pdf', 'FontMode', 'fixed', 'FontSize', 6, 'color',
%'cmyk', 'height', 4, 'width',6, 'Resolution', 1024);
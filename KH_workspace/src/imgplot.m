% Purpose: To use imagesc to plot the data with a reasnable appearance!
%          This may need further cosmetics
% Note: in the Argo data: lat = y and  lon = x

function [hndl]= imgplot(lon,lat,data, cmin, cmax)

data2   = [data(:,lon>180)  , data(:,lon<180)];
lon2    = [lon(lon>180)-360 , lon(lon<180)];

h=imagesc(lon2,lat,data2);
set(h,'alphadata',~isnan(data2))
axis xy;
xlim([-180 180]);   ylim([-80 80]);
caxis([cmin,cmax]);
hndl = colorbar;
geoshow('landareas.shp','FaceColor',[0.5 0.5 0.5]);

set(gca,'YMinorTick','on',...
      'XMinorTick','on',...
      'XTick',-180:30:180,...
      'YTick',[-60 -30 0 30 60],...
      'LineWidth',2);

% Append Degree Symbolts to X and Y Axes of a Plot
[hx,hy] = format_ticks(gca,'^{\circ}','^{\circ}');

set(gca,'XTickMode','manual', 'YTickMode','manual',...
    'XTickLabel','','YTickLabel','');

end

% set(h,'Location','NorthOutside','XAxisLocation','top');

% set(gca,...
%     'XTickMode','manual', 'YTickMode','manual',...
%     'XTick',-180:30:180,...
%     'YTick',[-90 -60 -30 0 30 60 90],...
%     'LineWidth',2)

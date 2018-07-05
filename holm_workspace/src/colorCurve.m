function [] = colorCurve()
% Purpose : color code the curves inside a figure

hf = get(gca,'Children');
ncurve = length(hf);
palette = distinguishable_colors(ncurve);
for i=1:ncurve
    ih = ncurve - i + 1;
    set(hf(ih),  'Color', palette(i,:));
end
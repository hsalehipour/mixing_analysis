function drawvline(xp,max_curve,min_curve)
% Purpose: for the given curve(s) specified with maxucrve, mincurve, 
% plot a straight vertical line at xp; length(xp)=1

maxValue = max_curve;
minValue = min_curve;

if (length(max_curve) > 1); maxValue = max(max_curve); end;
if (length(min_curve) > 1); minValue = min(min_curve); end;

fac = 0.25;
fac_min = 1 - sign(minValue)*fac;
fac_max = 1 + sign(maxValue)*fac;

y = [fac_min*minValue,fac_max*maxValue];
hold on;
plot([xp xp],y,'k--','LineWidth',0.5)
end

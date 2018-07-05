function []=shadedpatch(x1,x2,y1,y2,fstr)
% Purpose: draws a patch of a rectangular shaped between x=x1 and x=x2
% and y=y1 and y=y2;
%

px=[x1 x2 x2 x1]'; % make closed patch
py=[y1 y1 y2 y2]';
patch(px,py,1,'FaceColor',fstr,'EdgeColor','none');
 
alpha(.2); % make patch transparent
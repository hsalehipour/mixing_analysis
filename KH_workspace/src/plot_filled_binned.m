function hndl = plot_filled_binned(xdata,ydata,errdata,ecolor,fcolor)

% Purpose: Plot data (xdata,ydata) as a curve in log scale with ecolor 
% (edge color) with two surrounding lines filled with fcolor (face color). 
% The plot handle is returned. 
% Written by H. Salehipour (May 17, 2016). 

loglog(xdata,ydata,'-','Color',ecolor,'LineWidth',2); 
hold on; 

% options for shaded region of errors
opts={'EdgeColor', ecolor,...       
      'FaceColor', fcolor};
err_p = ydata + errdata;
err_m = ydata - errdata;
[y1handle, y2handle, hndl] = fill_between(xdata,err_p, err_m, [], opts{:});  
set(y1handle,'Visible','off');
set(y2handle,'Visible','off');  

end


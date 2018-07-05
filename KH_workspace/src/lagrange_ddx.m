function ddx_f = lagrange_ddx(f,x)
% Take second order first derivative of a function f
% NOTE: For general non-uniform grids (x)

nx = length(x);
ddx_f= x*0;

% find the surrounding three points;
for i=1:nx;
    if i==1;
        x0  = x(1);   x1=x(2);     x2=x(3);
        fx0 = f(1);  fx1=f(2);    fx2=f(3);
        xj = x0;
    elseif i == nx
        x0 = x(nx-2);   x1=x(nx-1);     x2=x(nx);
       fx0 = f(nx-2);  fx1=f(nx-1);    fx2=f(nx);
        xj = x2;
    else
        x0 = x(i-1);   x1=x(i);     x2=x(i+1);
       fx0 = f(i-1);  fx1=f(i);    fx2=f(i+1);
        xj = x1;
    end
    ddx_f(i)  = (2*xj-x1-x2)./(x0-x1)./(x0-x2).*fx0 + ...
                (2*xj-x0-x2)./(x1-x0)./(x1-x2).*fx1 + ...
                (2*xj-x0-x1)./(x2-x0)./(x2-x1).*fx2;
end
     

return
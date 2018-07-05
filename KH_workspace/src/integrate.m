function intF = integrate(x,fx)
nx = length(x);
intF = fx*0;
for i=2:nx
     intF(i)  = trapz(x(1:i),fx(1:i));
end
function ddx_f = ddx(f,x)
% Take fourth order FD derivative of a function f
% NOTE: only for uniform grid (x)

% 2nd order lagrange interpolation, 1st derivative, general grid
ddx_f = f*0;
dim_f = size(f);
dim_x = size(x);
xx = x(:,1);
if dim_f(1)==dim_x(1)    
    for i=1:dim_f(2)
        ff = f(:,i);
        if (dim_x(2)~=1); xx = x(:,i); end
        ddx_f(:,i) = lagrange_ddx(ff,xx);
    end
end

% 4rth order uniform spacing, finite diff
% Nx = length(x);
% A=zeros(Nx,Nx);
% %order = 4;
% 
% for i=3:Nx-3
%     A(i,i-2) = 1/12;
%     A(i,i-1) = -2/3;
%     A(i,i+1) = 2/3;
%     A(i,i+2) = -1/12;
% end
% for i=1:2
%     A(i,i) = -25/12;
%     A(i,i+1) = 4;
%     A(i,i+2) = -3;
%     A(i,i+3) = 4/3;
%     A(i,i+4) = -1/4;
% end
% for i=Nx-2:Nx
%     A(i,i) = 25/12;
%     A(i,i-1) = -4;
%     A(i,i-2) = 3;
%     A(i,i-3) = -4/3;
%     A(i,i-4) = 1/4;
% end
% 
% % dx = x(2)- x(1);
% dx = x(end)- x(end-1);
% 
% if size(f,1)==1;    f=f';   end
% ddx_f = A*f/dx;

return

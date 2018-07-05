clear all; clc; close all;
nz= 5e3;
z=linspace(-15,15,nz);
rho=1-tanh(z);
%rho = 2*rand(1,nz);


% Sort using bin pdf approach
np=250;
%nbins = 3*(np+1)-2;
nbins = 3*np-1;

rmax=max(rho);
rmin=min(rho);
dr=rmax-rmin;

theta=zeros(1,np+1);
rhob =zeros(1,np+1);
for i=1:np+1
    theta(i)=pi*(i-1)/np;
    rhob(i)=rmin+dr/2*(1+cos(theta(i)));
end
theta = theta/2;
rhobmax = rhob(2) + (rhob(1)-rhob(2))*cos(theta);
rhobmin = rhob(np)-(rhob(np)-rhob(np+1))*cos(theta(np+1:-1:1));
%rhobmid = 0.5*(rhob(1:np)+rhob(2:np+1));
rhobmid = rhob(1:np);

rho_base=zeros(1,nbins);
rpdf = zeros(1,nbins);
dzi = z(2)-z(1);

rho_base(np+1:2*np-2) = rhobmid(2:np-1);
for i=1:np
    rho_base(i) = (rhobmax(i)+rhobmax(i+1))/2;
    rho_base(i+2*np-2) = (rhobmin(i)+rhobmin(i+1))/2;
    for j=1:nz
        if (rho(j)>rhobmax(i+1) && rho(j)<=rhobmax(i))
            rpdf(i) = rpdf(i) + dzi;
        elseif (rho(j)>rhobmin(i+1) && rho(j)<=rhobmin(i))          
            rpdf(i+2*np-2) = rpdf(i+2*np-2) + dzi;
        end
    end
end
for i=2:np-1
    for j=1:nz
        if (rho(j)>rhob(i+1) && rho(j)<=rhob(i))
            rpdf(i+np-1) = rpdf(i+np-1) + dzi;
        end
    end
end


zb=zeros(1,nbins);
for i=2:nbins      % Running sum of z - eliminate voids
    zb(i)=zb(i-1)+rpdf(i-1);
end

plot(rho_base,zb,'bo-')
hold on;
plot(rho,z+15,'r')

% zb=zeros(1,nbins);
% k=1;
% for i=1:nbins      % Running sum of z - eliminate voids
%     if (rpdf(i)>0) 
%         k=k+1;
%         zb(k)=zb(k-1)+rpdf(i);
%     end
% end
% 
% plot(rho_base(1:k),zb(1:k),'bo-')
% hold on;
% plot(rho,z+15,'r')
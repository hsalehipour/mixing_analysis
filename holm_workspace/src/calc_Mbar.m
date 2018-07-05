% Vertical profiles of irreversible mixing,
% Ref: Salehipour & Peltier (2015)
% Note: chi_ref=<|grad rho|>_{S*};

Dpbar = -kappa0*g/rho0*gradrhob(:,1);
Mbar  = -kappa0*g/rho0*chi_ref./gradrhob - Dpbar*ones(1,ntime);
Mbar(isnan(Mbar))=0.0;
for i=1:ntime      
      res = abs(mean(Mbar(:,i))-M(i));
      tol = 1e-1;
      itr = 0;
      time(i)
%       while(res/M(i)>0.1 && tol<=2/Irho(i))
%           tol = tol+1e-2;
          indrho_nonzero = abs(gradrhob(:,i))>=tol;
          if(~isempty(find(indrho_nonzero)==1))
          tzmax = max(z(indrho_nonzero),[],1)';
          tzmin = min(z(indrho_nonzero),[],1)';
          indz_outlier = z>=tzmax | z<=tzmin;
          Mbar(indz_outlier,i)=0.0;
          end
%           Mbar(:,i) = smooth(Mbar(:,i),17);
          res = abs(mean(Mbar(:,i))-M(i));
          itr = itr+1;
%           figure(1); hold all;
%           plot(itr,res/M(i),'o');
%       end
%       pause;
end


% Dpbar = -kappa0*g/rho0*gradrhob(:,1);
% Mbar  = -kappa0*g/rho0*chi_ref./gradrhob - Dpbar*ones(1,ntime);
% Mbar(isnan(Mbar))=0.0;
% Mbar_avg = zeros(nz/2,1);
% for i=1:ntime      
%       for nn=1:nz/2;
%           itr=0.;
%           Mbar_avg(nn)=sum(Mbar(nz/2-nn+1:nz/2+nn,i))/nz;
%       end
%       [~, indmatch]=min(abs(Mbar_avg-M(i)));
%       indz_inside = z<=indmatch/2*dz & z>=-indmatch/2*dz;
%       Mbar(~indz_inside,i) = 0.0;
% 
% end

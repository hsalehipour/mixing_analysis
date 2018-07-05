% Purpose: Write the initial field based on the results of the stability
% analysis

function output_initial_eigv(fname, k_fgm,z,uz_eigv,rho_eigv)

ij = sqrt(-1);
nz = length(z);

%1st derivative matrix with 1-sided boundary terms
D1=ddz(z); 

% Continuity: ux_eigv based on uz_eigv 
ux_eigv = ij/k_fgm*(D1*uz_eigv);    

% Write the initial fields
A = [z';    real(ux_eigv)'; imag(ux_eigv)';...
            real(uz_eigv)'; imag(uz_eigv)';...
            real(rho_eigv)';imag(rho_eigv)'];

fid=fopen(fname,'w');
fprintf(fid,'%g \t %g\n',k_fgm, nz);
fprintf(fid,'Z\t ux\t uz \t rho\n');

formatSpec = '%12.8e (%12.8e,%12.8e) (%12.8e,%12.8e) (%12.8e,%12.8e)\n';
fprintf(fid,formatSpec,A);
end


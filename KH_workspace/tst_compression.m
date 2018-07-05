close all; clc;

% fname = 'kh_t168.okc';
% nelx = 47;
% nelz = 53;

fname = 'holm_t184.okc';
% fname = 'holm_t274.okc';
nelx = 71;
nelz = 116;

% Read data
A = loadtxt(fname, 6, 13);
x_in  = A(1,:)';
z_in  = A(2,:)';
vortx_in  = A(4,:)';
vorty_in  = A(5,:)';
vortz_in  = A(6,:)';

[vorty_out, x_out,z_out] = xmdv_reader(vorty_in, x_in, z_in, nelx, nelz);
vortx_out = xmdv_reader(vortx_in, x_in, z_in, nelx, nelz);
vortz_out = xmdv_reader(vortz_in, x_in, z_in, nelx, nelz);

% calculate the magnitude of the total vorticity field
vort_tot = sqrt(vortx_out.^2+vorty_out.^2+vortz_out.^2);

N = 9;                               % Decomposition Level 
wname = 'coif5';                     % Near symmetric wavelet
[c, s] = wavedec2(vort_tot,N,wname);    % Multilevel 2-D wavelet decomposition.
opt = 'gbl';                    % Global threshold
%thr = 2;                       % Threshold
Z = 0.5*mean(vort_tot(:).^2);
thr = sqrt(4/3*Z*(log(numel(vort_tot))));
% thr = var(vorty_out(:))*sqrt(2*log(numel(vorty_out)));
sorh = 'h';                     % Hard thresholding
keepapp = 0;                    % Approximation coefficients cannot be thresholded
[vort_cmp,cden,sden,perf0,perfl2] = wdencmp(opt,c,s,wname,N,thr,sorh,keepapp);

figure; mesh(x_out,z_out,vort_tot); 
view(0,90); axis equal;
title('Original Image')
caxis([-2 2]); colorbar

figure; mesh(x_out,z_out,vort_cmp); 
view(0,90); axis equal;
title(['Compressed Image - Global Threshold = ',num2str(thr)])
caxis([-2 2]); colorbar

figure; mesh(x_out,z_out,vort_tot-vort_cmp); 
view(0,90); axis equal;
title('Filtered "Noise"')
caxis([-2 2]); colorbar

perf0
perfl2


% PSD = zeros(N,1);
% for i=1:N
%     A1 = appcoef2(c,s,wname,i);
%     figure;
%     imagesc(flipud(A1));
%     caxis([0 12]); colorbar
%     PSD(i) = mean(A1(:).^2);
% end
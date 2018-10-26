
% relevant parameters of a given conv-layer that has a size 
% [ntime,nprofile,nz,nfilter]
% ntime is given by number of time step in one DNS simulation for example;
nz          = 128;
ntime       = 267;
nprofile    = 2;
nfilter     = 64;

% note: conv_output.dat file has been created by adding proper hooks
% in the eval mode of my TensorFlow implementation of CNN;
fname = [fadrs,'conv3_output.dat'];
A = loadtxt(fname,2,1);
conv_output = reshape(A(2,:), [ntime, nprofile, nz ,nfilter]);

% Read original DNS data that could be fed into CNN network
B = dlmread([fadrs,'dl_input.dat'],',',1,0);
zz = reshape(B(:,1), [512, ntime]);
zz = imresize(zz,[nz, ntime]);
tt = ones(nz,1)*time';

for i=1:2:nfilter; 
    disp_filtered = squeeze(conv_output(:,1,:,i))';
    N2_filtered   = squeeze(conv_output(:,2,:,i))';
    figure; 
    subplot(2,1,1); 
    contourf(tt,zz, N2_filtered  ,128,'LineStyle','none'); 
    colorbar; 
    subplot(2,1,2); 
    contourf(tt,zz, disp_filtered,128,'LineStyle','none'); 
    colorbar; 
end;

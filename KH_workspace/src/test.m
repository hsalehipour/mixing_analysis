% An idea script to parameterize Ri_g based on Reb thus having a peicewise
% parameterization for epsilon based on N2 and S2! 
% would that be useful?

%y= Ri
%x= Reb

% clear all; clc; close all;

x = 3:0.1:1e4;
ind1 = x<=100;          x1=x(ind1);
ind2 = x>100 & x<1e3;   x2=x(ind2);
ind3 = x>1e3;           x3=x(ind3);

y(ind1) = (-0.4*x1.^(-0.9)+0.2).*(167.4*x1.^(-2.9) + 3.6*x1.^(-0.2));
y(ind2) = (2*x2.^(-0.5)).*(167.4*x2.^(-2.9)+3.6*x2.^(-0.2));
y(ind3) = 2*x3.^(-0.5);

figure(1)
semilogy(6*y,x)
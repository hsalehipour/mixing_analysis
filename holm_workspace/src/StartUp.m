%% 0) Simulation startup

% Problem Setup
fname = [fadrs,'siminfo.dat'];
fid = fopen(fname);     A=fscanf(fid, '%f');    fclose(fid); 

% Flow parameters
Re = A(1);      Ri = A(2);      Pr = A(3);      R = A(4);
nu = 1/Re;
kappa0 = nu./Pr;

% Geometrical setup
% NOTE: Lz is nonunifrom, use the elements inside middle box
xmin = A(5);    xmax=A(6);      Lx = xmax-xmin; 
ymin = A(7);    ymax=A(8);      Ly = ymax-ymin;
zmin = A(9);    zmax=A(10);     Lz = zmax-zmin;

nelx = A(11);   nely = A(12);   nelz = A(13);
norder = A(14);
nzp = A(15);


% defining constants
% g = 9.81;
% rho0 = 1000.0;
g = Ri/R;
rho0 = 1;
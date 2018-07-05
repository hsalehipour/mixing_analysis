% A = loadtxt('../data/test/Re4000/s.0000002.dat', 11, 2);
% A = loadtxt('../data/test/Re500-hr/s.0000002.dat', 11, 2);
% A = loadtxt('../data/test/Re500/s.0000002.dat', 11, 2);
% A = loadtxt('../data/test/s.0000002.dat', 11, 2);

nkx = A(2,1);
nky = A(3,1);
nkz = A(4,1);
kx  = A(6,:);
ky  = A(7,:);
kz  = A(8,:);
Ekx = A(9,:);
Eky = A(10,:);
Ekz = A(11,:);

figure(1); hold all;
loglog(kx,Ekx,'r-o')
loglog(ky,Eky,'b-o')
loglog(kz,Ekz,'g-o')

% loglog(kx,Ekx,'r-*')
% loglog(ky,Eky,'b-*')
% loglog(kz,Ekz,'g-*')

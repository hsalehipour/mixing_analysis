import math

# Base resolution Ref: Mashayek & Peltier JFM, 2013
# Kolmogorov length scale: L_k = L/(Re)^{3/4}
# Batchelor  length scale: L_b = L_k/sqrt(Pr) = L/(Re)^{3/4}/Pr^{0.5}
#Nx0 = 76;      Ny0 = 10;       Nz0 = 27;
#Nx0 = 129;      Ny0 = 28;       Nz0 = 91;
#Nx0 = 90;      Ny0 = 18;       Nz0 = 62;   # for Re=6000, Pr=4
Nx0 = 64;      Ny0 = 13;       Nz0 = 44;    # for Re=6000, Pr=1


Ntot0= Nx0*Ny0*Nz0;
Re0 = 6000.0
Pr0 = 2.0
Lx0 = 14.27

ndim = 3
Re = 48000.0
Pr = 2.0
Lx = 2.5040

fac = (Re/Re0)**0.75*(Pr/Pr0)**0.5
Nx = Nx0*fac*Lx/Lx0
Ny = Ny0*fac
Nz = Nz0*fac
if ndim == 2: Nz = 1;
Ntot = Nx*Ny*Nz;

print "For Re = ", Re, " and Pr = ", Pr, " then comp cost is ", fac**(ndim+1), "more."
print "Nx = ", Nx
print "Ny = ", Ny
print "Nz = ", Nz
print 'Ntot = ', Ntot/1e6,  'million'
print 'Ntot : ', fac**ndim, 'x more'
print 'DT   : ', fac,       'x less'
print 'NP   : ', fac**(ndim+1), 'x more.'
print 'For best Nek performance of 50 elements per core:'
print 'Np <=', Ntot/50, 'processors are required'

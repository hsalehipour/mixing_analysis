import math

# Finding the number of elements in the turbulent region (Ri-based)
# Only the following 4 lines would possibly change
Ri = 0.12;
Nelybox2 = 125   # num of y elements inside box2
r = 1.25         # expansion ratio
Hbox1 = 10.      # depends on Ri



#===================================================
# DO NOT TOUCH HERE

Htot = 30.       # height of our domain (fixed)
Hbox3 = Hbox1
Hbox2 = Htot-Hbox1-Hbox3
DY = Hbox2/Nelybox2

# Geometric series \sum ar^k = a(1-r^n)/(1-r)
# Looking for n here
nom = math.log10(Hbox1/DY*(r-1.) + 1.)/math.log10(r)
Nelbox1 = round(nom)


print "Nely = ", Nelbox1

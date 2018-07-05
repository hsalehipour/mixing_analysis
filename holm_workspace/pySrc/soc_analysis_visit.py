#================================
# Written by Hesam Salehipour
# The purpose of this code is to drive visit in order to analyze the HWI data for the SOC paper
# Date: 18 Dec 2017
#================================

visit_path           = "/usr/local/visit/bin"
visit_python_package = "/usr/local/visit/current/linux-x86_64/lib/site-packages"

import sys
sys.path.append(visit_python_package)
import numpy as np
from visit import *
Launch(vdir=visit_path)


#================================
# User changes this line
simname ="Re-6000-Ri-016-Pr-8"
shear_layer   = False
density_layer = True
#================================


# Initial setup
fname = "viz01:/bg01/homescinet/p/peltier/hsalehip/bgq-scratch/holm/"+simname+"/pp/socppholm.nek5000"
Mesh = "mesh"

# Read the position of upper and lower flanks of shear and density layers
siminfo     = np.loadtxt(simname+"/siminfo.dat")
data_flanks = np.loadtxt(simname+"/flanks.dat")
Irho = data_flanks[:,1]
Iu   = data_flanks[:,2]
xmin, xmax = siminfo[4:6]
if shear_layer:
    zl = Iu/2.
elif density_layer:
    zl = Irho/2.

p0 = (xmin,  zl[0], 0.0)
p1 = (xmax,  zl[0], 0.0)
p2 = (xmin, -zl[0], 0.0)
p3 = (xmax, -zl[0], 0.0)
p4 = (xmin, 0.0, 0.0)
p5 = (xmax, 0.0, 0.0)

# Define the scalar variables
DefineScalarExpression("xc", "coord(%s)[0]" % Mesh)
DefineScalarExpression("zc", "coord(%s)[2]" % Mesh)
DefineScalarExpression("vort3d_mag_sqr", "velocity_mag^2")
DefineScalarExpression("vort2d_mag_sqr", "temperature^2")

# Read database
OpenDatabase(fname)

# add the main plot (Window 1)
AddPlot("Pseudocolor", "vort3d_mag_sqr", 1, 1)
AddOperator("Slice", 1)
DrawPlots()

# add the lineout curve plots (Window 2)
Lineout(p0, p1, ("xc", "zc", "vort2d_mag_sqr", "default"))
Lineout(p2, p3, ("vort2d_mag_sqr", "default"))
Lineout(p4, p5, ("vort2d_mag_sqr", "default"))

# loop over all time steps and store the line-out data
for i in range(TimeSliderGetNStates()):

    # Change time slider
    SetTimeSliderState(i)

    # activate curve plots
    SetActiveWindow(2)

    # This is necessary when the line-out points are time-dependent
    # top flank
    SetActivePlots((0,1,2,3))
    a1 = LineoutAttributes()
    a1.point1 = (xmin, zl[i], 0.0)
    a1.point2 = (xmax, zl[i], 0.0)
    SetOperatorOptions(a1)

    # bottom flank
    SetActivePlots((4,5))
    a2 = LineoutAttributes()
    a2.point1 = (xmin, -zl[i], 0.0)
    a2.point2 = (xmax, -zl[i], 0.0)
    SetOperatorOptions(a2)

    # Get the (x-z) coord locations of the line-out operation
    SetActivePlots(0)
    xc = GetPlotInformation()["Curve"]
    SetActivePlots(1)
    zc = GetPlotInformation()["Curve"]

    # get the values of the plotted field at the line-out operation
    SetActivePlots(2)
    vort2d_top = GetPlotInformation()["Curve"]
    SetActivePlots(3)
    vort3d_top = GetPlotInformation()["Curve"]

    SetActivePlots(4)
    vort2d_bot = GetPlotInformation()["Curve"]
    SetActivePlots(5)
    vort3d_bot = GetPlotInformation()["Curve"]

    # middle interface
    SetActivePlots(6)
    vort2d_z0 = GetPlotInformation()["Curve"]
    SetActivePlots(7)
    vort3d_z0 = GetPlotInformation()["Curve"]

    filename = simname + '/soc.'+str(i).zfill(4) + '.dat'
    headertxt = 'x \t z \t vort2d(top) \t vort2d(bot) \t vort3d(top) \t vort3d(bot) \t vort2d(z=0) \t vort3d(z=0)'
    data =[        [xc[2 * idx + 1],
                    zc[2 * idx + 1],
            vort2d_top[2 * idx + 1],
            vort2d_bot[2 * idx + 1],
            vort3d_top[2 * idx + 1],
            vort3d_bot[2 * idx + 1],
            vort2d_z0 [2 * idx + 1],
            vort3d_z0 [2 * idx + 1]] for idx in range(len(xc) / 2)]

    data = np.array(data)
    np.savetxt(filename, data, header=headertxt.expandtabs(16))








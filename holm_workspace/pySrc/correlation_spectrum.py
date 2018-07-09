#================================
# Written by Hesam Salehipour
# The purpose of this code is to drive visit in order to analyze the HWI data for the SOC paper
# Date: 18 Dec 2017
#================================
from visitOps import *

#================================
# User changes this line
simname ="Re-6000-Ri-016-Pr-8"
shear_layer   = False
density_layer = True
#================================

# Initial setup
path_remote = 'niagara.scinet.utoronto.ca:/gpfs/fs1/home/p/peltier/hsalehip/scratch/soc_paper/holm/'
path_local  = '/home/hesam/REPOs/mixing_analysis/holm_workspace/data/'
dbname = path_remote + 'holm.nek5000'

# Read the position of upper and lower flanks of shear and density layers
# siminfo     = np.loadtxt(path_local+simname+"/siminfo.dat")
data_flanks = np.loadtxt(path_local+simname+"/flanks.dat")
Irho = data_flanks[:,1]
Iu   = data_flanks[:,2]
if shear_layer:
    zl = Iu/2.
elif density_layer:
    zl = Irho/2.

# Read database
rho = 'temperature'
zvel= 'z_velocity'
nbins = 1000
launch_visit(visit_path)
visit.OpenDatabase(dbname)

# Apply initial driver to get the points for the linout operation
driver = VisitSetupBase()
driver.plotSlice(rho)
xmin, xmax, zmin, zmax = driver.mesh_spatial_extents()

# caculate rho' = rho - <rho>_{xy}
rho_bar   = calc_ubar(rho, nbins)
rho_fluct = calc_ufluct(rho, rho_bar)
visit.ChangeActivePlotsVar(rho_fluct)

# calculate w' = w - <w>_{xy}
zvel_bar   = calc_ubar(zvel, nbins)
zvel_fluct = calc_ufluct(zvel, zvel_bar)

# define the initial points
p0 = (xmin, zl[0], 0.0)
p1 = (xmax, zl[0], 0.0)
p2 = (xmin, -zl[0], 0.0)
p3 = (xmax, -zl[0], 0.0)
p4 = (xmin, 0.0, 0.0)
p5 = (xmax, 0.0, 0.0)

# add the lineout curve plots (Window 2)
var_name_list = [rho_fluct, zvel_fluct]
line_dix = {'top_flank': (p0, p1),
            'bot_flank': (p2, p3),
            'interface': (p4, p5)}
curves = LinoutOps(line_dix, var_name_list)

# loop over all time steps and store the line-out data
for i in range(visit.TimeSliderGetNStates()):

    # Change time slider
    visit.SetTimeSliderState(i)

    if i == 0:
        curves.create()

    # extract the information at the top flank
    if i > 0:
        p0 = (xmin,  zl[i], 0.0)
        p1 = (xmax,  zl[i], 0.0)
        p2 = (xmin, -zl[i], 0.0)
        p3 = (xmax, -zl[i], 0.0)
        curves.update('top_flank', p0, p1, window_id=2)
        curves.update('bot_flank', p2, p3, window_id=2)

    # extract the information at the top flank
    rho_fluct_top = curves.extract('top_flank', rho_fluct , window_id=2)
    zvel_fluct_top= curves.extract('top_flank', zvel_fluct, window_id=2)
    rho_fluct_bot = curves.extract('bot_flank', rho_fluct , window_id=2)
    zvel_fluct_bot= curves.extract('bot_flank', zvel_fluct, window_id=2)

    # perform cross-correlation
    buoy_flux_top = np.correlate(rho_fluct_top, zvel_fluct_top, 'same')
    buoy_flux_bot = np.correlate(rho_fluct_bot, zvel_fluct_bot, 'same')


    # filename = simname + '/soc.'+str(i).zfill(4) + '.dat'
    # headertxt = 'x \t z \t vort2d(top) \t vort2d(bot) \t vort3d(top) \t vort3d(bot) \t vort2d(z=0) \t vort3d(z=0)'
    # data =[        [xc[2 * idx + 1],
    #                 zc[2 * idx + 1],
    #         vort2d_top[2 * idx + 1],
    #         vort2d_bot[2 * idx + 1],
    #         vort3d_top[2 * idx + 1],
    #         vort3d_bot[2 * idx + 1],
    #         vort2d_z0 [2 * idx + 1],
    #         vort3d_z0 [2 * idx + 1]] for idx in range(len(xc) / 2)]
    #
    # data = np.array(data)
    # np.savetxt(filename, data, header=headertxt.expandtabs(16))









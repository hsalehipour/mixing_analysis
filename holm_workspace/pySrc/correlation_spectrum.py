#================================
# Written by Hesam Salehipour
# The purpose of this code is to drive visit in order to analyze the HWI data for the SOC paper
# Date: 18 Dec 2017
#================================
from visitOps import *
import matplotlib.pyplot as plt

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
time = data_flanks[:,0]
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


# loop over all time steps and store the line-out data
nts = visit.TimeSliderGetNStates()
for ts in range(nts):
    dbtime = driver.get_times()[ts]
    print "[analyzing database at time = %d]" % dbtime

    # Change time slider
    visit.SetActiveWindow(1)
    visit.SetTimeSliderState(ts)

    # caculate rho' = rho - <rho>_{xy}
    rho_bar = calc_ubar(rho, nbins)
    rho_fluct = calc_ufluct(rho, rho_bar)
    visit.ChangeActivePlotsVar(rho_fluct)

    # calculate w' = w - <w>_{xy}
    zvel_bar = calc_ubar(zvel, nbins)
    zvel_fluct = calc_ufluct(zvel, zvel_bar)

    # define all the points for the lineout operation
    dt = time[1]-time[0]
    it = np.where(np.abs(time - dbtime) < 0.2*dt)[0][0]
    p0 = (xmin,  zl[it], 0.0)
    p1 = (xmax,  zl[it], 0.0)
    p2 = (xmin, -zl[it], 0.0)
    p3 = (xmax, -zl[it], 0.0)
    p4 = (xmin, 0.0, 0.0)
    p5 = (xmax, 0.0, 0.0)

    # create all the curves
    if ts == 0:
        # add the lineout curve plots (Window 2)
        var_name_list = [rho_fluct, zvel_fluct]
        line_dix = {'top_flank': (p0, p1),
                    'bot_flank': (p2, p3),
                    'interface': (p4, p5)}
        curves = LinoutOps(line_dix, var_name_list)
        curves.create()

    # update the curves based on time-dependent lines
    if ts > 0:
        visit.SetActiveWindow(2)
        visit.SetTimeSliderState(ts)
        curves.update('top_flank', p0, p1, window_id=2)
        curves.update('bot_flank', p2, p3, window_id=2)

    # extract the curve data
    rho_fluct_top = curves.extract('top_flank', rho_fluct , window_id=2)
    zvel_fluct_top= curves.extract('top_flank', zvel_fluct, window_id=2)
    rho_fluct_bot = curves.extract('bot_flank', rho_fluct , window_id=2)
    zvel_fluct_bot= curves.extract('bot_flank', zvel_fluct, window_id=2)

    # perform cross-correlation
    cc_buoy_flux_top = np.correlate(rho_fluct_top, zvel_fluct_top, 'same')
    cc_buoy_flux_bot = np.correlate(rho_fluct_bot, zvel_fluct_bot, 'same')

    aa = np.fft.fft(cc_buoy_flux_top)
    plt.loglog(aa * aa.conj())
    plt.show()










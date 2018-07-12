#================================
# Written by Hesam Salehipour
# The purpose of this code is to drive visit in order to analyze the HWI data for the SOC paper
# Date: 18 Dec 2017
#================================
from visitOps import *



#================================
# User changes this line
# simname ="Re-6000-Ri-016-Pr-8"
simname ="Re-6000-Ri-016-Pr-8"

# Initial setup
path_remote = 'niagara.scinet.utoronto.ca:/gpfs/fs1/home/p/peltier/hsalehip/scratch/soc_paper/holm/r3/'
path_local  = '/home/hesam/REPOs/mixing_analysis/holm_workspace/data/'
dbname = path_remote + 'holm.nek5000'


# Read the position of upper and lower flanks of shear and density layers
data_flanks = np.loadtxt(path_local+simname+"/flanks.dat")
time = data_flanks[:, 0]
Irho = data_flanks[:, 1]
Iu   = data_flanks[:, 2]

# Read database
rho = 'temperature'
zvel= 'z_velocity'
xmesh  = 'mesh_x_nodal'
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
    visit.SetTimeSliderState(ts)

    # caculate rho' = rho - <rho>_{xy}
    visit.SetActiveWindow(1)
    rho_bar = calc_ubar(rho, nbins)
    rho_fluct = calc_ufluct(rho, rho_bar)
    visit.ChangeActivePlotsVar(rho_fluct)

    # calculate w' = w - <w>_{xy}
    zvel_bar = calc_ubar(zvel, nbins)
    zvel_fluct = calc_ufluct(zvel, zvel_bar)

    # define all the points for the lineout operation
    dt = time[1]-time[0]
    it = np.where(np.abs(time - dbtime) < 0.2*dt)[0][0]
    p0 = (xmin,  Irho[it]/2., 0.0)
    p1 = (xmax,  Irho[it]/2., 0.0)
    p2 = (xmin, -Irho[it]/2 , 0.0)
    p3 = (xmax, -Irho[it]/2 , 0.0)
    p4 = (xmin,  Iu[it]/2.  , 0.0)
    p5 = (xmax,  Iu[it]/2.  , 0.0)
    p6 = (xmin, -Iu[it]/2   , 0.0)
    p7 = (xmax, -Iu[it]/2   , 0.0)
    p8 = (xmin, 0.0, 0.0)
    p9 = (xmax, 0.0, 0.0)
    line_dix = {'z=Irho/2'  : (p0, p1),
                'z=-Irho/2' : (p2, p3),
                'z=Iu/2'    : (p4, p5),
                'z=-Iu/2'   : (p6, p7),
                'z=0'       : (p8, p9)}

    # create all the curves
    if ts == 0:
        # add the lineout curve plots (Window 2)
        var_name_list = [rho_fluct, zvel_fluct]
        curves = LinoutOps(line_dix, var_name_list)
        curves.create()

        #add lineout to get x-coordinate
        centerline_dic = {'z=0': line_dix['z=0']}
        centerline = LinoutOps(centerline_dic, [xmesh])
        centerline.create()

    # update the curves based on time-dependent lines
    if ts > 0:
        visit.SetActiveWindow(2)
        visit.SetTimeSliderState(ts)
        curves.update(line_dix, window_id=2)

    # extract all the curve data in DataFrame form
    centerline_plotid = len(var_name_list) * len(line_dix)
    xmesh_values = centerline.extract('z=0', xmesh, window_id=2, plot_id=centerline_plotid)
    df = curves.extract_all(window_id=2)
    df['time'] = dbtime
    df['xmesh']= xmesh_values

    # save data for offline processing
    filename = path_local+simname + '/cc.' + str(ts).zfill(4) + '.dat'
    save_data(df, filename)



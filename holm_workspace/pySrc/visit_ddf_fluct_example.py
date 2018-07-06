##################################
#
# file: visit_ddf_fluct_example.py
#
# Example showing how to use DDFs and the "apply_ddf" expression
# to create a fluctuation field.
#
# This example uses the dataset "varying.visit" from:
#   http://visitusers.org/index.php?title=Tutorial_Data
#
# usage:
#  visit -cli -s visit_ddf_fluct_example.py varying.visit
#
#
##################################

visit_path           = "/usr/local/visit/bin"
# visit_python_package = "/usr/local/visit/current/linux-x86_64/lib/site-packages"

# import sys
# sys.path.append(visit_python_package)

from visit import *
from visit_utils import *
Launch(vdir=visit_path)

dbname = "~/Downloads/tutorial_data/varying.visit"

def ddf(atts, var_name, ddf_op):
    ddf_op_map = {"avg": atts.Average,
                  "min": atts.Minimum,
                  "max": atts.Maximum,
                  "stddev": atts.StandardDeviation,
                  "var": atts.Variance,
                  "sum": atts.Sum,
                  "count": atts.Count,
                  "rms": atts.RMS,
                  "pdf": atts.PDF}
    atts.statisticalOperator = ddf_op_map[ddf_op]
    visit.ConstructDDF(atts)
    # ndims = len(atts.numSamples)
    # ddf_oname = "%s_%s_%dd" % (var_name, ddf_op, ndims)
    # if len(atts.numSamples) == 1:
    #     src_fname = "%s.ultra" % atts.ddfName
    #     des_fname = "%s.ult" % (atts.ddfName)
    #     common.sexe("mv %s %s" % (src_fname, des_fname))
    #     lines = open(des_fname).readlines()
    #     f = open(des_fname, "w")
    #     f.write("# %s\n" % (ddf_oname))
    #     for l in lines[1:]:
    #         f.write(l)
    #     f.close()
    # else:
    #     src_fname = "%s.vtk" % atts.ddfName
    #     orig_vtk_var = "SCALARS %s float" % var_name
    #     ddf_vtk_var = "SCALARS %s float" % ddf_oname
    #     des_fname = "%s%s_%04d.vtk" % (Params.output.file_base, atts.ddfName, ts())
    #     common.sexe("mv %s %s" % (src_fname, des_fname))
    #     data = open(des_fname).read()
    #     f = open(des_fname, "w")
    #     data = data.replace(orig_vtk_var, ddf_vtk_var)
    #     f.write(data)
    # print "[ddf output: %s]" % des_fname
    # return des_fname


def active_db():
    """
    Returns the path to the active database.
    """
    return visit.GetWindowInformation().activeSource


def mesh_md():
    """
    Fetches the metadata about the active database.
    """
    return visit.GetMetaData(active_db())


def mesh_spatial_extents():
    """
    Queries the spatial extents of active plot.
    """
    return query("SpatialExtents")


def mesh_name():
    """
    Fetches the first mesh name in the active database.
    """
    return mesh_md().GetMeshes(0).name


def setup_exprs():
    """
    Sets up spatial expressions to use with DDF.
    """
    mname = mesh_name()
    exprs.define("mesh_x_nodal", "coord(%s)[0]" % mname)
    exprs.define("mesh_y_nodal", "coord(%s)[1]" % mname)
    exprs.define("mesh_z_nodal", "coord(%s)[2]" % mname)
    exprs.define("mesh_x_zonal", "recenter(coord(%s)[0])" % mname)
    exprs.define("mesh_y_zonal", "recenter(coord(%s)[1])" % mname)
    exprs.define("mesh_z_zonal", "recenter(coord(%s)[2])" % mname)


def plot_ddf_fluct(var_name, num_samples, ts=None):
    """
    Uses DDF machinery to create a spatial scalar fluctuation field.
    """
    setup_exprs()
    AddPlot("Pseudocolor", var_name)
    mname = mesh_name()
    DrawPlots()
    sext = mesh_spatial_extents()
    if not ts is None:
        ddf_name = "ddf_%s_%04d" % (var_name, ts)
    else:
        ddf_name = "ddf_%s" % (var_name)
    atts = visit.ConstructDDFAttributes()
    atts.ddfName = ddf_name
    atts.codomainName = var_name
    atts.varnames = ("mesh_z_nodal",)
    atts.ranges = (sext[4], sext[5])
    atts.numSamples = (num_samples,)
    ddf(atts, var_name, "avg")
    # After we have the ddf, we can map it back onto the mesh
    # using the apply_ddf expression
    avg_ename = "%s_avg" % var_name
    fluct_ename = "%s_fluct" % var_name
    exprs.define(avg_ename, "apply_ddf(%s,%s)" % (mname, ddf_name))
    # from here we can create an expression with the fluct
    exprs.define(fluct_ename, "%s - %s" % (var_name, avg_ename))
    ChangeActivePlotsVar(fluct_ename)
    ChangeActivePlotsVar(avg_ename)
    pass


def main():
    #
    # Open the database passed on the command line
    #
    # dbname = Argv()[0]
    OpenDatabase(dbname)
    ## if you need to use a parallel engine, visit_utils.engine
    ## provides a simple interface:
    # engine.open(nprocs=2)
    #
    # setup the plot we want to analyze
    #
    plot_ddf_fluct(var_name="temp", num_samples=50)
    #
    # Or Analyze all timesteps
    #
    # nts = TimeSliderGetNStates()
    # for ts in xrange(nts):
    #   print "[analyzing ts = %d]" % ts
    #   TimeSliderSetState(ts)
    #   DeleteAllPlots()
    #   plot_ddf_fluct(var_name = "d", num_samples = 20, ts = ts)
    #   SaveWindow()
    # sys.exit(0)


# if __visit_script_file__ == __visit_source_file__:
if __name__ == "__main__":
    main()
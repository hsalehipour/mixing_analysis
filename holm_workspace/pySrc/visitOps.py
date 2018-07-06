from visit import *
from visit_utils import exprs, query, common


visit_path = "/usr/local/visit/bin"
# visit_python_package = "/usr/local/visit/current/linux-x86_64/lib/site-packages"

# helper function to lunch visit
def launch_visit(visit_path):
    return visit.Launch(vdir=visit_path)



class VisitSetupBase(object):

    def __init__(self):
        self.mesh_name = self.get_mesh_name()
        self.setup_exprs()

    def active_db(self):
        """
        Returns the path to the active database.
        """
        return visit.GetWindowInformation().activeSource

    def mesh_md(self):
        """
        Fetches the metadata about the active database.
        """
        return visit.GetMetaData(self.active_db())

    def mesh_spatial_extents(self):
        """
        Queries the spatial extents of active plot.
        """
        return query("SpatialExtents")

    def get_mesh_name(self):
        """
        Fetches the first mesh name in the active database.
        """
        return self.mesh_md().GetMeshes(0).name

    def setup_exprs(self):
        """
        Sets up spatial expressions to use with DDF.
        """
        mesh_name = self.mesh_name
        exprs.define("mesh_x_nodal", "coord(%s)[0]" % mesh_name)
        exprs.define("mesh_y_nodal", "coord(%s)[1]" % mesh_name)
        exprs.define("mesh_z_nodal", "coord(%s)[2]" % mesh_name)
        exprs.define("mesh_x_zonal", "recenter(coord(%s)[0])" % mesh_name)
        exprs.define("mesh_y_zonal", "recenter(coord(%s)[1])" % mesh_name)
        exprs.define("mesh_z_zonal", "recenter(coord(%s)[2])" % mesh_name)

    def plotPseudocolor(self):
        visit.AddPlot("Pseudocolor", self.fld_name)
        visit.DrawPlots()


class ReductionOps(VisitSetupBase):

    def __init__(self, op, fld_name):
        super(ReductionOps, self).__init__()
        self.ddf_op = op
        self.fld_name = fld_name
        self.atts = visit.ConstructDDFAttributes()
        self.ddf_op_map = { "avg"   : self.atts.Average,
                            "min"   : self.atts.Minimum,
                            "max"   : self.atts.Maximum,
                            "stddev": self.atts.StandardDeviation,
                            "var"   : self.atts.Variance,
                            "sum"   : self.atts.Sum,
                            "count" : self.atts.Count,
                            "rms"   : self.atts.RMS,
                            "pdf"   : self.atts.PDF}
        self.atts.statisticalOperator = self.ddf_op_map[op]
        self.atts.codomainName = fld_name


    def __set_attributes__(self):
        return

    def ddf(self):
        visit.ConstructDDF(self.atts)
        return

    def output(self, ddf_op, var_name):
        ndims = len(self.atts.numSamples)
        ddf_oname = "%s_%s_%dd" % (var_name, ddf_op, ndims)
        if len(self.atts.numSamples) == 1:
            src_fname = "%s.ultra" % self.atts.ddfName
            des_fname = "%s.ult" % (self.atts.ddfName)
            common.sexe("mv %s %s" % (src_fname, des_fname))
            lines = open(des_fname).readlines()
            f = open(des_fname, "w")
            f.write("# %s\n" % (ddf_oname))
            for l in lines[1:]:
                f.write(l)
            f.close()
        else:
            src_fname = "%s.vtk" % self.atts.ddfName
            orig_vtk_var = "SCALARS %s float" % var_name
            ddf_vtk_var = "SCALARS %s float" % ddf_oname
            des_fname = "%s%s_%04d.vtk" % (Params.output.file_base, self.atts.ddfName, ts())
            common.sexe("mv %s %s" % (src_fname, des_fname))
            data = open(des_fname).read()
            f = open(des_fname, "w")
            data = data.replace(orig_vtk_var, ddf_vtk_var)
            f.write(data)
        print "[ddf output: %s]" % des_fname
        return des_fname




def average_xy(fld_name, num_samples, ts = None):
    """
    averages a given field in the xy plane
    """

    #  ts : time slider
    if ts is not None:
        ddf_name = "%s_avg_xy_%04d" % (fld_name, ts)
    else:
        ddf_name = "%s_avg_xy" % (fld_name)

    # create the reduction object with the right operator
    fld_bar = ReductionOps(op='avg', fld_name=fld_name)

    # plot Pseudo-color
    fld_bar.plotPseudocolor()

    # find the spatial extents of the comp. domain
    sext =fld_bar.mesh_spatial_extents()

    # set the remainder of the attributes
    # average in xy by binning based on z
    fld_bar.atts.varnames = ("mesh_z_nodal",)
    fld_bar.atts.ranges = (sext[4], sext[5])
    fld_bar.atts.ddfName = ddf_name
    fld_bar.atts.numSamples = (num_samples,)

    # construct the ddf once all attributes are set
    fld_bar.ddf()

    mesh_name = fld_bar.mesh_name
    avg_ename = "%s_avg" % fld_name
    exprs.define(avg_ename, "apply_ddf(%s,%s)" % (mesh_name, ddf_name))

    return avg_ename


def calc_ubar(fld_name, num_samples, time_slider=None):
    return average_xy(fld_name, num_samples, ts=time_slider)


def main():
    dbname = "~/Downloads/tutorial_data/varying.visit"
    launch_visit(visit_path)
    OpenDatabase(dbname)
    ubar = calc_ubar('temp', 500)
    ChangeActivePlotsVar(ubar)
    pass

# if __visit_script_file__ == __visit_source_file__:
if __name__ == "__main__":
    main()

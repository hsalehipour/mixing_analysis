import visit
from visit_utils import exprs, query, common

class VisitSetupBase(object):

    visit_path = "/usr/local/visit/bin"
    # visit_python_package = "/usr/local/visit/current/linux-x86_64/lib/site-packages"

    def __init__(self):
        self.mname = self.mesh_name()
        self.launch_visit(self.visit_path)

    def launch_visit(visit_path):
        return visit.Launch(vdir=visit_path)

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

    def mesh_name(self):
        """
        Fetches the first mesh name in the active database.
        """
        return self.mesh_md().GetMeshes(0).name

    def setup_exprs(self):
        """
        Sets up spatial expressions to use with DDF.
        """
        mname = self.mesh_name()
        exprs.define("mesh_x_nodal", "coord(%s)[0]" % mname)
        exprs.define("mesh_y_nodal", "coord(%s)[1]" % mname)
        exprs.define("mesh_z_nodal", "coord(%s)[2]" % mname)
        exprs.define("mesh_x_zonal", "recenter(coord(%s)[0])" % mname)
        exprs.define("mesh_y_zonal", "recenter(coord(%s)[1])" % mname)
        exprs.define("mesh_z_zonal", "recenter(coord(%s)[2])" % mname)

    def plotPseudocolor(self, var_name):
        visit.AddPlot("Pseudocolor", var_name)
        visit.DrawPlots()


class ReductionOps(VisitSetupBase):

    def __init__(self,):
        super(ReductionOps, self).__init__()
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

    def ddf(self, operation = "avg"):
        self.atts.statisticalOperator = self.ddf_op_map[operation]
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


    def calc_ubar(self, var_name, num_samples, time_slider=None):

        if time_slider is not None:
            ddf_name = "%s_avg_xy_%04d" % (var_name, time_slider)
        else:
            ddf_name = "%s_avg_xy" % (var_name)
        sext = self.mesh_spatial_extents()
        self.atts.ddfName = ddf_name
        self.atts.codomainName = var_name
        self.atts.varnames = ("mesh_z_nodal",)
        self.atts.ranges = (sext[4], sext[5])
        self.atts.numSamples = (num_samples,)
        self.ddf(operation="avg")

        mname = self.mesh_name()
        avg_ename = "%s_avg" % var_name
        exprs.define(avg_ename, "apply_ddf(%s,%s)" % (mname, ddf_name))

        return
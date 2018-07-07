from visit import *
from visit_utils import exprs, query, common
import numpy as np

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

    def plotSlice(self, var_name):
        # plot Pseudo-color
        visit.AddPlot("Pseudocolor", var_name)
        visit.AddOperator("Slice")
        visit.DrawPlots()




class ReductionOps(VisitSetupBase):

    def __init__(self, op, var_name):
        super(ReductionOps, self).__init__()
        self.ddf_op = op
        self.var_name = var_name
        self.__set_attributes__()

    def __set_attributes__(self):
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
        self.atts.statisticalOperator = self.ddf_op_map[self.ddf_op]
        self.atts.codomainName = self.var_name
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



class LinoutOps(VisitSetupBase):

    def __init__(self, point1, point2, var_name_list):
        super(LinoutOps, self).__init__()
        self.p1 = point1
        self.p2 = point2
        self.var_name_list = var_name_list
        self.plot_id_dic = {key: value for (key, value) in enumerate(var_name_list)}
        self.__create__()

    def __create__(self):
        visit.Lineout(self.p1, self.p2, self.var_name_list)
        return

    def update_curve(self, var_name_list, p1, p2, window_id=2):
        """
        Updates the curve by adjusting points of the lineout operator
        :param var_name_list: a list of variable names to be modified
        :param p1: point 1
        :param p2: point 2
        :param window_id  : the active window id (usually=2)
        """
        plot_id_list = self.plot_id_dic[var_name_list]
        visit.SetActiveWindow(window_id)
        visit.SetActivePlots(plot_id_list)
        atts = visit.LineoutAttributes()
        atts.point1 = p1
        atts.point2 = p2
        visit.SetOperatorOptions(atts)
        return

    @staticmethod
    def extract_curve(window_id, plot_id):
        """
        Extracts curve data in the form of numpy array
        :param window_id    = a scalar for the active window id where Lineout curves are plotted
        :param plot_id_list = a list of plots to
        """
        visit.SetActiveWindow(window_id)
        visit.SetActivePlots(plot_id)
        data = visit.GetPlotInformation()["Curve"]
        return np.array(data).reshape((len(data) / 2, -1))[:, 1]



def average_xy(var_name, num_samples, ts = None):
    """
    averages a given field in the xy plane
    """

    #  ts : time slider
    if ts is not None:
        ddf_name = "%s_avg_xy_%04d" % (var_name, ts)
    else:
        ddf_name = "%s_avg_xy" % (var_name)

    # create the reduction object with the right operator
    fld_bar = ReductionOps(op='avg', var_name=var_name)

    # find the spatial extents of the comp. domain
    sext =fld_bar.mesh_spatial_extents()

    # set the remainder of the attributes
    # average in xy by binning based on z
    fld_bar.atts.varnames = ("mesh_z_nodal",)
    fld_bar.atts.ranges = (sext[2], sext[3])
    fld_bar.atts.ddfName = ddf_name
    fld_bar.atts.numSamples = (num_samples,)

    # construct the ddf once all attributes are set
    fld_bar.ddf()

    mesh_name = fld_bar.mesh_name
    avg_ename = "%s_avg" % var_name
    exprs.define(avg_ename, "apply_ddf(%s,%s)" % (mesh_name, ddf_name))

    return avg_ename


def calc_ubar(ufld, num_samples, time_slider=None):
    ubar = average_xy(ufld, num_samples, ts=time_slider)
    return ubar

def calc_ufluct(ufld, ubar):
    """
    create the expression for u' = u - ubar where ubar = <u>_{xy}
    """
    u_fluct = "%s_fluct" % ufld
    exprs.define(u_fluct, "%s - %s" % (ufld, ubar))
    return u_fluct

def main():
    # dbname = "~/Downloads/tutorial_data/varying.visit"
    dbname  ="~/REPOs/mixing_analysis/holm_workspace/data/holm.nek5000"

    rho = 'temperature'
    launch_visit(visit_path)
    OpenDatabase(dbname)
    rho_bar   = calc_ubar(rho, 500)
    rho_fluct = calc_ufluct(rho, rho_bar)
    ChangeActivePlotsVar(rho_bar)
    ChangeActivePlotsVar(rho_fluct)

    pass

# if __visit_script_file__ == __visit_source_file__:
if __name__ == "__main__":
    main()

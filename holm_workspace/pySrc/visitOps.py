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

    def get_times(self):
        """
        Fetches the times in the active database
        """
        return self.mesh_md().times

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

    def __init__(self, line_dic, var_name_list):
        """
        :param line_list: a list of input lines including N couple of tuples, denoting the two points of a line.
        :param var_name_list: a list of variables for which curves are plotted.
        """
        super(LinoutOps, self).__init__()
        self.line_dic = line_dic
        self.var_name_list = var_name_list
        self.nvar   = len(var_name_list)
        self.nlines = len(line_dic)
        self.nplots = self.nvar * self.nlines
        self.__set_plotID__()


    def __set_plotID__(self):
        """
        Sets a mapping dictionary between line_id and a list of associated plot_ids.
        e.g. {'line01' : (plot_id1, plot_id2), 'line02': (plot_id3, plot_id4)}
        """
        line_names = self.line_dic.keys()
        plot_ids   = [ii for ii in range(self.nplots)]
        plot_ids  = [ tuple(plot_ids[ii:ii+self.nvar]) for ii in range(0, self.nplots, self.nvar) ]
        self.plotID_dic = dict(zip(line_names, plot_ids))

    def create(self):
        for ii in range(self.nlines):
            p1, p2 = self.line_dic.values()[ii]
            visit.Lineout(p1, p2, self.var_name_list)
        return

    def update(self, line_name, p1, p2, window_id=2):
        """
        Updates the curve by adjusting points of the lineout operator
        :param line_name: the "name" of the line whose points are to be updated.
        :param p1: point 1
        :param p2: point 2
        :param window_id  : the active window id (usually=2)
        """
        self.line_dic[line_name] = (p1, p2)
        plot_ids = self.plotID_dic[line_name]
        visit.SetActiveWindow(window_id)
        visit.SetActivePlots(plot_ids)
        atts = visit.LineoutAttributes()
        atts.point1 = p1
        atts.point2 = p2
        visit.SetOperatorOptions(atts)
        return

    def extract(self, line_name, var_name, window_id=2):
        """
        Extracts curve data in the form of numpy array
        :param line_name: "name" of the line where Lineout curves are plotted
        :param var_name : "name" of the variable to be extracted
        :param window_id : a scalar for the active window id where Lineout curves are plotted
        """
        plot_id = self.plotID_dic[line_name][self.var_name_list.index(var_name)]
        visit.SetActiveWindow(window_id)
        visit.SetActivePlots(plot_id)
        data = visit.GetPlotInformation()["Curve"]
        return np.array(data).reshape((-1, 2))[:, 1]



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

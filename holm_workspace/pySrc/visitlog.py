# Visit 2.10.0 log file
ScriptVersion = "2.10.0"
if ScriptVersion != Version():
    print "This script is for VisIt %s. It may not work with version %s" % (ScriptVersion, Version())
visit.ShowAllWindows()
OpenDatabase("~/Downloads/tutorial_data/varying.visit", 0)
# The UpdateDBPluginInfo RPC is not supported in the VisIt module so it will not be logged.
metadata = GetMetaData("localhost:/home/hesam/Downloads/tutorial_data/varying.visit", -1)
DefineScalarExpression("mesh_x_nodal", "coord(mesh)[0]")
DefineScalarExpression("mesh_x_nodal", "coord(mesh)[0]")
DefineScalarExpression("mesh_y_nodal", "coord(mesh)[1]")
DefineScalarExpression("mesh_x_nodal", "coord(mesh)[0]")
DefineScalarExpression("mesh_y_nodal", "coord(mesh)[1]")
DefineScalarExpression("mesh_z_nodal", "coord(mesh)[2]")
DefineScalarExpression("mesh_x_nodal", "coord(mesh)[0]")
DefineScalarExpression("mesh_y_nodal", "coord(mesh)[1]")
DefineScalarExpression("mesh_z_nodal", "coord(mesh)[2]")
DefineScalarExpression("mesh_x_zonal", "recenter(coord(mesh)[0])")
DefineScalarExpression("mesh_x_nodal", "coord(mesh)[0]")
DefineScalarExpression("mesh_y_nodal", "coord(mesh)[1]")
DefineScalarExpression("mesh_z_nodal", "coord(mesh)[2]")
DefineScalarExpression("mesh_x_zonal", "recenter(coord(mesh)[0])")
DefineScalarExpression("mesh_y_zonal", "recenter(coord(mesh)[1])")
DefineScalarExpression("mesh_x_nodal", "coord(mesh)[0]")
DefineScalarExpression("mesh_y_nodal", "coord(mesh)[1]")
DefineScalarExpression("mesh_z_nodal", "coord(mesh)[2]")
DefineScalarExpression("mesh_x_zonal", "recenter(coord(mesh)[0])")
DefineScalarExpression("mesh_y_zonal", "recenter(coord(mesh)[1])")
DefineScalarExpression("mesh_z_zonal", "recenter(coord(mesh)[2])")
AddPlot("Pseudocolor", "temp", 1, 1)
metadata = GetMetaData("localhost:/home/hesam/Downloads/tutorial_data/varying.visit", -1)
DrawPlots()
Query("SpatialExtents")
ConstructDataBinningAtts = ConstructDataBinningAttributes()
ConstructDataBinningAtts.name = "ddf_temp"
ConstructDataBinningAtts.varnames = ("mesh_z_nodal")
ConstructDataBinningAtts.binType = ()
ConstructDataBinningAtts.binBoundaries = (-10, 10)
ConstructDataBinningAtts.reductionOperator = ConstructDataBinningAtts.Average  # Average, Minimum, Maximum, StandardDeviation, Variance, Sum, Count, RMS, PDF
ConstructDataBinningAtts.varForReductionOperator = "temp"
ConstructDataBinningAtts.undefinedValue = 0
ConstructDataBinningAtts.binningScheme = ConstructDataBinningAtts.Uniform  # Uniform, Unknown
ConstructDataBinningAtts.numBins = (50)
ConstructDataBinningAtts.overTime = 0
ConstructDataBinningAtts.timeStart = 0
ConstructDataBinningAtts.timeEnd = 1
ConstructDataBinningAtts.timeStride = 1
ConstructDataBinningAtts.outOfBoundsBehavior = ConstructDataBinningAtts.Clamp  # Clamp, Discard
ConstructDataBinning(ConstructDataBinningAtts)
DefineScalarExpression("mesh_x_nodal", "coord(mesh)[0]")
DefineScalarExpression("mesh_y_nodal", "coord(mesh)[1]")
DefineScalarExpression("mesh_z_nodal", "coord(mesh)[2]")
DefineScalarExpression("mesh_x_zonal", "recenter(coord(mesh)[0])")
DefineScalarExpression("mesh_y_zonal", "recenter(coord(mesh)[1])")
DefineScalarExpression("mesh_z_zonal", "recenter(coord(mesh)[2])")
DefineScalarExpression("temp_avg", "apply_ddf(mesh,ddf_temp)")
DefineScalarExpression("mesh_x_nodal", "coord(mesh)[0]")
DefineScalarExpression("mesh_y_nodal", "coord(mesh)[1]")
DefineScalarExpression("mesh_z_nodal", "coord(mesh)[2]")
DefineScalarExpression("mesh_x_zonal", "recenter(coord(mesh)[0])")
DefineScalarExpression("mesh_y_zonal", "recenter(coord(mesh)[1])")
DefineScalarExpression("mesh_z_zonal", "recenter(coord(mesh)[2])")
DefineScalarExpression("temp_avg", "apply_ddf(mesh,ddf_temp)")
DefineScalarExpression("temp_fluct", "temp - temp_avg")
ChangeActivePlotsVar("temp_fluct")
# Begin spontaneous state
View3DAtts = View3DAttributes()
View3DAtts.viewNormal = (0.424362, 0.403224, 0.810757)
View3DAtts.focus = (0, 0, 0)
View3DAtts.viewUp = (-0.0330915, 0.901685, -0.431126)
View3DAtts.viewAngle = 30
View3DAtts.parallelScale = 17.3205
View3DAtts.nearPlane = -34.641
View3DAtts.farPlane = 34.641
View3DAtts.imagePan = (0, 0)
View3DAtts.imageZoom = 1
View3DAtts.perspective = 1
View3DAtts.eyeAngle = 2
View3DAtts.centerOfRotationSet = 0
View3DAtts.centerOfRotation = (0, 0, 0)
View3DAtts.axis3DScaleFlag = 0
View3DAtts.axis3DScales = (1, 1, 1)
View3DAtts.shear = (0, 0, 1)
View3DAtts.windowValid = 1
SetView3D(View3DAtts)
# End spontaneous state

ChangeActivePlotsVar("temp_avg")
# Begin spontaneous state
View3DAtts = View3DAttributes()
View3DAtts.viewNormal = (0.141788, 0.98217, 0.12344)
View3DAtts.focus = (0, 0, 0)
View3DAtts.viewUp = (0.381317, 0.0608849, -0.922437)
View3DAtts.viewAngle = 30
View3DAtts.parallelScale = 17.3205
View3DAtts.nearPlane = -34.641
View3DAtts.farPlane = 34.641
View3DAtts.imagePan = (0, 0)
View3DAtts.imageZoom = 1
View3DAtts.perspective = 1
View3DAtts.eyeAngle = 2
View3DAtts.centerOfRotationSet = 0
View3DAtts.centerOfRotation = (0, 0, 0)
View3DAtts.axis3DScaleFlag = 0
View3DAtts.axis3DScales = (1, 1, 1)
View3DAtts.shear = (0, 0, 1)
View3DAtts.windowValid = 1
SetView3D(View3DAtts)
# End spontaneous state

# Begin spontaneous state
View3DAtts = View3DAttributes()
View3DAtts.viewNormal = (-0.754331, -0.259737, 0.602927)
View3DAtts.focus = (0, 0, 0)
View3DAtts.viewUp = (-0.503658, 0.818045, -0.277726)
View3DAtts.viewAngle = 30
View3DAtts.parallelScale = 17.3205
View3DAtts.nearPlane = -34.641
View3DAtts.farPlane = 34.641
View3DAtts.imagePan = (0, 0)
View3DAtts.imageZoom = 1
View3DAtts.perspective = 1
View3DAtts.eyeAngle = 2
View3DAtts.centerOfRotationSet = 0
View3DAtts.centerOfRotation = (0, 0, 0)
View3DAtts.axis3DScaleFlag = 0
View3DAtts.axis3DScales = (1, 1, 1)
View3DAtts.shear = (0, 0, 1)
View3DAtts.windowValid = 1
SetView3D(View3DAtts)
# End spontaneous state

# Begin spontaneous state
View3DAtts = View3DAttributes()
View3DAtts.viewNormal = (0.0559778, 0.562183, 0.825116)
View3DAtts.focus = (0, 0, 0)
View3DAtts.viewUp = (-0.975406, 0.207252, -0.0750347)
View3DAtts.viewAngle = 30
View3DAtts.parallelScale = 17.3205
View3DAtts.nearPlane = -34.641
View3DAtts.farPlane = 34.641
View3DAtts.imagePan = (0, 0)
View3DAtts.imageZoom = 1
View3DAtts.perspective = 1
View3DAtts.eyeAngle = 2
View3DAtts.centerOfRotationSet = 0
View3DAtts.centerOfRotation = (0, 0, 0)
View3DAtts.axis3DScaleFlag = 0
View3DAtts.axis3DScales = (1, 1, 1)
View3DAtts.shear = (0, 0, 1)
View3DAtts.windowValid = 1
SetView3D(View3DAtts)
# End spontaneous state

# Begin spontaneous state
View3DAtts = View3DAttributes()
View3DAtts.viewNormal = (0.932607, -0.0344767, 0.359244)
View3DAtts.focus = (0, 0, 0)
View3DAtts.viewUp = (-0.241705, 0.679532, 0.692687)
View3DAtts.viewAngle = 30
View3DAtts.parallelScale = 17.3205
View3DAtts.nearPlane = -34.641
View3DAtts.farPlane = 34.641
View3DAtts.imagePan = (0, 0)
View3DAtts.imageZoom = 1
View3DAtts.perspective = 1
View3DAtts.eyeAngle = 2
View3DAtts.centerOfRotationSet = 0
View3DAtts.centerOfRotation = (0, 0, 0)
View3DAtts.axis3DScaleFlag = 0
View3DAtts.axis3DScales = (1, 1, 1)
View3DAtts.shear = (0, 0, 1)
View3DAtts.windowValid = 1
SetView3D(View3DAtts)
# End spontaneous state

# Begin spontaneous state
View3DAtts = View3DAttributes()
View3DAtts.viewNormal = (0.716695, -0.697345, 0.0076491)
View3DAtts.focus = (0, 0, 0)
View3DAtts.viewUp = (0.10287, 0.11656, 0.987842)
View3DAtts.viewAngle = 30
View3DAtts.parallelScale = 17.3205
View3DAtts.nearPlane = -34.641
View3DAtts.farPlane = 34.641
View3DAtts.imagePan = (0, 0)
View3DAtts.imageZoom = 1
View3DAtts.perspective = 1
View3DAtts.eyeAngle = 2
View3DAtts.centerOfRotationSet = 0
View3DAtts.centerOfRotation = (0, 0, 0)
View3DAtts.axis3DScaleFlag = 0
View3DAtts.axis3DScales = (1, 1, 1)
View3DAtts.shear = (0, 0, 1)
View3DAtts.windowValid = 1
SetView3D(View3DAtts)
# End spontaneous state


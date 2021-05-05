from visit_utils import *

OpenDatabase("localhost:/Users/mnataraj/Desktop/Research/PeleC_repos/PeleC_Moving/Submodules/AMReX/Tutorials/EB/CNS_EBoft_ShockAndExpansion/Exec/AcceleratingPiston/movie.visit", 0)
AddPlot("Pseudocolor", "x_velocity", 1, 0)

for i in range(0,78,1):
	DrawPlots()
	Query("Lineout", end_point=(0.9, 0, 0), num_samples=200, start_point=(-0.1, 0, 0), use_sampling=0)
	SetActiveWindow(2)
	SetActiveWindow(2)
	SaveWindowAtts = SaveWindowAttributes()
	SaveWindowAtts.outputToCurrentDirectory = 0
	SaveWindowAtts.outputDirectory = "/Users/mnataraj/Desktop/Research/PeleC_repos/PeleC_Moving/Submodules/AMReX/Tutorials/EB/CNS_EBoft_ShockAndExpansion/Exec/AcceleratingPiston/Comparison/Shock"
	filename="temp%04d"%i
	SaveWindowAtts.fileName = filename
	SaveWindowAtts.family = 0
	SaveWindowAtts.format = SaveWindowAtts.CURVE  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
	SaveWindowAtts.width = 1024
	SaveWindowAtts.height = 1024
	SaveWindowAtts.screenCapture = 0
	SaveWindowAtts.saveTiled = 0
	SaveWindowAtts.quality = 80
	SaveWindowAtts.progressive = 0
	SaveWindowAtts.binary = 0
	SaveWindowAtts.stereo = 0
	SaveWindowAtts.compression = SaveWindowAtts.None  # None, PackBits, Jpeg, Deflate
	SaveWindowAtts.forceMerge = 0
	SaveWindowAtts.resConstraint = SaveWindowAtts.ScreenProportions  # NoConstraint, EqualWidthHeight, ScreenProportions
	SaveWindowAtts.advancedMultiWindowSave = 0
	SaveWindowAtts.subWindowAtts.win1.position = (0, 0)
	SaveWindowAtts.subWindowAtts.win1.size = (128, 128)
	SaveWindowAtts.subWindowAtts.win1.layer = 0
	SaveWindowAtts.subWindowAtts.win1.transparency = 0
	SaveWindowAtts.subWindowAtts.win1.omitWindow = 0
	SaveWindowAtts.subWindowAtts.win2.position = (0, 0)
	SaveWindowAtts.subWindowAtts.win2.size = (128, 128)
	SaveWindowAtts.subWindowAtts.win2.layer = 0
	SaveWindowAtts.subWindowAtts.win2.transparency = 0
	SaveWindowAtts.subWindowAtts.win2.omitWindow = 0
	SaveWindowAtts.subWindowAtts.win3.position = (0, 0)
	SaveWindowAtts.subWindowAtts.win3.size = (128, 128)
	SaveWindowAtts.subWindowAtts.win3.layer = 0
	SaveWindowAtts.subWindowAtts.win3.transparency = 0
	SaveWindowAtts.subWindowAtts.win3.omitWindow = 0
	SaveWindowAtts.subWindowAtts.win4.position = (0, 0)
	SaveWindowAtts.subWindowAtts.win4.size = (128, 128)
	SaveWindowAtts.subWindowAtts.win4.layer = 0
	SaveWindowAtts.subWindowAtts.win4.transparency = 0
	SaveWindowAtts.subWindowAtts.win4.omitWindow = 0
	SaveWindowAtts.subWindowAtts.win5.position = (0, 0)
	SaveWindowAtts.subWindowAtts.win5.size = (128, 128)
	SaveWindowAtts.subWindowAtts.win5.layer = 0
	SaveWindowAtts.subWindowAtts.win5.transparency = 0
	SaveWindowAtts.subWindowAtts.win5.omitWindow = 0
	SaveWindowAtts.subWindowAtts.win6.position = (0, 0)
	SaveWindowAtts.subWindowAtts.win6.size = (128, 128)
	SaveWindowAtts.subWindowAtts.win6.layer = 0
	SaveWindowAtts.subWindowAtts.win6.transparency = 0
	SaveWindowAtts.subWindowAtts.win6.omitWindow = 0
	SaveWindowAtts.subWindowAtts.win7.position = (0, 0)
	SaveWindowAtts.subWindowAtts.win7.size = (128, 128)
	SaveWindowAtts.subWindowAtts.win7.layer = 0
	SaveWindowAtts.subWindowAtts.win7.transparency = 0
	SaveWindowAtts.subWindowAtts.win7.omitWindow = 0
	SaveWindowAtts.subWindowAtts.win8.position = (0, 0)
	SaveWindowAtts.subWindowAtts.win8.size = (128, 128)
	SaveWindowAtts.subWindowAtts.win8.layer = 0
	SaveWindowAtts.subWindowAtts.win8.transparency = 0
	SaveWindowAtts.subWindowAtts.win8.omitWindow = 0
	SaveWindowAtts.subWindowAtts.win9.position = (0, 0)
	SaveWindowAtts.subWindowAtts.win9.size = (128, 128)
	SaveWindowAtts.subWindowAtts.win9.layer = 0
	SaveWindowAtts.subWindowAtts.win9.transparency = 0
	SaveWindowAtts.subWindowAtts.win9.omitWindow = 0
	SaveWindowAtts.subWindowAtts.win10.position = (0, 0)
	SaveWindowAtts.subWindowAtts.win10.size = (128, 128)
	SaveWindowAtts.subWindowAtts.win10.layer = 0
	SaveWindowAtts.subWindowAtts.win10.transparency = 0
	SaveWindowAtts.subWindowAtts.win10.omitWindow = 0
	SaveWindowAtts.subWindowAtts.win11.position = (0, 0)
	SaveWindowAtts.subWindowAtts.win11.size = (128, 128)
	SaveWindowAtts.subWindowAtts.win11.layer = 0
	SaveWindowAtts.subWindowAtts.win11.transparency = 0
	SaveWindowAtts.subWindowAtts.win11.omitWindow = 0
	SaveWindowAtts.subWindowAtts.win12.position = (0, 0)
	SaveWindowAtts.subWindowAtts.win12.size = (128, 128)
	SaveWindowAtts.subWindowAtts.win12.layer = 0
	SaveWindowAtts.subWindowAtts.win12.transparency = 0
	SaveWindowAtts.subWindowAtts.win12.omitWindow = 0
	SaveWindowAtts.subWindowAtts.win13.position = (0, 0)
	SaveWindowAtts.subWindowAtts.win13.size = (128, 128)
	SaveWindowAtts.subWindowAtts.win13.layer = 0
	SaveWindowAtts.subWindowAtts.win13.transparency = 0
	SaveWindowAtts.subWindowAtts.win13.omitWindow = 0
	SaveWindowAtts.subWindowAtts.win14.position = (0, 0)
	SaveWindowAtts.subWindowAtts.win14.size = (128, 128)
	SaveWindowAtts.subWindowAtts.win14.layer = 0
	SaveWindowAtts.subWindowAtts.win14.transparency = 0
	SaveWindowAtts.subWindowAtts.win14.omitWindow = 0
	SaveWindowAtts.subWindowAtts.win15.position = (0, 0)
	SaveWindowAtts.subWindowAtts.win15.size = (128, 128)
	SaveWindowAtts.subWindowAtts.win15.layer = 0
	SaveWindowAtts.subWindowAtts.win15.transparency = 0
	SaveWindowAtts.subWindowAtts.win15.omitWindow = 0
	SaveWindowAtts.subWindowAtts.win16.position = (0, 0)
	SaveWindowAtts.subWindowAtts.win16.size = (128, 128)
	SaveWindowAtts.subWindowAtts.win16.layer = 0
	SaveWindowAtts.subWindowAtts.win16.transparency = 0
	SaveWindowAtts.subWindowAtts.win16.omitWindow = 0
	SetSaveWindowAttributes(SaveWindowAtts)
	SaveWindow()
	TimeSliderNextState()
	
